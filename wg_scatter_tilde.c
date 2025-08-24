/* wg_scatter_tilde.c
   Implements Pd external: wg_scatter_tilde~  (file: wg_scatter_tilde~.pd_darwin)
   Setup symbol: wg_scatter_tilde_tilde_setup
   Optional alias: wg_scatter~

   Purpose:
   - Safely couple N waveguide lines (e.g., wg_line_tilde~ "couple_send" outs)
     back to their "couple_in" inlets using physically-inspired diffusive scattering.
   - NO DSP LOOPS: uses 1-sample history per branch (z^-1) so Pd won't complain.
   - Global signal inlets: coupling_all~ (g 0..1), couple_loss~ (L 0..1).
   - Per-branch weights & gains via messages: weight <i> <v>, gain <i> <v>.
     (i is 0-based). Also: set_all_weights <v>, set_all_gains <v>.

   Math (per sample, using previous-sample inputs x_prev[i]):
     S  = sum_j (w[j] * x_prev[j]),   W = sum_j w[j]
     mean_all = (W > 0 ? S/W : 0)
     mean_ex_i = (W - w[i] > 1e-9 ? (S - w[i]*x_prev[i]) / (W - w[i]) : 0)
     y[i] = (1 - L) * [ (1 - g) * x_prev[i] + g * mean_ex_i ] * gain[i]

   This conserves energy for L=0 and uniform weights/gains when 0<=g<=1.

   Creation:
     [wg_scatter_tilde~ N]
     N = number of coupled branches (2..16). Defaults to 2.

   Signal inlets (total N + 2):
     1..N: input signals from each branch's couple_send~
     N+1 : coupling_all~ (g 0..1)
     N+2 : couple_loss~  (L 0..1)

   Signal outlets (N):
     1..N: output signals to feed each branch's couple_in~
*/
#include "m_pd.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define WGSC_MAXN 16

static inline t_float clampf(t_float x, t_float lo, t_float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
static inline t_float zapgremlins(t_float x) {
    if (!isfinite(x) || fabsf(x) < 1e-30f) return 0.f;
    return x;
}

typedef struct _wg_scatter_tilde {
    t_object x_obj;

    int N;
    t_inlet  **sig_in;
    t_outlet **sig_out;

    /* state: previous-sample inputs to break DSP loop */
    t_float *x_prev;

    /* per-branch controls via messages */
    t_float *weight;  /* contribution weight into the shared mean */
    t_float *gain;    /* per-output gain after mixing */

} t_wg_scatter_tilde;

static t_class *wg_scatter_tilde_tilde_class;

static void *wg_scatter_tilde_new(t_floatarg fN) {
    int N = (int)fN;
    if (N < 2) N = 2;
    if (N > WGSC_MAXN) N = WGSC_MAXN;

    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)pd_new(wg_scatter_tilde_tilde_class);
    x->N = N;
    x->sig_in  = (t_inlet  **)getbytes((N + 2) * sizeof(t_inlet *));  /* N inputs + 2 globals */
    x->sig_out = (t_outlet **)getbytes(N * sizeof(t_outlet *));
    x->x_prev  = (t_float *)getbytes(N * sizeof(t_float));
    x->weight  = (t_float *)getbytes(N * sizeof(t_float));
    x->gain    = (t_float *)getbytes(N * sizeof(t_float));

    if (!x->sig_in || !x->sig_out || !x->x_prev || !x->weight || !x->gain) {
        post("wg_scatter_tilde~: memory alloc failed");
        return x;
    }

    for (int i = 0; i < N; ++i) {
        x->sig_in[i] = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
        x->sig_out[i] = outlet_new(&x->x_obj, &s_signal);
        x->x_prev[i] = 0.f;
        x->weight[i] = 1.f;
        x->gain[i] = 1.f;
    }
    /* global signal inlets: coupling_all~ and couple_loss~ */
    x->sig_in[N]   = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* g */
    x->sig_in[N+1] = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* L */

    return x;
}

static void wg_scatter_tilde_free(t_wg_scatter_tilde *x) {
    if (x->sig_in)  freebytes(x->sig_in,  (x->N + 2) * sizeof(t_inlet *));
    if (x->sig_out) freebytes(x->sig_out, x->N * sizeof(t_outlet *));
    if (x->x_prev)  freebytes(x->x_prev,  x->N * sizeof(t_float));
    if (x->weight)  freebytes(x->weight,  x->N * sizeof(t_float));
    if (x->gain)    freebytes(x->gain,    x->N * sizeof(t_float));
}

/* messages: weight <i> <v>, gain <i> <v>, set_all_weights <v>, set_all_gains <v> */
static void wg_scatter_tilde_weight(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg v) {
    int i = (int)idx;
    if (i < 0 || i >= x->N) return;
    x->weight[i] = (t_float)(v);
}
static void wg_scatter_tilde_gain(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg v) {
    int i = (int)idx;
    if (i < 0 || i >= x->N) return;
    x->gain[i] = (t_float)(v);
}
static void wg_scatter_tilde_set_all_weights(t_wg_scatter_tilde *x, t_floatarg v) {
    for (int i = 0; i < x->N; ++i) x->weight[i] = (t_float)v;
}
static void wg_scatter_tilde_set_all_gains(t_wg_scatter_tilde *x, t_floatarg v) {
    for (int i = 0; i < x->N; ++i) x->gain[i] = (t_float)v;
}

static t_int *wg_scatter_tilde_perform(t_int *w) {
    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)(w[1]);
    int n = (int)(w[2]);
    int N = x->N;

    /* Build pointers to input vectors: first N inputs, then g, then L */
    t_sample **in = (t_sample **)(w + 3);
    t_sample *in_g = in[N];     /* coupling_all~ */
    t_sample *in_L = in[N + 1]; /* couple_loss~ */

    /* Output vectors (N of them) follow after inputs (N+2) */
    t_sample **out = (t_sample **)(w + 3 + (N + 2));

    for (int i = 0; i < n; ++i) {
        /* read current sample and compute previous-sample mean */
        t_float S = 0.f, W = 0.f;
        for (int k = 0; k < N; ++k) {
            S += x->weight[k] * x->x_prev[k];
            W += x->weight[k];
        }

        t_float g = clampf(in_g[i], 0.f, 1.f);
        t_float L = clampf(in_L[i], 0.f, 1.f);

        for (int k = 0; k < N; ++k) {
            t_float wk = x->weight[k];
            t_float denom = W - wk;
            t_float mean_ex = (denom > 1e-9f) ? (S - wk * x->x_prev[k]) / denom : 0.f;
            t_float yk = (1.f - L) * ( (1.f - g) * x->x_prev[k] + g * mean_ex );
            out[k][i] = zapgremlins( yk * x->gain[k] );
        }

        /* update history with CURRENT input sample AFTER producing outputs */
        for (int k = 0; k < N; ++k) {
            x->x_prev[k] = in[k][i];
        }
    }

    return (w + 3 + (N + 2) + N);
}

static void wg_scatter_tilde_dsp(t_wg_scatter_tilde *x, t_signal **sp) {
    /* We need to build a dsp_addv call with variable (N) arguments:
       Args: x, n, N input vecs, g vec, L vec, N output vecs
    */
    int N = x->N;
    int n_sig_in  = N + 2;
    int n_sig_out = N;
    int nargs = 3 + n_sig_in + n_sig_out;
    t_int *vec = (t_int *)getbytes(nargs * sizeof(t_int));
    if (!vec) return;

    vec[0] = (t_int)x;
    vec[1] = (t_int)sp[0]->s_n;
    int w = 2;

    /* inputs */
    for (int i = 0; i < N; ++i) {
        vec[++w] = (t_int)sp[i]->s_vec;
    }
    vec[++w] = (t_int)sp[N]->s_vec;     /* g */
    vec[++w] = (t_int)sp[N+1]->s_vec;   /* L */

    /* outputs */
    for (int o = 0; o < N; ++o) {
        vec[++w] = (t_int)sp[N+2+o]->s_vec;
    }

    dsp_addv(wg_scatter_tilde_perform, nargs, vec);
    freebytes(vec, nargs * sizeof(t_int));
}

void wg_scatter_tilde_tilde_setup(void) {
    wg_scatter_tilde_tilde_class = class_new(gensym("wg_scatter_tilde~"),
                                             (t_newmethod)wg_scatter_tilde_new,
                                             (t_method)wg_scatter_tilde_free,
                                             sizeof(t_wg_scatter_tilde),
                                             CLASS_DEFAULT,
                                             A_DEFFLOAT, 0);

    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_weight, gensym("weight"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_gain, gensym("gain"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_set_all_weights, gensym("set_all_weights"), A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_set_all_gains, gensym("set_all_gains"), A_FLOAT, 0);

    CLASS_MAINSIGNALIN(wg_scatter_tilde_tilde_class, t_wg_scatter_tilde, N); /* dummy float for left inlet */

    /* optional alias */
    class_addcreator((t_newmethod)wg_scatter_tilde_new, gensym("wg_scatter~"), A_DEFFLOAT, 0);
}
