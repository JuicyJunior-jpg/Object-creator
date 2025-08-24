/* wg_scatter_tilde.c (SAFE)
   Pd external: wg_scatter_tilde~  (alias: wg_scatter~)
   Setup: wg_scatter_tilde_tilde_setup
   - Fixed-dsp_add version (no var-arg packing), caches vectors in object.
   - 1-sample history per input to break feedback loops cleanly.
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
    t_float  x_f; /* dummy for mainsignalin float */

    int N;
    t_inlet  **sig_in;
    t_outlet **sig_out;

    t_sample *invec[WGSC_MAXN + 2];
    t_sample *outvec[WGSC_MAXN];

    t_float *x_prev;
    t_float *weight;
    t_float *gain;
} t_wg_scatter_tilde;

static t_class *wg_scatter_tilde_tilde_class;

static void *wg_scatter_tilde_new(t_floatarg fN) {
    int N = (int)fN; if (N < 2) N = 2; if (N > WGSC_MAXN) N = WGSC_MAXN;

    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)pd_new(wg_scatter_tilde_tilde_class);
    x->x_f = 0.f;
    x->N = N;

    x->sig_in  = (t_inlet  **)getbytes((N + 2) * sizeof(t_inlet *));
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
        x->invec[i] = 0;
        x->outvec[i] = 0;
    }
    x->sig_in[N]   = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* g */
    x->sig_in[N+1] = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* L */
    x->invec[N] = x->invec[N+1] = 0;

    return x;
}

static void wg_scatter_tilde_free(t_wg_scatter_tilde *x) {
    if (x->sig_in)  freebytes(x->sig_in,  (x->N + 2) * sizeof(t_inlet *));
    if (x->sig_out) freebytes(x->sig_out, x->N * sizeof(t_outlet *));
    if (x->x_prev)  freebytes(x->x_prev,  x->N * sizeof(t_float));
    if (x->weight)  freebytes(x->weight,  x->N * sizeof(t_float));
    if (x->gain)    freebytes(x->gain,    x->N * sizeof(t_float));
}

static void wg_scatter_tilde_weight(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg v) {
    int i = (int)idx; if (i < 0 || i >= x->N) return; x->weight[i] = (t_float)v;
}
static void wg_scatter_tilde_gain(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg v) {
    int i = (int)idx; if (i < 0 || i >= x->N) return; x->gain[i] = (t_float)v;
}
static void wg_scatter_tilde_set_all_weights(t_wg_scatter_tilde *x, t_floatarg v) {
    for (int i=0;i<x->N;++i) x->weight[i]=(t_float)v;
}
static void wg_scatter_tilde_set_all_gains(t_wg_scatter_tilde *x, t_floatarg v) {
    for (int i=0;i<x->N;++i) x->gain[i]=(t_float)v;
}

static t_int *wg_scatter_tilde_perform(t_int *w) {
    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)(w[1]);
    int n = (int)(w[2]);
    int N = x->N;

    t_sample **in = x->invec;
    t_sample **out = x->outvec;
    t_sample *in_g = in[N];
    t_sample *in_L = in[N + 1];

    for (int i = 0; i < n; ++i) {
        t_float S = 0.f, W = 0.f;
        for (int k = 0; k < N; ++k) { S += x->weight[k] * x->x_prev[k]; W += x->weight[k]; }

        t_float g = clampf(in_g ? in_g[i] : 0.f, 0.f, 1.f);
        t_float L = clampf(in_L ? in_L[i] : 0.f, 0.f, 1.f);

        for (int k = 0; k < N; ++k) {
            t_float wk = x->weight[k];
            t_float denom = W - wk;
            t_float mean_ex = (denom > 1e-9f) ? (S - wk * x->x_prev[k]) / denom : 0.f;
            t_float yk = (1.f - L) * ( (1.f - g) * x->x_prev[k] + g * mean_ex );
            if (out[k]) out[k][i] = zapgremlins( yk * x->gain[k] );
        }
        for (int k = 0; k < N; ++k) {
            t_sample *ink = in[k]; x->x_prev[k] = ink ? ink[i] : 0.f;
        }
    }
    return (w + 3);
}

static void wg_scatter_tilde_dsp(t_wg_scatter_tilde *x, t_signal **sp) {
    int N = x->N;
    for (int i=0;i<N;++i) { x->invec[i] = sp[i]->s_vec; x->outvec[i] = sp[N+2+i]->s_vec; }
    x->invec[N]   = sp[N]->s_vec;
    x->invec[N+1] = sp[N+1]->s_vec;
    dsp_add(wg_scatter_tilde_perform, 2, x, sp[0]->s_n);
}

void wg_scatter_tilde_tilde_setup(void) {
    wg_scatter_tilde_tilde_class = class_new(gensym("wg_scatter_tilde~"),
                                             (t_newmethod)wg_scatter_tilde_new,
                                             (t_method)wg_scatter_tilde_free,
                                             sizeof(t_wg_scatter_tilde),
                                             CLASS_DEFAULT, A_DEFFLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_weight, gensym("weight"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_gain, gensym("gain"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_set_all_weights, gensym("set_all_weights"), A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_set_all_gains, gensym("set_all_gains"), A_FLOAT, 0);
    CLASS_MAINSIGNALIN(wg_scatter_tilde_tilde_class, t_wg_scatter_tilde, x_f);
    class_addcreator((t_newmethod)wg_scatter_tilde_new, gensym("wg_scatter~"), A_DEFFLOAT, 0);
}
