// wg_scatter_tilde.c
// Implements Pd external: wg_scatter_tilde~  (file: wg_scatter_tilde~.pd_darwin)
// Setup symbol: wg_scatter_tilde_tilde_setup
// Optional alias: wg_scatter~
#include "m_pd.h"
#include <math.h>
#include <string.h>

#define WG_SCATTER_MAXPORTS 32

typedef struct _wg_scatter_tilde {
    t_object x_obj;

    int N; /* number of ports */

    /* DSP vectors: stored per spec (do not pass malloc'd temps to perform) */
    t_sample **invecs;
    t_sample **outvecs;

    /* per-port state */
    unsigned char enabled[WG_SCATTER_MAXPORTS]; /* 0/1 */
    t_float       Y[WG_SCATTER_MAXPORTS];       /* admittance (1/Z) */
    t_float       c_override[WG_SCATTER_MAXPORTS]; /* per-port coupling [0..1], <0 = none */

    /* global controls */
    t_float c_all;   /* [0..1] */
    t_float loss;    /* l in [0..0.99] */

    /* built-in one-sample output delay (z^-1) to break algebraic loop in Pd graph */
    t_sample *y_prev; /* length N, previous-sample outputs */
} t_wg_scatter_tilde;

static t_class *wg_scatter_tilde_tilde_class;

/* clamp helper */
static inline t_float clampf(t_float x, t_float lo, t_float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

/* perform routine: energy-conserving N-port junction with open-circuit disables.
   Outputs are delayed by one sample (y_prev) to ensure graph causality in Pd. */
static t_int *wg_scatter_tilde_perform(t_int *w) {
    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)(w[1]);
    int n = (int)(w[2]);

    int N = x->N;
    t_sample **invecs = x->invecs;
    t_sample **outvecs = x->outvecs;

    t_float loss = x->loss;
    if (loss < 0.f) loss = 0.f;
    if (loss > 0.99f) loss = 0.99f;

    for (int i = 0; i < n; ++i) {
        /* sums for enabled ports */
        t_float sumY = 0.f;
        t_float sumYa = 0.f;

        for (int k = 0; k < N; ++k) {
            if (x->enabled[k]) {
                t_float a = invecs[k][i];
                t_float Yk = x->Y[k];
                sumY  += Yk;
                sumYa += Yk * a;
            }
        }

        t_float p = 0.f;
        if (sumY > 1e-20f) p = sumYa / sumY;

        /* compute current-sample outputs, but write the PREVIOUS sample to the graph */
        for (int k = 0; k < N; ++k) {
            t_float y_curr = 0.f;
            if (x->enabled[k]) {
                t_float a = invecs[k][i];
                t_float b = 2.f * p - a; /* physical outgoing */
                t_float c_eff = (x->c_override[k] >= 0.f)
                                  ? clampf(x->c_override[k], 0.f, 1.f)
                                  : clampf(x->c_all, 0.f, 1.f);
                /* y = (1 - l) * ((1 - c)a + c b) */
                y_curr = (1.f - loss) * ( (1.f - c_eff) * a + c_eff * b );
            } else {
                y_curr = 0.f; /* open circuit */
            }

            /* write previous sample to output to break algebraic loop */
            outvecs[k][i] = x->y_prev[k];
            /* store current for next sample */
            x->y_prev[k] = y_curr;
        }
    }

    return (w + 3);
}

static void wg_scatter_tilde_dsp(t_wg_scatter_tilde *x, t_signal **sp) {
    /* store vectors as required: first N inlets, then N outlets */
    for (int k = 0; k < x->N; ++k) {
        x->invecs[k]  = sp[k]->s_vec;
        x->outvecs[k] = sp[k + x->N]->s_vec;
    }
    dsp_add(wg_scatter_tilde_perform, 2, x, sp[0]->s_n);
}

/* messages */

static void wg_scatter_tilde_enable(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg on) {
    int i = (int)idx;
    if (i < 1 || i > x->N) return;
    x->enabled[i - 1] = (on != 0.f) ? 1 : 0;
}

static void wg_scatter_tilde_coupling_all(t_wg_scatter_tilde *x, t_floatarg c) {
    x->c_all = clampf((t_float)c, 0.f, 1.f);
}

static void wg_scatter_tilde_coupling(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg c) {
    int i = (int)idx;
    if (i < 1 || i > x->N) return;
    x->c_override[i - 1] = clampf((t_float)c, 0.f, 1.f);
}

static void wg_scatter_tilde_impedance_ratio(t_wg_scatter_tilde *x, t_floatarg idx, t_floatarg rho) {
    int i = (int)idx;
    if (i < 1 || i > x->N) return;
    t_float r = (t_float)rho;
    if (r < 0.01f) r = 0.01f;
    x->Y[i - 1] = 1.f / r; /* admittance */
}

static void wg_scatter_tilde_couple_loss(t_wg_scatter_tilde *x, t_floatarg l) {
    x->loss = clampf((t_float)l, 0.f, 0.99f); /* spec sweet-spot 0..0.3; clamp hard to 0.99 */
}

static void *wg_scatter_tilde_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    int N = 2;
    if (argc >= 1 && argv[0].a_type == A_FLOAT) {
        N = (int)atom_getfloat(argv);
    }
    if (N < 2) N = 2;
    if (N > WG_SCATTER_MAXPORTS) N = WG_SCATTER_MAXPORTS;

    t_wg_scatter_tilde *x = (t_wg_scatter_tilde *)pd_new(wg_scatter_tilde_tilde_class);
    x->N = N;

    /* allocate pointer tables using Pd memory */
    x->invecs  = (t_sample **)getbytes(N * sizeof(t_sample *));
    x->outvecs = (t_sample **)getbytes(N * sizeof(t_sample *));
    for (int i = 0; i < N; ++i) { x->invecs[i] = NULL; x->outvecs[i] = NULL; }

    /* per-port defaults */
    for (int i = 0; i < N; ++i) {
        x->enabled[i] = 1;
        x->Y[i] = 1.f;             /* default Zref ratio = 1 */
        x->c_override[i] = -1.f;   /* <0 means "use global" */
    }

    /* GLOBAL DEFAULTS — hearable by default */
    x->c_all = 1.f;   /* coupling ON by default */
    x->loss  = 0.05f; /* gentle junction loss */

    /* built-in 1-sample delay state */
    x->y_prev = (t_sample *)getbytes(N * sizeof(t_sample));
    for (int i = 0; i < N; ++i) x->y_prev[i] = 0.f;

    /* create signal inlets/outlets
       Leftmost inlet becomes signal via CLASS_MAINSIGNALIN.
       We therefore add N-1 additional signal inlets here.
    */
    for (int i = 1; i < N; ++i) {
        inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    }
    for (int i = 0; i < N; ++i) {
        outlet_new(&x->x_obj, &s_signal);
    }

    return (x);
}

static void wg_scatter_tilde_free(t_wg_scatter_tilde *x) {
    if (x->invecs)  freebytes(x->invecs,  x->N * sizeof(t_sample *));
    if (x->outvecs) freebytes(x->outvecs, x->N * sizeof(t_sample *));
    if (x->y_prev)  freebytes(x->y_prev,  x->N * sizeof(t_sample));
}

void wg_scatter_tilde_tilde_setup(void) {
    wg_scatter_tilde_tilde_class = class_new(gensym("wg_scatter_tilde~"),
                                             (t_newmethod)wg_scatter_tilde_new,
                                             (t_method)wg_scatter_tilde_free,
                                             sizeof(t_wg_scatter_tilde),
                                             CLASS_DEFAULT,
                                             A_GIMME, 0);

    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_dsp, gensym("dsp"), A_CANT, 0);

    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_enable, gensym("enable"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_coupling_all, gensym("coupling_all"), A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_coupling, gensym("coupling"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_impedance_ratio, gensym("impedance_ratio"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(wg_scatter_tilde_tilde_class, (t_method)wg_scatter_tilde_couple_loss, gensym("couple_loss"), A_FLOAT, 0);

    /* leftmost inlet is a signal inlet; we still keep c_all as the CLASS_MAINSIGNALIN float */
    CLASS_MAINSIGNALIN(wg_scatter_tilde_tilde_class, t_wg_scatter_tilde, c_all);

    /* optional alias */
    class_addcreator((t_newmethod)wg_scatter_tilde_new, gensym("wg_scatter~"), A_GIMME, 0);
}
