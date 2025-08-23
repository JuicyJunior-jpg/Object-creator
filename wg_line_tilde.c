/* wg_line_tilde.c
   Implements Pd external: wg_line_tilde~  (file outputs: wg_line_tilde~.pd_darwin)
   Setup symbol: wg_line_tilde_tilde_setup
   Optional alias: wg_line~
*/
#include "m_pd.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ----- helpers ----- */
static inline t_float clampf(t_float x, t_float lo, t_float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
static inline t_float zapgremlins(t_float x) {
    /* denormal & NaN protection */
    if (!isfinite(x) || fabsf(x) < 1e-30f) return 0.f;
    return x;
}

/* 1st-order allpass in canonical form: H(z) = (a + z^-1) / (1 + a z^-1) */
static inline t_float ap1(t_float x, t_float *z1, t_float a) {
    t_float y = -a * x + *z1;
    *z1 = x + a * y;
    return y;
}

/* 4-point (3rd-order) Lagrange interpolation
   Uses x[-1], x[0], x[1], x[2] around integer index with fractional mu in [0,1).
*/
static inline t_float lagrange4(t_float xm1, t_float x0, t_float x1, t_float x2, t_float mu) {
    t_float c_m1 = -(mu) * (mu - 1.f) * (mu - 2.f) / 6.f;
    t_float c_0  =  (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f;
    t_float c_1  = -(mu + 1.f) * (mu) * (mu - 2.f) / 2.f;
    t_float c_2  =  (mu + 1.f) * (mu) * (mu - 1.f) / 6.f;
    return (c_m1 * xm1) + (c_0 * x0) + (c_1 * x1) + (c_2 * x2);
}

/* ----- object ----- */
typedef struct _wg_line_tilde {
    t_object  x_obj;

    /* signal inlets provided via dsp vectors; polarity is a message */
    t_inlet  *in_freq;
    t_inlet  *in_couple;
    t_inlet  *in_feedback;
    t_inlet  *in_dampAmt;
    t_inlet  *in_dampFollow;
    t_inlet  *in_disp;
    t_inlet  *in_nl;
    t_inlet  *in_gain;
    t_inlet  *in_coupleFollow;

    t_outlet *out_main;
    t_outlet *out_couple;
    t_outlet *out_report; /* float outlet */

    /* runtime */
    t_float  sr;
    t_float  inv_sr;

    /* delay buffer */
    t_sample *buf;
    int       bufsize;    /* samples */
    int       widx;

    /* states */
    t_float ap_z1[4];     /* dispersion cascade */
    t_float damp_lp;      /* damping lowpass state */
    t_float couple_lp;    /* low-shelf tilt state for couple_send */

    /* delay read smoothing */
    t_float d_read_smooth;
    t_float d_smooth_coef;

    /* status/reporting */
    t_float last_f0_reported;
    t_float polarity_sign;  /* -1 (negative), +1 (positive); default -1 */

} t_wg_line_tilde;

/* forward decls */
static t_class *wg_line_tilde_tilde_class;

static void wg_line_tilde_report(t_wg_line_tilde *x, t_floatarg f);

/* resize/prepare buffer (>=120 ms) */
static void wg_line_tilde_allocbuf(t_wg_line_tilde *x, t_float sr) {
    int need = (int)ceilf(sr * 0.12f) + 4096; /* headroom to avoid wrap hazards */
    if (need < 2048) need = 2048;
    if (need != x->bufsize) {
        if (x->buf) freebytes(x->buf, x->bufsize * sizeof(t_sample));
        x->buf = (t_sample *)getbytes(need * sizeof(t_sample));
        x->bufsize = x->buf ? need : 0;
        x->widx = 0;
        if (x->buf) memset(x->buf, 0, x->bufsize * sizeof(t_sample));
    }
}

/* wrap helper */
static inline int wrapi(int i, int n) {
    while (i < 0) i += n;
    while (i >= n) i -= n;
    return i;
}

static t_int *wg_line_tilde_perform(t_int *w) {
    t_wg_line_tilde *x = (t_wg_line_tilde *)(w[1]);
    int n = (int)(w[2]);

    /* in/out vectors (left-to-right signals) */
    t_sample *in_exciter      = (t_sample *)(w[3]);
    t_sample *in_freqHz       = (t_sample *)(w[4]);
    t_sample *in_couple       = (t_sample *)(w[5]);
    t_sample *in_feedback     = (t_sample *)(w[6]);
    t_sample *in_dampAmt      = (t_sample *)(w[7]);
    t_sample *in_dampFollow   = (t_sample *)(w[8]);
    t_sample *in_disp         = (t_sample *)(w[9]);
    t_sample *in_nl           = (t_sample *)(w[10]);
    t_sample *in_gain         = (t_sample *)(w[11]);
    t_sample *in_coupleFollow = (t_sample *)(w[12]);

    t_sample *out_main        = (t_sample *)(w[13]);
    t_sample *out_couple      = (t_sample *)(w[14]);

    /* locals */
    t_float sr = x->sr;
    t_float inv_sr = x->inv_sr;
    int     N = x->bufsize;
    if (!x->buf || N < 8) {
        /* safety: zero outputs */
        memset(out_main, 0, n * sizeof(t_sample));
        memset(out_couple, 0, n * sizeof(t_sample));
        return (w + 15);
    }

    int widx = x->widx;
    t_float z1_0 = x->ap_z1[0], z1_1 = x->ap_z1[1], z1_2 = x->ap_z1[2], z1_3 = x->ap_z1[3];
    t_float damp_lp = x->damp_lp;
    t_float couple_lp = x->couple_lp;
    t_float d_smooth = x->d_read_smooth;
    const t_float d_coef = x->d_smooth_coef;
    const t_float pol = x->polarity_sign;

    for (int i = 0; i < n; ++i) {
        /* controls per-sample */
        t_float f0 = in_freqHz[i];
        if (!isfinite(f0) || f0 <= 0.f) f0 = 1.f;
        f0 = clampf(f0, 0.1f, sr * 0.45f); /* avoid silly extremes */

        /* dispersion mapping (quadratic into stable a<1) */
        t_float dknob = clampf(in_disp[i], 0.f, 1.f);
        t_float a = 0.8f * dknob * dknob;

        /* group delay of single 1st-order allpass @ omega0 */
        t_float omega0 = 2.f * (t_float)M_PI * f0 * inv_sr;
        t_float cosw = cosf(omega0);
        t_float one_minus_a2 = 1.f - a * a;
        t_float denom = 1.f + a * a - 2.f * a * cosw;
        if (denom < 1e-12f) denom = 1e-12f;
        t_float tau_ap = one_minus_a2 / denom; /* samples */
        t_float tau_disp = 4.f * tau_ap;

        /* target delay and compensated read */
        t_float D_target = sr / f0;
        t_float D_read_target = D_target - tau_disp;
        if (D_read_target < 1.f) D_read_target = 1.f;

        /* smooth D_read to avoid zipper */
        d_smooth += d_coef * (D_read_target - d_smooth);
        if (d_smooth < 1.f) d_smooth = 1.f;

        /* fractional read index relative to current write index */
        t_float rpos = (t_float)widx - d_smooth;
        while (rpos < 0.f) rpos += (t_float)N;
        while (rpos >= (t_float)N) rpos -= (t_float)N;

        int i0 = (int)floorf(rpos);
        t_float mu = rpos - (t_float)i0;

        /* gather 4 taps with wrapping */
        int i_m1 = wrapi(i0 - 1, N);
        int i_p1 = wrapi(i0 + 1, N);
        int i_p2 = wrapi(i0 + 2, N);

        t_float xm1 = x->buf[i_m1];
        t_float x0  = x->buf[i0];
        t_float x1  = x->buf[i_p1];
        t_float x2  = x->buf[i_p2];

        /* interpolated output before processing */
        t_float y = lagrange4(xm1, x0, x1, x2, mu);

        /* dispersion cascade: 4x allpass */
        y = ap1(y, &z1_0, a);
        y = ap1(y, &z1_1, a);
        y = ap1(y, &z1_2, a);
        y = ap1(y, &z1_3, a);

        /* damping inside the loop via LP/HP split with freq-following cutoff */
        t_float dampFollow = clampf(in_dampFollow[i], 0.f, 1.f);
        t_float base_fc = 1000.f;
        t_float fc = base_fc * (1.f - dampFollow) + f0 * dampFollow;
        fc = clampf(fc, 20.f, sr * 0.45f);
        t_float a1 = expf(-2.f * (t_float)M_PI * fc * inv_sr); /* one-pole */
        damp_lp = (1.f - a1) * y + a1 * damp_lp;
        t_float y_hp = y - damp_lp;

        t_float dampAmt = clampf(in_dampAmt[i], 0.f, 1.f);
        t_float s = 2.f * dampAmt - 1.f; /* -1..+1 */
        t_float strength = fabsf(s);
        t_float y_lp = damp_lp;

        if (s <= 0.f) {
            /* LP loss (highs die faster): attenuate HP component */
            y_hp *= (1.f - strength);
        } else {
            /* HP loss (lows die faster): attenuate LP component */
            y_lp *= (1.f - strength);
        }
        y = y_lp + y_hp;

        /* cubic nonlinearity with blend */
        t_float nl = clampf(in_nl[i], 0.f, 1.f);
        t_float alpha = 0.3f * nl * nl;
        t_float y_nl = y + alpha * y * y * y;
        y = (1.f - nl) * y + nl * y_nl;

        /* feedback branch */
        t_float fb = clampf(in_feedback[i], 0.f, 0.99f);
        t_float fb_sample = x->polarity_sign * fb * y;

        /* write to buffer: exciter + couple_in + feedbackBranch */
        t_float in_sum = in_exciter[i] + in_couple[i] + fb_sample;
        in_sum = zapgremlins(in_sum);
        x->buf[widx] = in_sum;

        /* advance write index */
        if (++widx >= N) widx = 0;

        /* outputs with gain and couple low-shelf tilt */
        t_float gain = clampf(in_gain[i], 0.f, 1.5f);
        t_float main_out = gain * y;

        /* low-shelf tilt for couple_send: emphasize lows by mixing LP(main) */
        t_float cf = clampf(in_coupleFollow[i], 0.f, 1.f);
        /* shelf cutoff roughly follows pitch for musical behavior */
        t_float shelf_fc = clampf(0.5f * f0 + 150.f, 30.f, sr * 0.25f);
        t_float a1s = expf(-2.f * (t_float)M_PI * shelf_fc * inv_sr);
        couple_lp = (1.f - a1s) * main_out + a1s * couple_lp;
        t_float couple_send = (1.f - cf) * main_out + cf * couple_lp;

        out_main[i] = zapgremlins(main_out);
        out_couple[i] = zapgremlins(couple_send);

        /* auto-report on ~1% pitch change */
        t_float last = x->last_f0_reported > 0.f ? x->last_f0_reported : f0;
        t_float rel = fabsf(f0 - last) / (last > 1e-6f ? last : 1.f);
        if (rel > 0.01f) {
            outlet_float(x->out_report, d_smooth); /* effective delay (samples) */
            outlet_float(x->out_report, f0);       /* current f0 (Hz) */
            x->last_f0_reported = f0;
        }
    }

    x->ap_z1[0] = z1_0; x->ap_z1[1] = z1_1; x->ap_z1[2] = z1_2; x->ap_z1[3] = z1_3;
    x->damp_lp = damp_lp;
    x->couple_lp = couple_lp;
    x->d_read_smooth = d_smooth;
    x->widx = widx;

    return (w + 15);
}

static void wg_line_tilde_dsp(t_wg_line_tilde *x, t_signal **sp) {
    t_float newsr = (t_float)sp[0]->s_sr;
    if (newsr <= 0) newsr = 48000.f;
    if (newsr != x->sr || !x->buf) {
        x->sr = newsr;
        x->inv_sr = 1.f / x->sr;
        wg_line_tilde_allocbuf(x, x->sr);
        memset(x->ap_z1, 0, sizeof(x->ap_z1));
        x->damp_lp = 0.f;
        x->couple_lp = 0.f;
        x->d_read_smooth = x->sr / 440.f; /* start around A4 */
        x->d_smooth_coef = 0.0015f; /* gentle, zipper-safe */
        x->last_f0_reported = 0.f;
    }

    /* Inlets/Outlets order in perform:
       1: x, 2: n,
       3: exciter~, 4: freqHz~, 5: couple_in~,
       6: feedback, 7: dampAmt, 8: dampFollow,
       9: dispersion, 10: nonlinearity,
       11: gain, 12: coupleFollow,
       13: main_out~, 14: couple_send~
    */
    dsp_add(wg_line_tilde_perform, 14,
            x, sp[0]->s_n,
            sp[0]->s_vec,  /* exciter~ */
            sp[1]->s_vec,  /* freqHz~ */
            sp[2]->s_vec,  /* couple_in~ */
            sp[3]->s_vec,  /* feedback */
            sp[4]->s_vec,  /* dampAmt */
            sp[5]->s_vec,  /* dampFollow */
            sp[6]->s_vec,  /* dispersion */
            sp[7]->s_vec,  /* nonlinearity */
            sp[8]->s_vec,  /* gain */
            sp[9]->s_vec,  /* coupleFollow */
            sp[10]->s_vec, /* main_out~ */
            sp[11]->s_vec  /* couple_send~ */
    );
}

/* messages */

static void wg_line_tilde_polarity(t_wg_line_tilde *x, t_floatarg f) {
    /* 1 = negative feedback, 2 = positive; default negative */
    int mode = (int)f;
    if (mode == 2) x->polarity_sign = +1.f;
    else x->polarity_sign = -1.f;
}

static void wg_line_tilde_report(t_wg_line_tilde *x, t_floatarg f) {
    (void)f;
    /* Emit current effective delay (samples) and estimated f0 from it */
    t_float D = (x->d_read_smooth > 1.f) ? x->d_read_smooth : 1.f;
    t_float f0 = (D > 0.f) ? (x->sr / D) : 0.f;
    outlet_float(x->out_report, D);
    outlet_float(x->out_report, f0);
    x->last_f0_reported = f0;
}

static void *wg_line_tilde_new(void) {
    t_wg_line_tilde *x = (t_wg_line_tilde *)pd_new(wg_line_tilde_tilde_class);

    x->buf = NULL;
    x->bufsize = 0;
    x->widx = 0;
    memset(x->ap_z1, 0, sizeof(x->ap_z1));
    x->damp_lp = 0.f;
    x->couple_lp = 0.f;
    x->sr = sys_getsr();
    if (x->sr <= 0) x->sr = 48000.f;
    x->inv_sr = 1.f / x->sr;
    x->d_read_smooth = x->sr / 440.f;
    x->d_smooth_coef = 0.0015f;
    x->last_f0_reported = 0.f;
    x->polarity_sign = -1.f;

    wg_line_tilde_allocbuf(x, x->sr);

    /* create signal inlets: total signals = 10
       Left inlet is exciter~ (implicit via CLASS_MAINSIGNALIN).
       We add: freqHz~, couple_in~, feedback, dampAmt, dampFollow, dispersion, nonlinearity, gain, coupleFollow
    */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* freqHz~       (2) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* couple_in~    (3) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* feedback      (4) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* dampAmt       (5) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* dampFollow    (6) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* dispersion    (7) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* nonlinearity  (8) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* gain          (9) */
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); /* coupleFollow  (10) */

    /* outlets */
    x->out_main   = outlet_new(&x->x_obj, &s_signal);
    x->out_couple = outlet_new(&x->x_obj, &s_signal);
    x->out_report = outlet_new(&x->x_obj, &s_float);

    return (x);
}

static void wg_line_tilde_free(t_wg_line_tilde *x) {
    if (x->buf) {
        freebytes(x->buf, x->bufsize * sizeof(t_sample));
        x->buf = NULL;
        x->bufsize = 0;
    }
}

void wg_line_tilde_tilde_setup(void) {
    wg_line_tilde_tilde_class = class_new(gensym("wg_line_tilde~"),
                                          (t_newmethod)wg_line_tilde_new,
                                          (t_method)wg_line_tilde_free,
                                          sizeof(t_wg_line_tilde),
                                          CLASS_DEFAULT,
                                          0);

    class_addmethod(wg_line_tilde_tilde_class, (t_method)wg_line_tilde_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(wg_line_tilde_tilde_class, (t_method)wg_line_tilde_polarity, gensym("polarity"), A_FLOAT, 0);
    class_addmethod(wg_line_tilde_tilde_class, (t_method)wg_line_tilde_report, gensym("report"), A_DEFFLOAT, 0);

    CLASS_MAINSIGNALIN(wg_line_tilde_tilde_class, t_wg_line_tilde, sr); /* dummy float field; left inlet is signal */

    /* optional alias */
    class_addcreator((t_newmethod)wg_line_tilde_new, gensym("wg_line~"), 0);
}
