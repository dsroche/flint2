/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "flint.h"
#include "fmpz_spoly.h"
#include "profiler.h"

#define MINWALL (1000)
#define DO_BOT (0)
#define DO_BPI (0)
#define DO_KAL (0)
#define DO_GS (1)
#define DO_SPI (1)

#include "fmpz_spoly/ptimer.h"
PTIMER_DECLARE(BOTTIME, 10)
PTIMER_DECLARE(KALTIME, 10)

PTIMER_EXTERN(BPITIME)

void fmpz_disclog(fmpz_t res, const fmpz_t a, const fmpz_t w, const fmpz_t winv, ulong k, const fmpz_t p)
{
    fmpz_t wpow, apow, arem, wipow;
    ulong shift;
    slong i;

    fmpz_init_set(wpow, w);
    fmpz_init_set(wipow, winv);
    fmpz_init_set(arem, a);
    fmpz_init(apow);

    fmpz_zero(res);
    for (shift = 0; shift < k; ++shift)
    {
        fmpz_set(apow, arem);
        for (i = 1; (ulong)i < k - shift; ++i)
        {
            fmpz_mul(apow, apow, apow);
            fmpz_mod(apow, apow, p);
        }
        if (fmpz_is_one(apow))
        {
        }
        else
        {
            fmpz_setbit(res, shift);
            fmpz_mul(arem, arem, wipow);
            fmpz_mod(arem, arem, p);
        }
        fmpz_mul(wpow, wpow, wpow);
        fmpz_mod(wpow, wpow, p);
        fmpz_mul(wipow, wipow, wipow);
        fmpz_mod(wipow, wipow, p);
    }

    fmpz_clear(wpow);
    fmpz_clear(wipow);
    fmpz_clear(arem);
    fmpz_clear(apow);
}

void fmpz_poly_simple_roots(fmpz* roots, flint_rand_t state, const fmpz_poly_t poly, const fmpz_t maxroot)
{
    slong i;
    nmod_poly_t fp;
    ulong p;
    nmod_poly_factor_t fact;
    fmpz_t pk;
    slong n = fmpz_poly_length(poly) - 1;
    ulong* cofacs;
    fmpz_t temp;

    if (n <= 0)
        return;

    FLINT_ASSERT(fmpz_is_one(fmpz_poly_lead(poly)));

    while(1) 
    {
        p = n_randprime(state, FLINT_BITS - 2, 0);
        nmod_poly_init(fp, p);
        fmpz_poly_get_nmod_poly(fp, poly);
        if (nmod_poly_is_squarefree(fp))
            break;
        else
            nmod_poly_clear(fp);
    }

    nmod_poly_factor_init(fact);
    nmod_poly_factor_fit_length(fact, n);

    nmod_poly_factor_equal_deg(fact, fp, 1);
    FLINT_ASSERT(fact->num == n);
    nmod_poly_derivative(fp, fp);

    cofacs = flint_malloc(n * sizeof *cofacs);

    for (i = 0; i < n; ++i)
    {
        FLINT_ASSERT(fact->exp[i] == 1);
        FLINT_ASSERT(nmod_poly_length(fact->p + i) == 2);
        cofacs[i] = nmod_neg(nmod_poly_get_coeff_ui(fact->p + i, 0), fp->mod);
        fmpz_set_ui(roots + i, cofacs[i]);
        cofacs[i] = nmod_poly_evaluate_nmod(fp, cofacs[i]);
        FLINT_ASSERT(cofacs[i] != 0);
        cofacs[i] = p - nmod_inv(cofacs[i], fp->mod);
    }

    nmod_poly_factor_clear(fact);

    fmpz_init(temp);
    fmpz_init_set_ui(pk, p);
    while (fmpz_cmp(pk, maxroot) <= 0)
    {
        for (i = 0; i < n; ++i)
        {
            fmpz_poly_evaluate_fmpz(temp, poly, roots + i);
            fmpz_divexact(temp, temp, pk);
            fmpz_addmul_ui(roots + i, pk,
                nmod_mul(fmpz_fdiv_ui(temp, p), cofacs[i], fp->mod));
        }
        fmpz_mul_ui(pk, pk, p);
    }
    fmpz_clear(temp);

    fmpz_clear(pk);
    flint_free(cofacs);
    nmod_poly_clear(fp);
}

typedef struct {
    fmpz* coeffs;
    fmpz** expons;
    slong nvars;
    slong len;
} mvp_t[1];

void mvp_init(mvp_t f, slong nvars, slong len) 
{
    slong i;
    f->nvars = nvars;
    f->len = len;
    f->coeffs = _fmpz_vec_init(len);
    f->expons = flint_malloc(len * sizeof *f->expons);
    f->expons[0] = _fmpz_vec_init(len * nvars);
    for (i = 1; i < len; ++i)
    {
        f->expons[i] = f->expons[i - 1] + nvars;
    }
}

void mvp_clear(mvp_t f)
{
    _fmpz_vec_clear(f->coeffs, f->len);
    _fmpz_vec_clear(f->expons[0], f->len * f->nvars);
    flint_free(f->expons);
}

void mvp_set_fmpz_spoly(mvp_t f, const fmpz_spoly_t g, ulong shift)
{
    slong i, j;
    fmpz_t e;

    FLINT_ASSERT(f->len == fmpz_spoly_terms(g));
    fmpz_init(e);

    for (i = 0; i < f->len; ++i)
    {
        fmpz_spoly_get_term_coeff(f->coeffs + i, g, i);
        fmpz_spoly_get_term_expon(e, g, i);
        for (j = 0; j < f->nvars; ++j)
        {
            fmpz_fdiv_r_2exp(f->expons[i] + j, e, shift);
            fmpz_fdiv_q_2exp(e, e, shift);
        }
        FLINT_ASSERT(fmpz_is_zero(e));
    }

    fmpz_clear(e);
}

void mvp_get_fmpz_spoly(fmpz_spoly_t g, const mvp_t f, ulong shift)
{
    slong i, j;
    fmpz_t e;

    fmpz_spoly_zero(g);
    _fmpz_spoly_reserve(g, f->len);
    fmpz_init(e);

    for (i = 0; i < f->len; ++i)
    {
        fmpz_set(g->coeffs + i, f->coeffs + i);
        fmpz_set_ui(g->expons + i, UWORD(0));
        for (j = 0; j < f->nvars; ++j)
        {
            fmpz_mul_2exp(e, f->expons[i] + j, shift * (ulong)j);
            fmpz_add(g->expons + i, g->expons + i, e);
        }
    }

    fmpz_clear(e);

    _fmpz_spoly_set_length(g, f->len);
    _fmpz_spoly_normalise(g);
}

void mvp_eval1(fmpz_t res, const mvp_t f, const fmpz* xx)
{
    slong i, j;
    fmpz_t t1, t2;
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_zero(res);
    for (i = 0; i < f->len; ++i)
    {
        fmpz_set(t2, f->coeffs + i);
        for(j = 0; j < f->nvars; ++j)
        {
            fmpz_pow_ui(t1, xx + j, fmpz_get_ui(f->expons[i] + j));
            fmpz_mul(t2, t2, t1);
        }
        fmpz_add(res, res, t2);
    }
    fmpz_clear(t1);
    fmpz_clear(t2);
}

void mvp_eval2(fmpz_t res, const mvp_t f, const fmpz* xx, const fmpz_t m)
{
    slong i, j;
    fmpz_t t1, t2;
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_zero(res);
    for (i = 0; i < f->len; ++i)
    {
        fmpz_set(t2, f->coeffs + i);
        for(j = 0; j < f->nvars; ++j)
        {
            fmpz_powm(t1, xx + j, f->expons[i] + j, m);
            fmpz_mul(t2, t2, t1);
            fmpz_mod(t2, t2, m);
        }
        fmpz_add(res, res, t2);
    }
    fmpz_mod(res, res, m);
    fmpz_clear(t1);
    fmpz_clear(t2);
}

void mvp_print(const mvp_t f)
{
    slong i, j;
    flint_printf("mvp");
    for (i = 0; i < f->len; ++i)
    {
        flint_printf(" ");
        fmpz_print(f->coeffs + i);
        flint_printf(" [");
        for (j = 0; j < f->nvars; ++j)
        {
            if (j) flint_printf(" ");
            fmpz_print(f->expons[i] + j);
        }
        flint_printf("]");
    }
    flint_printf("\n");
}

/* ben-or & tiwari */

void bot(fmpz_spoly_t out, flint_rand_t state, const mvp_t in, const fmpz_t deg, ulong cbits)
{
    slong t, n, i, j, k;
    ulong* primes;
    fmpz* xx;
    fmpz* evals;
    fmpz_mat_t mat, vec, bb;
    fmpz_t tempint;
    fmpz_poly_t prony;
    fmpz* roots;
    mvp_t res;
    fmpz_factor_t ifac;
    int check;

    fmpz_init(tempint);

    t = in->len;
    n = in->nvars;
    
    xx = _fmpz_vec_init(n);
    for (i = 0; i < n; ++i)
    {
        fmpz_one(xx + i);
    }
    
    primes = flint_malloc(n * sizeof *primes);

    primes[0] = UWORD(2);
    for (i = 1; i < n; ++i)
    {
        primes[i] = n_nextprime(primes[i - 1], 1);
    }

    /* evals */

    evals = _fmpz_vec_init(2 * t);
    for (i = 0; i < 2*t; ++i)
    {
        mvp_eval1(evals + i, in, xx);
        for (j = 0; j < n; ++j)
        {
            fmpz_mul_ui(xx + j, xx + j, primes[j]);
        }
    }

    PTIMER_BEGIN(BOTTIME, "Solve Hankel system");

    /* berlekamp massey */

    fmpz_mat_init(mat, t, t);
    fmpz_mat_init(vec, t, 1);
    fmpz_mat_init(bb, t, 1);

    for (i = 0; i < t; ++i)
    {
        for (j = 0; j < t; ++j)
        {
            fmpz_set(fmpz_mat_entry(mat, i, j), evals + (i + j));
        }
        fmpz_set(fmpz_mat_entry(bb, i, 0), evals + (i + t));
    }

    check = fmpz_mat_solve(vec, tempint, mat, bb);
    FLINT_ASSERT(check == 1);
    fmpz_mul_si(tempint, tempint, WORD(-1));
    
    fmpz_poly_init(prony);
    fmpz_poly_set_coeff_si(prony, t, WORD(1));
    for (i = 0; i < t; ++i)
    {
        fmpz_divexact(fmpz_poly_get_coeff_ptr(prony, i), fmpz_mat_entry(vec, i, 0), tempint);
    }

    PTIMER_NEXT(BOTTIME, "Root finding");

    /* root finding */

    roots = _fmpz_vec_init(t);
    fmpz_set_ui(tempint, 1);
    for (j = 0; j < n; ++j)
        fmpz_mul_ui(tempint, tempint, primes[j]);
    fmpz_pow_ui(tempint, tempint, fmpz_get_ui(deg));
    fmpz_poly_simple_roots(roots, state, prony, tempint);

    fmpz_poly_clear(prony);

    PTIMER_NEXT(BOTTIME, "Factoring the roots");

    /* discrete logs */

    mvp_init(res, n, t);
    fmpz_factor_init(ifac);

    for (i = 0; i < t; ++i)
    {
        fmpz_factor_trial_range(ifac, roots + i, 0, n);
        FLINT_ASSERT(ifac->sign == 1);
        FLINT_ASSERT(fmpz_cmp_ui(ifac->p + (ifac->num - 1), primes[n - 1]) <= 0);
        for (j = 0, k = 0; k < ifac->num; ++k)
        {
            for (; fmpz_cmp_ui(ifac->p + k, primes[j]) > 0; ++j);
                fmpz_zero(res->expons[i] + j);
            FLINT_ASSERT(j < n);
            FLINT_ASSERT(fmpz_equal_ui(ifac->p + k, primes[j]));
            fmpz_set_ui(res->expons[i] + j, ifac->exp[k]);
        }
        for (++j; j < n; ++j)
            fmpz_zero(res->expons[i] + j);
    }

    fmpz_factor_clear(ifac);

    PTIMER_NEXT(BOTTIME, "Solving Vandermonde system");

    /* get coeffs */

    for (j = 0; j < t; ++j)
        fmpz_set_si(fmpz_mat_entry(mat, 0, j), WORD(1));

    for (i = 1; i < t; ++i)
    {
        for (j = 0; j < t; ++j)
        {
            fmpz_mul(fmpz_mat_entry(mat, i, j),
                     fmpz_mat_entry(mat, i - 1, j),
                     roots + j);
        }
    }

    for (i = 0; i < t; ++i)
        fmpz_set(fmpz_mat_entry(bb, i, 0), evals + i);

    check = fmpz_mat_solve(vec, tempint, mat, bb);
    FLINT_ASSERT(check == 1);

    for (i = 0; i < t; ++i)
    {
        fmpz_divexact(res->coeffs + i, fmpz_mat_entry(vec, i, 0), tempint);
    }

    PTIMER_END(BOTTIME);

    /* cleanup */

    _fmpz_vec_clear(roots, t);

    fmpz_mat_clear(mat);
    fmpz_mat_clear(vec);
    fmpz_mat_clear(bb);

    mvp_get_fmpz_spoly(out, res, fmpz_bits(deg));
    mvp_clear(res);

    _fmpz_vec_clear(evals, 2*t);
    _fmpz_vec_clear(xx, n);
    flint_free(primes);
    fmpz_clear(tempint);
}

/* Kaltofen */

void kal(fmpz_spoly_t out, flint_rand_t state, const fmpz_spoly_t in, ulong cbits)
{
    fmpz_mod_poly_t G;
    fmpz_mod_poly_factor_t gfac;
    fmpz * roots;
    slong i, t;
    fmpz_t temp, winv;

    fmpz_spoly_bp_interp_basis_t basis;
    fmpz_spoly_bp_interp_eval_t eval;

    fmpz_spoly_bp_interp_basis_init(basis, state, fmpz_spoly_terms(in), fmpz_bits(fmpz_spoly_degree_ptr(in)), cbits);
    fmpz_spoly_bp_interp_eval_init(eval, basis);
    
    fmpz_spoly_bp_interp_eval(eval, in);

    fmpz_spoly_zero(out);

    if (eval->basis->length == 0) 
    {
        return;
    }

    /* Berlekamp-Massey to discover Prony polynomial */

    fmpz_mod_poly_init(G, eval->basis->q);
    PTIMER_BEGIN(KALTIME, "minpoly");
    fmpz_mod_poly_minpoly(G, eval->evals, eval->basis->length);

    t = fmpz_mod_poly_degree(G);
    if (t > eval->basis->length / 2)
    {
        /* sparsity estimate was too low */
        fmpz_mod_poly_clear(G);
        return;
    }

    _fmpz_spoly_reserve(out, t);
    for (i=out->length; i<t; ++i)
    {
        fmpz_init(out->coeffs + i);
        fmpz_init(out->expons + i);
    }
    _fmpz_spoly_set_length(out, t);

    /* find roots of Prony polynomial, and their orders */

    PTIMER_NEXT(KALTIME, "root finding");
    fmpz_mod_poly_factor_init(gfac);
    fmpz_mod_poly_factor_fit_length(gfac, t);
    fmpz_mod_poly_factor_equal_deg(gfac, G, 1);

    FLINT_ASSERT(gfac->num == t);
    roots = _fmpz_vec_init(t);
    for (i = 0; i < t; ++i)
    {
        FLINT_ASSERT(gfac->exp[i] == 1);
        FLINT_ASSERT(fmpz_mod_poly_length(gfac->poly + i) == 2);
        if (fmpz_is_zero(gfac->poly[i].coeffs + 0))
            fmpz_zero(roots + i);
        else
            fmpz_sub(roots + i, &G->p, gfac->poly[i].coeffs + 0);
    }

    fmpz_mod_poly_factor_clear(gfac);

    PTIMER_NEXT(KALTIME, "discrete logs");

    fmpz_init(temp);
    fmpz_init(winv);
    fmpz_setbit(temp, eval->basis->log2_order);
    fmpz_sub_ui(temp, temp, UWORD(1));
    fmpz_powm(winv, eval->basis->points + 1, temp, eval->basis->q);

    for (i = 0; i < t; ++i)
    {
        fmpz_disclog(out->expons + i, 
            roots + i, eval->basis->points + 1, winv, eval->basis->log2_order, eval->basis->q);
    }

    fmpz_clear(temp);
    fmpz_clear(winv);

    /* solve transposed Vandermode to get coeffs */
    /* Varndermonde(roots)^T * x = evals, truncated to length t */

    PTIMER_NEXT(KALTIME, "transposed vandermonde");
    _fmpz_spoly_transp_vandermonde_inv(out->coeffs,
            roots, eval->evals, t, eval->basis->q);

    /* sort terms and remove zero coeffs */
    PTIMER_END(KALTIME);
    _fmpz_spoly_normalise(out);

    /* clean-up */
    _fmpz_vec_clear(roots, t);
    fmpz_mod_poly_clear(G);

    fmpz_spoly_bp_interp_eval_clear(eval);
    fmpz_spoly_bp_interp_basis_clear(basis);
}

/* bp_interp */

void bpi(fmpz_spoly_t out, flint_rand_t state, const fmpz_spoly_t in, ulong cbits)
{
    fmpz_spoly_bp_interp_basis_t basis;
    fmpz_spoly_bp_interp_eval_t eval;

    fmpz_spoly_bp_interp_basis_init(basis, state, fmpz_spoly_terms(in), fmpz_bits(fmpz_spoly_degree_ptr(in)), cbits);
    fmpz_spoly_bp_interp_eval_init(eval, basis);
    
    fmpz_spoly_bp_interp_eval(eval, in);
    fmpz_spoly_bp_interp(out, eval);

    fmpz_spoly_bp_interp_eval_clear(eval);
    fmpz_spoly_bp_interp_basis_clear(basis);
}

/* garg & schost */

void gs(fmpz_spoly_t out, flint_rand_t state, const fmpz_spoly_t in, ulong cbits)
{
    fmpz_t Q;
    slong T = fmpz_spoly_terms(in);
    fmpz_t D;
    fmpz_t temp;
    slong i, j;
    ulong pbits;
    fmpz_poly_t eval, X;
    nmod_poly_t Xi;
    fmpz_t M;
    ulong p;
    slong maxt = 0;
    mp_ptr roots;

    fmpz_init(D);
    fmpz_spoly_degree(D, in);
    fmpz_init_set_ui(Q, UWORD(2));
    fmpz_init(temp);

    for (i = -T + 2; i <= 1; ++i)
    {
        fmpz_sub_ui(temp, D, T);
        fmpz_add_ui(temp, temp, (i + T));
        fmpz_mul(Q, Q, temp);
    }

    pbits = FLINT_BIT_COUNT((ulong)((10.0 * log(2) / 6) * (double)T * (double)(T-1) * fmpz_bits(D)));

    fmpz_init_set_ui(M, UWORD(1));
    fmpz_poly_init(X);
    fmpz_poly_init(eval);
    roots = flint_malloc(T * sizeof *roots);

    while (fmpz_cmp(M, Q) <= 0)
    {
        do
        {
            p = n_randprime(state, pbits, 0);
        } while (fmpz_divisible_si(M, (slong)p));
        fmpz_spoly_rem_cyc_dense(eval, in, p);
        for (i = 0, j = 0; i < fmpz_poly_length(eval); ++i)
        {
            if (! fmpz_is_zero(fmpz_poly_get_coeff_ptr(eval, i)))
            {
                FLINT_ASSERT(j < T);
                roots[j] = (ulong) i;
                ++j;
            }
        }
        if (j < maxt)
        {
            continue;
        }
        else if (j > maxt)
        {
            maxt = j;
            fmpz_poly_zero(X);
            fmpz_set_ui(M, UWORD(1));
        }

        nmod_poly_init(Xi, p);
        nmod_poly_product_roots_nmod_vec(Xi, roots, j);
        fmpz_poly_CRT_ui(X, X, M, Xi, 1);
        fmpz_mul_ui(M, M, p);
        nmod_poly_clear(Xi);
    }

    flint_free(roots);

    _fmpz_spoly_reserve(out, maxt);
    fmpz_poly_simple_roots(out->expons, state, X, D);
    
    for (i = 0; i < maxt; ++i)
    {
        fmpz_poly_get_coeff_fmpz(out->coeffs + i,
            eval, fmpz_fdiv_ui(out->expons + i, p));
        FLINT_ASSERT(! fmpz_is_zero(out->coeffs + i));
    }

    _fmpz_spoly_set_length(out, maxt);
    _fmpz_spoly_normalise(out);

    fmpz_poly_clear(X);
    fmpz_poly_clear(eval);
    fmpz_clear(M);
    fmpz_clear(temp);
    fmpz_clear(Q);
    fmpz_clear(D);
}

/* arnold, giesbrecht, roche */

/* sp_interp */

void spi(fmpz_spoly_t out, flint_rand_t state, const fmpz_spoly_t in, ulong cbits)
{
    fmpz_spoly_sp_interp_basis_t basis;
    fmpz_spoly_sp_interp_eval_t eval;

    fmpz_spoly_sp_interp_basis_init(basis, state, fmpz_spoly_terms(in), fmpz_bits(fmpz_spoly_degree_ptr(in)), cbits);
    fmpz_spoly_sp_interp_eval_init(eval, basis);
    
    fmpz_spoly_sp_interp_eval(eval, in);
    fmpz_spoly_sp_interp(out, eval);

    fmpz_spoly_sp_interp_eval_clear(eval);
    fmpz_spoly_sp_interp_basis_clear(basis);
}

int main(int argc, char** argv)
{
    fmpz_spoly_t orig;
    mvp_t origm;
    fmpz_spoly_t res;
    fmpz_t D;
    ulong dbits, hbits, nvars;
    timeit_t timer;
    slong T, l, loops;
    double ctime, wtime;
    int retval = 0;

    FLINT_TEST_INIT(state);

    if (argc != 5)
    {
        flint_printf("usage: %s terms degree nvars cbits\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    T = (slong) strtoul(argv[1], NULL, 10);
    fmpz_init(D);
    fmpz_set_str(D, argv[2], 10);
    dbits = fmpz_bits(D);
    nvars = strtoul(argv[3], NULL, 10);
    hbits = strtoul(argv[4], NULL, 10);

    flint_printf("Testing sparse interpolation with %wu variables, degree ", nvars);
    fmpz_print(D);
    flint_printf(", %wd terms, and %wu bits in each coefficient.\n", T, hbits);

    fmpz_spoly_init(orig);
    fmpz_spoly_randtest_kron(orig, state, T, D, hbits, dbits, nvars);
    mvp_init(origm, nvars, T);
    mvp_set_fmpz_spoly(origm, orig, dbits);
    fmpz_spoly_init(res);

    if (fmpz_spoly_terms(orig) != T)
    {
        flint_printf("\nERROR: only room for %wd terms\n",
            fmpz_spoly_terms(orig));
        abort();
    }

    if (DO_BOT) {
        flint_printf("\n========== BEN-OR & TIWARI =======================================\n");
        
        loops = 1;
        PTIMER_ENABLE(BOTTIME);

        flint_printf("\nrunning...");
        fflush(stdout);

        while (1)
        {
            PTIMER_CLEAR(BOTTIME);
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                bot(res, state, origm, D, hbits);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) 
                break;
            else
                loops *= 2;
        }
        PTIMER_DISABLE(BOTTIME);

        flint_printf("complete...");
        fflush(stdout);

        if (fmpz_spoly_equal(res, orig))
        {
            flint_printf("check passed.\n");
        }
        else
        {
            flint_printf("check FAILED.\n");
            retval = 1;
        }

        flint_printf("\n");
        flint_printf("loops: %wd\n", loops);
        ctime = ((double)timer->cpu) / loops;
        flint_printf("  cpu: %lf ms avg\n", ctime);
        wtime = ((double)timer->wall) / loops;
        flint_printf(" wall: %lf ms avg\n", wtime);

        PTIMER_PRINT(BOTTIME, loops);
    }

    if (DO_KAL) {
        flint_printf("\n========== BIG PRIME WITH KALTOFEN'S IMPROVEMENTS ================\n");
        
        loops = 1;
        PTIMER_ENABLE(KALTIME);

        flint_printf("\nrunning...");
        fflush(stdout);

        while (1)
        {
            PTIMER_CLEAR(KALTIME);
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                kal(res, state, orig, hbits);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) 
                break;
            else
                loops *= 2;
        }
        PTIMER_DISABLE(KALTIME);

        flint_printf("complete...");
        fflush(stdout);

        if (fmpz_spoly_equal(res, orig))
        {
            flint_printf("check passed.\n");
        }
        else
        {
            flint_printf("check FAILED.\n");
            retval = 1;
        }

        flint_printf("\n");
        flint_printf("loops: %wd\n", loops);
        ctime = ((double)timer->cpu) / loops;
        flint_printf("  cpu: %lf ms avg\n", ctime);
        wtime = ((double)timer->wall) / loops;
        flint_printf(" wall: %lf ms avg\n", wtime);

        PTIMER_PRINT(KALTIME, loops);
    }
        
    if (DO_BPI) {
        flint_printf("\n========== BIG PRIME WITH COMBINED ROOT FINDING AND LOGS =========\n");
        
        loops = 1;
        PTIMER_ENABLE(BPITIME);

        flint_printf("\nrunning...");
        fflush(stdout);

        while (1)
        {
            PTIMER_CLEAR(BPITIME);
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                bpi(res, state, orig, hbits);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) 
                break;
            else
                loops *= 2;
        }
        PTIMER_DISABLE(BPITIME);

        flint_printf("complete...");
        fflush(stdout);

        if (fmpz_spoly_equal(res, orig))
        {
            flint_printf("check passed.\n");
        }
        else
        {
            flint_printf("check FAILED.\n");
            retval = 1;
        }

        flint_printf("\n");
        flint_printf("loops: %wd\n", loops);
        ctime = ((double)timer->cpu) / loops;
        flint_printf("  cpu: %lf ms avg\n", ctime);
        wtime = ((double)timer->wall) / loops;
        flint_printf(" wall: %lf ms avg\n", wtime);

        PTIMER_PRINT(BPITIME, loops);
    }
        
    if (DO_GS) {
        flint_printf("\n========== GARG & SCHOST (SMALL PRIMES) ==========================\n");
        
        loops = 1;
        PTIMER_ENABLE(BPITIME);

        flint_printf("\nrunning...");
        fflush(stdout);

        while (1)
        {
            PTIMER_CLEAR(BPITIME);
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                gs(res, state, orig, hbits);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) 
                break;
            else
                loops *= 2;
        }
        PTIMER_DISABLE(BPITIME);

        flint_printf("complete...");
        fflush(stdout);

        if (fmpz_spoly_equal(res, orig))
        {
            flint_printf("check passed.\n");
        }
        else
        {
            flint_printf("check FAILED.\n");
            retval = 1;
        }

        flint_printf("\n");
        flint_printf("loops: %wd\n", loops);
        ctime = ((double)timer->cpu) / loops;
        flint_printf("  cpu: %lf ms avg\n", ctime);
        wtime = ((double)timer->wall) / loops;
        flint_printf(" wall: %lf ms avg\n", wtime);

        /* 
        PTIMER_PRINT(BPITIME, loops);
        */
    }
        
    if (DO_SPI) {
        flint_printf("\n========== DOUBLE SMALL PRIMES WITH CHEATING =====================\n");
        
        loops = 1;
        PTIMER_ENABLE(BPITIME);

        flint_printf("\nrunning...");
        fflush(stdout);

        while (1)
        {
            PTIMER_CLEAR(BPITIME);
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                spi(res, state, orig, hbits);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) 
                break;
            else
                loops *= 2;
        }
        PTIMER_DISABLE(BPITIME);

        flint_printf("complete...");
        fflush(stdout);

        if (fmpz_spoly_equal(res, orig))
        {
            flint_printf("check passed.\n");
        }
        else
        {
            flint_printf("check FAILED.\n");
            retval = 1;
        }

        flint_printf("\n");
        flint_printf("loops: %wd\n", loops);
        ctime = ((double)timer->cpu) / loops;
        flint_printf("  cpu: %lf ms avg\n", ctime);
        wtime = ((double)timer->wall) / loops;
        flint_printf(" wall: %lf ms avg\n", wtime);

        /* 
        PTIMER_PRINT(BPITIME, loops);
        */
    }
        
    /* clean up */
    FLINT_TEST_CLEANUP(state);
    fmpz_spoly_clear(orig);
    mvp_clear(origm);
    fmpz_spoly_clear(res);
    fmpz_clear(D);

    return retval;
}
