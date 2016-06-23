/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "flint.h"
#include "fmpz_spoly.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "profiler.h"

/* imitation of bp_interp or sp_interp, for testing only */
typedef struct 
{
    mp_ptr primes;
    slong len;
} poly_interp_basis_struct;
typedef poly_interp_basis_struct poly_interp_basis_t[1];

void poly_interp_basis_init(poly_interp_basis_t basis, 
    flint_rand_t state, slong deg, mp_bitcnt_t hbits)
{
    slong i;
    ulong pbits = FLINT_BITS - 2;
    basis->len = hbits / (pbits - 1) + 1;
    basis->primes = flint_malloc(basis->len * sizeof *basis->primes);
    for (i = 0; i < basis->len; ++i)
        basis->primes[i] = n_randprime(state, pbits, 0);
}

void poly_interp_basis_clear(poly_interp_basis_t basis)
{
    flint_free(basis->primes);
}

typedef struct 
{
    nmod_poly_struct* evals;
    const poly_interp_basis_struct* basis;
} poly_interp_eval_struct;
typedef poly_interp_eval_struct poly_interp_eval_t[1];

void poly_interp_eval_init(poly_interp_eval_t res, const poly_interp_basis_t basis)
{
    slong i;
    res->basis = basis;
    res->evals = flint_malloc(basis->len * sizeof *res->evals);
    for (i = 0; i < basis->len; ++i)
        nmod_poly_init(res->evals + i, basis->primes[i]);
}

void poly_interp_eval_clear(poly_interp_eval_t res)
{
    slong i;
    for (i = 0; i < res->basis->len; ++i)
        nmod_poly_clear(res->evals + i);
    flint_free(res->evals);
}

void poly_interp_eval(poly_interp_eval_t res, const fmpz_poly_t poly)
{
    slong i, j;
    fmpz_comb_t comb;
    fmpz_comb_temp_t ctemp;
    mp_ptr mods;
    mods = flint_malloc(res->basis->len * sizeof *mods);
    fmpz_comb_init(comb, res->basis->primes, res->basis->len);
    fmpz_comb_temp_init(ctemp, comb);
    for (i = 0; i < res->basis->len; ++i)
    {
        nmod_poly_zero(res->evals + i);
        nmod_poly_fit_length(res->evals + i, fmpz_poly_length(poly));
    }
    for (i = 0; i < fmpz_poly_length(poly); ++i)
    {
        fmpz_multi_mod_ui(mods, fmpz_poly_get_coeff_ptr(poly, i), comb, ctemp);
        for (j = 0; j < res->basis->len; ++j)
            nmod_poly_set_coeff_ui(res->evals + j, i, mods[j]);
    }
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(ctemp);
    flint_free(mods);
}

void poly_interp(fmpz_poly_t poly, const poly_interp_eval_t eval)
{
    slong i, j, plen = 0;
    fmpz_comb_t comb;
    fmpz_comb_temp_t ctemp;
    mp_ptr mods;
    fmpz_t coeff;
    mods = flint_malloc(eval->basis->len * sizeof *mods);
    fmpz_comb_init(comb, eval->basis->primes, eval->basis->len);
    fmpz_comb_temp_init(ctemp, comb);
    fmpz_init(coeff);
    for (i = 0; i < eval->basis->len; ++i)
        plen = FLINT_MAX(plen, nmod_poly_length(eval->evals + i));
    fmpz_poly_zero(poly);
    fmpz_poly_fit_length(poly, plen);
    for (i = 0; i < plen; ++i)
    {
        for (j = 0; j < eval->basis->len; ++j)
            mods[j] = nmod_poly_get_coeff_ui(eval->evals + j, i);
        fmpz_multi_CRT_ui(coeff, mods, comb, ctemp, 1);
        fmpz_poly_set_coeff_fmpz(poly, i, coeff);
    }
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(ctemp);
    flint_free(mods);
    fmpz_clear(coeff);
}

#define NUMEX (10)
#define MINWALL (1000)

int main(int argc, char** argv)
{
    fmpz_spoly_struct orig[NUMEX];
    fmpz_poly_struct dense[NUMEX];
    poly_interp_basis_struct bases[NUMEX];
    poly_interp_eval_struct evals[NUMEX];
    fmpz_poly_t res;
    fmpz_t D;
    slong wD;
    ulong dbits, hbits;
    timeit_t timer;
    slong T, i, l, loops;
    double ctime, wtime;
    int retval = 0;

    FLINT_TEST_INIT(state);

    if (argc != 4)
    {
        flint_printf("usage: %s terms log2_degree log2_height\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    T = (slong) strtoul(argv[1], NULL, 10);
    dbits = strtoul(argv[2], NULL, 10);
    hbits = strtoul(argv[3], NULL, 10);

    if (dbits > FLINT_BITS - 2)
    {
        flint_printf("ERROR: degree is way too big for dense.\n");
        abort();
    }
    wD = (WORD(1) << dbits) - 1;
    fmpz_init_set_si(D, wD);

    flint_printf("Testing time for sp_interp with %wd terms, degree ", T);
    fmpz_print(D);
    flint_printf(" (%wu bits), and %wd-bit coefficients.\n", dbits, hbits);

    flint_printf("Generating examples"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        fmpz_spoly_init(orig + i);
        fmpz_spoly_randtest(orig + i, state, T, D, hbits);
        fmpz_poly_init(dense + i);
        fmpz_spoly_get_fmpz_poly(dense + i, orig + i);
        poly_interp_basis_init(bases + i, state, wD, hbits);
        poly_interp_eval_init(evals + i, bases + i);
        if (fmpz_spoly_terms(orig + i) != T)
        {
            flint_printf("\nERROR: only room for %wd terms\n",
                fmpz_spoly_terms(orig + i));
            abort();
        }
        putchar('.'); fflush(stdout);
    }
    fmpz_poly_init(res);
    putchar('\n'); fflush(stdout);

    /* warm up and check time */
    flint_printf("Initial time estimate"); fflush(stdout);

    timeit_start(timer);
    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        poly_interp_eval(evals + i, dense + i);
        poly_interp(res, evals + i);
    }
    timeit_stop(timer);
    flint_printf(" %lf ms cpu, %lf ms wall\n",
        ((double) timer->cpu) / NUMEX, ((double) timer->wall) / NUMEX);

    loops = 2 * MINWALL / timer->wall + 1;

    flint_printf("Critical timing"); fflush(stdout);
    while (1)
    {
        putchar('.'); fflush(stdout);
        timeit_start(timer);
        for (l = 0; l < loops; ++l)
        {
            for (i = 0; i < NUMEX; ++i)
            {
                poly_interp_eval(evals + i, dense + i);
                poly_interp(res, evals + i);
            }
        }
        timeit_stop(timer);

        if (timer->wall >= MINWALL) break;
        else loops *= 2;
    }
    putchar('\n'); fflush(stdout);

    /* show time */
    flint_printf("\n");
    flint_printf("loops: %wd\n", loops);
    ctime = ((double)timer->cpu) / (NUMEX * loops);
    flint_printf("  cpu: %lf ms avg\n", ctime);
    wtime = ((double)timer->wall) / (NUMEX * loops);
    flint_printf(" wall: %lf ms avg\n", wtime);

    /* cool down and check results */
    flint_printf("\n");
    flint_printf("Final checks"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        poly_interp_eval(evals + i, dense + i);
        poly_interp(res, evals + i);
        if (!fmpz_poly_equal(res, dense + i))
        {
            flint_printf("FAIL (iteration %wd)\n", i);
            retval = 1;
        }
    }
    putchar('\n'); fflush(stdout);

    /* clean up */
    FLINT_TEST_CLEANUP(state);
    for (i = 0; i < NUMEX; ++i)
    {
        fmpz_spoly_clear(orig + i);
        fmpz_poly_clear(dense + i);
        poly_interp_eval_clear(evals + i);
        poly_interp_basis_clear(bases + i);
    }
    fmpz_poly_clear(res);
    fmpz_clear(D);

    return retval;
}
