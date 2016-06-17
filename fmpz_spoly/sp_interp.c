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
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_spoly.h"

typedef struct {
    ulong dcoeff;
    ulong res;
    nmod_t mod;
} image;

static int image_cmp(const void* a, const void* b)
{
    ulong ca = ((image*) a)->dcoeff;
    ulong cb = ((image*) b)->dcoeff;
    if (ca < cb) return -1;
    else if (ca > cb) return 1;
    else return 0;
}


/* NOTE: tight dependency with fmpz/comb_init.c; mostly copied from there */
static void _fmpz_comb_reusable_init(fmpz_comb_t comb, slong num_primes)
{
    slong i, j;
    slong n;

    comb->num_primes = num_primes;

    n = FLINT_BIT_COUNT(num_primes);
    comb->n = n;

    /* Create nmod_poly modulus information */
    comb->mod = (nmod_t *) flint_malloc(sizeof(nmod_t) * num_primes);

    /* Nothing to do */
    if (n == 0) return;

    /* Allocate space for comb and res */
    comb->comb = (fmpz **) flint_malloc(n * sizeof(fmpz *));
    comb->res = (fmpz **) flint_malloc(n * sizeof(fmpz *));

    /* Size of top level */
    j = (WORD(1) << (n - 1));

    /* Initialise arrays at each level */
    for (i = 0; i < n; i++)
    {
        comb->comb[i] = _fmpz_vec_init(j);
        comb->res[i] = _fmpz_vec_init(j);
        j /= 2;
    }
}

/* NOTE: tight dependency with fmpz/comb_init.c; mostly copied from there */
static void _fmpz_comb_reusable_set(fmpz_comb_t comb, mp_srcptr primes, const image* images)
{
    slong i, j;
    slong num, log_res;
    ulong log_comb;
    fmpz_t temp, temp2;

    comb->primes = primes;

    for (i = 0; i < comb->num_primes; i++) comb->mod[i] = images[i].mod;

    /* Compute products of pairs of primes and place in comb */
    for (i = 0, j = 0; i + 2 <= comb->num_primes; i += 2, j++)
    {
        fmpz_set_ui(comb->comb[0] + j, primes[i]);
        fmpz_mul_ui(comb->comb[0] + j, comb->comb[0] + j, primes[i+1]);
    }

    /* In case number of primes is odd */
    if (i < comb->num_primes)
    {
        fmpz_set_ui(comb->comb[0] + j, primes[i]);
        i += 2;
        j++;
    }

    /* Set the rest of the entries on that row of the comb to 1 */
    num = (WORD(1) << comb->n);
    for (; i < num; i += 2, j++)
    {
        fmpz_one(comb->comb[0] + j);
    }

    /* Compute rest of comb by multiplying in pairs */
    log_comb = 1;
    num /= 2;
    while (num >= 2)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            fmpz_mul(comb->comb[log_comb] + j, comb->comb[log_comb-1] + i,
                comb->comb[log_comb-1] + i + 1);
        }
        log_comb++;
        num /= 2;
    }

    /* Compute inverses from pairs of primes */
    fmpz_init(temp);
    fmpz_init(temp2);

    for (i = 0, j = 0; i + 2 <= comb->num_primes; i += 2, j++)
    {
        fmpz_set_ui(temp, primes[i]);
        fmpz_set_ui(temp2, primes[i+1]);
        fmpz_invmod(comb->res[0] + j, temp, temp2);
    }

    fmpz_clear(temp);
    fmpz_clear(temp2);

    /* Compute remaining inverses, each level
       combining pairs from the level below */
    log_res = 1;
    num = (WORD(1) << (comb->n - 1));

    while (log_res < comb->n)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            fmpz_invmod(comb->res[log_res] + j, comb->comb[log_res-1] + i,
                comb->comb[log_res-1] + i + 1);
        }
        log_res++;
        num /= 2;
    }
}

int fmpz_spoly_sp_interp(fmpz_spoly_t res, 
        const fmpz_spoly_sp_interp_eval_t eval)
{
    fmpz_spoly_sp_interp_eval_t roundeval, remeval;
    fmpz_spoly_t round_poly;
    image *eimgs, *cimgs;
    fmpz_comb_t ecomb, ccomb;
    fmpz_comb_temp_t ectemp, cctemp;
    fmpz_t full_coeff, full_expon;
    mp_ptr crt_res, crt_primes;
    slong i, j, k, l, gstart = 0;
    slong eimg_len = 0, cimg_len = 0;
    ulong coeff;
    const fmpz_spoly_sp_interp_basis_struct* basis = eval->basis;
    int result = (basis->length == 0);

    fmpz_spoly_zero(res);
    fmpz_spoly_sp_interp_eval_init(roundeval, basis);
    fmpz_spoly_sp_interp_eval_init(remeval, basis);
    fmpz_spoly_sp_interp_eval_set(remeval, eval);

    fmpz_spoly_init(round_poly);

    eimgs = flint_malloc(basis->eimg_per_round * sizeof *eimgs);
    cimgs = flint_malloc(basis->cimg_per_round * sizeof *cimgs);

    _fmpz_comb_reusable_init(ecomb, basis->eimg_needed);
    _fmpz_comb_reusable_init(ccomb, basis->cimg_needed);
    fmpz_comb_temp_init(ectemp, ecomb);
    fmpz_comb_temp_init(cctemp, ccomb);

    l = FLINT_MAX(basis->eimg_needed, basis->cimg_needed);
    crt_res = flint_malloc(l * sizeof *crt_res);
    crt_primes = flint_malloc(l * sizeof *crt_primes);

    fmpz_init(full_coeff);
    fmpz_init(full_expon);

    for (i = 0; i < basis->length; ++i)
    {
        if (basis->shifts[i] > 1)
        {
            /* start of a group */
            gstart = i;

            /* collect exponent images */
            for (j = 0; j < nmod_poly_length(remeval->evals + gstart); ++j)
            {
                ulong coeff = nmod_poly_get_coeff_ui(remeval->evals + gstart, j);
                if (coeff != 0 && eimg_len < basis->eimg_per_round)
                {
                    eimgs[eimg_len].dcoeff = coeff;
                    eimgs[eimg_len].res = j;
                    eimgs[eimg_len].mod = basis->emods[i];
                    ++eimg_len;
                }
            }
        }
        else
        {
            /* collect coefficient images */
            for (j = 0; j < nmod_poly_length(remeval->evals + gstart); ++j)
            {
                coeff = nmod_poly_get_coeff_ui(remeval->evals + gstart, j);
                if (coeff != 0 && cimg_len < basis->cimg_per_round)
                {
                    cimgs[cimg_len].dcoeff = coeff;
                    cimgs[cimg_len].res = 
                        nmod_poly_get_coeff_ui(remeval->evals + i, j);
                    cimgs[cimg_len].mod = basis->cmods[i];
                    ++cimg_len;
                }
            }
        }

        if (i + 1 == basis->length 
            || (basis->shifts[i + 1] > 1 
                && (basis->cmods[gstart].n != basis->cmods[i + 1].n
                    || basis->shifts[gstart] != basis->shifts[i + 1])))
        {
            /* end of round */
            qsort(eimgs, eimg_len, sizeof *eimgs, image_cmp);
            qsort(cimgs, cimg_len, sizeof *cimgs, image_cmp);

            j = k = 0;
            while (j + basis->eimg_needed <= eimg_len 
                   && k + basis->cimg_needed <= cimg_len)
            {
                if (eimgs[j].dcoeff < cimgs[k].dcoeff)
                {
                    /* skip past eimg */
                    for (++j; j < eimg_len && eimgs[j].dcoeff < cimgs[k].dcoeff; ++j);
                }
                else if (cimgs[k].dcoeff < eimgs[j].dcoeff)
                {
                    /* skip past cimg */
                    for (++k; k < cimg_len && cimgs[k].dcoeff < eimgs[j].dcoeff; ++k);
                }
                else
                {
                    /* check sufficiently many images for coeff and expon */
                    if (eimgs[j].dcoeff == eimgs[j + basis->eimg_needed - 1].dcoeff
                        && cimgs[k].dcoeff == cimgs[k + basis->cimg_needed - 1].dcoeff)
                    {
                        /* recover a new term */
                        for (l = 0; l < basis->eimg_needed; ++l) 
                        {
                            crt_res[l] = eimgs[j + l].res;
                            crt_primes[l] = eimgs[j + l].mod.n;
                        }
                        _fmpz_comb_reusable_set(ecomb, crt_primes, eimgs + j);
                        fmpz_multi_CRT_ui(full_expon, crt_res, ecomb, ectemp, 0);

                        for (l = 0; l < basis->cimg_needed; ++l) 
                        {
                            crt_res[l] = cimgs[k + l].res;
                            crt_primes[l] = cimgs[k + l].mod.n;
                        }
                        _fmpz_comb_reusable_set(ccomb, crt_primes, cimgs + k);
                        fmpz_multi_CRT_ui(full_coeff, crt_res, ccomb, cctemp, 1);

                        /* add new term to round_poly */
                        fmpz_spoly_set_coeff(round_poly, full_coeff, full_expon);
                    }

                    /* skip to next images */
                    for (++j; j < eimg_len && eimgs[j].dcoeff == eimgs[j-1].dcoeff; ++j);
                    for (++k; k < cimg_len && cimgs[k].dcoeff == cimgs[k-1].dcoeff; ++k);
                }
            }

            /* update evaluations */
            fmpz_spoly_sp_interp_eval(roundeval, round_poly);
            fmpz_spoly_sp_interp_addmul_si(remeval, WORD(-1), roundeval);

            /* update res */
            if (fmpz_spoly_is_zero(res)) fmpz_spoly_swap(res, round_poly);
            else
            {
                fmpz_spoly_add(res, res, round_poly);
                fmpz_spoly_zero(round_poly);
            }

            /* check if done */
            if (fmpz_spoly_sp_interp_eval_is_zero(remeval))
            {
                result = 1;
                break;
            }
            
            eimg_len = cimg_len = 0;
        }
    }

    flint_free(eimgs);
    flint_free(cimgs);
    fmpz_spoly_clear(round_poly);
    fmpz_spoly_sp_interp_eval_clear(roundeval);
    fmpz_spoly_sp_interp_eval_clear(remeval);
    fmpz_comb_temp_clear(ectemp);
    fmpz_comb_temp_clear(cctemp);
    fmpz_comb_clear(ecomb);
    fmpz_comb_clear(ccomb);
    flint_free(crt_res);
    flint_free(crt_primes);
    fmpz_clear(full_coeff);
    fmpz_clear(full_expon);

    return result;
}
