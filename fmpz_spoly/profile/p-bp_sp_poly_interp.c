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

    Authored 2016 by A. Whitman Groves; US Government work in the public domain

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
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

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Step size for the length
   bitslo   Minimum bit size
   bitshi   Maximum bit size
   bitsh    Step size for the bit size
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define bits     100
#define lenlo    100
#define lenhi    10000
#define lenh     495
#define deglo    10000
#define degrat   2.0
#define rows     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define cols     20 /* ((slong) ceil((log(deghi) - log(deglo)) / log(degrat))) */
#define cpumin   1000
#define ncases   1
#define nalgs    3
#define img      1
#define imgname  "bp_sp_poly_interp.ppm"

/*
   Write a binary 24-bit ppm image.
 */
int write_rgb_ppm(const char* file_name, unsigned char* pixels, 
                   unsigned int width, unsigned int height)
{
    FILE* file = fopen(file_name, "wb");
    if (file == NULL)
        return -1;
    flint_fprintf(file, "P6\n%d %d\n255\n", width, height);
    fwrite(pixels, sizeof(unsigned char), width * height * 3, file);
    fclose(file);
    return 0;
}

int
main(void)
{
    slong i, j, len;
    ulong deg;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_spoly_t f, h;
    fmpz_poly_t fdense, hdense;
    fmpz_t degree;

    FLINT_TEST_INIT(state);
       
    fmpz_spoly_init(f);
    fmpz_spoly_init(h);

    fmpz_init(degree);
    
    fmpz_poly_init(fdense);
    fmpz_poly_init(hdense);

    for (len = lenhi, i = 0; len >= lenlo; len -= lenh, i++)
    {
        slong times[nalgs];
        int c;

        j = 0;

        for (c = 0; c < nalgs; c++)
            times[c] = WORD(0);
        
        for (deg = deglo, j = 0; j < cols; deg = (slong)(deg * degrat), j++)
        {
            int n;
            int reps[nalgs];
            int doit[nalgs];
            fmpz_spoly_bp_interp_basis_t bp_basis;
            fmpz_spoly_bp_interp_eval_t bp_eval;
            fmpz_spoly_sp_interp_basis_t sp_basis;
            fmpz_spoly_sp_interp_eval_t sp_eval;
            poly_interp_basis_t poly_basis;
            poly_interp_eval_t poly_eval;

            fmpz_set_si(degree, deg);

            /* whether to actually do each algorithm or skip it */
            doit[0] = len <= 5000;
            doit[1] = 1;
            doit[2] = fmpz_bits(degree) <= 25;
            
            if (doit[0])
            {
                fmpz_spoly_bp_interp_basis_init(bp_basis, state, len, fmpz_bits(degree), bits);
                fmpz_spoly_bp_interp_eval_init(bp_eval, bp_basis);
            }
            if (doit[1])
            {
                fmpz_spoly_sp_interp_basis_init(sp_basis, state, len, fmpz_bits(degree), bits);
                fmpz_spoly_sp_interp_eval_init(sp_eval, sp_basis);
            }
            if (doit[2])
            {
                poly_interp_basis_init(poly_basis, state, deg, bits);
                poly_interp_eval_init(poly_eval, poly_basis);
            }

            for (c = 0; c < nalgs; c++)
            {
                times[c] = WORD(0);
                reps[c] = 0;
            }
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t timer;
                int l, loops;

                /*
                   Construct random sparse polynomial f
                 */
                {
                    fmpz_spoly_randtest(f, state, len, degree, bits);
                    if (doit[2])
                        fmpz_spoly_get_fmpz_poly(fdense, f);
                    if (fmpz_spoly_terms(f) != len)
                    {
                        flint_printf("ERROR: wrong number of terms\n");
                        abort();
                    }
                }

                /* bp_interp */
                c = 0;
                if (doit[c])
                {
                    loops = 1;
                    do
                    {
                        timeit_start(timer);
                        for (l = 0; l < loops; ++l)
                        {
                            fmpz_spoly_bp_interp_eval(bp_eval, f);
                            fmpz_spoly_bp_interp(h, bp_eval);
                        }
                        timeit_stop(timer);
                        times[c] += timer->cpu;
                        reps[c] += loops;
                        loops *= 2;
                    }
                    while (timer->cpu < cpumin);
                }
                else
                {
                    reps[c] = 1;
                    times[c] = WORD_MAX;
                }

                /* sp_interp */
                c = 1;
                if (doit[c])
                {
                    loops = 1;
                    do
                    {
                        timeit_start(timer);
                        for (l = 0; l < loops; ++l)
                        {
                            fmpz_spoly_sp_interp_eval(sp_eval, f);
                            fmpz_spoly_sp_interp(h, sp_eval);
                        }
                        timeit_stop(timer);
                        times[c] += timer->cpu;
                        reps[c] += loops;
                        loops *= 2;
                    }
                    while (timer->cpu < cpumin);
                }
                else
                {
                    reps[c] = 1;
                    times[c] = WORD_MAX;
                }

                /* poly_interp */
                c = 2;
                if (doit[c])
                {
                    loops = 1;
                    do
                    {
                        timeit_start(timer);
                        for (l = 0; l < loops; ++l)
                        {
                            poly_interp_eval(poly_eval, fdense);
                            poly_interp(hdense, poly_eval);
                        }
                        timeit_stop(timer);
                        times[c] += timer->cpu;
                        reps[c] += loops;
                        loops *= 2;
                    }
                    while (timer->cpu < cpumin);
                }
                else
                {
                    reps[c] = 1;
                    times[c] = WORD_MAX;
                }

            }
            
            for (c = 0; c < nalgs; c++)
                T[i][j][c] = times[c] / (double) reps[c];
            
            if (T[i][j][1] <= FLINT_MIN(T[i][j][0], T[i][j][2]))
                X[i][j] = 1;
            else if (T[i][j][0] <= T[i][j][2])
                X[i][j] = 0;
            else
                X[i][j] = 2;
           flint_printf("len = %wd, deg = %wd, winner = %d\n", len, deg, X[i][j]), fflush(stdout);
           flint_printf("bptime: %lf, sptime: %lf, polytime: %lf\n",
                   T[i][j][0], T[i][j][1], T[i][j][2]);

           if (doit[0])
           {
               fmpz_spoly_bp_interp_eval_clear(bp_eval);
               fmpz_spoly_bp_interp_basis_clear(bp_basis);
           }
           if (doit[1])
           {
               fmpz_spoly_sp_interp_eval_clear(sp_eval);
               fmpz_spoly_sp_interp_basis_clear(sp_basis);
           }
           if (doit[2])
           {
               poly_interp_eval_clear(poly_eval);
               poly_interp_basis_clear(poly_basis);
           }
        }
    }
    
    /* 
       Print 2-D ASCII image of the winning algorithms
     */
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            flint_printf("%d", X[i][j]);
        flint_printf("\n");
    }

    fmpz_spoly_clear(f);
    fmpz_spoly_clear(h);
    fmpz_poly_clear(fdense);
    fmpz_poly_clear(hdense);
    fmpz_clear(degree);

    /*
       Print 2-D coloured image to file
     */
    if (img)
    {
        unsigned char * PIXELS;
        int k;
        
        PIXELS = (unsigned char *) flint_malloc(3 * rows * cols * sizeof(unsigned char));
        k = 0;
        for (i = 0; i < rows; i++)
        {
            for (j = 0; j < cols; j++)
            {
                double v[nalgs];
                int m;
                
                for (m = 0; m < FLINT_MIN(3, nalgs); m++)
                {
                    v[m] = T[i][j][X[i][j]] / T[i][j][m];
                    PIXELS[k++] = (unsigned char) (v[m] * 255);
                }
                for(; m < 3; m++)
                    PIXELS[k++] = (unsigned char) 0;
            }
        }

        k = write_rgb_ppm(imgname, PIXELS, cols, rows);
        flint_free(PIXELS);
        if (k)
        {
            flint_printf("Exception:  writing ppm image failed\n");
        }
    }

    flint_randclear(state);

    return 0;
}
