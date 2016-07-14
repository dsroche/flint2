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
#include "profiler.h"

typedef struct 
{
    slong vars;
    ulong kshift;
    fmpz_t maxdeg;
    fmpz_spoly_t poly;
} fmpz_spoly_mvar_struct;

typedef fmpz_spoly_mvar_struct fmpz_spoly_mvar_t[1];

void fmpz_spoly_mvar_init2(fmpz_spoly_mvar_t poly, slong vars, const fmpz_t maxdeg)
{
    poly->vars = vars;
    poly->kshift = fmpz_bits(maxdeg);
    fmpz_init(poly->maxdeg);
    fmpz_spoly_init(poly->poly);
}

void fmpz_spoly_mvar_clear(fmpz_spoly_mvar_t poly)
{
    fmpz_clear(poly->maxdeg);
    fmpz_spoly_clear(poly->poly);
}

FMPZ_SPOLY_INLINE
void _fmpz_spoly_mvar_rekron(const fmpz_spoly_mvar_t poly, ulong kshift)
{
    if (kshift > poly->kshift)
    {
        flint_printf("ERROR: Kronecker re-shift not implemented!\n");
        abort();
    }
}

void fmpz_spoly_mvar_get_coeff(fmpz_t c, const fmpz_spoly_mvar_t poly, const fmpz* expons)
{
    fmpz_t e;
    slong i;

    i = _fmpz_vec_height_index(expons, poly->vars);
    if (fmpz_cmp(expons + i, poly->maxdeg) > 0)
    {
        fmpz_zero(c);
        return;
    }
    
    fmpz_init_set(e, expons + 0);
    for (i = 1; i < poly->vars; ++i)
    {
        fmpz_mul_2exp(e, e, poly->kshift);
        fmpz_add(e, e, expons + i);
    }
    fmpz_spoly_get_coeff(c, poly->poly, e);
    fmpz_clear(e);
}

void fmpz_spoly_mvar_set_coeff(fmpz_spoly_mvar_t poly, const fmpz_t c, const fmpz* expons)
{
    fmpz_t e;
    slong i;

    i = _fmpz_vec_height_index(expons, poly->vars);
    if (fmpz_cmp(expons + i, poly->maxdeg) > 0)
    {
        fmpz_set(poly->maxdeg, expons + i);
        _fmpz_spoly_mvar_rekron(poly, fmpz_bits(poly->maxdeg));
    }
    
    fmpz_init_set(e, expons + 0);
    for (i = 1; i < poly->vars; ++i)
    {
        fmpz_mul_2exp(e, e, poly->kshift);
        fmpz_add(e, e, expons + i);
    }
    fmpz_spoly_set_coeff(poly->poly, c, e);
    fmpz_clear(e);
}

FMPZ_SPOLY_INLINE
slong fmpz_spoly_mvar_terms(fmpz_spoly_mvar_t poly)
{
    return fmpz_spoly_terms(poly->poly);
}

void fmpz_spoly_mvar_mul_heaps(fmpz_spoly_mvar_t res, 
    const fmpz_spoly_mvar_t poly1, const fmpz_spoly_mvar_t poly2)
{
    ulong shift;

    fmpz_add(res->maxdeg, poly1->maxdeg, poly2->maxdeg);
    shift = FLINT_MAX(FLINT_MAX(poly1->kshift, poly2->kshift), 
                      fmpz_bits(res->maxdeg));
    if (shift != poly1->kshift)
        _fmpz_spoly_mvar_rekron(poly1, shift);
    if (shift != shift)
        _fmpz_spoly_mvar_rekron(poly2, shift);
    res->kshift = shift;

    fmpz_spoly_mul_heaps(res->poly, poly1->poly, poly2->poly);
}

void fmpz_spoly_mvar_mul_OS(fmpz_spoly_mvar_t res, flint_rand_t state, 
    const fmpz_spoly_mvar_t poly1, const fmpz_spoly_mvar_t poly2)
{
    ulong shift;

    fmpz_add(res->maxdeg, poly1->maxdeg, poly2->maxdeg);
    shift = FLINT_MAX(FLINT_MAX(poly1->kshift, poly2->kshift), 
                      fmpz_bits(res->maxdeg));
    if (shift != poly1->kshift)
        _fmpz_spoly_mvar_rekron(poly1, shift);
    if (shift != shift)
        _fmpz_spoly_mvar_rekron(poly2, shift);
    res->kshift = shift;

    fmpz_spoly_mul_OS(res->poly, state, poly1->poly, poly2->poly);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_mvar_set(fmpz_spoly_mvar_t res, const fmpz_spoly_mvar_t poly)
{
    res->vars = poly->vars;
    res->kshift = poly->kshift;
    fmpz_set(res->maxdeg, poly->maxdeg);
    fmpz_spoly_set(res->poly, poly->poly);
}

#define NUMEX 6
#define MINWALL 1000

int main(int argc, char** argv)
{
    fmpz_spoly_mvar_t tpoly;
    fmpz_spoly_mvar_struct f[NUMEX];
    fmpz_spoly_mvar_struct g[NUMEX];
    fmpz_spoly_mvar_struct res_heaps[NUMEX];
    fmpz_spoly_mvar_struct res_OS[NUMEX];
    double heapT[NUMEX], OST[NUMEX];
    const char* names[NUMEX];
    slong interms[NUMEX];
    slong outerms[NUMEX];
    timeit_t timer;
    slong i, j, l, loops, vars;
    fmpz_t degree, c;
    fmpz* expons;

    FLINT_TEST_INIT(state);

    fmpz_init(degree);
    fmpz_init(c);
    expons = _fmpz_vec_init(5);

    flint_printf("Generating inputs..."); fflush(stdout);

    i = 0;


    /***************/ names[i] = "$f_6 (f_6 + 1)$";
    vars = 2;
    fmpz_set_ui(degree, UWORD(3000));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_set_si(c, WORD(3));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 0, UWORD(11));
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 1, UWORD(13));
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 100; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_set(g + i, f + i);
    _fmpz_vec_zero(expons, vars);
    fmpz_spoly_mvar_get_coeff(c, f + i, expons);
    fmpz_add_ui(c, c, UWORD(1));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;

    /***************/ names[i] = "$f_1 (f_1 + 1)$";
    vars = 3;
    fmpz_set_ui(degree, UWORD(40));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 0);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 1);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 2);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 20; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_set(g + i, f + i);
    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(3));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;


    /***************/ names[i] = "$f_2 (f_2 + 1)$";
    vars = 3;
    fmpz_set_ui(degree, UWORD(80));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 0, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 1, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 2, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 20; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_set(g + i, f + i);
    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(3));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;


    /***************/ names[i] = "$f_3 (f_3 + 1)$";
    vars = 3;
    fmpz_set_ui(degree, UWORD(60));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 0);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 1);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 2);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 30; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_set(g + i, f + i);
    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(3));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;


    /***************/ names[i] = "$f_4 (f_4 + 1)$";
    vars = 4;
    fmpz_set_ui(degree, UWORD(40));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 0);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 1);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 2);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_one(expons + 3);
    fmpz_one(c);
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 20; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_set(g + i, f + i);
    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(3));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;


    /***************/ names[i] = "$f_5 g_5$";
    vars = 5;
    fmpz_set_ui(degree, UWORD(40));
    fmpz_spoly_mvar_init2(tpoly, vars, degree);
    fmpz_spoly_mvar_init2(f + i, vars, degree);
    fmpz_spoly_mvar_init2(g + i, vars, degree);
    fmpz_spoly_mvar_init2(res_heaps + i, vars, degree);
    fmpz_spoly_mvar_init2(res_OS + i, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 0, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 1, UWORD(1));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 2, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 3, UWORD(1));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 4, UWORD(1));
    fmpz_set_si(c, WORD(-1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(f + i, tpoly);
    for (j = 1; j < 10; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(f + i, f + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(f + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    fmpz_spoly_mvar_init2(tpoly, vars, degree);

    _fmpz_vec_zero(expons, vars);
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 0, UWORD(1));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 1, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 2, UWORD(1));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 3, UWORD(2));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    _fmpz_vec_zero(expons, vars); fmpz_set_ui(expons + 4, UWORD(1));
    fmpz_set_si(c, WORD(1));
    fmpz_spoly_mvar_set_coeff(tpoly, c, expons);

    fmpz_spoly_mvar_set(g + i, tpoly);
    for (j = 1; j < 10; ++j)
    {
        fmpz_spoly_mvar_mul_heaps(g + i, g + i, tpoly);
    }

    _fmpz_vec_zero(expons, vars);
    fmpz_set_ui(c, UWORD(2));
    fmpz_spoly_mvar_set_coeff(g + i, c, expons);

    fmpz_spoly_mvar_clear(tpoly);
    ++i;


    flint_printf("\n");

    for (i = 0; i < NUMEX; ++i)
    {
        flint_printf("Testing %s...", names[i]); fflush(stdout);

        /* heaps */
        timeit_start(timer);
        fmpz_spoly_mvar_mul_heaps(res_heaps + i, f + i, g + i);
        timeit_stop(timer);

        loops = 2 * MINWALL / timer->wall + 1;

        while (1)
        {
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                fmpz_spoly_mvar_mul_heaps(res_heaps + i, f + i, g + i);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) break;
            else loops *= 2;
        }

        fmpz_spoly_mvar_mul_heaps(res_heaps + i, f + i, g + i);
        heapT[i] = ((double)timer->wall) / loops / 1000;

        /* OS */
        timeit_start(timer);
        fmpz_spoly_mvar_mul_OS(res_OS + i, state, f + i, g + i);
        timeit_stop(timer);

        loops = 2 * MINWALL / timer->wall + 1;

        while (1)
        {
            timeit_start(timer);
            for (l = 0; l < loops; ++l)
            {
                fmpz_spoly_mvar_mul_OS(res_OS + i, state, f + i, g + i);
            }
            timeit_stop(timer);

            if (timer->wall >= MINWALL) break;
            else loops *= 2;
        }

        fmpz_spoly_mvar_mul_OS(res_OS + i, state, f + i, g + i);
        OST[i] = ((double)timer->wall) / loops / 1000;

        /* show times */
        interms[i] = fmpz_spoly_mvar_terms(f + i);
        outerms[i] = fmpz_spoly_mvar_terms(res_heaps + i);
        flint_printf(" interms: %wd, outerms: %wd, heaps: %.3lf, OS: %.3lf\n", 
            interms[i], outerms[i], heapT[i], OST[i]);
    }

    flint_printf("\n\\begin{tabular}{l|r|r}\nExample & Input terms & Output terms & Heaps & Output-sensitive \\\\ \\hline\n");
    for (i = 0; i < NUMEX; ++i)
    {
        flint_printf("  %s & %wd & %wd & %.3lf & %.3lf \\\\\n",
            names[i], interms[i], outerms[i], heapT[i], OST[i]);
    }
    flint_printf("\\end{tabular}\n");


    for (i = 0; i < NUMEX; ++i)
    {
        fmpz_spoly_mvar_clear(f + i);
        fmpz_spoly_mvar_clear(g + i);
        fmpz_spoly_mvar_clear(res_heaps + i);
        fmpz_spoly_mvar_clear(res_OS + i);
    }
    fmpz_clear(degree);
    fmpz_clear(c);
    _fmpz_vec_clear(expons, 5);
    FLINT_TEST_CLEANUP(state);

    return 0;
}


#if 0
    if (argc != 6)
    {
        flint_printf("usage: %s terms degree log2_height kron_shift nvars\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    fmpz_init(D);
    T = (slong) strtoul(argv[1], NULL, 10);
    fmpz_set_str(D, argv[2], 10);
    hbits = strtoul(argv[3], NULL, 10);
    kshift = strtoul(argv[4], NULL, 10);
    nvars = strtoul(argv[5], NULL, 10);

    flint_printf("Testing time for mul_OS with %wd terms, degree ", T);
    fmpz_print(D);
    flint_printf(",\n  %wu bits per coefficient, %wu-bit Kronecker shift, and %wu variables.\n\n",
            hbits, kshift, nvars);

    flint_printf("Generating examples"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        fmpz_spoly_init(f + i);
        fmpz_spoly_randtest_kron(f + i, state, T, D, hbits, kshift, nvars);
        fmpz_spoly_init(g + i);
        fmpz_spoly_randtest_kron(g + i, state, T, D, hbits, kshift, nvars);
        if (FLINT_MIN(fmpz_spoly_terms(f + i), fmpz_spoly_terms(g + i)) < T)
        {
            flint_printf("\nERROR: only room for %wd terms\n",
                FLINT_MIN(fmpz_spoly_terms(f + i), fmpz_spoly_terms(g + i)));
            abort();
        }
        putchar('.'); fflush(stdout);
    }
    fmpz_spoly_init(res);
    fmpz_spoly_init(check);
    putchar('\n'); fflush(stdout);

    /* warm up and check time */
    flint_printf("Initial time estimate"); fflush(stdout);

    timeit_start(timer);
    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        fmpz_spoly_mul_OS(res, state, f + i, g + i);
        if (i == 0)
        {
            minterms = maxterms = fmpz_spoly_terms(res);
            mindbits = maxdbits = fmpz_bits(fmpz_spoly_degree_ptr(res));
            minhbits = maxhbits = fmpz_spoly_height_bits(res);
        }
        else
        {
            minterms = FLINT_MIN(minterms, fmpz_spoly_terms(res));
            maxterms = FLINT_MAX(maxterms, fmpz_spoly_terms(res));
            mindbits = FLINT_MIN(mindbits, fmpz_bits(fmpz_spoly_degree_ptr(res)));
            maxdbits = FLINT_MAX(maxdbits, fmpz_bits(fmpz_spoly_degree_ptr(res)));
            minhbits = FLINT_MIN(minhbits, fmpz_spoly_height_bits(res));
            maxhbits = FLINT_MAX(maxhbits, fmpz_spoly_height_bits(res));
        }
    }
    timeit_stop(timer);
    flint_printf(" %lf ms cpu, %lf ms wall\n",
        ((double) timer->cpu) / NUMEX, ((double) timer->wall) / NUMEX);

    /* show params */
    flint_printf("\n");
    flint_printf("terms: %wd - %wd\n", minterms, maxterms);
    flint_printf("degree bits: %wu - %wu\n", mindbits, maxdbits);
    flint_printf("height bits: %wu - %wu\n", minhbits, maxhbits);
    flint_printf("\n");

    loops = 2 * MINWALL / timer->wall + 1;

    flint_printf("Critical timing"); fflush(stdout);
    while (1)
    {
        putchar('.'); fflush(stdout);
        PTIMER_CLEAR(MITIME);
        timeit_start(timer);
        for (l = 0; l < loops; ++l)
        {
            for (i = 0; i < NUMEX; ++i)
            {
                fmpz_spoly_mul_OS(res, state, f + i, g + i);
            }
        }
        timeit_stop(timer);

        if (timer->wall >= MINWALL) break;
        else loops *= 2;
    }
    PTIMER_DISABLE(MITIME);
    putchar('\n'); fflush(stdout);

    /* show time */
    flint_printf("\n");
    flint_printf("loops: %wd\n", loops);
    ctime = ((double)timer->cpu) / (NUMEX * loops);
    flint_printf("  cpu: %lf ms avg\n", ctime);
    wtime = ((double)timer->wall) / (NUMEX * loops);
    flint_printf(" wall: %lf ms avg\n", wtime);

    /* 
    PTIMER_PRINT(MITIME, NUMEX * loops);
    */

    /* cool down and check results */
    flint_printf("\n");
    flint_printf("Final checks"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        fmpz_spoly_mul_OS(res, state, f + i, g + i);
        fmpz_spoly_mul_heaps(check, f + i, g + i);
        if (!fmpz_spoly_equal(res, check))
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
        fmpz_spoly_clear(f + i);
        fmpz_spoly_clear(g + i);
    }
    fmpz_spoly_clear(res);
    fmpz_spoly_clear(check);
    fmpz_clear(D);

    return retval;
}

#endif
