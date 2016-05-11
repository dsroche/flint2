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

    Authored 2016 by Daniel S. Roche and A. Whitman Groves;
    US Government work in the public domain. 

******************************************************************************/

#include "flint.h"
#include "fmpz_spoly.h"
#include "profiler.h"

#define NUMEX (20)
#define MINCPU (1000)

int main(int argc, char** argv)
{
    fmpz_spoly_struct orig[NUMEX];
    fmpz * res, * check;
    fmpz_t D, H;
    timeit_t timer;
    ulong T, hbits, dbits;
    slong i, l, loops, len1, len2;
    double ctime, wtime;
    int retval = 0;

    FLINT_TEST_INIT(state);

    if (argc != 4)
    {
        flint_printf("usage: %s terms log2_degree log2_height\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    T = strtoul(argv[1], NULL, 10);
    dbits = strtoul(argv[2], NULL, 10);
    hbits = strtoul(argv[3], NULL, 10);
    fmpz_init(D);
    fmpz_randbits(D, state, dbits);
    fmpz_abs(D, D);
    fmpz_init(H);
    fmpz_randbits(H, state, hbits);

    for (i=0; i<NUMEX; i++)
    {
        fmpz_spoly_init(orig+i);
        fmpz_spoly_randtest(orig+i, state, T, D, hbits-1);
    }

    res = NULL;
    check = NULL;

    /* get sumset time */
    timeit_start(timer);
    for (i=0; i<NUMEX-1; i += 2)
    {
        fmpz_spoly_sumset(&res, state, orig+i, orig+i+1);
    }
    timeit_stop(timer);

    loops = 2*MINCPU / timer->cpu + 1;

    while (1)
    {
        timeit_start(timer);
        for (l=0; l<loops; ++l)
        {
            for (i=0; i<NUMEX - 1; i += 2)
            {
                fmpz_spoly_sumset(&res, state, orig+i, orig+i+1);
            }
        }
        timeit_stop(timer);

        if (timer->cpu >= MINCPU) break;
        else loops *= 2;
    }

    /* cool down and check results */
    for (i=0; i<NUMEX-1; i+=2)
    {
        len1 = fmpz_spoly_sumset(&res, state, orig+i, orig+i+1);
        len2 = fmpz_spoly_sumcheck(&check, orig+i, orig+i+1);
        if (!_fmpz_vec_equal(res, check, len1) || len1 != len2)
        {
            flint_printf("FAIL\n");
            /*fmpz_spoly_print(orig+i); flint_printf("\n\n");
            fmpz_spoly_print(res); flint_printf("\n\n");
            fmpz_print(D); flint_printf("\n"); *//* FIXME */
            /*{ 
                fmpz_t order, t; 
                fmpz_init(order);  fmpz_init(t);
                fmpz_setbit(order, interps[i].log2_order); 
                fmpz_print(order); flint_printf("\n");
                fmpz_powm(t, interps[i].sample_points + 1, order, interps[i].q);
                if (fmpz_is_one(t)) flint_printf("good\n");
                fmpz_divexact_ui(order, order, UWORD(2));
                fmpz_powm(t, interps[i].sample_points + 1, order, interps[i].q);
                if (!fmpz_is_one(t)) flint_printf("great\n");
                fmpz_clear(order); fmpz_clear(t);
            }
            fmpz_print(interps[i].q); flint_printf("\n");*/
            retval = 1;
        }
    }

    /* show time */
    flint_printf("loops: %wd\n", loops);
    ctime = ((double)timer->cpu) / (NUMEX * loops);
    flint_printf("  cpu: %lf ms avg\n", ctime);
    wtime = ((double)timer->wall) / (NUMEX * loops);
    flint_printf(" wall: %lf ms avg\n", wtime);

    /* clean up */
    FLINT_TEST_CLEANUP(state);
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_clear(orig+i);
    }
    _fmpz_vec_clear(res, len1);
    _fmpz_vec_clear(check, len2);
    fmpz_clear(H);
    fmpz_clear(D);

    return retval;
}
