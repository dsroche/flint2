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
#include "fmpz_sparse.h"
#include "profiler.h"

#define NUMEX (10)
#define MINCPU (1000)

int main(int argc, char** argv)
{
    fmpz_sparse_struct orig[NUMEX];
    fmpz_sparse_bp_interp_struct interps[NUMEX];
    fmpz_sparse_t res;
    fmpz_t D, H;
    timeit_t timer;
    ulong T, hbits, dbits;
    slong i, l, loops;
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

    for (i=0; i<NUMEX; ++i)
    {
        fmpz_sparse_init(orig+i);
        fmpz_sparse_randtest(orig+i, state, T, D, hbits-1);
        fmpz_sparse_bp_interp_init(interps+i, T, H, D, state);
    }
    fmpz_sparse_init(res);

    /* warm up and check time */
    timeit_start(timer);
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_sparse_bp_interp_eval(interps+i, orig+i);
        fmpz_sparse_bp_interp(res, interps+i);
    }
    timeit_stop(timer);

    loops = 2*MINCPU / timer->cpu + 1;

    while (1)
    {
        timeit_start(timer);
        for (l=0; l<loops; ++l)
        {
            for (i=0; i<NUMEX; ++i)
            {
                fmpz_sparse_bp_interp_eval(interps+i, orig+i);
                fmpz_sparse_bp_interp(res, interps+i);
            }
        }
        timeit_stop(timer);

        if (timer->cpu >= MINCPU) break;
        else loops *= 2;
    }

    /* cool down and check results */
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_sparse_bp_interp_eval(interps+i, orig+i);
        fmpz_sparse_bp_interp(res, interps+i);
        if (!fmpz_sparse_equal(res, orig+i))
        {
            flint_printf("FAIL\n");
            fmpz_sparse_print(orig+i); flint_printf("\n\n");
            fmpz_sparse_print(res); flint_printf("\n\n");
            fmpz_print(D); flint_printf("\n"); /* FIXME */
            { 
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
            fmpz_print(interps[i].q); flint_printf("\n");
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
        fmpz_sparse_clear(orig+i);
        fmpz_sparse_bp_interp_clear(interps+i);
    }
    fmpz_sparse_clear(res);
    fmpz_clear(H);
    fmpz_clear(D);

    return retval;
}
