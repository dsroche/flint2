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

#define NUMEX (10)
#define MINCPU (1000)

int main(int argc, char** argv)
{
    fmpz_spoly_struct orig[NUMEX];
    fmpz_spoly_sp_interp_basis_struct bases[NUMEX];
    fmpz_spoly_sp_interp_eval_struct evals[NUMEX];
    fmpz_spoly_t res;
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
        fmpz_spoly_init(orig+i);
        fmpz_spoly_randtest(orig+i, state, T, D, hbits-1);
        fmpz_spoly_sp_interp_basis_init(bases + i, state, T, dbits, hbits);
        fmpz_spoly_sp_interp_eval_init(evals + i, bases + i);
    }
    fmpz_spoly_init(res);

    /* warm up and check time */
    timeit_start(timer);
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_sp_interp_eval(evals+i, orig+i);
        fmpz_spoly_sp_interp(res, evals+i);
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
                fmpz_spoly_sp_interp_eval(evals+i, orig+i);
                fmpz_spoly_sp_interp(res, evals+i);
            }
        }
        timeit_stop(timer);

        if (timer->cpu >= MINCPU) break;
        else loops *= 2;
    }

    /* cool down and check results */
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_sp_interp_eval(evals+i, orig+i);
        fmpz_spoly_sp_interp(res, evals+i);
        if (!fmpz_spoly_equal(res, orig+i))
        {
            flint_printf("FAIL\n");
            fmpz_spoly_print(orig+i); flint_printf("\n\n");
            fmpz_spoly_print(res); flint_printf("\n\n");
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
        fmpz_spoly_sp_interp_eval_clear(evals + i);
        fmpz_spoly_sp_interp_basis_clear(bases + i);
    }
    fmpz_spoly_clear(res);
    fmpz_clear(H);
    fmpz_clear(D);

    return retval;
}
