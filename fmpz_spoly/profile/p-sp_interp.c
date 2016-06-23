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
#define MINWALL (1000)

#include "fmpz_spoly/ptimer.h"
PTIMER_EXTERN(BPITIME)

int main(int argc, char** argv)
{
    fmpz_spoly_struct orig[NUMEX];
    fmpz_spoly_sp_interp_basis_struct bases[NUMEX];
    fmpz_spoly_sp_interp_eval_struct evals[NUMEX];
    fmpz_spoly_t res;
    fmpz_t D;
    ulong dbits, hbits;
    timeit_t timer;
    slong T, i, l, loops;
    double ctime, wtime;
    int retval = 0;

    FLINT_TEST_INIT(state);

    PTIMER_ENABLE(BPITIME);

    if (argc != 4)
    {
        flint_printf("usage: %s terms log2_degree log2_height\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    T = (slong) strtoul(argv[1], NULL, 10);
    dbits = strtoul(argv[2], NULL, 10);
    hbits = strtoul(argv[3], NULL, 10);

    fmpz_init(D);
    fmpz_setbit(D, dbits);
    fmpz_sub_ui(D, D, UWORD(1));

    flint_printf("Testing time for sp_interp with %wd terms, degree ", T);
    fmpz_print(D);
    flint_printf(" (%wu bits), and %wd-bit coefficients.\n", dbits, hbits);

    flint_printf("Generating examples"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        fmpz_spoly_init(orig + i);
        fmpz_spoly_randtest(orig + i, state, T, D, hbits);
        fmpz_spoly_sp_interp_basis_init(bases + i, state, T, dbits, hbits);
        fmpz_spoly_sp_interp_eval_init(evals + i, bases + i);
        if (fmpz_spoly_terms(orig + i) != T)
        {
            flint_printf("\nERROR: only room for %wd terms\n",
                fmpz_spoly_terms(orig + i));
            abort();
        }
        putchar('.'); fflush(stdout);
    }
    fmpz_spoly_init(res);
    putchar('\n'); fflush(stdout);

    /* warm up and check time */
    flint_printf("Initial time estimate"); fflush(stdout);

    timeit_start(timer);
    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        fmpz_spoly_sp_interp_eval(evals + i, orig + i);
        fmpz_spoly_sp_interp(res, evals + i);
    }
    timeit_stop(timer);
    flint_printf(" %lf ms cpu, %lf ms wall\n",
        ((double) timer->cpu) / NUMEX, ((double) timer->wall) / NUMEX);

    loops = 2 * MINWALL / timer->wall + 1;

    flint_printf("Critical timing"); fflush(stdout);
    while (1)
    {
        putchar('.'); fflush(stdout);
        PTIMER_CLEAR(BPITIME);
        timeit_start(timer);
        for (l = 0; l < loops; ++l)
        {
            for (i = 0; i < NUMEX; ++i)
            {
                fmpz_spoly_sp_interp_eval(evals + i, orig + i);
                fmpz_spoly_sp_interp(res, evals + i);
            }
        }
        timeit_stop(timer);

        if (timer->wall >= MINWALL) break;
        else loops *= 2;
    }
    PTIMER_DISABLE(BPITIME);
    putchar('\n'); fflush(stdout);

    /* show time */
    flint_printf("\n");
    flint_printf("loops: %wd\n", loops);
    ctime = ((double)timer->cpu) / (NUMEX * loops);
    flint_printf("  cpu: %lf ms avg\n", ctime);
    wtime = ((double)timer->wall) / (NUMEX * loops);
    flint_printf(" wall: %lf ms avg\n", wtime);

    /* 
    PTIMER_PRINT(BPITIME, NUMEX * loops);
    */

    /* cool down and check results */
    flint_printf("\n");
    flint_printf("Final checks"); fflush(stdout);

    for (i = 0; i < NUMEX; ++i)
    {
        putchar('.'); fflush(stdout);
        fmpz_spoly_sp_interp_eval(evals + i, orig + i);
        fmpz_spoly_sp_interp(res, evals + i);
        if (!fmpz_spoly_equal(res, orig + i))
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
        fmpz_spoly_clear(orig+i);
        fmpz_spoly_sp_interp_basis_clear(bases + i);
        fmpz_spoly_sp_interp_eval_clear(evals + i);
    }
    fmpz_spoly_clear(res);
    fmpz_clear(D);

    return retval;
}
