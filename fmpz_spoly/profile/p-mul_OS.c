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
PTIMER_EXTERN(MITIME)

int main(int argc, char** argv)
{
    fmpz_spoly_struct f[NUMEX];
    fmpz_spoly_struct g[NUMEX];
    fmpz_spoly_t res;
    fmpz_spoly_t check;
    fmpz_t D;
    timeit_t timer;
    ulong hbits, nvars, kshift;
    slong T, i, l, loops;
    double ctime, wtime;
    int retval = 0;
    slong minterms = 0, maxterms = 0;
    mp_bitcnt_t mindbits = 0, maxdbits = 0, minhbits = 0, maxhbits = 0;

    FLINT_TEST_INIT(state);

    PTIMER_ENABLE(MITIME);

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
