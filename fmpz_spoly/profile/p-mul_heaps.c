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
    fmpz_spoly_struct f[NUMEX];
    fmpz_spoly_struct g[NUMEX];
    fmpz_spoly_t res;
    fmpz_spoly_t check;
    fmpz_t D, H;
    timeit_t timer;
    ulong T, hbits, dbits, nvars;
    slong i, l, loops;
    double ctime, wtime;
    int retval = 0;
    slong minterms = 0, maxterms = 0;
    mp_bitcnt_t mindbits = 0, maxdbits = 0, minhbits = 0, maxhbits = 0;

    FLINT_TEST_INIT(state);

    if (argc != 5)
    {
        flint_printf("usage: %s terms log2_degree nvars log2_height\n", argv[0]);
        FLINT_TEST_CLEANUP(state);
        return 2;
    }

    T = strtoul(argv[1], NULL, 10);
    dbits = strtoul(argv[2], NULL, 10);
    nvars = strtoul(argv[3], NULL, 10);
    hbits = strtoul(argv[4], NULL, 10);
    fmpz_init(D);
    fmpz_randbits(D, state, dbits);
    fmpz_abs(D, D);
    fmpz_init(H);
    fmpz_randbits(H, state, hbits);

    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_init(f + i);
        fmpz_spoly_randtest_kron(f + i, state, T, D, hbits - 1, dbits + 1, nvars);
        fmpz_spoly_init(g + i);
        fmpz_spoly_randtest_kron(g + i, state, T, D, hbits - 1, dbits + 1, nvars);
        if ((ulong)FLINT_MIN(fmpz_spoly_terms(f + i), fmpz_spoly_terms(g + i)) < T)
        {
            flint_printf("WARNING: only got %wd and %wd terms\n",
                    fmpz_spoly_terms(f + i), fmpz_spoly_terms(g + i));
        }
    }
    fmpz_spoly_init(res);
    fmpz_spoly_init(check);

    /* warm up and check time */
    timeit_start(timer);
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_mul_heaps(res, f + i, g + i);
    }
    timeit_stop(timer);

    loops = 2 * MINCPU / timer->cpu + 1;

    while (1)
    {
        timeit_start(timer);
        for (l=0; l<loops; ++l)
        {
            for (i=0; i<NUMEX; ++i)
            {
                fmpz_spoly_mul_heaps(res, f + i, g + i);
            }
        }
        timeit_stop(timer);

        if (timer->cpu >= MINCPU) break;
        else loops *= 2;
    }

    /* cool down and check results */
    for (i=0; i<NUMEX; ++i)
    {
        fmpz_spoly_mul_heaps(res, f + i, g + i);
        fmpz_spoly_mul_classical(check, f + i, g + i);
        if (!fmpz_spoly_equal(res, check))
        {
            flint_printf("FAIL\n");
            fmpz_spoly_print(f + i); flint_printf("\n\n");
            fmpz_spoly_print(g + i); flint_printf("\n\n");
            fmpz_spoly_print(res); flint_printf("\n\n");
            fmpz_spoly_print(check); flint_printf("\n\n");
            retval = 1;
        }
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

    /* show params */
    flint_printf("terms: %wd - %wd\n", minterms, maxterms);
    flint_printf("degree bits: %wu - %wu\n", mindbits, maxdbits);
    flint_printf("height bits: %wu - %wu\n", minhbits, maxhbits);
    flint_printf("\n");

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
        fmpz_spoly_clear(f + i);
        fmpz_spoly_clear(g + i);
    }
    fmpz_spoly_clear(res);
    fmpz_spoly_clear(check);
    fmpz_clear(H);
    fmpz_clear(D);

    return retval;
}
