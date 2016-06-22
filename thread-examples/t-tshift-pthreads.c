/*
    Authored 2016 by Daniel S. Roche; U.S. Government work product in the public domain

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"
#include "profiler.h"

#define BITS (1000)
#define LEN (400)

void usage(char** argv)
{
    flint_printf("usage: %s [num_threads]\n", argv[0]);
    abort();
}

int main(int argc, char** argv)
{
    timeit_t timer;
    fmpz_poly_t f, g, check;
    fmpz_t c;
    int nthreads = 4;
    int result = 0;

    FLINT_TEST_INIT(state);

    if (argc > 2)
        usage(argv);
    else if (argc == 2)
    {
        nthreads = atoi(argv[1]);
        if (nthreads <= 0) usage(argv);
    }

    flint_set_num_threads(nthreads);

    flint_printf("PTHREADS taylor shift (%d bits, length %d, %d threads)...", 
            BITS, LEN, nthreads);
    fflush(stdout);

    /* generate input */

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(check);
    fmpz_init(c);

    do 
    {
        fmpz_poly_randtest(f, state, LEN, BITS);
    } 
    while (fmpz_poly_length(f) != LEN);

    fmpz_randbits(c, state, BITS);

    /* warm up */

    fmpz_poly_taylor_shift_multi_mod(g, f, c);

    /* actual timing */

    timeit_start(timer);
    fmpz_poly_taylor_shift_multi_mod(g, f, c);
    timeit_stop(timer);

    /* cool down and check results */

    fmpz_poly_taylor_shift_multi_mod(g, f, c);
    fmpz_poly_taylor_shift(check, f, c);

    if (!fmpz_poly_equal(g, check))
    {
        flint_printf("FAIL\n");
        result = 1;
    }

    if (result == 0)
        flint_printf("PASS\n");

    /* show time */
    flint_printf("\n");
    flint_printf("  cpu: %lf ms\n", ((double)timer->cpu));
    flint_printf(" wall: %lf ms\n", ((double)timer->wall));

    /* clean up */
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(check);
    fmpz_clear(c);
    FLINT_TEST_CLEANUP(state);

    return result;
}
