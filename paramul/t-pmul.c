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
#include "omp.h"

#define NUMEX (10)
#define MINWALL (1000)
#define NUMALG (2)

void fmpz_poly_pmul(fmpz_poly_t, const fmpz_poly_t, const fmpz_poly_t);

void usage(char** argv)
{
    flint_printf("usage: %s len cbits [num_threads]\n", argv[0]);
    abort();
}

#define START_ALG(NLOOPS) \
{ \
    slong loopshere = (NLOOPS); \
    timeit_start(timer); \
    for (loop = 0; loop < loopshere; ++loop) \
    { \
        for (ex = 0; ex < NUMEX; ++ex) \
        { 

#define END_ALG(CPU, WALL) \
        } \
    } \
    timeit_stop(timer); \
    (CPU) = (double)timer->cpu / (loopshere * NUMEX); \
    (WALL) = (double)timer->wall / (loopshere * NUMEX); \
}

int main(int argc, char** argv)
{
    timeit_t timer;
    fmpz_poly_struct in1[NUMEX];
    fmpz_poly_struct in2[NUMEX];
    fmpz_poly_struct res[NUMALG];
    slong len;
    ulong cbits;
    slong ex, alg, loop;
    slong nloops[NUMALG];
    double cputimes[NUMALG + 1];
    double walltimes[NUMALG + 1];
    const char* algnames[] = {"parallel mul", "regular mul"};
    
    FLINT_TEST_INIT(state);

    if (argc < 3 || argc > 4)
        usage(argv);

    len = strtoul(argv[1], NULL, 10);
    cbits = strtoul(argv[2], NULL, 10);

    if (argc >= 4)
    {
        int nthreads = atoi(argv[3]);
        omp_set_num_threads(nthreads);
    }

    /* initialize random examples */
    for (ex = 0; ex < NUMEX; ++ex)
    {
        fmpz_poly_init(in1 + ex);
        fmpz_poly_init(in2 + ex);
        fmpz_poly_randtest(in1 + ex, state, len, cbits);
        fmpz_poly_randtest(in2 + ex, state, len, cbits);
    }

    for (alg = 0; alg < NUMALG; ++alg)
    {
        fmpz_poly_init(res + alg);
    }

    /* warm up and check times */

    alg = 0;

    START_ALG(1)
        fmpz_poly_pmul(res + alg, in1 + ex, in2 + ex);
    END_ALG(cputimes[alg], walltimes[alg])
    nloops[alg] = MINWALL / timer->wall + 1;
    flint_printf("Prelim time for %s: %.2f msec\n", algnames[alg], walltimes[alg]);
    ++alg;

    START_ALG(1)
        fmpz_poly_mul(res + alg, in1 + ex, in2 + ex);
    END_ALG(cputimes[alg], walltimes[alg])
    nloops[alg] = MINWALL / timer->wall + 1;
    flint_printf("Prelim time for %s: %.2f msec\n", algnames[alg], walltimes[alg]);
    ++alg;

    flint_printf("\n");

    /* run critical timing */
    alg = 0;

    do
    {
        START_ALG(nloops[alg])
            fmpz_poly_pmul(res + alg, in1 + ex, in2 + ex);
        END_ALG(cputimes[alg], walltimes[alg])
        nloops[alg] *= 2;
    } while (timer->wall < MINWALL);
    flint_printf("Final time for %s: %.2f msec wall, %.2f msec cpu\n", algnames[alg], walltimes[alg], cputimes[alg]);
    ++alg;

    do
    {
        START_ALG(nloops[alg])
            fmpz_poly_mul(res + alg, in1 + ex, in2 + ex);
        END_ALG(cputimes[alg], walltimes[alg])
        nloops[alg] *= 2;
    } while (timer->wall < MINWALL);
    flint_printf("Final time for %s: %.2f msec wall, %.2f msec cpu\n", algnames[alg], walltimes[alg], cputimes[alg]);
    ++alg;

    /* cool down and check results */
    START_ALG(1)
        fmpz_poly_pmul(res + 0, in1 + ex, in2 + ex);
        fmpz_poly_mul(res + 1, in1 + ex, in2 + ex);
        if (!fmpz_poly_equal(res + 0, res + 1))
        {
            flint_printf("ERROR in example %wd:\n", ex);
            fmpz_poly_print(in1 + ex); flint_printf("\n\n");
            fmpz_poly_print(in2 + ex); flint_printf("\n\n");
            fmpz_poly_print(res + 0); flint_printf("\n\n");
            fmpz_poly_print(res + 1); flint_printf("\n\n");
            return 1;
        }
    END_ALG(cputimes[alg], walltimes[alg])

    /* clean up */
    for (ex = 0; ex < NUMEX; ++ex)
    {
        fmpz_poly_clear(in1 + ex);
        fmpz_poly_clear(in2 + ex);
    }

    for (alg = 0; alg < NUMALG; ++alg)
    {
        fmpz_poly_clear(res + alg);
    }

    FLINT_TEST_CLEANUP(state);

    return 0;
}
