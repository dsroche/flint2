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

#ifndef FMPZ_SPOLY_PTIMER
#define FMPZ_SPOLY_PTIMER

#include "profiler.h"

#define PTIMER_DECLARE(TNAME, NUM) \
    const int num_##TNAME = (NUM) + 1; \
    timeit_struct timers_##TNAME[(NUM) + 1]; \
    timeit_t ptimer1_##TNAME, ptimer2_##TNAME; \
    const char* titles_##TNAME[(NUM) + 1]; \
    int counter_##TNAME = 0; \
    int do_##TNAME = 0;

#define PTIMER_BEGIN(TNAME, TITLE) do { \
    titles_##TNAME[0] = "TOTAL"; \
    counter_##TNAME = 0; \
    titles_##TNAME[++counter_##TNAME] = (TITLE); \
    timeit_start(ptimer1_##TNAME); \
    timeit_start(ptimer2_##TNAME); \
} while(0)

#define PTIMER_NEXT(TNAME, TITLE) do { \
    timeit_stop(ptimer2_##TNAME); \
    if (do_##TNAME) \
    { \
        timers_##TNAME[counter_##TNAME].cpu += ptimer2_##TNAME->cpu; \
        timers_##TNAME[counter_##TNAME].wall += ptimer2_##TNAME->wall; \
    } \
    titles_##TNAME[++counter_##TNAME] = (TITLE); \
    timeit_start(ptimer2_##TNAME); \
} while(0)

#define PTIMER_END(TNAME) do { \
    timeit_stop(ptimer2_##TNAME); \
    timeit_stop(ptimer1_##TNAME); \
    if (do_##TNAME) \
    { \
        timers_##TNAME[counter_##TNAME].cpu += ptimer2_##TNAME->cpu; \
        timers_##TNAME[counter_##TNAME].wall += ptimer2_##TNAME->wall; \
        timers_##TNAME[0].cpu += ptimer1_##TNAME->cpu; \
        timers_##TNAME[0].wall += ptimer1_##TNAME->wall; \
    } \
} while(0)


#define PTIMER_EXTERN(TNAME) \
    extern const int num_##TNAME; \
    extern timeit_struct timers_##TNAME[]; \
    extern int do_##TNAME; \
    extern int counter_##TNAME; \
    extern const char* titles_##TNAME[];

#define PTIMER_CLEAR(TNAME) do { \
    int pti; \
    for (pti = 0; pti < num_##TNAME; ++pti) \
    { \
        timers_##TNAME[pti].cpu = 0; \
        timers_##TNAME[pti].wall = 0; \
    } \
} while(0)

#define PTIMER_ENABLE(TNAME) do { \
    do_##TNAME = 1; \
} while(0)

#define PTIMER_DISABLE(TNAME) do { \
    do_##TNAME = 0; \
} while(0)

#define PTIMER_PRINT(TNAME, NRUNS) do { \
    int pti; \
    double totc = ((double) timers_##TNAME[0].cpu) / (NRUNS); \
    double totw = ((double) timers_##TNAME[0].wall) / (NRUNS); \
    flint_printf("\n==== TIMES FOR " #TNAME " ====\n"); \
    for (pti = 0; pti <= counter_##TNAME; ++pti) \
    { \
        double tc = ((double) timers_##TNAME[pti].cpu) / (NRUNS); \
        double tw = ((double) timers_##TNAME[pti].wall) / (NRUNS); \
        flint_printf("%s\n", titles_##TNAME[pti]); \
        flint_printf("  cpu: %lf ms = %.2lf%%\n", tc, tc * 100 / totc); \
        flint_printf(" wall: %lf ms = %.2lf%%\n", tw, tw * 100 / totw); \
    } \
} while(0)

#endif
