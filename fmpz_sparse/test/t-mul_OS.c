/*============================================================================

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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1101 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 109 William Hart
    Copyright (C) 1010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result, count, total;
    FLINT_TEST_INIT(state);

    flint_printf("mul_OS....");
    fflush(stdout);

    count = 0;
    total = 0;

    /* Check aliasing of a and b */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b, c, f;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(d, state, 30);
        fmpz_randtest(e, state, 30);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_init(f);
        fmpz_sparse_randtest(b, state, n_randint(state, 15), d, 30);
        fmpz_sparse_randtest(c, state, n_randint(state, 15), e, 30);
        
        fmpz_sparse_mul_OS(a, state, b, c);
        fmpz_sparse_mul_classical(f, b, c);

        result = (fmpz_sparse_equal(a, f));
        if (!result)
        {
          flint_printf("FAIL PHASE 1:\n");
          flint_printf("\non the %w try\n", i);
          fmpz_sparse_print(a), flint_printf("\n\n");
          fmpz_sparse_print(f), flint_printf("\n\n");
          count++;
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_sparse_clear(f);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    total += i;

    /* Check aliasing of a and c */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b, c;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest(d, state, 30);
        fmpz_randtest(e, state, 30);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_randtest(b, state, n_randint(state, 15), d, 30);
        fmpz_sparse_randtest(c, state, n_randint(state, 15), e, 30);

        fmpz_sparse_mul_classical(a, b, c);
        fmpz_sparse_mul_OS(b, state, b, c);

        result = (fmpz_sparse_equal(a, b));
        if (!result)
        {
          flint_printf("FAIL PHASE 2:\n");
          flint_printf("\non the %w try\n", i + total);
          fmpz_sparse_print(a), flint_printf("\n\n");
          fmpz_sparse_print(c), flint_printf("\n\n");
          count++;
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    total += i;

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a1, a2, b, c, d;
        fmpz_t e, f, g;

        fmpz_init(e);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest(e, state, 30);
        fmpz_randtest(f, state, 30);
        fmpz_randtest(g, state, 30);

        fmpz_sparse_init(a1);
        fmpz_sparse_init(a2);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_init(d);
        fmpz_sparse_randtest(b, state, n_randint(state, 15), e, 30);
        fmpz_sparse_randtest(c, state, n_randint(state, 15), f, 30);
        fmpz_sparse_randtest(d, state, n_randint(state, 15), g, 30);

        fmpz_sparse_mul_OS(a1, state, b, c);
        fmpz_sparse_mul_classical(a2, b, d);
        fmpz_sparse_add(a1, a1, a2);

        fmpz_sparse_add(c, c, d);
        fmpz_sparse_mul_OS(a2, state, b, c);

        result = (fmpz_sparse_equal(a1, a2));
        if (!result)
        {
          flint_printf("FAIL PHASE 3:\n");
          flint_printf("\non the %w try\n", i + total);
          fmpz_sparse_print(a1), flint_printf("\n\n");
          fmpz_sparse_print(a2), flint_printf("\n\n");
          count++;
        }

        fmpz_sparse_clear(a1);
        fmpz_sparse_clear(a2);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_sparse_clear(d);
        fmpz_clear(e);
        fmpz_clear(f);
        fmpz_clear(g);
    }

    total += i;

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    flint_printf("failed %w times out of %w\n", count, total);
    return 0;
}
