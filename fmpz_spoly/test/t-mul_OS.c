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

  Authored 2016 by A. Whitman Groves; US Government work in the public domain.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mul_OS....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c, f;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(d, state, 20);
        fmpz_randtest(e, state, 20);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_init(f);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), d, 100);
        fmpz_spoly_randtest(c, state, n_randint(state, 30), e, 100);
        
        fmpz_spoly_mul_OS(a, state, b, c);
        fmpz_spoly_mul_classical(f, b, c);

        result = (fmpz_spoly_equal(a, f));
        if (!result)
        {
          flint_printf("FAIL PHASE 1:\n");
          flint_printf("\non the %w try\n", i);
          fmpz_spoly_print(a), flint_printf("\n\n");
          fmpz_spoly_print(f), flint_printf("\n\n");
          abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_spoly_clear(f);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest(d, state, 20);
        fmpz_randtest(e, state, 20);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), d, 100);
        fmpz_spoly_randtest(c, state, n_randint(state, 30), e, 100);

        fmpz_spoly_mul_classical(a, b, c);
        fmpz_spoly_mul_OS(b, state, b, c);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
          flint_printf("FAIL PHASE 2:\n");
          flint_printf("\non the %w try\n", i);
          fmpz_spoly_print(a), flint_printf("\n\n");
          fmpz_spoly_print(c), flint_printf("\n\n");
          abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a1, a2, b, c, d;
        fmpz_t e, f, g;

        fmpz_init(e);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest(e, state, 20);
        fmpz_randtest(f, state, 20);
        fmpz_randtest(g, state, 20);

        fmpz_spoly_init(a1);
        fmpz_spoly_init(a2);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_init(d);
        fmpz_spoly_randtest(b, state, n_randint(state, 30), e, 100);
        fmpz_spoly_randtest(c, state, n_randint(state, 30), f, 100);
        fmpz_spoly_randtest(d, state, n_randint(state, 30), g, 100);

        fmpz_spoly_mul_OS(a1, state, b, c);
        fmpz_spoly_mul_classical(a2, b, d);
        fmpz_spoly_add(a1, a1, a2);

        fmpz_spoly_add(c, c, d);
        fmpz_spoly_mul_OS(a2, state, b, c);

        result = (fmpz_spoly_equal(a1, a2));
        if (!result)
        {
          flint_printf("FAIL PHASE 3:\n");
          flint_printf("\non the %w try\n", i);
          fmpz_spoly_print(a1), flint_printf("\n\n");
          fmpz_spoly_print(a2), flint_printf("\n\n");
          abort();
        }

        fmpz_spoly_clear(a1);
        fmpz_spoly_clear(a2);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_spoly_clear(d);
        fmpz_clear(e);
        fmpz_clear(f);
        fmpz_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
