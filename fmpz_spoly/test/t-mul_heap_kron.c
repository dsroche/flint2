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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

  Authored 2015 by A. Whitman Groves; US Government in the public domain.

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

    flint_printf("mul_heap_kron....");
    fflush(stdout);

    
    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e;
        ulong limit, vars;
        slong terms1, terms2;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(d, state, 50);
        fmpz_randtest_unsigned(e, state, 50);

        limit = fmpz_bits(d) + fmpz_bits(e)+1;

        vars = n_randint(state, 6) + 1;
        
        terms1 = n_randint(state, 50);
        terms2 = n_randint(state, 50);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest_kron(b, state, terms1, d, 200, limit, vars);
        fmpz_spoly_randtest_kron(c, state, terms2, e, 200, limit, vars);
        
        fmpz_spoly_mul_heaps(a, b, c);
        fmpz_spoly_mul_heaps(b, b, c);
        
        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
          flint_printf("FAIL PHASE 1:\n");
          fmpz_spoly_print(a), flint_printf("\n\n");
          fmpz_spoly_print(b), flint_printf("\n\n");
          abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b, c;
        fmpz_t d, e;
        ulong limit, vars;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(d, state, 50);
        fmpz_randtest_unsigned(e, state, 50);

        limit = fmpz_bits(d) + fmpz_bits(e);

        vars = 2;

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_randtest_kron(b, state, n_randint(state, 50), d, 200, limit, vars);
        fmpz_spoly_randtest_kron(c, state, n_randint(state, 50), e, 200, limit, vars);

        fmpz_spoly_mul_heaps(a, b, c);
        fmpz_spoly_mul_heaps(c, b, c);

        result = (fmpz_spoly_equal(a, c));
        if (!result)
        {
          flint_printf("FAIL PHASE 2:\n");
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
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a1, a2, b, c, d;
        fmpz_t e, f, g;
        ulong limit, vars;

        fmpz_init(e);
        fmpz_init(f);
        fmpz_init(g);
        fmpz_randtest_unsigned(e, state, 50);
        fmpz_randtest_unsigned(f, state, 50);
        fmpz_randtest_unsigned(g, state, 50);

        limit = fmpz_bits(e) + fmpz_bits(f) + fmpz_bits(g);

        vars = 2;

        fmpz_spoly_init(a1);
        fmpz_spoly_init(a2);
        fmpz_spoly_init(b);
        fmpz_spoly_init(c);
        fmpz_spoly_init(d);
        fmpz_spoly_randtest_kron(b, state, n_randint(state, 50), e, 200, limit, vars);
        fmpz_spoly_randtest_kron(c, state, n_randint(state, 50), f, 200, limit, vars);
        fmpz_spoly_randtest_kron(d, state, n_randint(state, 50), g, 200, limit, vars);

        fmpz_spoly_mul_heaps(a1, b, c);
        fmpz_spoly_mul_heaps(a2, b, d);
        fmpz_spoly_add(a1, a1, a2);

        fmpz_spoly_add(c, c, d);
        fmpz_spoly_mul_heaps(a2, b, c);

        result = (fmpz_spoly_equal(a1, a2));
        if (!result)
        {
          flint_printf("FAIL PHASE 3:\n");
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
