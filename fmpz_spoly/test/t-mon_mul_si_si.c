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

    flint_printf("mon_mul_si_si....");
    fflush(stdout);

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, m, n;
        slong i = 0, j = 0;

        fmpz_init(m);
        fmpz_init(n);
        fmpz_randtest_unsigned(m, state, 20);
        fmpz_randtest_unsigned(n, state, 20);

        j = fmpz_get_si(m);
        i = fmpz_get_si(n);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 20);
        
        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 10), c, 20);

        fmpz_spoly_mon_mul_si_si(b, a, j, i);
        fmpz_spoly_mon_mul_si_si(a, a, j, i);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL: 1\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_clear(m);
        fmpz_clear(n);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, d, e;
        slong f, g, i;

        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest_unsigned(c, state, 10);
        fmpz_randtest_unsigned(d, state, 10);
        fmpz_randtest_unsigned(e, state, 10);

        f = fmpz_get_si(e);
        g = fmpz_get_si(d);

        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 10), c, 20);

        fmpz_spoly_mon_mul_si_si(b, a, g, f);
        fmpz_spoly_scalar_mul_si(a, a, g);

        for(i = 0; i < a->length; i++)
        {
            fmpz_add_ui(a->expons + i, a->expons + i, f);
        }

        _fmpz_spoly_normalise(a);

        result = (fmpz_spoly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL: 2\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            fmpz_print(e), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
