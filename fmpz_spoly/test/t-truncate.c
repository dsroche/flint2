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

#include <stdlib.h>
#include <stdio.h>
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

    flint_printf("truncate....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_t a, b;
        fmpz_t c, d, temp;

        fmpz_init(temp);

        fmpz_init(c);
        fmpz_randtest_unsigned(c, state, 10);

        fmpz_spoly_init(a);
        fmpz_spoly_randtest(a, state, n_randint(state, 10), c, 20);

        fmpz_init(d);
        fmpz_randtest_unsigned(d, state, 8);

        while(fmpz_sgn(d) < 0)
        {
            fmpz_randtest_unsigned(d, state, 8);
        }

        fmpz_spoly_init(b);
        fmpz_spoly_set(b, a);

        fmpz_spoly_truncate(b, d);
        fmpz_spoly_degree(temp, b);

        result = ((fmpz_cmp(d, temp) >= 0 || fmpz_sgn(d) < 0) && (_fmpz_spoly_index(b,d) == 0 || _fmpz_spoly_index(b,d) == -1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_spoly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_spoly_print(b), flint_printf("\n\n");
            flint_printf("deg = "), fmpz_print(d), flint_printf("\n\n");
            flint_printf("index = %w", _fmpz_spoly_index(b,d)), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
