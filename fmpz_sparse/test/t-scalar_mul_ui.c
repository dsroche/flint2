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

    Copyright (C) 2009 William Hart

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
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_ui....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b;
        fmpz_t c;
        ulong n = n_randtest(state);

        fmpz_init(c);
        fmpz_randtest(c, state, 100);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_randtest(a, state, n_randint(state, 100), c, 200);

        fmpz_sparse_scalar_mul_ui(b, a, n);
        fmpz_sparse_scalar_mul_ui(a, a, n);

        result = (fmpz_sparse_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_sparse_print(a), flint_printf("\n\n");
            fmpz_sparse_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_clear(c);
    }

    /* Check (a*n1)*n2 = a*(n1*n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t a, b, c;
        fmpz_t d;
        ulong n1 = n_randbits(state, FLINT_BITS / 2);
        ulong n2 = n_randbits(state, FLINT_BITS / 2);

        fmpz_init(d);
        fmpz_randtest(d, state, 100);

        fmpz_sparse_init(a);
        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_randtest(a, state, n_randint(state, 100), d, 200);

        fmpz_sparse_scalar_mul_ui(b, a, n1);
        fmpz_sparse_scalar_mul_ui(c, b, n2);
        fmpz_sparse_scalar_mul_ui(b, a, n1 * n2);

        result = (fmpz_sparse_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL n1 = %wu, n2 = %wu:\n", n1, n2);
            fmpz_sparse_print(a), flint_printf("\n\n");
            fmpz_sparse_print(b), flint_printf("\n\n");
            fmpz_sparse_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_sparse_clear(a);
        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        fmpz_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}