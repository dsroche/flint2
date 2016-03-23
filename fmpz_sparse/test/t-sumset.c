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
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("sumset....");
    fflush(stdout);

    
    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t b, c;
        fmpz_t d, e;
        fmpz * vec1, * vec2;
        slong len1, len2;

        fmpz_init(d);
        fmpz_init(e);
        fmpz_randtest(d, state, 10);
        fmpz_randtest(e, state, 10);

        fmpz_sparse_init(b);
        fmpz_sparse_init(c);
        fmpz_sparse_randtest(b, state, n_randint(state, 10), d, 10);
        fmpz_sparse_randtest(c, state, n_randint(state, 10), e, 10);

        vec1 = NULL;
        vec2 = NULL;

        flint_printf("poly1: ");
        fmpz_sparse_print(b);
        flint_printf("\npoly2: ");
        fmpz_sparse_print(c);
        flint_printf("\n");

        len1 = fmpz_sparse_sumcheck(&vec1, b, c);
        len2 = fmpz_sparse_sumset(&vec2, state, b, c);
        
        result = (len1 == len2 && _fmpz_vec_equal(vec1, vec2, len1));
        if (!result)
        {
          flint_printf("FAIL:\nvector_1: ");
          _fmpz_vec_print(vec1, len1), flint_printf("\n\nvector_2: ");
          _fmpz_vec_print(vec2, len2), flint_printf("\n\n");
          abort();
        }

        fmpz_sparse_clear(b);
        fmpz_sparse_clear(c);
        if(len1 > 1)
          _fmpz_vec_clear(vec1, len1);
        if(len2 > 1)
          _fmpz_vec_clear(vec2, len2);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
