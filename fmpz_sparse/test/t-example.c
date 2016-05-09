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

    Copyright (C) 2010, 2012 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_sparse.h"
#include "ulong_extras.h"

int
main(void)
{
  
    int i, result = 1;
    FLINT_TEST_INIT(state);
    
    
    flint_printf("evaluate_fmpz_vec_fast....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_sparse_t poly1, poly2, res1, res2;

        fmpz_sparse_init(poly1);
        fmpz_sparse_init(poly2);
        fmpz_sparse_init(res1);
        fmpz_sparse_init(res2);

        _fmpz_sparse_reserve(poly1, 4);
        _fmpz_sparse_reserve(poly2, 3);
        _fmpz_sparse_reserve(res1, 5);
        _fmpz_sparse_reserve(res2, 5);

        fmpz_set_si(poly1->coeffs + 0, 20);
        fmpz_set_si(poly1->coeffs + 1, 65);
        fmpz_set_si(poly1->coeffs + 2, 16);
        fmpz_set_si(poly1->coeffs + 3, 26);
        fmpz_set_si(poly1->expons + 0, 4913);
        fmpz_set_si(poly1->expons + 1, 3631);
        fmpz_set_si(poly1->expons + 2, 2520);
        fmpz_set_si(poly1->expons + 3, 1238);
        
        fmpz_set_si(poly2->coeffs + 0, 60);
        fmpz_set_si(poly2->coeffs + 1, -48);
        fmpz_set_si(poly2->coeffs + 2, 78);
        fmpz_set_si(poly2->expons + 0, 4316);
        fmpz_set_si(poly2->expons + 1, 1923);
        fmpz_set_si(poly2->expons + 2, 641);

        poly1->length = 4;
        poly2->length = 3;

        _fmpz_sparse_normalise(poly1);
        _fmpz_sparse_normalise(poly2);

        fmpz_set_ui(res1->coeffs + 0, 3900);
        fmpz_set_ui(res1->coeffs + 1, 1200);
        fmpz_set_ui(res1->coeffs + 2, 5070);
        fmpz_set_ui(res1->coeffs + 3, 2028);
        fmpz_set_si(res1->coeffs + 4, -768);
        fmpz_set_ui(res1->expons + 0, 7947);
        fmpz_set_ui(res1->expons + 1, 9229);
        fmpz_set_ui(res1->expons + 2, 4272);
        fmpz_set_ui(res1->expons + 3, 1879);
        fmpz_set_ui(res1->expons + 4, 4443);

        res1->length = 5;

        _fmpz_sparse_normalise(res1);

        fmpz_sparse_mul_OS(res2, state, poly1, poly2);
        
        flint_printf("\npoly1: "), fmpz_sparse_print(poly1);
        flint_printf("\npoly2: "), fmpz_sparse_print(poly2);
        flint_printf("\nres1: "), fmpz_sparse_print(res1);
        flint_printf("\nres2: "), fmpz_sparse_print(res2);
        flint_printf("\n");


        result = fmpz_sparse_equal(res1, res2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            abort();
        }

        fmpz_sparse_clear(poly1);
        fmpz_sparse_clear(poly2);
        fmpz_sparse_clear(res1);
        fmpz_sparse_clear(res2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
