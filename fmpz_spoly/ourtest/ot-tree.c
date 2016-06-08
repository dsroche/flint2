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
#include "fmpz_spoly.h"
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
        fmpz_poly_struct ** tree;
        fmpz_mod_poly_t res1, res2;
        fmpz * roots;
        fmpz_t mod;
        slong j, i, npoints;

        fmpz_init(mod);

        /*TODO play with this and change the number of points*/

        npoints = 19;

        roots = _fmpz_vec_init(npoints);

        tree = _fmpz_mod_poly_tree_alloc(npoints);

        fmpz_set_ui(mod, 29);
        
        fmpz_mod_poly_init(res1, mod);
        fmpz_mod_poly_init(res2, mod);

        for(j = 0; j < npoints; j++)
        {
          fmpz_set_ui(roots + j, j);
        }

        _fmpz_mod_poly_tree_build(tree, roots, npoints, mod);

        flint_printf("\n\n");

        /*for(i = 0; i < FLINT_CLOG2(npoints); i++)
        {
          for(j = 0; j < (npoints >> i); j++)
          {
            fmpz_poly_print((tree[i]) + j), flint_printf("\n\n");
          }
        }*/

        for(i = 0, j = WORD(2); j < npoints; i++, j*=2);

        j = j/2 + 1;

        fmpz_mod_poly_fit_length(res1, npoints + 1);
        _fmpz_mod_poly_mul(res1->coeffs, tree[i][0].coeffs, j, tree[i][1].coeffs, npoints + 2 - j, &(res1->p));

        _fmpz_mod_poly_set_length(res1, npoints+1);


        fmpz_mod_poly_print(res1), flint_printf("\n\n");

        _fmpz_mod_poly_build_roots(res2, roots, npoints);

        result = fmpz_mod_poly_equal(res1, res2);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("mod=");
            fmpz_print(mod);
            flint_printf(", npoints=%wd\n\n", npoints);
            flint_printf("\n\n");
            abort();
        }

        fmpz_clear(mod);
        _fmpz_vec_clear(roots, npoints);
        fmpz_mod_poly_clear(res1);
        fmpz_mod_poly_clear(res2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
