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

    flint_printf("interp....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_spoly_sp_interp_basis_t basis;
        fmpz_spoly_sp_interp_eval_t eval;
        fmpz_spoly_t a, b;
        fmpz_t d, h;
        int ret;

        /*random fmpz*/
        fmpz_init(d);
        fmpz_init(h);
        fmpz_randtest(d, state, 20);

        /*random fmpz_spoly*/
        fmpz_spoly_init(a);
        fmpz_spoly_init(b);
        fmpz_spoly_randtest(a, state, n_randint(state, 100), d, 1 + n_randint(state, 200));

        /*initialize interpolation struct and eval random fmpz_spoly*/
        fmpz_spoly_sp_interp_basis_init(basis, state, a->length, 
                fmpz_bits(fmpz_spoly_degree_ptr(a)), fmpz_spoly_height_bits(a));
        fmpz_spoly_sp_interp_eval_init(eval, basis);
        fmpz_spoly_sp_interp_eval(eval, a);
        
        /* TODO FIXME REMOVE */ flint_printf("params: %wd %wd %wd %wd %wd\n", basis->length, basis->cimg_per_round, basis->cimg_needed, basis->eimg_per_round, basis->eimg_needed);
        /*b gets result*/
        ret = fmpz_spoly_sp_interp(b, eval);

        /*check that a == b*/
        result = ret > 0 && fmpz_spoly_equal(a, b);

        if (!result)
        {
            flint_printf("FAIL: %d %wd %wu %wu\n", ret, fmpz_spoly_terms(a), 
                fmpz_bits(fmpz_spoly_degree_ptr(a)), fmpz_spoly_height_bits(a));
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }
        else flint_printf("pass %wd %wu %wu\n", fmpz_spoly_terms(a), fmpz_bits(fmpz_spoly_degree_ptr(a)), fmpz_spoly_height_bits(a));
        /* TODO FIXME REMOVE */

        fmpz_clear(d);
        fmpz_clear(h);
        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
        fmpz_spoly_sp_interp_eval_clear(eval);
        fmpz_spoly_sp_interp_basis_clear(basis);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
