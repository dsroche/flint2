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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("rem_cyc....");
    fflush(stdout);

    /* Check correctness and aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t poly, divisor, res_reg, res_alias, res_check;
        ulong deg, n;
        mp_limb_t mod;
        
        /* ensure modulus is at least 2 */
        do
        { mod = n_randtest_not_zero(state); }
        while (mod < 2);

        nmod_poly_init(poly, mod);
        nmod_poly_init(divisor, mod);
        nmod_poly_init(res_reg, mod);
        nmod_poly_init(res_alias, mod);
        nmod_poly_init(res_check, mod);

        deg = n_randint(state, 100);
        nmod_poly_randtest(poly, state, deg);
        n = n_randint(state, deg+10) + 1; /* allow mod to be larger than degree */

        /* divisor = x^n - 1 */
        nmod_poly_set_coeff_ui(divisor, n, UWORD(1));
        nmod_poly_set_coeff_ui(divisor, WORD(0), mod-1);

        nmod_poly_set(res_alias, poly);

        nmod_poly_rem(res_check, poly, divisor);
        nmod_poly_rem_cyc(res_reg, poly, n);
        nmod_poly_rem_cyc(res_alias, res_alias, n);

        if (! nmod_poly_equal(res_check, res_reg))
        {
            flint_printf("FAIL:\n");
            flint_printf("regular %wu %wu\n\n", n, mod);
            nmod_poly_print(poly); flint_printf("\n\n");
            nmod_poly_print(res_reg); flint_printf("\n\n");
            nmod_poly_print(res_check); flint_printf("\n\n");
            abort();
        }

        if (! nmod_poly_equal(res_check, res_alias))
        {
            flint_printf("FAIL:\n");
            flint_printf("alias %wu %wu\n\n", n, mod);
            nmod_poly_print(poly); flint_printf("\n\n");
            nmod_poly_print(res_alias); flint_printf("\n\n");
            nmod_poly_print(res_check); flint_printf("\n\n");
            abort();
        }

        nmod_poly_clear(poly);
        nmod_poly_clear(divisor);
        nmod_poly_clear(res_reg);
        nmod_poly_clear(res_alias);
        nmod_poly_clear(res_check);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
