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

    flint_printf("randtest_kron....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int j;
        fmpz_spoly_t a;
        fmpz_t degree, abs_degree, bits, lim, h;
        slong height, terms, limit, vars;
                
        fmpz_init(h);
        fmpz_init(lim);
        fmpz_init(bits);
        fmpz_init(degree);
        fmpz_init(abs_degree);
        fmpz_randtest(degree, state, 100);
        fmpz_abs(abs_degree, degree);
        
        terms = n_randint(state, 100);
        height = 200;
        limit = fmpz_bits(abs_degree) + 1;
        vars = 2;

        fmpz_spoly_init(a);
        fmpz_spoly_randtest_kron(a, state, terms, abs_degree, limit, height, vars);

        result = (terms == a->length || terms - fmpz_get_si(abs_degree) >= 1);
        if(!result)
        {
            flint_printf("FAIL (undesired length):\n");
            flint_printf("DESIRED: %lu RECEIVED: %lu and %d\n", terms, a->length, fmpz_cmp_si(abs_degree, terms));
        }

        /*
        result = (fmpz_equal(abs_degree, a->expons) || terms == 0);
        if (!result)
        {
            flint_printf("FAIL (undesired degree):\n");
            flint_printf("DESIRED: ");
            fmpz_print(abs_degree);
            flint_printf(" RECEIVED: ");
            fmpz_print(a->expons);
            flint_printf("\n");
        }
        */

        _fmpz_vec_height(bits, a->coeffs, a->length);
        result = (abs(height - fmpz_bits(bits)) < 15 || terms == 0);
        if(!result)
        {
            flint_printf("FAIL (undesired height):\n");
            flint_printf("DESIRED: %lu RECEIVED: %lu\n", height, fmpz_bits(bits));
        }

        result = 1;
        for(j = 0; j < a->length -1; ++j)
        {
            if(fmpz_cmp(a->expons + j, a->expons + j + 1) <= 0)
                result = 0;
        }
        
        if(!result)
        {
            flint_printf("FAIL (unsorted fmpz_spoly):\n");
            fmpz_spoly_print(a), flint_printf("\n\n");
            flint_printf("terms: %wd\n", terms);
            flint_printf("degree: "); fmpz_print(degree); flint_printf("\n");
            flint_printf("height: %wd\n", height);
            {
                slong j;
                for (j = 0; j < fmpz_spoly_terms(a); ++j)
                {
                    fmpz_print(a->expons + j); flint_printf("\n");
                }
            }
            abort();
        }

        fmpz_set_si(h, 2);
        fmpz_set_si(lim, 2);
        fmpz_pow_ui(h, h, height);
        fmpz_pow_ui(lim, lim, limit);

        flint_printf("\n");
        fmpz_spoly_print(a), flint_printf("\n");
        flint_printf("Degree: "), fmpz_print(degree), flint_printf("\n"); 
        flint_printf("Height: %w, ", height), fmpz_print(h), flint_printf("\n"); 
        flint_printf("Limit: %w, ", limit), fmpz_print(lim), flint_printf("\n");
        flint_printf("Terms: %w", terms), flint_printf("\n"); 
        flint_printf("Vars: %w", vars), flint_printf("\n"); 

        fmpz_spoly_clear(a);
        fmpz_clear(degree);
        fmpz_clear(abs_degree);
        fmpz_clear(bits);
    }



    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
