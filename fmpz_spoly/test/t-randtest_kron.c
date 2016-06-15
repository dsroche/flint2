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
        fmpz_t degree, lim, h, maxterms;
        ulong limit, vars;
        slong height, terms;

        fmpz_init(h);
        fmpz_init(lim);
        fmpz_init(maxterms);
        fmpz_init(degree);
        fmpz_randtest(degree, state, 100);
        fmpz_abs(degree, degree);
        
        terms = n_randint(state, 100);
        height = 200;
        limit = fmpz_bits(degree) + 1;
        vars = 2;

        fmpz_spoly_init(a);
        fmpz_spoly_randtest_kron(a, state, terms, degree, height, limit, vars);

        fmpz_set(maxterms, degree);
        fmpz_add_ui(maxterms, maxterms, UWORD(1));
        fmpz_pow_ui(maxterms, maxterms, vars);

        if (fmpz_cmp_ui(maxterms, (ulong)terms) < 0)
        {
            result = fmpz_cmp_ui(maxterms, fmpz_spoly_terms(a)) == 0;
        }
        else
        {
            result = terms == fmpz_spoly_terms(a);
        }

        if(!result)
        {
            flint_printf("FAIL (undesired length):\n");
            flint_printf("DESIRED: min of %wd and ", terms); fmpz_print(maxterms);
            flint_printf(" RECEIVED: %wd\n", fmpz_spoly_terms(a));
            abort();
        }

        fmpz_set(lim, degree);
        for (j = 1; (ulong)j < vars; ++j)
        {
            fmpz_mul_2exp(lim, lim, limit);
            fmpz_add(lim, lim, degree);
        }

        if (terms == 0)
        {
            result = fmpz_equal_si(fmpz_spoly_degree_ptr(a), WORD(-1));
        }
        else
        {
            result = fmpz_equal(lim, fmpz_spoly_degree_ptr(a));
        }
        if (!result)
        {
            flint_printf("FAIL (undesired degree):\n");
            flint_printf("DESIRED: ");
            fmpz_print(lim);
            flint_printf(" RECEIVED: ");
            fmpz_print(fmpz_spoly_degree_ptr(a));
            flint_printf("\n");
            abort();
        }

        if (terms == 0)
        {
            result = fmpz_spoly_height_bits(a) == 0;
        }
        else
        {
            result = fmpz_spoly_height_bits(a) == (ulong)height;
        }
        if(!result)
        {
            flint_printf("FAIL (undesired height):\n");
            flint_printf("DESIRED: %lu RECEIVED: %lu\n", fmpz_spoly_height_bits(a), height);
            abort();
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

        fmpz_spoly_clear(a);
        fmpz_clear(degree);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
