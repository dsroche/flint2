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

  Authored 2015 by A. Whitman Groves; US Government work in the public domain.

******************************************************************************/

#include "flint.h"
#include "fmpz_sparse.h"
#include "long_extras.h"

void check_degree_si(const fmpz_sparse_t poly, slong deg)
{
    fmpz_t res;
    fmpz_init(res);

    if (! fmpz_equal_si(fmpz_sparse_degree_ptr(poly), deg)) 
    {
        flint_printf("FAIL (small degree_ptr): %wd\n", deg);
        fmpz_sparse_print(poly), flint_printf("\n\n");
        abort();
    }

    fmpz_sparse_degree(res, poly);
    if (! fmpz_equal_si(res, deg)) 
    {
        flint_printf("FAIL (small degree): %wd\n", deg);
        fmpz_sparse_print(poly), flint_printf("\n\n");
        abort();
    }

    if (fmpz_sparse_degree_si(poly) != deg)
    {
        flint_printf("FAIL (small degree_si): %wd\n", deg);
        fmpz_sparse_print(poly), flint_printf("\n\n");
        abort();
    }

    fmpz_clear(res);
}

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("degree....");
    fflush(stdout);

    {
        fmpz_sparse_t poly;
        fmpz_sparse_init(poly);

        check_degree_si(poly, WORD(-1));

        fmpz_sparse_clear(poly);
    }

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        slong deg;
      
        fmpz_sparse_t poly;
        fmpz_sparse_init(poly);

        deg = z_randtest(state);
        fmpz_sparse_set_coeff_si_si(poly, z_randtest_not_zero(state), deg);

        check_degree_si(poly, deg);

        fmpz_sparse_clear(poly);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
