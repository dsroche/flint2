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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_spoly.h"

int
fmpz_spoly_equal_fmpz_poly(const fmpz_spoly_t spoly,
                          const fmpz_poly_t dpoly)
{
    fmpz_t i, j, temp, foo;
    slong terms;

    fmpz_init(i);
    fmpz_init(j);
    fmpz_init(temp);
    fmpz_init(foo);

    fmpz_zero(j);

    if (fmpz_spoly_is_zero(spoly) ^ (fmpz_poly_length(dpoly) == 0))
    {
        return 0;
    }

    fmpz_set_si(i, fmpz_poly_degree(dpoly));
    fmpz_spoly_degree(temp, spoly);

    if(fmpz_cmp(temp, i) != 0)
    {
        return 0;
    }

    terms = fmpz_spoly_terms(spoly);

    while(fmpz_cmp_si(i, 0) >= 0)
    {
        fmpz_poly_get_coeff_fmpz(temp, dpoly, fmpz_get_si(i));

        if(!fmpz_is_zero(temp))
        {
            if(fmpz_cmp_si(j, terms) > 0)
            {
                return 0;
            }

            fmpz_spoly_get_coeff(foo, spoly, i);

            if(fmpz_cmp(foo, temp) != 0)
            {
                return 0;
            }

            fmpz_add_ui(j, j, 1);
        }

        fmpz_sub_ui(i, i, 1);
    }
    return 1;
}
