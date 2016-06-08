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
    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_spoly.h"

int
fmpz_spoly_equal_fmpz_poly(const fmpz_spoly_t spoly, const fmpz_poly_t dpoly)
{
    if (fmpz_spoly_is_zero(spoly) && fmpz_poly_is_zero(dpoly))
    {
        return 1;
    }
    else if (! fmpz_equal_si(fmpz_spoly_degree_ptr(spoly), fmpz_poly_degree(dpoly)))
    {
        return 0;
    }
    else
    {
        slong si = 0, di = fmpz_poly_degree(dpoly);
        while (1)
        {
            if (! fmpz_equal_si(fmpz_spoly_get_term_expon_ptr(spoly, si), di))
            {
                return 0;
            }

            if (! fmpz_equal(fmpz_spoly_get_term_coeff_ptr(spoly, si), 
                             fmpz_poly_get_coeff_ptr(dpoly, di)))
            {
                return 0;
            }

            ++si;
            do
            {
                --di;
            }
            while (di >= 0 && fmpz_is_zero(fmpz_poly_get_coeff_ptr(dpoly, di)));

            if (di < 0 && si == fmpz_spoly_terms(spoly)) return 1;
            else if (di < 0 || si == fmpz_spoly_terms(spoly)) return 0;
        }
    }
}
