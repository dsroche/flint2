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
#include "fmpz_spoly.h"

void
fmpz_spoly_scalar_mul(fmpz_spoly_t poly1, const fmpz_spoly_t poly2,
                          const fmpz_t x)
{
    /* Either scalar or input poly is zero */
    if (fmpz_is_zero(x) || (poly2->length == 0))
    {
        fmpz_spoly_zero(poly1);
    }
    else if (fmpz_is_one(x))
    {
        fmpz_spoly_set(poly1, poly2);
    }
    else
    {
        slong t2 = fmpz_spoly_terms(poly2);

        if (poly1 != poly2)
        {
            _fmpz_spoly_reserve(poly1, t2);
            _fmpz_vec_set(poly1->expons, poly2->expons, t2);
        }

        _fmpz_vec_scalar_mul_fmpz(poly1->coeffs, poly2->coeffs, t2, x);

        _fmpz_spoly_set_length(poly1, t2);
    }
}
