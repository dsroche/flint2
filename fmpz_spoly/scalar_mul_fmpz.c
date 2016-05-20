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
#include "fmpz_spoly.h"

void
fmpz_spoly_scalar_mul(fmpz_spoly_t poly1, const fmpz_spoly_t poly2,
                          const fmpz_t x)
{
    /* Either scalar or input poly is zero */
    if (fmpz_is_zero(x) || (poly2->length == 0))
    {
        fmpz_spoly_zero(poly1);
        return;
    }

    if (fmpz_is_one(x))
    {
      fmpz_spoly_set(poly1, poly2);
      return;
    }

    if (poly1 == poly2)
    {
      fmpz_spoly_t temp;
      fmpz_spoly_init(temp);
      fmpz_spoly_scalar_mul(temp, poly2, x);
      poly1->length = poly2->length;
      fmpz_spoly_set(poly1, temp);
      fmpz_spoly_clear(temp);
    }
    else
    {
      fmpz_spoly_init2(poly1, poly2->length);
      poly1->length = poly2->length;
    
      _fmpz_vec_scalar_mul_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
      _fmpz_vec_set(poly1->expons, poly2->expons, poly2->length);
    }
}
