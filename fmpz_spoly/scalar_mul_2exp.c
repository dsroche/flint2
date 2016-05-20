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
fmpz_spoly_scalar_mul_2exp(fmpz_spoly_t res, const fmpz_spoly_t poly,
                          ulong exp)
{
    if (poly->length == 0)
    {
        fmpz_spoly_zero(res);
        return;
    }

    if (res == poly)
    {
      fmpz_spoly_t temp;
      fmpz_spoly_init(temp);
      fmpz_spoly_scalar_mul_2exp(temp, poly, exp);
      res->length = poly->length;
      fmpz_spoly_set(res, temp);
      fmpz_spoly_clear(temp);
    }
    else
    {
      fmpz_spoly_init2(res, poly->length);
      res->length = poly->length;
    
      _fmpz_vec_scalar_mul_2exp(res->coeffs, poly->coeffs, poly->length, exp);
      _fmpz_vec_set(res->expons, poly->expons, poly->length);
      _fmpz_spoly_normalise(res);
    }
}
