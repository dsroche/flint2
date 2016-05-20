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

#include "fmpz_spoly.h"
#include "fmpz_vec.h"
void 
fmpz_spoly_new_mul_classical(fmpz_spoly_t res, const fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2)
{
  slong i;

  if((poly1->length == 0) || (poly2->length == 0))
  {
    fmpz_spoly_zero(res);
    return;
  }

  if (res == poly1 || res == poly2) 
  {
    fmpz_spoly_t temp;
    fmpz_spoly_init(temp);
    fmpz_spoly_new_mul_classical(temp, poly1, poly2);
    fmpz_spoly_init2(res, temp->length);
    fmpz_spoly_set(res, temp);
    fmpz_spoly_clear(temp);
  }
  else
  {
    fmpz_spoly_zero(res);
    for (i = 0; i < poly2->length; ++i) 
    {
      fmpz_spoly_t temp;
      fmpz_spoly_init2(temp, poly1->length);
      temp->length = poly1->length;

      _fmpz_vec_scalar_mul_fmpz(temp->coeffs, poly1->coeffs, poly1->length, poly2->coeffs + i);
      fmpz_spoly_shift_left(temp, poly1, poly2->expons + i);
    
      fmpz_spoly_add(res, temp, res);
      
      fmpz_spoly_clear(temp);
    }
  }
  _fmpz_spoly_normalise(res);
}
