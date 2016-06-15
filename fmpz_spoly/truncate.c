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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_spoly.h"

void
fmpz_spoly_truncate(fmpz_spoly_t poly, const fmpz_t deg)
{
    slong index = 0, len;

    index = _fmpz_spoly_index(poly, deg);
    if(index < 0)
      index = index*-1 - 1;

    len = poly->length - index;
    if(len <= 0 || index == poly->length)
    {
      fmpz_spoly_zero(poly);
      return;
    }

    _fmpz_spoly_vec_shift(poly, index, poly->length, -1*index);

    poly->length = len;
}
