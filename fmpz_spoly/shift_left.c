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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"

void fmpz_spoly_shift_left(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t n)
{
    int i;
    slong newlen = poly->length;
    
    if (res != poly)
    {
        _fmpz_spoly_reserve(res, newlen);
        _fmpz_vec_set(res->coeffs, poly->coeffs, newlen);
        _fmpz_spoly_set_length(res, newlen);
    }

    for (i = 0; i < newlen; ++i)
    {
        fmpz_add(res->expons+i, poly->expons+i, n);
    }
}
