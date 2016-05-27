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
#include "fmpz_poly.h"

void fmpz_spoly_get_fmpz_poly(fmpz_poly_t out, const fmpz_spoly_t in)
{
    if(fmpz_spoly_is_zero(in))
    {
        fmpz_poly_zero(out);
        return;
    }
    else if(fmpz_sgn(in->expons) == -1)
    {
        fmpz_poly_zero(out);
        return;
    }
    else
    {
        slong i;
        slong deg = fmpz_spoly_degree_si(in);
        
        fmpz_poly_fit_length(out, deg + 1); 
        
        for (i = 0; i < fmpz_spoly_terms(in); ++i)
        {
            fmpz_set(out->coeffs + fmpz_get_si(in->expons + i), in->coeffs + i);
        }

        _fmpz_poly_set_length(out, deg + 1);
    }
}
