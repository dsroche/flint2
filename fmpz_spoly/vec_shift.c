/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and / or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
     
    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110 - 1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_spoly.h"
#include <string.h>

void _fmpz_spoly_vec_shift_arr(fmpz* coeffs, fmpz* expons, 
    slong start, slong end, slong dist)
{
    memmove(coeffs + start + dist, coeffs + start, 
            (end - start) * sizeof *coeffs);
    memmove(expons + start + dist, expons + start, 
            (end - start) * sizeof *expons);
    
    if (dist > 0) 
    {
        memset(coeffs + start, 0, dist * sizeof *coeffs);
        memset(expons + start, 0, dist * sizeof *expons);
    }
    else if (dist < 0) 
    {
        memset(coeffs + end + dist, 0, (-dist) * sizeof *coeffs);
        memset(expons + end + dist, 0, (-dist) * sizeof *expons);
    }
}
