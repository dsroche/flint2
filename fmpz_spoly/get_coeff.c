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

slong fmpz_spoly_get_coeff_si_si(const fmpz_spoly_t poly, slong e)
{
    slong ind = _fmpz_spoly_index_si(poly, e);
    if (ind < 0) return WORD(0);
    else return fmpz_get_si(poly->coeffs + ind);
}

slong fmpz_spoly_get_coeff_si_fmpz(const fmpz_spoly_t poly, const fmpz_t e)
{
    slong ind = _fmpz_spoly_index(poly, e);
    if (ind < 0) return WORD(0);
    else return fmpz_get_si(poly->coeffs + ind);
}

void fmpz_spoly_get_coeff_fmpz_si(fmpz_t res, const fmpz_spoly_t poly, slong e)
{
    slong ind = _fmpz_spoly_index_si(poly, e);
    if (ind < 0) fmpz_zero(res);
    else fmpz_set(res, poly->coeffs + ind);
}

void fmpz_spoly_get_coeff(fmpz_t res, const fmpz_spoly_t poly, const fmpz_t e)
{
    slong ind = _fmpz_spoly_index(poly, e);
    if (ind < 0) fmpz_zero(res);
    else fmpz_set(res, poly->coeffs + ind);
}

