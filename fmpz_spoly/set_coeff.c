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

#define FMPZ_SPOLY_SET_COEFF_MACRO(COEFF_IS_ZERO, SET_COEFF, SET_EXPON) \
    if (ind < 0) { \
        if (! (COEFF_IS_ZERO)) { \
            /* nonzero coeff, not in the polynomial; we have to insert it. */ \
            _fmpz_spoly_reserve(poly, poly->length+1); \
            ind = WORD(-1) - ind; \
            _fmpz_spoly_vec_shift(poly, ind, poly->length, 1); \
            SET_COEFF(poly->coeffs+ind, c); \
            SET_EXPON(poly->expons+ind, e); \
            ++ poly->length; \
        } \
    } \
    else if (COEFF_IS_ZERO) { \
        /* zero coeff, in the polynomial; we have to remove it. */ \
        fmpz_clear(poly->coeffs+ind); \
        fmpz_clear(poly->expons+ind); \
        _fmpz_spoly_vec_shift(poly, ind+1, poly->length, WORD(-1)); \
        poly->length -= 1; \
    } \
    else { \
        /* nonzero coeff, in the polynomial; just change it. */ \
        SET_COEFF(poly->coeffs + ind, c); \
    }

void fmpz_spoly_set_coeff_si_si(fmpz_spoly_t poly, slong c, slong e)
{
    slong ind = _fmpz_spoly_index_si(poly, e);
    FMPZ_SPOLY_SET_COEFF_MACRO((c == 0), fmpz_set_si, fmpz_set_si)
}

void fmpz_spoly_set_coeff_si_fmpz(fmpz_spoly_t poly, slong c, const fmpz_t e)
{
    slong ind = _fmpz_spoly_index(poly, e);
    FMPZ_SPOLY_SET_COEFF_MACRO((c == 0), fmpz_set_si, fmpz_set)
}

void fmpz_spoly_set_coeff_fmpz_si(fmpz_spoly_t poly, const fmpz_t c, slong e)
{
    slong ind = _fmpz_spoly_index_si(poly, e);
    FMPZ_SPOLY_SET_COEFF_MACRO(fmpz_is_zero(c), fmpz_set, fmpz_set_si)
}

void fmpz_spoly_set_coeff(fmpz_spoly_t poly, const fmpz_t c, const fmpz_t e)
{
    slong ind = _fmpz_spoly_index(poly, e);
    FMPZ_SPOLY_SET_COEFF_MACRO(fmpz_is_zero(c), fmpz_set, fmpz_set)
}

#undef FMPZ_SPOLY_SET_COEFF_MACRO
