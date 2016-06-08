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

static const int EVAL_XOVER = 32;

void fmpz_spoly_evaluate_powers(fmpz* res, slong len,
    const fmpz_spoly_t poly, const fmpz_t w, const fmpz_t p)
{
    if (len == 0)
    {
    }
    else if (fmpz_spoly_terms(poly) < EVAL_XOVER)
    {
        slong i;
        fmpz_t wpow;

        fmpz_init_set_ui(wpow, UWORD(1));
        for (i = 0; ; ++i)
        {
            fmpz_spoly_evaluate_mod(res + i, poly, wpow, p);
            if (i + 1 == len) break;
            fmpz_mul(wpow, wpow, w);
            fmpz_mod(wpow, wpow, p);
        }

        fmpz_clear(wpow);
    }
    else if (len < poly->length)
    {
        fmpz* temp = _fmpz_vec_init(len);
        fmpz* wpows = _fmpz_vec_init(len);
        slong i, j;

        for (i = 0; i + len < poly->length; i += len)
        {
            for (j = 0; j < len; ++j)
            {
                fmpz_powm(wpows + j, w, poly->expons + (i + j), p);
            }
            _fmpz_spoly_transp_vandermonde(temp, len, 
                    wpows, poly->coeffs + i, len, p);
            _fmpz_vec_add(res, res, temp, len);
        }

        for (j = 0; i + j < poly->length; ++j)
        {
            fmpz_powm(wpows + j, w, poly->expons + (i + j), p);
        }
        _fmpz_spoly_transp_vandermonde(temp, len,
                wpows, poly->coeffs + i, j, p);
        _fmpz_vec_add(res, res, temp, len);
        _fmpz_vec_scalar_mod_fmpz(res, res, len, p);
        
        _fmpz_vec_clear(temp, len);
        _fmpz_vec_clear(wpows, len);
    }
    else
    {
        fmpz* wpows = _fmpz_vec_init(poly->length);
        slong i;

        for (i = 0; i < poly->length; ++i)
        {
            fmpz_powm(wpows + i, w, poly->expons + i, p);
        }

        _fmpz_spoly_transp_vandermonde(res, len,
                wpows, poly->coeffs, poly->length, p);

        _fmpz_vec_clear(wpows, poly->length);
    }
}
