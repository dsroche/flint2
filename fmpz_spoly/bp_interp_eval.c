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

void fmpz_spoly_bp_interp_eval(fmpz_spoly_bp_interp_t res,
    const fmpz_spoly_t poly)
{
    fmpz_spoly_evaluate_powers(res->evaluations, res->length, 
            poly, res->sample_points + 1, res->q);
}

void fmpz_spoly_bp_interp_mul(fmpz_spoly_bp_interp_t res,
        const fmpz_spoly_t poly)
{
    if (res->length < EVAL_XOVER)
    {
        slong i;
        fmpz_t temp;
        
        fmpz_init2(temp, fmpz_size(res->q));

        for (i = 0; i < res->length; ++i)
        {
            fmpz_spoly_evaluate_mod(temp, poly, res->sample_points + i, res->q);
            fmpz_mul(res->evaluations + i, res->evaluations + i, temp);
            fmpz_mod(res->evaluations + i, res->evaluations + i, res->q);
        }

        fmpz_clear(temp);
    }
    else
    {
        slong i;
        fmpz* temp = _fmpz_vec_init(res->length);

        fmpz_spoly_evaluate_powers(temp, res->length,
                poly, res->sample_points + 1, res->q);

        for (i = 0; i < res->length; ++i)
        {
            fmpz_mul(res->evaluations + i, res->evaluations + i, temp + i);
            fmpz_mod(res->evaluations + i, res->evaluations + i, res->q);
        }

        _fmpz_vec_clear(temp, res->length);
    }
}


void fmpz_spoly_bp_interp_add(fmpz_spoly_bp_interp_t res,
        const fmpz_t c, const fmpz_spoly_t poly)
{
    if (res->length < EVAL_XOVER)
    {
        slong i;
        fmpz_t temp;
        
        fmpz_init2(temp, fmpz_size(res->q));

        for (i = 0; i<res->length; ++i)
        {
            fmpz_spoly_evaluate_mod(temp, poly, res->sample_points + i, res->q);
            fmpz_addmul(res->evaluations + i, c, temp);
            fmpz_mod(res->evaluations + i, res->evaluations + i, res->q);
        }

        fmpz_clear(temp);
    }
    else
    {
        fmpz* temp = _fmpz_vec_init(res->length);

        fmpz_spoly_evaluate_powers(temp, res->length,
                poly, res->sample_points + 1, res->q);

        _fmpz_vec_add(res->evaluations, res->evaluations, temp, res->length);
        _fmpz_vec_scalar_mod_fmpz(res->evaluations, res->evaluations, res->length, res->q);

        _fmpz_vec_clear(temp, res->length);
    }
}

void fmpz_spoly_bp_interp_pow(fmpz_spoly_bp_interp_t res, ulong pow)
{
    slong i;
    
    for (i=0; i<res->length; ++i)
    {
        fmpz_powm_ui(res->evaluations+i, res->evaluations+i, pow, res->q);
    }
}

