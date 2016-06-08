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

    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "nmod_poly.h"
#include "fmpz_spoly.h"

void fmpz_spoly_sp_interp_addmul_si(fmpz_spoly_sp_interp_eval_t res,
    slong c, const fmpz_spoly_sp_interp_eval_t poly2)
{
    slong i;
    ulong uc;

    FLINT_ASSERT(res->basis == poly2->basis);

    for (i = 0; i < res->basis->length; ++i)
    {
        nmod_poly_fit_length(res->evals + i, nmod_poly_length(poly2->evals + i));
        /* XXX: does this work for any slong value c? */
        NMOD_RED(uc, c + res->basis->cmods[i].n, res->basis->cmods[i]);
        _nmod_vec_scalar_addmul_nmod(res->evals[i].coeffs, 
                poly2->evals[i].coeffs, nmod_poly_length(poly2->evals + i),
                uc, res->basis->cmods[i]);
        _nmod_poly_normalise(res->evals + i);
    }
}
