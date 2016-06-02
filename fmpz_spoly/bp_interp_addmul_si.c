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

#include "fmpz_spoly.h"

void fmpz_spoly_bp_interp_addmul_si(fmpz_spoly_bp_interp_eval_t res,
    const fmpz_spoly_bp_interp_eval_t poly1,
    slong c,
    const fmpz_spoly_bp_interp_eval_t poly2,
    const fmpz_spoly_bp_interp_basis_t basis)
{
    FLINT_ASSERT(res->length == basis->length);
    FLINT_ASSERT(poly1->length == basis->length);
    FLINT_ASSERT(poly2->length == basis->length);

    if (res->length == 0) return;

    if (res == poly2)
    {
        fmpz_spoly_bp_interp_eval_t temp;
        temp->length = poly2->length;
        temp->evals = _fmpz_vec_init(temp->length);
        _fmpz_vec_set(temp->evals, poly2->evals, temp->length);
        fmpz_spoly_bp_interp_addmul_si(res, poly1, c, temp, basis);
        _fmpz_vec_clear(temp->evals, temp->length);
        return;
    }

    if (res != poly1) _fmpz_vec_set(res->evals, poly1->evals, res->length);

    _fmpz_vec_scalar_addmul_si(res->evals, poly2->evals, res->length, c);
    _fmpz_vec_scalar_mod_fmpz(res->evals, res->evals, res->length, basis->q);
}
