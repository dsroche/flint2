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
#include "fmpz_vec.h"

void 
fmpz_spoly_mul_interp(fmpz_spoly_t res, flint_rand_t state, 
        const fmpz_spoly_t poly1, const fmpz_spoly_t poly2, slong terms)
{
    fmpz_spoly_bp_interp_basis_t basis;
    fmpz_spoly_bp_interp_eval_t p1e, p2e;
    mp_bitcnt_t d, h;
    
    d = FLINT_MAX(fmpz_bits(fmpz_spoly_degree_ptr(poly1)), 
                  fmpz_bits(fmpz_spoly_degree_ptr(poly2))) + 1;
    h = fmpz_spoly_height_bits(poly1) + fmpz_spoly_height_bits(poly2) 
        + FLINT_BIT_COUNT(terms);

    fmpz_spoly_bp_interp_basis_init(basis, state, terms, d, h);

    fmpz_spoly_bp_interp_eval_init(p1e, basis);
    fmpz_spoly_bp_interp_eval(p1e, poly1);

    fmpz_spoly_bp_interp_eval_init(p2e, basis);
    fmpz_spoly_bp_interp_eval(p2e, poly2);

    fmpz_spoly_bp_interp_mul(p1e, p1e, p2e);
    fmpz_spoly_bp_interp(res, p1e);

    fmpz_spoly_bp_interp_eval_clear(p1e);
    fmpz_spoly_bp_interp_eval_clear(p2e);
    fmpz_spoly_bp_interp_basis_clear(basis);
}
