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

#include "flint.h"
#include "nmod_poly.h"
#include "fmpz_spoly.h"

void fmpz_spoly_sp_interp_eval_init(fmpz_spoly_sp_interp_eval_t res,
        const fmpz_spoly_sp_interp_basis_t basis)
{
    slong i;

    res->basis = basis;
    
    if (basis->length == 0) return;

    res->evals = flint_malloc(basis->length * sizeof *res->evals);
    
    for (i = 0; i < basis->length; ++i)
    {
        nmod_poly_init2_preinv(res->evals + i, 
                basis->cmods[i].n, basis->cmods[i].ninv, basis->emods[i].n);
    }
}
