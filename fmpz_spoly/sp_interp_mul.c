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

#ifdef _OPENMP
#include "omp.h"
#endif

void fmpz_spoly_sp_interp_mul(fmpz_spoly_sp_interp_eval_t res,
    const fmpz_spoly_sp_interp_eval_t poly1,
    const fmpz_spoly_sp_interp_eval_t poly2)
{
    slong i;

    FLINT_ASSERT(res->basis == poly1->basis);
    FLINT_ASSERT(poly1->basis == poly2->basis);

#ifdef _OPENMP
#pragma omp parallel 
{
#pragma omp for
#endif
    for (i = 0; i < res->basis->length; ++i)
    {
        nmod_poly_mul(res->evals + i, poly1->evals + i, poly2->evals + i);
        nmod_poly_rem_cyc(res->evals + i, res->evals + i, res->basis->emods[i].n);
    }
#ifdef _OPENMP
if (omp_get_thread_num()) flint_cleanup();
}
#endif
}
