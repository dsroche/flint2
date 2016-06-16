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

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz_spoly.h"

void fmpz_spoly_sp_interp_eval(fmpz_spoly_sp_interp_eval_t res,
        const fmpz_spoly_t poly)
{
    slong i = 0, j, poly_len;
    ulong round_cmod = 0, round_shift = 0, group_emod;
    slong group_len;
    ulong *round_coeffs, *round_expons;
    const nmod_t *cmods = res->basis->cmods;
    const ulong *shifts = res->basis->shifts;
    const nmod_t *emods = res->basis->emods;
    mp_ptr eval_coeffs;

    poly_len = fmpz_spoly_terms(poly);
    round_coeffs = flint_malloc(2 * poly_len * sizeof *round_coeffs);
    round_expons = round_coeffs + poly_len;
    
    while (i < res->basis->length)
    {
        if (cmods[i].n != round_cmod || shifts[i] != round_shift)
        {
            /* start new round */
            round_cmod = cmods[i].n;
            round_shift = shifts[i];

            for (j = 0; j < poly_len; ++j)
            {
                round_coeffs[j] = n_mulmod2_preinv(
                    n_powmod2_ui_preinv(shifts[i], 
                        fmpz_fdiv_ui(poly->expons + j, cmods[i].n - 1),
                        cmods[i].n, cmods[i].ninv),
                    fmpz_fdiv_ui(poly->coeffs + j, cmods[i].n),
                    cmods[i].n, cmods[i].ninv);
            }
        }

        /* start new group */
        group_emod = emods[i].n;
        group_len = 0;

        for (j = 0; j < poly_len; ++j)
        {
            round_expons[j] = fmpz_fdiv_ui(poly->expons + j, group_emod);
            if (round_expons[j] >= (ulong) group_len) 
            {
                group_len = (slong) round_expons[j] + 1;
            }
        }

        /* evaluate group leader */
        nmod_poly_fit_length(res->evals + i, group_len);
        eval_coeffs = res->evals[i].coeffs;
        memset(eval_coeffs, 0, group_len * sizeof *eval_coeffs);

        for (j = 0; j < poly_len; ++j)
        {
            ulong rese = round_expons[j];
            if (eval_coeffs[rese])
            {
                eval_coeffs[rese] = n_addmod(eval_coeffs[rese], 
                        round_coeffs[j], round_cmod);
            }
            else eval_coeffs[rese] = round_coeffs[j];
        }

        _nmod_poly_set_length(res->evals + i, group_len);
        _nmod_poly_normalise(res->evals + i);

        /* evaluate rest of the group */
        while (++i < res->basis->length 
               && shifts[i] == UWORD(1) && emods[i].n == group_emod)
        {
            nmod_poly_fit_length(res->evals + i, group_len);
            eval_coeffs = res->evals[i].coeffs;
            memset(eval_coeffs, 0, group_len * sizeof *eval_coeffs);

            for (j = 0; j < poly_len; ++j)
            {
                ulong rese = round_expons[j];
                if (eval_coeffs[rese])
                {
                    eval_coeffs[rese] = n_addmod(eval_coeffs[rese], 
                        fmpz_fdiv_ui(poly->coeffs + j, cmods[i].n), cmods[i].n);
                }
                else 
                {
                    eval_coeffs[rese] = 
                        fmpz_fdiv_ui(poly->coeffs + j, cmods[i].n);
                }
            }

            _nmod_poly_set_length(res->evals + i, group_len);
            _nmod_poly_normalise(res->evals + i);
        }
    }

    flint_free(round_coeffs);
}
