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
    Authored 2016 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_spoly.h"

void
fmpz_spoly_scalar_addmul(fmpz_spoly_t poly1, const fmpz_spoly_t poly2,
                          const fmpz_t x)
{
    if (fmpz_is_zero(x) || fmpz_spoly_is_zero(poly2))
    {
    }
    else if (fmpz_spoly_is_zero(poly1))
    {
        fmpz_spoly_scalar_mul(poly1, poly2, x);
    }
    else if (fmpz_cmp(fmpz_spoly_degree_ptr(poly2),
                      fmpz_spoly_lowdeg_ptr(poly1)) < 0)
    {
        /* no overlap case #1 */
        slong t1 = fmpz_spoly_terms(poly1);
        slong t2 = fmpz_spoly_terms(poly2);

        _fmpz_spoly_reserve(poly1, t1 + t2);

        /* copy poly2*x to the end of poly1 */
        _fmpz_vec_set(poly1->expons + t1, poly2->expons, t2);
        _fmpz_vec_scalar_mul_fmpz(poly1->coeffs + t1, poly2->coeffs, t2, x);

        _fmpz_spoly_set_length(poly1, t1 + t2);
    }
    else if (fmpz_cmp(fmpz_spoly_degree_ptr(poly1), 
                      fmpz_spoly_lowdeg_ptr(poly2)) < 0)
    {
        /* no overlap case #2 */
        slong t1 = fmpz_spoly_terms(poly1);
        slong t2 = fmpz_spoly_terms(poly2);
        slong i;

        _fmpz_spoly_reserve(poly1, t1 + t2);

        /* move everything in poly1 down t2 spots */
        for (i = t1 - 1; i >= 0; --i)
        {
            fmpz_swap(poly1->coeffs + (i + t2), poly1->coeffs + i);
            fmpz_swap(poly1->expons + (i + t2), poly1->expons + i);
        }

        /* copy poly2*x to the beginning of poly1 */
        _fmpz_vec_set(poly1->expons, poly2->expons, t2);
        _fmpz_vec_scalar_mul_fmpz(poly1->coeffs, poly2->coeffs, t2, x);

        _fmpz_spoly_set_length(poly1, t1 + t2);
    }
    else
    {
        /* general case */
        fmpz_spoly_t res;
        slong p1 = 0, p2 = 0, rp = 0;
        slong t1 = fmpz_spoly_terms(poly1);
        slong t2 = fmpz_spoly_terms(poly2);

        fmpz_spoly_init2(res, t1 + t2);

        while (p1 < t1 && p2 < t2)
        {
            if (fmpz_cmp(poly1->expons + p1, poly2->expons + p2) > 0)
            {
                fmpz_swap(poly1->expons + p1, res->expons + rp);
                fmpz_swap(poly1->coeffs + p1, res->coeffs + rp);
                ++p1;
                ++rp;
            }
            else if (fmpz_cmp(poly2->expons + p2, poly1->expons + p1) > 0)
            {
                fmpz_set(res->expons + rp, poly2->expons + p2);
                fmpz_mul(res->coeffs + rp, poly2->coeffs + p2, x);
                ++p2;
                ++rp;
            }
            else
            {
                fmpz_swap(poly1->expons + p1, res->expons + rp);
                fmpz_mul(res->coeffs + rp, poly2->coeffs + p2, x);
                fmpz_add(res->coeffs + rp, res->coeffs + rp, poly1->coeffs + p1);
                ++p1;
                ++p2;
                if (! fmpz_is_zero(res->coeffs + rp)) ++rp;
            }
        }

        if (p1 < t1)
        {
            _fmpz_vec_swap(res->expons + rp, poly1->expons + p1, t1 - p1);
            _fmpz_vec_swap(res->coeffs + rp, poly1->coeffs + p1, t1 - p1);
            rp += t1 - p1;
        }
        else if (p2 < t2)
        {
            _fmpz_vec_set(res->expons + rp, poly2->expons + p2, t2 - p2);
            _fmpz_vec_scalar_mul_fmpz(res->coeffs + rp, 
                    poly2->coeffs + p2, t2 - p2, x);
            rp += t2 - p2;
        }

        _fmpz_spoly_set_length(res, rp);
        fmpz_spoly_swap(poly1, res);
        fmpz_spoly_clear(res);
    }
}
