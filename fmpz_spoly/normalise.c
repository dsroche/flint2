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

const slong FMPZ_SPOLY_QSORT_XOVER = 128;

#define FMPZ_SPOLY_NORM_SWAP(coeffs, expons, pos1, pos2) do { \
    if ((pos1) != (pos2)) \
    { \
        fmpz t = (coeffs)[pos1]; \
        (coeffs)[pos1] = (coeffs)[pos2]; \
        (coeffs)[pos2] = t; \
        t = (expons)[pos1]; \
        (expons)[pos1] = (expons)[pos2]; \
        (expons)[pos2] = t; \
    } \
} while(0)

slong _fmpz_spoly_normalise_rec(fmpz* coeffs, fmpz* expons, slong len)
{
    if (len <= FMPZ_SPOLY_QSORT_XOVER)
    {
        /* A modified insertion sort collects terms with the same exponent. */

        /* nfree will hold the number of emptied terms not yet removed */
        slong cur = 1, free1 = 0, free2 = 0, i; 
        fmpz tcoeff, texpon;

        /* 1st pass: sort and combine like terms */
        for (; cur < len; ++cur) 
        {
            i = cur - free1 - 1;

            /* find the position where cur belongs */
            for (; i >= 0 && fmpz_cmp(expons + i, expons + cur) < 0; --i);

            if (i >= 0 && fmpz_equal(expons + i, expons + cur))
            {
                /* collision */
                fmpz_add(coeffs + i, coeffs + i, coeffs + cur);
                fmpz_zero(coeffs + cur);
                fmpz_zero(expons + cur);
                ++free1;
            }
            else if (i + 1 < cur) 
            {
                if (free1 == 0)
                {
                    tcoeff = coeffs[cur];
                    texpon = expons[cur];
                    _fmpz_spoly_vec_shift_arr(coeffs, expons, i + 1, cur, 1);
                    coeffs[i + 1] = tcoeff;
                    expons[i + 1] = texpon;
                }
                else 
                {
                    _fmpz_spoly_vec_shift_arr(coeffs, expons, i + 1, cur - free1, 1);
                    FMPZ_SPOLY_NORM_SWAP(coeffs, expons, i + 1, cur);
                }
            }
        }

        /* 2nd pass: remove zero terms */
        for (cur = 0; cur < len - free1 && ! fmpz_is_zero(coeffs + cur); 
             ++cur);

        while (cur < len - free1) 
        {
            i = cur;
            for (; cur < len - free1 && fmpz_is_zero(coeffs + cur); ++cur);
            free2 += (cur - i);
            i = cur;
            for (; cur < len - free1 && ! fmpz_is_zero(coeffs + cur); ++cur);
            _fmpz_spoly_vec_shift_arr(coeffs, expons, i, cur, -free2);
        }

        return len - free1 - free2;
    }
    else
    {
        /* quicksort that also combines terms equal to the pivot */

        /* initially a = front, b = mid, c = rear */
        slong a = 0, b = len / 2, c = len - 1;
        int frontCpiv = fmpz_cmp(expons + a, expons + b);
        int rearCpiv = fmpz_cmp(expons + c, expons + b);

        if ((frontCpiv ^ rearCpiv) < 0)
        {
            /* make midpoint the pivot */
            FMPZ_SPOLY_NORM_SWAP(coeffs, expons, b, c);
        }
        else
        {
            int frontCrear = fmpz_cmp(expons + a, expons + c);
            if ((frontCpiv ^ frontCrear) < 0)
            {
                /* make front the pivot */
                FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, c);
                rearCpiv = -frontCpiv;
                frontCpiv = -frontCrear;
            }
            else
            {
                /* make rear the pivot */
                frontCpiv = frontCrear;
                rearCpiv = -rearCpiv;
            }
        }

        FMPZ_SPOLY_NORM_SWAP(coeffs, expons, b, c - 1);
        b = c - 1;
        /* now c is the pivot, a is the front, b is the rear,
         * and frontCpiv and rearCpiv are correct */

        while (1)
        {
            FLINT_ASSERT(a < b);

            while (frontCpiv >= 0)
            {
                if (frontCpiv == 0)
                {
                    /* add front coeff to pivot and move bubble to end */
                    fmpz_add(coeffs + c, coeffs + c, coeffs + a);
                    FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, b);
                    frontCpiv = rearCpiv;
                    --len;
                    FMPZ_SPOLY_NORM_SWAP(coeffs, expons, b, len - 1);
                    if (--b == a)
                        goto abequal;
                    rearCpiv = fmpz_cmp(expons + b, expons + c);
                }
                else
                {
                    if (++a == b)
                    {
                        frontCpiv = rearCpiv;
                        goto abequal;
                    }
                    frontCpiv = fmpz_cmp(expons + a, expons + c);
                }
            }

            while (rearCpiv <= 0)
            {
                if (rearCpiv == 0)
                {
                    /* add rear coeff to pivot and move bubble to end */
                    fmpz_add(coeffs + c, coeffs + c, coeffs + b);
                    --len;
                    FMPZ_SPOLY_NORM_SWAP(coeffs, expons, b, len - 1);
                }
                if (--b == a)
                    goto abequal;
                rearCpiv = fmpz_cmp(expons + b, expons + c);
            }

            /* swap front and rear */
            FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, b);
            if (++a == --b)
                goto abequal;
            else if (a > b)
                goto finish;

            frontCpiv = fmpz_cmp(expons + a, expons + c);
            rearCpiv = fmpz_cmp(expons + b, expons + c);
        }

abequal:
        if (frontCpiv > 0)
            ++a;
        else if (frontCpiv == 0)
        {
            fmpz_add(coeffs + c, coeffs + c, coeffs + a);
            --len;
            FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, len - 1);
        }

finish:
        /* a is the index where the pivot belongs. */
        if (fmpz_is_zero(coeffs + c))
        {
            --len;
            b = a;
        }
        else
        {
            FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, len - 1);
            FMPZ_SPOLY_NORM_SWAP(coeffs, expons, a, c);
            b = a + 1;
        }

        /* recursive calls */
        c = _fmpz_spoly_normalise_rec(coeffs, expons, a);
        if (c < a)
        {
            /* size shrank in first recursive call */
            _fmpz_spoly_vec_shift_arr(coeffs, expons, b, len, c - a);
            len -= a - c;
            b -= a - c;
        }

        return c + _fmpz_spoly_normalise_rec(coeffs + b, expons + b, len - b);
    }
}


void _fmpz_spoly_normalise(fmpz_spoly_t poly)
{
    _fmpz_spoly_set_length(poly,
        _fmpz_spoly_normalise_rec(poly->coeffs, poly->expons, poly->length));
}
