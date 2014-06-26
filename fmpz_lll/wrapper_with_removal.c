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

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

int
fmpz_lll_wrapper_with_removal(fmpz_mat_t B, const fmpz_t gs_B,
                              const fmpz_lll_t fl)
{
    int res = fmpz_lll_d_with_removal(B, gs_B, fl);

    if (res == -1)
    {
        if (fl->gt == APPROX)
        {
            res = fmpz_lll_d_heuristic_with_removal(B, gs_B, fl);
            if (res == -1)
            {
                res = fmpz_lll_mpf_with_removal(B, gs_B, fl);
            }
        }
        else
        {
            res = fmpz_lll_mpf_with_removal(B, gs_B, fl);
        }
    }

    return res;
}
