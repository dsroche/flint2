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

#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

void nmod_poly_rem_cyc(nmod_poly_t res, const nmod_poly_t poly, ulong n)
{
    slong i;
    ulong rese=0;
    ulong q = nmod_poly_modulus(poly);

    FLINT_ASSERT(nmod_poly_modulus(res) == q);

    if (n == 0)
    {
        flint_printf("Exception (nmod_poly_rem_cyc). Division by zero.\n");
        abort();
    }
    
    if (res != poly) { /* set low-order coeffs of res to those of poly*/
        nmod_poly_truncate(res, n);
        for (i=n-1; i>=0; --i) {
            nmod_poly_set_coeff_ui(res, i, nmod_poly_get_coeff_ui(poly, i));
        }
    }
    
    for (i=n; i<nmod_poly_length(poly); ++i) {
        nmod_poly_set_coeff_ui(res, rese, n_addmod(nmod_poly_get_coeff_ui(res, rese), 
                    nmod_poly_get_coeff_ui(poly, i), q));
        if (++rese == n) 
            rese = 0;
    }
    
    if (res == poly) { /* aliasing; time to truncate the high-order coeffs*/
        nmod_poly_truncate(res, n);
    }
}
