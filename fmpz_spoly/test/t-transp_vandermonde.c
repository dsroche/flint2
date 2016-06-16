
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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("transp_vandermonde....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz* vv;
        fmpz* xx;
        fmpz* bb;
        fmpz* bb_check;
        fmpz* xx_check;
        fmpz_t temp;
        slong len, blen, j, k;

        fmpz_init(p);
        fmpz_init(temp);

        len = n_randbits(state, n_randint(state, UWORD(10)));
        if (i % 3 == 0) blen = len;
        else if (i % 3 == 1)
        {
            if (len == 0) len = 1;
            blen = len + n_randint(state, len);
        }
        else
        {
            if (len == 0) len = 1;
            blen = 2 * len + n_randint(state, 2 * len);
        }

        vv = _fmpz_vec_init(len);
        xx = _fmpz_vec_init(len);
        bb = _fmpz_vec_init(blen);
        bb_check = _fmpz_vec_init(blen);
        xx_check = _fmpz_vec_init(len);

        do
        {
            fmpz_randprime(p, state, n_randint(state, UWORD(100)) + 2, 0);
        }
        while (fmpz_cmp_si(p, 2*len) < 0);

        for (j = 0; j < len; ++j)
        {
            int check;
            do 
            {
                fmpz_randm(vv + j, state, p);

                check = 1;
                for (k = 0; k < j; ++k)
                {
                    if (fmpz_equal(vv + k, vv + j)) check = 0;
                }
            }
            while (fmpz_is_zero(vv + j) || ! check);

            fmpz_randtest_unsigned_mod_signed(xx + j, state, p);
        }

        /* quadratic computation of result */
        for (j = 0; j < blen; ++j)
        {
            fmpz_zero(bb + j);
            for (k = 0; k < len; ++k) {
                fmpz_powm_ui(temp, vv + k, j, p);
                fmpz_mul(temp, temp, xx + k);
                fmpz_add(bb + j, bb + j, temp);
            }
            fmpz_mod(bb + j, bb + j, p);
        }

        _fmpz_spoly_transp_vandermonde(bb_check, blen, vv, xx, len, p);

        result = _fmpz_vec_equal(bb_check, bb, blen);

        if (!result)
        {
            flint_printf("FAIL forward:\n");
            flint_printf("len: %wd\n", len);
            flint_printf("p: "); fmpz_print(p); flint_printf("\n");
            _fmpz_vec_print(vv, len); flint_printf("\n\n");
            _fmpz_vec_print(xx, len); flint_printf("\n\n");
            _fmpz_vec_print(bb_check, blen); flint_printf("\n\n");
            _fmpz_vec_print(bb, blen); flint_printf("\n\n");
            abort();
        }

        _fmpz_spoly_transp_vandermonde_inv(xx_check, vv, bb, len, p);

        result = _fmpz_vec_equal(xx_check, xx, len);

        if (!result)
        {
            flint_printf("FAIL inverse:\n");
            flint_printf("len: %wd\n", len);
            flint_printf("p: "); fmpz_print(p); flint_printf("\n");
            _fmpz_vec_print(vv, len); flint_printf("\n\n");
            _fmpz_vec_print(bb, len); flint_printf("\n\n");
            _fmpz_vec_print(xx_check, len); flint_printf("\n\n");
            _fmpz_vec_print(xx, len); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(p);
        fmpz_clear(temp);
        _fmpz_vec_clear(vv, len);
        _fmpz_vec_clear(xx, len);
        _fmpz_vec_clear(bb, blen);
        _fmpz_vec_clear(bb_check, blen);
        _fmpz_vec_clear(xx_check, len);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
