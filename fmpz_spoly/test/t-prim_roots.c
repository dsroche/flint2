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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

ulong my_randtest(flint_rand_t state, ulong maxbits)
{
    return n_randbits(state, n_randint(state, maxbits+1));
}

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("prim_roots....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_t p1, p2;
        fmpz *qq1, *qq2, *ww1, *ww2;
        slong len1, len2, res1, res2;
        ulong pbits, qbits; 
        fmpz_t qprod, a, b, c;
        slong j;

        pbits = my_randtest(state, 6)+2;
        qbits = pbits * (my_randtest(state, 4) + 1);
        len1 = (qbits / pbits) + 1; /* conservative over-estimate */
        len2 = my_randtest(state, 2); /* frequent under-estimate */
        fmpz_init(qprod);
        fmpz_init(a); fmpz_init(b); fmpz_init(c);

        fmpz_init(p1);
        qq1 = _fmpz_vec_init(len1);
        ww1 = _fmpz_vec_init(len1);
        res1 = _fmpz_spoly_prim_roots(p1, qq1, ww1, state, len1, pbits, qbits);

        if (res1 < 0 || res1 > len1)
        {
            flint_printf("FAIL:\n");
            flint_printf("result length %wd\n", res1);
            flint_printf("%wu\n%wu\n%wd\n", pbits, qbits, len1);
            abort();
        }
        if (fmpz_bits(p1) < pbits || ! fmpz_is_prime(p1))
        {
            flint_printf("FAIL:\n");
            flint_printf("bad p\n");
            flint_printf("%wu\n%wu\n%wd\n", pbits, qbits, len1);
            fmpz_print(p1); flint_printf("\n");
            abort();
        }

        fmpz_set_si(qprod, 1);
        for (j=0; j<res1; ++j)
        {
            fmpz_gcd(a, qprod, qq1+j);
            fmpz_fdiv_r(b, qq1+j, p1);
            if (! fmpz_is_prime(qq1+j) || ! fmpz_is_one(a) || ! fmpz_is_one(b))
            {
                flint_printf("FAIL:\n");
                flint_printf("bad q\n");
                flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len1, res1);
                fmpz_print(p1); flint_printf("\n\n");
                fmpz_print(qq1+j); flint_printf("\n\n");
                _fmpz_vec_print(qq1, len1); flint_printf("\n\n");
                abort();
            }
            fmpz_mul(qprod, qprod, qq1+j);

            fmpz_powm(c, ww1+j, p1, qq1+j);
            if (fmpz_is_one(ww1+j) || ! fmpz_is_one(c))
            {
                flint_printf("FAIL:\n");
                flint_printf("bad w\n");
                flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len1, res1);
                fmpz_print(p1); flint_printf("\n\n");
                fmpz_print(qq1+j); flint_printf("\n\n");
                fmpz_print(ww1+j); flint_printf("\n\n");
                abort();
            }
        }

        if (fmpz_bits(qprod) < qbits)
        {
            flint_printf("FAIL:\n");
            flint_printf("q prod too small\n");
            flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len1, res1);
            fmpz_print(p1); flint_printf("\n\n");
            _fmpz_vec_print(qq1, len1); flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(qq1, len1);
        _fmpz_vec_clear(ww1, len1);

        fmpz_init(p2);
        qq2 = _fmpz_vec_init(len2);
        ww2 = _fmpz_vec_init(len2);
        res2 = _fmpz_spoly_prim_roots(p2, qq2, ww2, state, len2, pbits, qbits);

        if (res2 > len2)
        {
            flint_printf("FAIL:\n");
            flint_printf("result length %wd\n", res2);
            flint_printf("%wu\n%wu\n%wd\n", pbits, qbits, len2);
            abort();
        }
        if ((pbits >= 10 && fmpz_equal(p1, p2)) || fmpz_bits(p2) < pbits || ! fmpz_is_prime(p2))
        {
            flint_printf("FAIL:\n");
            flint_printf("bad p2\n");
            flint_printf("%wu\n%wu\n%wd\n", pbits, qbits, len2);
            fmpz_print(p1); flint_printf("\n");
            fmpz_print(p2); flint_printf("\n");
            abort();
        }

        fmpz_set_si(qprod, 1);
        for (j=0; j<(res2<0 ? len2 : res2); ++j)
        {
            fmpz_gcd(a, qprod, qq2+j);
            fmpz_fdiv_r(b, qq2+j, p2);
            if (! fmpz_is_prime(qq2+j) || ! fmpz_is_one(a) || ! fmpz_is_one(b))
            {
                flint_printf("FAIL:\n");
                flint_printf("bad q\n");
                flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len2, res2);
                fmpz_print(p2); flint_printf("\n\n");
                fmpz_print(qq2+j); flint_printf("\n\n");
                _fmpz_vec_print(qq2, len2); flint_printf("\n\n");
                abort();
            }
            fmpz_mul(qprod, qprod, qq2+j);

            fmpz_powm(c, ww2+j, p2, qq2+j);
            if (fmpz_is_one(ww2+j) || ! fmpz_is_one(c))
            {
                flint_printf("FAIL:\n");
                flint_printf("bad w\n");
                flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len2, res2);
                fmpz_print(p2); flint_printf("\n\n");
                fmpz_print(qq2+j); flint_printf("\n\n");
                fmpz_print(ww2+j); flint_printf("\n\n");
                abort();
            }
        }

        if (res2 >= 0 && fmpz_bits(qprod) < qbits)
        {
            flint_printf("FAIL:\n");
            flint_printf("q prod too small\n");
            flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len2, res2);
            fmpz_print(p2); flint_printf("\n\n");
            _fmpz_vec_print(qq2, len2); flint_printf("\n\n");
            abort();
        }
        else if (res2 < 0 && fmpz_bits(qprod) >= qbits)
        {
            flint_printf("FAIL:\n");
            flint_printf("q prod too large\n");
            flint_printf("%wu\n%wu\n%wd\n%wd\n\n", pbits, qbits, len2, res2);
            fmpz_print(p2); flint_printf("\n\n");
            _fmpz_vec_print(qq2, len2); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(p1);
        fmpz_clear(p2);
        _fmpz_vec_clear(qq2, len2);
        _fmpz_vec_clear(ww2, len2);

        fmpz_clear(qprod);
        fmpz_clear(a); fmpz_clear(b); fmpz_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
