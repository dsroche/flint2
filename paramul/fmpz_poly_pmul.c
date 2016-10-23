/*
    Authored 2016 by Daniel S. Roche; U.S. Government work product in the public domain

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"

void _fmpz_poly_pmul(fmpz* res,
        const fmpz* p1, slong len1, const fmpz* p2, slong len2, slong cbits)
{
    slong primebits = FLINT_BITS - 2;
    slong numprimes = 1 + cbits / (primebits - 1);
    mp_ptr primes;
    mp_ptr p1reds, p2reds, resreds;
    slong reslen = len1 + len2 - 1;
    slong i;

    primes = flint_malloc(numprimes * sizeof *primes);
    primes[0] = n_nextprime(UWORD(1) << (primebits - 1), 0);
    for (i = 1; i < numprimes; ++i)
    {
        primes[i] = n_nextprime(primes[i - 1], 0);
    }

    p1reds = flint_malloc(numprimes * len1 * sizeof *p1reds);
    p2reds = flint_malloc(numprimes * len2 * sizeof *p2reds);
    resreds = flint_malloc(numprimes * reslen * sizeof *resreds);

#pragma omp parallel private(i)
    {
        fmpz_comb_t comb;
        fmpz_comb_temp_t tcomb;
        mp_ptr imgs;
	nmod_t mod;
        slong j, k;

        fmpz_comb_init(comb, primes, numprimes);
        fmpz_comb_temp_init(tcomb, comb);
        imgs = flint_malloc(numprimes * sizeof *imgs);

#pragma omp for nowait
        for (i = 0; i < len1; ++i)
        {
            fmpz_multi_mod_ui(imgs, p1 + i, comb, tcomb);
	    for (j = 0, k = i; j < numprimes; ++j, k += len1)
		p1reds[k] = imgs[j];
        }

#pragma omp for
        for (i = 0; i < len2; ++i)
        {
            fmpz_multi_mod_ui(imgs, p2 + i, comb, tcomb);
	    for (j = 0, k = i; j < numprimes; ++j, k += len2)
		p2reds[k] = imgs[j];
        }
        
#pragma omp for
	for (i = 0; i < numprimes; ++i)
        {
	    nmod_init(&mod, primes[i]);
	    _nmod_poly_mul(resreds + i * reslen,
		p1reds + i * len1, len1, 
		p2reds + i * len2, len2, 
		mod);
        }

#pragma omp for nowait
        for (i = 0; i < len1 + len2 - 1; ++i)
        {
	    for (j = 0, k = i; j < numprimes; ++j, k += reslen)
		imgs[j] = resreds[k];
            fmpz_multi_CRT_ui(res + i, imgs, comb, tcomb, 1);
        }

        flint_free(imgs);
        fmpz_comb_temp_clear(tcomb);
        fmpz_comb_clear(comb);
    }

    flint_free(p1reds);
    flint_free(p2reds);
    flint_free(resreds);
    flint_free(primes);
}

void fmpz_poly_pmul(fmpz_poly_t res,
        const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong rlen;
    slong cbits;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    rlen = len1 + len2 - 1;
    fmpz_poly_fit_length(res, rlen);


    if (len1 >= len2)
    {
        cbits = FLINT_ABS(fmpz_poly_max_bits(poly1)) 
            + FLINT_ABS(fmpz_poly_max_bits(poly2))
            + FLINT_BIT_COUNT(len2);
        _fmpz_poly_pmul(res->coeffs, poly1->coeffs, len1,
                       poly2->coeffs, len2, cbits);
    }
    else
    {
        cbits = FLINT_ABS(fmpz_poly_max_bits(poly1)) 
            + FLINT_ABS(fmpz_poly_max_bits(poly2))
            + FLINT_BIT_COUNT(len1);
        _fmpz_poly_pmul(res->coeffs, poly2->coeffs, len2,
                       poly1->coeffs, len1, cbits);
    }

    _fmpz_poly_set_length(res, rlen);
}
