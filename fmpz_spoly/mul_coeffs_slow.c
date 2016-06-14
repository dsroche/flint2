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

#include "fmpz_spoly.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

static const double LN_2 = 0.693147180559945309417232121458;

slong _fmpz_spoly_prim_roots(fmpz_t p, fmpz * qq, fmpz * ww, flint_rand_t state,
    slong len, mp_bitcnt_t p_bits, mp_bitcnt_t q_prod_bits)
{
    /* This function computes single random prime p, with the proper number
     * of bits, and corresponding primes qq, as well as an array
     * ww of elements in GF(q).
     * For each index i, it holds that p divides (qq[i]-1), and
     * that ww[i] is a p-th primitive root mod qq[i].
     * Furthermore, the product of primes in qq has at least q_prod_bits bits.
     * All arrays must be allocated by the caller, and the fmpz values
     * initialized beforehand. 
     * The len parameter is the length of qq and ww.
     * The returned value gives the actual number of primes q computed, 
     * which will be at most len.
     * In case len is too small to compute sufficiently many primes, -1 is
     * returned.
     */
    slong count = 0;
    fmpz_t qprod;
    ulong a = 0;
    mp_bitcnt_t qpbits;

    fmpz_randprime(p, state, p_bits, 0);

    if (len == 0) 
    {
        return (q_prod_bits > 0) ? -1 : 0;
    }

    fmpz_init_set_si(qprod, 1);
    fmpz_set_si(qq+0, 1);

    while (count < len && fmpz_bits(qprod) < q_prod_bits)
    {
        /* qq[count] = a*p + 1 and a is even */
        fmpz_addmul_ui(qq+count, p, 2); 
        a += 2;

        if (fmpz_is_prime(qq+count))
        {
            do
            {
                fmpz_randm(ww+count, state, qq+count);
                fmpz_powm_ui(ww+count, ww+count, a, qq+count);
            } while (fmpz_cmp_si(ww+count, 1) <= 0);
            fmpz_mul(qprod, qprod, qq+count);
            if (++count < len) {
                fmpz_set(qq+count, qq+(count-1));
            }
        }
    }

    qpbits = fmpz_bits(qprod);
    fmpz_clear(qprod);

    if (qpbits >= q_prod_bits) return count;
    else return -1;
}

void _fmpz_spoly_mul_coeffs_slow(fmpz_spoly_t res, flint_rand_t state,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2, const fmpz * expons,
    slong len)
{
  fmpz * qq, * ww, * vv, * coeffs, * mod_expons, * eval1, * eval2, * coeffs_mod_q;
  fmpz_t p, C, H, D, temp, q_total;
  slong p_bits, q_prod_bits, n, num_primes, i, j, k;

  fmpz_init(p);
  fmpz_init(C);
  fmpz_init(H);
  fmpz_init(D);
  fmpz_init(temp);
  fmpz_init(q_total);

  fmpz_spoly_degree(D, poly1);
  fmpz_spoly_degree(temp, poly2);

  if(fmpz_cmp(temp, D) > 0)
    fmpz_set(D, temp);

  fmpz_spoly_height(H, poly1);
  fmpz_spoly_height(temp, poly2);
  fmpz_mul(C, temp, H);

  p_bits = 2*FLINT_CLOG2(len) + FLINT_CLOG2(fmpz_bits(D)) + 1;

  q_prod_bits = fmpz_bits(C) + FLINT_CLOG2(FLINT_MAX(poly1->length, poly2->length)) + 1;

  n = (q_prod_bits + p_bits - 1)/p_bits + 1;

  qq = _fmpz_vec_init(n);
  ww = _fmpz_vec_init(n);

  vv = _fmpz_vec_init(len);
  coeffs = _fmpz_vec_init(len);
  coeffs_mod_q = _fmpz_vec_init(len);
  mod_expons = _fmpz_vec_init(len);
  eval1 = _fmpz_vec_init(len);
  eval2 = _fmpz_vec_init(len);

  num_primes = _fmpz_spoly_prim_roots(p, qq, ww, state, n, p_bits, q_prod_bits);
  
  _fmpz_vec_scalar_mod_fmpz(mod_expons, expons, len, p);
  fmpz_one(q_total);
  
  for(i = 0; i < num_primes; i++)
  {
    for(j = 0; j < len; j++)
    {
      fmpz_powm(vv + j, ww + i, mod_expons + j, qq + i);
    }

    fmpz_spoly_evaluate_powers(eval1, len, poly1, ww + i, qq + i);
    fmpz_spoly_evaluate_powers(eval2, len, poly2, ww + i, qq + i);

    for(j = 0; j < len; j++)
    {
      fmpz_mul(eval1 + j, eval1 + j, eval2 + j);
      fmpz_mod(eval1 + j, eval1 + j, qq + i);
    }

    _fmpz_spoly_transp_vandermonde_inv(coeffs_mod_q, vv, eval1, len, qq + i);

    for(k = 0; k < len; k++)
    {
      fmpz_CRT(coeffs + k, coeffs + k, q_total, coeffs_mod_q + k, qq + i, 0);
    }
    fmpz_mul(q_total, q_total, qq + i);
  }

  for(i = 0; i < len; i++)
  {
    fmpz_cdiv_q_ui(C, q_total, 2);
    if(fmpz_cmp(coeffs + i, C) > 0)
    {
      fmpz_sub(coeffs + i, coeffs + i, q_total);
    }
  }

  _fmpz_spoly_reserve(res, len);

  for(i = 0; i < len; i++)
  {
    fmpz_set(res->coeffs + i, coeffs + len - i - 1);
    fmpz_set(res->expons + i, expons + len - i - 1);
  }
  res->length = len;

  _fmpz_spoly_normalise(res);

  _fmpz_vec_clear(qq, n);
  _fmpz_vec_clear(ww, n);
  _fmpz_vec_clear(vv, len);
  _fmpz_vec_clear(coeffs, len);
  _fmpz_vec_clear(coeffs_mod_q, len);
  _fmpz_vec_clear(eval1, len);
  _fmpz_vec_clear(eval2, len);

  fmpz_clear(p);
  fmpz_clear(C);
  fmpz_clear(H);
  fmpz_clear(D);
  fmpz_clear(temp);
  fmpz_clear(q_total);
}
