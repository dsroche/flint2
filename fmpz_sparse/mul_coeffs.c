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

******************************************************************************/

/*TODO
 * p should be larger than T^2log(D)
 * lg(P) > 2lg(T) + lglg(D)
 *
 * n > lg(H) / lg(P) = q_prod_bits / p_bits + 1
 *
 * CRT(coeff[i], q_total, coeff[i] % q[j], q[j])*/

#include "fmpz_sparse.h"
#include "fmpz_vec.h"

static const double LN_2 = 0.693147180559945309417232121458;

FLINT_DLL void 
_fmpz_sparse_mul_coeffs(fmpz_sparse_t res, flint_rand_t state,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2, const fmpz * expons,
    slong len)
{
  fmpz * qq, * ww, * vv, * coeffs, * mod_expons, * eval1, * eval2;
  fmpz_t p, C, H, T, D, temp, q_total;
  slong p_bits, q_prob_bits, n, num_primes, i, j, k;

  fmpz_init(p);
  fmpz_init(C);
  fmpz_init(H);
  fmpz_init(T);
  fmpz_init(D);
  fmpz_init(temp);
  fmpz_init(q_total);

  fmpz_sparse_get_degree(D, poly1);
  fmpz_sparse_get_degree(temp, poly2);

  if(fmpz_cmp(temp, D) > 0)
    fmpz_set(D, temp);

  fmpz_sparse_get_height(H, poly1);
  fmpz_sparse_get_height(temp, poly2);
  fmpz_mul(C, temp, H);

  p_bits = 2*fmpz_get_bits(T) + log(fmpz_get_bits(D))/LN_2 + 1;

  q_prod_bits = fmpz_get_bits(C) + log(FLINT_MAX(poly1->length, poly2->length))/LN_2 + 1;

  n = (q_prod_bits + p_bits - 1)/p_bits + 1;

  qq = _fmpz_vec_init(n);
  ww = _fmpz_vec_init(n);

  vv = _fmpz_vec_init(len);
  coeffs = _fmpz_vec_init(len);
  mod_expons = _fmpz_vec_init(len);
  eval1 = _fmpz_vec_init(len);
  eval2 = _fmpz_vec_init(len);

  num_primes = _fmpz_sparse_prim_roots(p, qq, ww, state, n, p_bits, q_prod_bits);
  
  _fmpz_vec_scalar_mod(mod_expons, expons, len, p);

  fmpz_one(q_total);

  for(i = 0, i < num_primes; i++)
  {
    for(j = 0; j < len; j++)
    {
      fmpz_powm(vv + j, ww + i, mod_expons + j, qq + i);
    }
  
    for(j = 0, j < len; j++)
    {
      /*
      * eval f at each vv + j
      * eval g at each vv + j
      * */
      fmpz_sparse_evaluate_mod(eval1 + j, poly1, vv + j, qq + i);
      fmpz_sparse_evaluate_mod(eval2 + j, poly2, vv + j, qq + i);
    }

    for(j = 0; j < len; j++)
    {
      fmpz_mul(eval1 + j, eval1 + j, eval2 + j);
    }

    /* 
     * transposed vandermonde calculates coeffs mod q
     *    go to bp_interp
     *      use void _fmpz_mod_poly_transposed_vandermonde(fmpz* coeffs_mod_q,
     *          const fmpz* vv, const fmpz* evals_prod, slong len, 
     *          const fmpz* poly, const fmpz_t qq + i)
     * */
    /*TODO use 35.28 and the test code to implement _fmpz_mod_poly_tree from vv*/
    _fmpz_mod_poly_transposed_vandermonde(coeffs_mod_q, vv, eval1, len, /*TODO*/poly, qq + i);

    for(k = 0; k < len; k++)
    {
      /*CRT
       * void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1,
       * fmpz_t r2, fmpz_t m2, int sign)
       * Use the Chinese Remainder Theorem to set \code{out} to the unique value
       * $0 \le x < M$ (if sign = 0) or $-M/2 < x \le M/2$ (if sign = 1)
       * congruent to $r_1$ modulo $m_1$ and $r_2$ modulo $m_2$,
       * where where $M = m_1 \times m_2$.
       * It is assumed that $m_1$ and $m_2$ are positive integers greater
       * than $1$ and coprime.
       * 
       * If sign = 0, it is assumed that $0 \le r_1 < m_1$ and $0 \le r_2 < m_2$.
       * Otherwise,it is assumed that $-m_1 \le r_1 < m_1$ and $0 \le r_2 < m_2$.
       * */
      /*TODO coeffs_mod_q or coeffs or res->coeffs? ??? TODO*/
      fmpz_CRT(coeffs + k, coeffs + k, q_total, coeffs_mod_q + k, qq + i, 1);
      fmpz_mul(q_total, q_total, qq + i);
    }
  }


  _fmpz_vec_clear(qq, n);
  _fmpz_vec_clear(ww, n);
  _fmpz_vec_clear(vv, len);
  _fmpz_vec_clear(coeffs, len);
  _fmpz_vec_clear(eval1, len);
  _fmpz_vec_clear(eval2, len);

  fmpz_clear(p);
  fmpz_clear(C);
  fmpz_clear(H);
  fmpz_clear(T);
  fmpz_clear(D);
  fmpz_clear(temp);
  fmpz_clear(q_total);
}
