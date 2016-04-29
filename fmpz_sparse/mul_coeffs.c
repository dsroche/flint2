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

#include "fmpz_sparse.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

static const double LN_2 = 0.693147180559945309417232121458;

FLINT_DLL void 
_fmpz_sparse_mul_coeffs(fmpz_sparse_t res, flint_rand_t state,
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2, const fmpz * expons,
    slong len)
{
  fmpz_mod_poly_t poly;
  fmpz * qq, * ww, * vv, * coeffs, * mod_expons, * eval1, * eval2, * coeffs_mod_q;
  fmpz_t p, C, H, T, D, temp, q_total;
  slong p_bits, q_prod_bits, n, num_primes, i, j, k;

  fmpz_init(p);
  fmpz_init(C);
  fmpz_init(H);
  fmpz_init(T);
  fmpz_init(D);
  fmpz_init(temp);
  fmpz_init(q_total);

  fmpz_sparse_degree(D, poly1);
  fmpz_sparse_degree(temp, poly2);

  if(fmpz_cmp(temp, D) > 0)
    fmpz_set(D, temp);

  fmpz_sparse_height(H, poly1);
  fmpz_sparse_height(temp, poly2);
  fmpz_mul(C, temp, H);

  p_bits = 2*fmpz_bits(T) + log(fmpz_bits(D))/LN_2 + 1;

  q_prod_bits = fmpz_bits(C) + log(FLINT_MAX(poly1->length, poly2->length))/LN_2 + 1;

  n = (q_prod_bits + p_bits - 1)/p_bits + 1;

  qq = _fmpz_vec_init(n);
  ww = _fmpz_vec_init(n);

  vv = _fmpz_vec_init(len);
  coeffs = _fmpz_vec_init(len);
  coeffs_mod_q = _fmpz_vec_init(len);
  mod_expons = _fmpz_vec_init(len);
  eval1 = _fmpz_vec_init(len);
  eval2 = _fmpz_vec_init(len);

  num_primes = _fmpz_sparse_prim_roots(p, qq, ww, state, n, p_bits, q_prod_bits);
  
  flint_printf("p: "), fmpz_print(p);
  flint_printf("\n");

  flint_printf("qq: "), _fmpz_vec_print(qq, num_primes);
  flint_printf("\n");

  flint_printf("ww: "), _fmpz_vec_print(ww, num_primes);
  flint_printf("\n");
  
  _fmpz_vec_scalar_mod_fmpz(mod_expons, expons, len, p);

  for(i = 0; i < num_primes; i++)
  {
    for(j = 0; j < len; j++)
    {
      fmpz_powm(vv + j, ww + i, mod_expons + j, qq + i);
    }
  
    for(j = 0; j < len; j++)
    {
      fmpz_sparse_evaluate_mod(eval1 + j, poly1, vv + j, qq + i);
      fmpz_sparse_evaluate_mod(eval2 + j, poly2, vv + j, qq + i);
    }

    for(j = 0; j < len; j++)
    {
      fmpz_mul(eval1 + j, eval1 + j, eval2 + j);
    }

    fmpz_mod_poly_init(poly, qq + i);
    _fmpz_mod_poly_build_roots(poly, vv, len);


    flint_printf("\ntree coeffs: "), _fmpz_vec_print(poly->coeffs, poly->length);
    flint_printf("\n");

    flint_printf("vv: "), _fmpz_vec_print(vv, len);
    flint_printf("\n");
    
    flint_printf("eval1: "), _fmpz_vec_print(eval1, len);
    flint_printf("\n");
    
    flint_printf("q: "), fmpz_print(qq+i);
    flint_printf("\n");
    
    _fmpz_mod_poly_transposed_vandermonde(coeffs_mod_q, vv, eval1, len, poly->coeffs, qq + i);

    fmpz_mod_poly_clear(poly);

    fmpz_one(q_total);

    flint_printf("\ncoeffs_mod_q: "), _fmpz_vec_print(coeffs_mod_q, len);
    flint_printf("\n");

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
      /*flint_printf("\ncoeffs + k: "), fmpz_print(coeffs + k);
      flint_printf("\nq_total: "), fmpz_print(q_total);
      flint_printf("\ncoeffs_mod_q + k: "), fmpz_print(coeffs_mod_q + k);
      flint_printf("\nqq + i: "), fmpz_print(qq + i);
      flint_printf("\n");*/
      fmpz_CRT(coeffs + k, coeffs + k, q_total, coeffs_mod_q + k, qq + i, 1);
    }
    fmpz_mul(q_total, q_total, qq + i);
  }

  for(i = 0; i < len; i++)
  {
    fmpz_set(res->coeffs + i, coeffs + i);
    fmpz_set(res->expons + i, expons + len - i - 1);
  }


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
  fmpz_clear(T);
  fmpz_clear(D);
  fmpz_clear(temp);
  fmpz_clear(q_total);
}
