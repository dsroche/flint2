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

#include "fmpz_spoly.h"
#include "fmpz_vec.h"
#include "math.h"

void
fmpz_van_prime(fmpz_t res, flint_rand_t state, slong support, mp_bitcnt_t deg_bits, double gamma)
{
  mp_bitcnt_t bits;
  double mu = .125;
  /*(10*ln(2)/3 = 2.31049060187)*/
  double lam = 2.31049060187*deg_bits/mu;

  if(support * (1 - gamma) <= 1)
    lam *= (double) support;
  else
    lam /= 1.0 - gamma;

  lam = ceil(log(lam)/0.693147180559945309417);

  bits = lam;
  if(bits < 5)
    bits = 5;

  fmpz_randprime(res, state, bits, 0);
}

void 
fmpz_diff_prime(fmpz_t res, flint_rand_t state, slong support, mp_bitcnt_t deg_bits, double gamma)
{
  mp_bitcnt_t bits;
  double mu = .75;
  /*(10*ln(2)/3 = 2.31049060187*/
  double lam = 2.31049060187*deg_bits*(support-1)/mu;

  if(support * (1 - gamma) <= 1)
    lam *= (double) support;
  else
    lam /= 1.0 - gamma;

  lam = ceil(log(lam)/0.693147180559945309417);

  bits = lam;
  if(lam < 5)
    bits = 5;
  else if(lam > deg_bits)
    bits = deg_bits;

  fmpz_randprime(res, state, bits, 0);
}

slong
fmpz_spoly_sumcheck(fmpz ** res, const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
  slong i, j, k, len;
  fmpz * test, * final;

  len = poly1->length * poly2->length;

  if(len == 0)
  {
    return 0;
  }
  test = _fmpz_vec_init(len);

  for(i = 0; i < poly1->length; i++)
  {
    for(j = 0; j < poly2->length; j++)
    {
      fmpz_add(test + j + i*poly2->length, poly1->expons + i, poly2->expons + j);
    }
  }
  
  _fmpz_vec_sort(test,len);

  j = 1;
  for(i = 0; i < len - 1; i++)
  {
    if(fmpz_cmp(test + i, test + i + 1) != 0)
      j++;
  }

  k = j;
  final = _fmpz_vec_init(k);

  for(i = len - 1; i > 0; i--)
  {
    if(fmpz_cmp(test + i, test + i - 1) != 0)
    {
      j--;
      FLINT_ASSERT(j >= 0 && j < k);
      FLINT_ASSERT(i - 1 >= 0 && i < len);
      fmpz_set(final + j, test + i);
    }
  }

  fmpz_set(final + 0, test + 0);

  *res = final;

  _fmpz_vec_clear(test, len);

  return k;
}

/* get sumset in a function */
slong
fmpz_spoly_sumset(fmpz ** res, flint_rand_t state, const fmpz_spoly_t poly1, const fmpz_spoly_t poly2) 
{
  fmpz_spoly_t f_1, g_1, f_2, g_2, h_1, h_2;
  nmod_poly_t f_nmod, g_nmod, h_nmod;
  slong len, Ss, R, i, j, q;
  fmpz_t L, p, temp, degree;
  fmpz * final;

  q = 0;
  if(fmpz_spoly_is_zero(poly1) || fmpz_spoly_is_zero(poly2))
  {
    return 0;
  }

  if(poly1->length == 1 && poly2->length == 1)
  {
    final = _fmpz_vec_init(1);
    fmpz_add(final, poly1->expons+0, poly2->expons+0);
    *res = final;
    return 1; 
  }

  if(poly1->length == 1)
  {
    final = _fmpz_vec_init(poly2->length);
    len = poly2->length;

    for(i = 0; i < len; i++)
      fmpz_add(final + i, poly2->expons + i, poly1->expons+0);
    _fmpz_vec_sort(final, len);
    *res = final; 
    return len;
  }
  else if(poly2->length == 1)
  {
    final = _fmpz_vec_init(poly1->length);
    len = poly1->length;
    
    for(i = 0; i < len; i++)
      fmpz_add(final + i, poly1->expons + i, poly2->expons+0);
    _fmpz_vec_sort(final, len);
    *res = final; 
    return len;
  }

  fmpz_spoly_init(f_1);
  fmpz_spoly_init(g_1);
  fmpz_spoly_init(f_2);
  fmpz_spoly_init(g_2);
  fmpz_spoly_init(h_1);
  fmpz_spoly_init(h_2);

  R = (poly1->length + poly2->length);
  
  nmod_poly_init(f_nmod, R*R);
  nmod_poly_init(g_nmod, R*R);
  nmod_poly_init(h_nmod, R*R);

  fmpz_init(L);
  fmpz_init(temp);
  fmpz_init(p);
  fmpz_init(degree);

  fmpz_spoly_set(f_1, poly1);
  fmpz_spoly_set(g_1, poly2);
 
  /*TODO 
   * Try to do some rough benchmarking (once base case multiply is implemented)
   * Do some FIXME's and replace some "poly->expon + i" with "fmpz_spoly_get_coeff"
   */

  if(fmpz_cmp(poly1->expons + 0, poly2->expons + 0) < 0)
    fmpz_set(degree, poly2->expons + 0);
  else
    fmpz_set(degree, poly1->expons + 0);

  if(fmpz_cmp(degree, poly2->expons + poly2->length - 1) < 0)
    fmpz_set(degree, poly2->expons + poly2->length - 1);
  
  if(fmpz_cmp(degree, poly1->expons + poly1->length - 1) < 0)
    fmpz_set(degree, poly1->expons + poly1->length - 1);

  fmpz_diff_prime(p, state, R * R, fmpz_bits(degree) + 2, 1.0);
  
  for(i = 0; i < poly1->length; i++)
    fmpz_one(f_1->coeffs + i);
 
  for(i = 0; i < poly2->length; i++)
    fmpz_one(g_1->coeffs + i);

  fmpz_spoly_rem_cyc(f_1, f_1, p);
  fmpz_spoly_rem_cyc(g_1, g_1, p);

  Ss = 2;
  for(i = 0; i < FLINT_MAX(22.1807097779182499013514, log(R)/0.693147180559945309417232121458 + 1); i++)
  {
    fmpz_diff_prime(temp, state, 2*Ss, fmpz_bits(p) + 1, .5);
  
    q = fmpz_get_si(temp);
    
    fmpz_spoly_rem_cyc_nmod(f_nmod, f_1, q, R*R);
    fmpz_spoly_rem_cyc_nmod(g_nmod, g_1, q, R*R);
    nmod_poly_mul(h_nmod, f_nmod, g_nmod);

    nmod_poly_rem_cyc(h_nmod, h_nmod, q);

    len = 0;
    for(j = 0; j < nmod_poly_length(h_nmod); ++j)
    {
      if(nmod_poly_get_coeff_ui(h_nmod, j) != 0)
        len++;
    }

    if(len > Ss)
      Ss *= 2;
  }
  
  /* h_1 and h_2 */
  fmpz_add(L, poly1->expons+0, poly2->expons+0);
  fmpz_mul_si(L, L, FLINT_MIN(poly1->length, poly2->length));
  fmpz_add_ui(L, L, 1);

  fmpz_spoly_set(f_2, poly1);
  fmpz_spoly_set(g_2, poly2);

  _fmpz_vec_scalar_mul_fmpz(f_2->coeffs, f_2->expons, f_2->length, L);
  _fmpz_vec_scalar_mul_fmpz(g_2->coeffs, g_2->expons, g_2->length, L);

  for(i = 0; i < f_2->length; i++)
    fmpz_add_ui(f_2->coeffs + i, f_2->coeffs + i, 1);

  for(i = 0; i < g_2->length; i++)
    fmpz_add_ui(g_2->coeffs + i, g_2->coeffs + i, 1);
  
  fmpz_spoly_rem_cyc(f_2, f_2, p);
  fmpz_spoly_rem_cyc(g_2, g_2, p);

  fmpz_spoly_mul_interp(h_1, state, f_1, g_1, Ss);
  _fmpz_spoly_reserve(h_2, h_1->length);
  _fmpz_vec_set(h_2->expons, h_1->expons, h_1->length);
  _fmpz_spoly_set_length(h_2, h_1->length);
  _fmpz_spoly_mul_coeffs(h_2, f_2, g_2);

  fmpz_spoly_rem_cyc(h_1, h_1, p);
  fmpz_spoly_rem_cyc(h_2, h_2, p);

  fmpz_mul(temp, L, L);
  
  _fmpz_vec_scalar_mod_fmpz(h_1->coeffs, h_1->coeffs, h_1->length, temp);
  _fmpz_vec_scalar_mod_fmpz(h_2->coeffs, h_2->coeffs, h_2->length, temp);
  
  final = _fmpz_vec_init(h_2->length);

  /* for every non_zero term of h_1 find the coefficient of h_2 and retrieve the exponent
   * corresponding to the non-zero term in h_1*/
  
  i = 0;
  R = 0;
  len = 0;

  while(i < h_1->length && R < h_2->length)
  {
    if(fmpz_equal(h_1->expons + i, h_2->expons + R))
    {
      fmpz_cdiv_q(final + len, h_2->coeffs + R, h_1->coeffs + i);
      fmpz_sub_ui(final + len, final + len, 1);
      fmpz_cdiv_q(final + len, final + len, L);
      i++;
      R++;
      len++;
    }
    else
    {
      if(fmpz_cmp(h_1->expons + i, h_2->expons + R) < 0)
        R++;
      else
        i++;
    }
  }

  _fmpz_vec_sort(final, len);
  *res = final;

  fmpz_clear(L);
  fmpz_clear(temp);
  fmpz_clear(p);
  fmpz_clear(degree);
  nmod_poly_clear(f_nmod);
  nmod_poly_clear(g_nmod);
  nmod_poly_clear(h_nmod);
  fmpz_spoly_clear(f_1);
  fmpz_spoly_clear(g_1);
  fmpz_spoly_clear(f_2);
  fmpz_spoly_clear(g_2);
  fmpz_spoly_clear(h_1);
  fmpz_spoly_clear(h_2);
  
  return len;
}
