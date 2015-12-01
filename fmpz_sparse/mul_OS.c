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
#include "math.h"

/* mu = .5 / 4 */
void
fmpz_van_prime(fmpz_t res, flint_rand_t state, slong support, slong degree, double gamma, double mu)
{
  if(support * (1 - gamma) <= 1)
    lam = max(21, 10 / (3.0 * mu) * support * eln(degree));
  else
    lam = max(21, 10 / (3.0 * mu) * 1 / (1.0 - mama) * eln(degree));

  fmpz_randprime(res, state, ceil(lam), 1);
}

/* mu = .75 */
slong 
diff_prime(flint_rand_t state, slong support, fmpz_t degree, double gamma, double mu)
{
  double lam;
  if(support * (1 - gamma) <= 1)
    lam = 10 / (3.0 * mu) * (support - 1) * support * fmpz_dlog(degree);
  else
    lam = 10 / (3.0 * mu) * (support - 1) / (1.0 * gamma) * fmpz_dlog(degree);

  if(lam < 21)
    return n_randbits(state, 21);
  else if(fmpz_cmp_si(degree, lam) == 0)
    return n_randbits(state, fmpz_get_si(degree));
  else
    return n_randbits(state, ceil(lam));
}

/* get prime roots to get p with subsequent q's for reduction and relaxation */
/* void get_prime_roots(fmpz * p, fmpz * q) */

/* get sumset in a function */
fmpz_sparse_sumset(fmpz * res, fmpz_sparse_t poly1, fmpz_sparse_t poly2) 
{
  fmpz_sparse_t f_1, g_1, f_2, g_2, h_1, h_2;
  nmod_poly_t f_nmod, g_nmod;
  slong Ss, p, q, degree, R, i;
  fmpz_t L, temp;

  fmpz_sparse_init(f_1);
  fmpz_sparse_init(g_1);
  fmpz_sparse_init(f_2);
  fmpz_sparse_init(g_2);
  fmpz_sparse_init(h_1);
  fmpz_sparse_init(h_2);

  nmod_poly_init(f_nmod, 1);
  nmod_poly_init(g_nmod, 1);

  fmpz_init(L);
  fmpz_init(temp);

  fmpz_sparse_set(f_1, poly1);
  fmpz_sparse_set(g_1, poly2);
  
  if(fmpz_cmp(poly1, poly2) < 0)
    degree = fmpz_get_si(poly1->expons);
  else
    degree = fmpz_get_si(poly2->expons);

  R = (poly1->length + poly2->length);

  for(i = 0; i < poly1->length; i++)
    fmpz_one(f_1->coeffs + i);
 
  for(i = 0; i < poly2->length; i++)
    fmpz_one(g_1->coeffs + i);

  p = diff_prime(state, R * R, 4*degree, 1.0, .125);

  fmpz_sparse_rem_cyc(f_1, f_1, p);
  fmpz_sparse_rem_cyc(g_1, g_1, p);

  Ss = 2;
  for(i = 0; i < max(8 * log(8.0/.5), log(R + 1, 2)); i++)
  {
    q = diff_prime(state, 2*Ss, 2*p, .5, .75);

    /*dense arithmetic to get h_1 = f_1 * g_1 mod R^2 ^ (mod q) */
    nmod_poly_init(h_nmod, q);

    fmpz_sparse_rem_cyc_nmod(f_nmod, f_1, q, R * R);
    fmpz_sparse_rem_cyc_nmod(g_nmod, g_1, q, R * R);
    nmod_poly_mul(h_nmod, f_nmod, g_nmod);

    /*need to somehow enforce the modulus of h_nmod or simply reduce the exponents of h_nmod by q*/

    len = 0;
    for(i = 0; i < h_nmod->length; ++i)
    {
      if(h_nmod->coeffs != 0)
        len++;
    }

    if(len > Ss)
      Ss *= 2;

    nmod_poly_clear(h_nmod);
  }

  /* h_1 and h_2 */
  fmpz_add(L, poly1->expons, poly2->expons);
  fmpz_mul_si(L, L, 8);

  _fmpz_vec_scalar_mul_fmpz(f_2->coeffs, poly1->expons, poly1->length, L);
  _fmpz_vec_scalar_mul_fmpz(g_2->coeffs, poly2->expons, poly2->length, L);

  for(i = 0; i < f_2->length; i++)
    fmpz_add_ui(f_2->coeffs + i, f_2->coeffs + i, 1);

  for(i = 0; i < g_2->length; i++)
    fmpz_add_ui(g_2->coeffs + i, g_2->coeffs + i, 1);

  fmpz_get_si(temp, p);
  fmpz_sparse_rem_cyc(f_2, f_2, temp);
  fmpz_sparse_rem_cyc(g_2, g_2, temp);

  fmpz_mul_interp(h_1, state, f_1, g_1, Ss);
  fmpz_mul_interp(h_2, state, f_2, g_2, Ss);

  fmpz_sparse_rem_cyc(h_1, h_1, temp);
  fmpz_sparse_rem_cyc(h_2, h_2, temp);

  fmpz_mul(temp, L, L);
  
  _fmpz_vec_scalar_mod_fmpz(h_1->coeffs, h_1->coeffs, h_1->length, temp);
  _fmpz_vec_scalar_mod_fmpz(h_2->coeffs, h_2->coeffs, h_2->length, temp);
  
  _fmpz_sparse_normalise(h_1);
  _fmpz_sparse_normalise(h_2);

  _fmpz_vec_init(res, h_2->length);

  /* for every non_zero term of h_1 find the coefficient of h_2 and retrieve the exponent
   * corresponding to the non-zero term in h_1*/
  
  for(i = 0; i < h_2->length; i++)
  {
    fmpz_cdiv_q(res + i, h_2->coeffs + i, h_1->coeffs + i);
    fmpz_sub_ui(res + i, h_2->coeffs + i, 1);
    fmpz_cdiv_q(res + i, h_2->coeffs + i, L);
  }

  fmpz_clear(L);
  fmpz_clear(temp);
  nmod_poly_clear(f_nmod);
  nmod_poly_clear(g_nmod);
  fmpz_sparse_clear(f_1);
  fmpz_sparse_clear(g_1);
  fmpz_sparse_clear(f_2);
  fmpz_sparse_clear(g_2);
  fmpz_sparse_clear(h_1);
  fmpz_sparse_clear(h_2);
}

void 
fmpz_sparse_mul_OS(fmpz_sparse_t res, flint_rand_t state, const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2)
{
  /*PRESSING PROBLEMS ARE THAT THERE ARE NO PRIME QUALITY EVALUATION IMPLEMENTED*/
  /*
   * (substitutions for kronecker substitution are ignored for this implementation)
   * make copies for f, g, f_k, g_k, f_s, and g_s
   * -----------f_s and g_s have coefficients of 1
   * estimate structural sparsity
   * -----------choose primes p and p' to modulate the polynomials by ((f_s*g_S)^mod p)mod p' until half dense
   * compute structural support
   * -----------choose same prime p and make h_1 = ((f_s*g_s)^mod p)mod p
   * -----------f_2 = the sum of (e*L + 1)*z^(e mod p)
   * -----------g_2 = the sum of (e*L + 1)*z^(e mod p)
   * -----------h_2 = ((f_2*g_2)^mod p)mod L^2
   * -----------compute exponents of S with ratio ((c_2/c_1)-1)/L
   * compute arithmetic support
   * -----------f_k and g_k mod p and q ((f_k*g_k)^mod p) mod q where p|(q-1)
   * compute the coefficients
   * -----------multiple snapshots of ((f_k*g_k)^mod p) mod q with ascending q and the same p
   * -----------group the like terms of the snapshots and then compute CRT
   */
  fmpz_sparse_t f, g, f_s, g_s, h_1, h_2, f_2, g_2, temp_1, temp_2;
  nmod_poly_t f_nmod, g_nmod, h_nmod;
  fmpz_t p, q;
  slong i, len;

  fmpz_sparse_init(f);
  fmpz_sparse_init(g);
  fmpz_sparse_init(f_s);
  fmpz_sparse_init(g_s);
  fmpz_sparse_init(f_2);
  fmpz_sparse_init(g_2);
  fmpz_sparse_init(h_1);
  fmpz_sparse_init(h_2);
  fmpz_sparse_init(temp_1);
  fmpz_sparse_init(temp_2);
 
  /*nmod must be initialized to some modulus*/
  nmod_poly_init(f_nmod, 1);
  nmod_poly_init(g_nmod, 1);

  fmpz_init(p);
  fmpz_init(q);

  fmpz_sparse_set(f, poly1);
  fmpz_sparse_set(g, poly2);

  flint_printf("BEFORE:\n");
  fmpz_sparse_print(f);
  flint_printf("\n\n");
  fmpz_sparse_print(g);
  flint_printf("\n\n");

  /* make f_s and g_s*/
  fmpz_sparse_set(f_s, poly1);
  fmpz_sparse_set(g_s, poly2);

  for(i = 0; i < poly1->length; i++)
  {
    fmpz_one(f_s->coeffs + i);
  }

  for(i = 0; i < poly2->length; i++)
  {
    fmpz_one(g_s->coeffs + i);
  }
  
  /* estimate structural sparsity */
  /* PRIMES NOTE: IMPORTANT TO DEFINE D and R
  fmpz_diff_prime(p, state, (2*R)**2, 4*D, 1, 1.0/32.0);*/
  
  fmpz_sparse_rem_cyc(f_s, f_s, p);
  fmpz_sparse_rem_cyc(g_s, g_s, p);
  
  /* PRIMES CAN ELIMINATE LOOP WITH A PROPER Q CANDIDATE PRODUCED VIA SAGE METHOD */
  do
  {

    /*need some form of get next prime (does this need to be an fmpz?)*/
    /* PRIMES fmpz_randprime(q, state, 6, 1);*/
    nmod_poly_init(h_nmod, fmpz_get_si(q));

    fmpz_sparse_rem_cyc_nmod(f_nmod, f_s, fmpz_get_si(q), (f_s->length + g_s->length)*(f_s->length + g_s->length));
    fmpz_sparse_rem_cyc_nmod(g_nmod, g_s, fmpz_get_si(q), (f_s->length + g_s->length)*(f_s->length + g_s->length));
    nmod_poly_mul(h_nmod, f_nmod, g_nmod);

    /*need to somehow enforce the modulus of h_nmod or simply reduce the exponents of h_nmod by q*/

    /*calculate how many terms are in h_nmod*/
    len = 0;
    for(i = 0; i < h_nmod->length; ++i)
    {
      if(h_nmod->coeffs != 0)
        len++;
    }

    nmod_poly_clear(h_nmod);
  } while(len >= (fmpz_get_si(p))/2);
  
  /*flint_printf("AFTER:\nLENGTH: %ld\n", len);
  fmpz_print(p);
  flint_printf("\n\n");
  fmpz_print(q);
  flint_printf("\n\n");
  fmpz_sparse_print(f_s);
  flint_printf("\n\n");
  fmpz_sparse_print(g_s);
  flint_printf("\n\n");
  nmod_poly_print(h_nmod);
  flint_printf("\n\n");*/


  /* compute structural support */

  /* h_1 */
  /* PRIMES f_s and g_s are modulated by p from previous calculations may need to be rebuilt if p is different */
  fmpz_sparse_mul_heaps(h_1, f_s, g_s);
  fmpz_sparse_rem_cyc(h_1, h_1, p);

  /* f_2 */
  
  fmpz_add(q, f->expons, g->expons);
  fmpz_mul_ui(q, q, 2);
 
  _fmpz_vec_set(f_s->coeffs, f_s->expons, f_s->length);
  _fmpz_vec_scalar_mul_fmpz(f_s->coeffs, f_s->coeffs, f_s->length, q);

  for(i = 0; i < f_s->length; i++)
  {
    fmpz_add_ui(f_s->coeffs + i, f_s->coeffs + i, 1);
  }

  fmpz_sparse_rem_cyc(f_s, f_s, p);

  /* g_2 */
 
  _fmpz_vec_set(g_s->coeffs, g_s->expons, g_s->length);
  _fmpz_vec_scalar_mul_fmpz(g_s->coeffs, g_s->coeffs, g_s->length, q);

  for(i = 0; i < g_s->length; i++)
  {
    fmpz_add_ui(g_s->coeffs + i, g_s->coeffs + i, 1);
  }

  fmpz_sparse_rem_cyc(g_s, g_s, p);

  /* h_2 MULINTERP (sparse polynomial interpolation is used here but how?)*/
  fmpz_sparse_mul_heaps(h_2, f_s, g_s);
  fmpz_sparse_rem_cyc(h_2, h_2, p);

  fmpz_mul(p, q, q);
  _fmpz_vec_scalar_mod_fmpz(h_2->coeffs, h_2->coeffs, h_2->length, p);

  _fmpz_sparse_normalise(h_2);

  /* make coefficient ratios */
  for(i = 0; i < h_2->length; i++)
  {
    fmpz_cdiv_q(h_2->coeffs + i, h_2->coeffs + i, h_1->coeffs + i);
    fmpz_sub_ui(h_2->coeffs + i, h_2->coeffs + i, 1);
    fmpz_cdiv_q(h_2->coeffs + i, h_2->coeffs + i, q);
  }

  /* NOTE: the coefficients of h_2 should be the structural support */
  
  /* compute arithmetic support */

  /* MULINTERP NEEDS SPARSE INTERPOLATION LIKE STEP ABOVE */
  /*fmpz_randprime(p, state, 5, 1);*/
  fmpz_mul_ui(q, p, 2);
  fmpz_add_ui(q, q, 1);
  fmpz_sparse_rem_cyc(f, f, p);
  fmpz_sparse_rem_cyc(g, g, p);

  flint_printf("AFTER:\n");
  fmpz_print(q);
  flint_printf("\n\n");
  fmpz_sparse_print(h_2);
  flint_printf("\n\n");

  fmpz_sparse_clear(f);
  fmpz_sparse_clear(g);
  fmpz_sparse_clear(f_s);
  fmpz_sparse_clear(g_s);
  fmpz_sparse_clear(f_2);
  fmpz_sparse_clear(g_2);
  fmpz_sparse_clear(h_1);
  fmpz_sparse_clear(h_2);
  fmpz_sparse_clear(temp_1);
  fmpz_sparse_clear(temp_2);

  /*nmod must be initialized to some modulus*/
  nmod_poly_clear(f_nmod, 1);
  nmod_poly_clear(g_nmod, 1);

  fmpz_clear(p);
  fmpz_clear(q);
}
