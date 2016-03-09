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

void
fmpz_sparse_mul_OS(fmpz_sparse_t res, flint_rand_t state, const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
  fmpz * test;
  slong length;

  test = NULL;

  
  flint_printf("\nthere is no god\n");
  length = fmpz_sparse_sumset(&test, state, poly1, poly2);

  if(length > 1)
  {
    _fmpz_vec_print(test, length);
    _fmpz_vec_clear(test, length);
  }
  else
    flint_printf("0");
  flint_printf("\n");
  
  /*

  length = fmpz_sparse_sumcheck(&temp, poly1, poly2);

  _fmpz_vec_print(temp, length);

  flint_printf("\n");

  _fmpz_vec_clear(temp, length);
  */

  fmpz_sparse_zero(res);
  
}

#if 0

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

#endif
