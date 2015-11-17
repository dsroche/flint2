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

    Authored 2015 by Daniel S. Roche and A. Whitman Groves;
    US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_vec.h"

void fmpz_sparse_new_randtest(fmpz_sparse_t res, flint_rand_t state, 
    slong terms, const fmpz_t degree, mp_bitcnt_t bits)
{
  slong i, j, k, length, unique, neg;
  fmpz * rands;
  fmpz_t abs, half, temp;

  unique = 0;
  
  neg = fmpz_sgn(degree);

  fmpz_sparse_zero(res);
  
  if(terms == 0)
  {
    return;
  }
  
  fmpz_abs(abs,degree);

  if(fmpz_cmp_si(abs, terms) < 0)
  {
    terms = fmpz_get_ui(abs);
    if(fmpz_is_zero(abs))
        terms = 1;
  }
  
  _fmpz_sparse_reserve(res, terms);
  
  fmpz_init(res->coeffs);
  fmpz_init(res->expons);
  fmpz_randbits(res->coeffs, state, bits);
  fmpz_set(res->expons, abs);
  
  if(fmpz_is_one(abs) || terms == 1)
  {
    res->length = 1;
    fmpz_clear(abs);
    return;
  }

  fmpz_sub_ui(abs, abs, 1);
  
  /*if terms > 1/2 * degree choose D - T unique exponents*/
  fmpz_fdiv_q_ui(half, abs, 2);
  fmpz_set_ui(temp, terms);

  if(fmpz_cmp(temp, half) > 0)
  {
    length = fmpz_get_ui(abs) + 1 - terms ;
  }
  /*if terms < 1/2 * degree choose T - 1 unique exponents*/
  else
  {
    length = terms - 1;
  }

  if(length == 0)
  {
    for(i=1; i<terms; ++i)
    {
      fmpz_init(res->coeffs + i);
      fmpz_init(res->expons + i);
      fmpz_randbits(res->coeffs + i, state, bits);
      fmpz_randm(res->expons + i, state, abs);
      if(neg == -1)
        fmpz_mul_si(res->expons + i, res->expons + i,  1 - 2*n_randint(state,2));
    }

    res->length = terms;
    _fmpz_sparse_normalise(res);

  }
  else
  {
    rands = _fmpz_vec_init(length);
    _fmpz_vec_randtest_unsigned(rands, state, length, fmpz_bits(abs));
    
    /*checks to make sure that all random exponents are unique 2*/
    while(unique == 0)
    {
      unique = 1;
      
      /*fmpz_vec sorts in ascending order*/
      _fmpz_vec_sort(rands, length);
      
      for(i = 0; i < length; i++)
      {
        if(fmpz_equal(rands + i, rands + i + 1))
        {
          fmpz_randm(rands + i, state, abs);
          unique = 0;
        }
      }
    }
    
    /*if terms <= .5 * degree*/
    if(length == terms - 1)
    {
      for(i=1; i < terms; ++i)
      {
        fmpz_set(res->expons + i, rands + i - 1);
      }
    }
    /*if terms > .5 * degree*/
    else
    {
      k = 0;
      j = 0;
      
      /*fmpz_vec rands is sorted in ascending order*/
      for(i=0; i < fmpz_get_si(abs); ++i)
      {
        if(fmpz_get_si(rands + j) == i)
          j++;
        else
        {
          fmpz_set_si(res->expons + k, i);
          k++;
        }
      }
    }
    
    for (i=1; i<terms; ++i) 
    {
      if(neg == -1)
        fmpz_mul_si(res->expons + i, res->expons + i,  1 - 2*n_randint(state,2));
      
      fmpz_init(res->coeffs+i);
      fmpz_randbits(res->coeffs+i, state, bits);
    }
    
    
    
    flint_free(rands);
  }
  
  res->length = terms;
  _fmpz_sparse_normalise(res);
  
  fmpz_clear(half);
  fmpz_clear(temp);
  fmpz_clear(abs);
}
