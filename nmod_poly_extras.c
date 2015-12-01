/* Extra utility functions for working with nmod_poly objects. */

#include <assert.h>
#include "nmod_poly_extras.h"

/* sets res = f mod (x^e - 1)*/
void nmod_poly_rem_cyc(nmod_poly_t res, const nmod_poly_t f, ulong e)
{
  slong i;
  ulong rese=0;
  assert (nmod_poly_modulus(res) == nmod_poly_modulus(f));
  ulong q = nmod_poly_modulus(f);
  
  if (res != f) { /* set low-order coeffs of res to those of f*/
    nmod_poly_truncate(res, e);
    for (i=e-1; i>=0; --i) {
      nmod_poly_set_coeff_ui(res, i, nmod_poly_get_coeff_ui(f, i));
    }
  }
  
  for (i=e; i<nmod_poly_length(f); ++i) {
    nmod_poly_set_coeff_ui(res, rese, n_addmod(nmod_poly_get_coeff_ui(res, rese), 
          nmod_poly_get_coeff_ui(f, i), q));
    if (++rese == e) 
      rese = 0;
  }
  
  if (res == f) { /* aliasing; time to truncate the high-order coeffs*/
    nmod_poly_truncate(res, e);
  }
}

/* convert to fmpz_poly*/
void nmod_poly_convert(fmpz_poly_t res, const nmod_poly_t f)
{
  slong i;
  ulong maxpos = (nmod_poly_modulus(f)-1)/2; /* the largest positive remainder*/
  ulong makeneg = (ulong) (-((slong)nmod_poly_modulus(f)));
  fmpz_poly_zero(res);
  for (i=nmod_poly_degree(f); i>=0; --i) {
    ulong c = nmod_poly_get_coeff_ui(f, i);
    if (c > maxpos) c += makeneg;
    fmpz_poly_set_coeff_si(res, i, (slong)c);
  }
}

/* pretty print functions just like in fmpz_poly*/
int nmod_poly_fprint_pretty(FILE* file, const nmod_poly_t f, const char* x)
{
  int ret;
  fmpz_poly_t temp;
  fmpz_poly_init2(temp, nmod_poly_length(f));
  nmod_poly_convert(temp, f);
  ret = fmpz_poly_fprint_pretty(file, temp, x);
  fmpz_poly_clear(temp);
  return ret;
}
