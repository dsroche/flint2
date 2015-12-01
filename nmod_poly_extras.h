/* Extra utility functions for working with nmod_poly objects. */

#ifndef NMOD_POLY_EXTRAS_H
#define NMOD_POLY_EXTRAS_H

#include <stdio.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>

/* sets res = f mod (x^e - 1)*/
void nmod_poly_rem_cyc(nmod_poly_t res, const nmod_poly_t f, ulong e);

/* convert to fmpz_poly*/
void nmod_poly_convert(fmpz_poly_t res, const nmod_poly_t f);

/* pretty print functions just like in fmpz_poly*/
int nmod_poly_fprint_pretty(FILE* file, const nmod_poly_t f, const char* x);

#endif /* NMOD_POLY_EXTRAS_H*/
