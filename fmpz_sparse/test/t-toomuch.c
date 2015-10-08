/* test program for fmpz_sparse */

#include "fmpz_sparse.h"

void nmod_poly_rem_cyc(nmod_poly_t res, const nmod_poly_t f, ulong e)
{
  slong i;
  ulong rese=0;
  ulong q = nmod_poly_modulus(f);
  FLINT_ASSERT (nmod_poly_modulus(res) == nmod_poly_modulus(f));

  if (res != f) { /* set low-order coeffs of res to those of f */
    nmod_poly_truncate(res, e);
    for (i=e-1; i>=0; --i) {
      nmod_poly_set_coeff_ui(res, i, nmod_poly_get_coeff_ui(f, i));
    }
  }

  for (i=e; i<nmod_poly_length(f); ++i) {
    nmod_poly_set_coeff_ui(res, rese, n_addmod(
      nmod_poly_get_coeff_ui(res, rese), nmod_poly_get_coeff_ui(f, i), q
    ));
    if (++rese == e) rese = 0;
  }

  if (res == f) { /* aliasing; time to truncate the high-order coeffs */
    nmod_poly_truncate(res, e);
  }
}

int main(void) {
  int retval = 0;

  flint_rand_t rs;
  flint_randinit(rs);

  { /* test 1 */
    fmpz_sparse_t f, g, h;
    fmpz_poly_t hmod;
    nmod_poly_t hmodq;
    nmod_poly_t hdiverse;
    fmpz_t D;
    slong e = 5;
    ulong q = 11;

    fmpz_init(D);
    fmpz_poly_init(hmod);
    nmod_poly_init(hmodq, q);
    nmod_poly_init(hdiverse, q);
    fmpz_sparse_init(f);
    fmpz_sparse_init(g);
    fmpz_sparse_init(h);

    fmpz_set_si(D, 50);
    fmpz_sparse_randtest(f, rs, 4, D, 7);
    printf("f is ");
    fmpz_sparse_print_pretty(f, "x");
    printf("\n\n");

    fmpz_set_si(D, 100);
    fmpz_sparse_randtest(g, rs, 3, D, 4);
    printf("g is ");
    fmpz_sparse_print_pretty(g, "x");
    printf("\n\n");

    fmpz_sparse_add(h, f, g);
    printf("f+g = ");
    fmpz_sparse_print_pretty(h, "x");
    printf("\n\n");

    fmpz_sparse_mul_classical(h, f, g);
    printf("f*g = ");
    fmpz_sparse_print_pretty(h, "x");
    printf("\n\n");

    fmpz_sparse_rem_cyc_dense(hmod, h, e);
    printf("f*g mod (x^%ld - 1) = ", e);
    fmpz_poly_print_pretty(hmod, "x");
    printf("\n\n");

    fmpz_sparse_rem_cyc_nmod(hmodq, h, e, q);
    printf("f*g mod (x^%ld - 1) mod %ld = ", e, q);
    nmod_poly_print_pretty(hmodq, "x");
    printf("\n\n");

    fmpz_sparse_rem_cyc_mod_diverse(hdiverse, h, 2, e, q);
    printf("f*g(2x) mod (x^%ld - 1) mod %ld = ", e, q);
    nmod_poly_print_pretty(hdiverse, "x");
    printf("\n\n");

    {
      fmpz_poly_t hdense;
      fmpz_poly_t shift;
      fmpz_poly_init(hdense);
      fmpz_poly_init(shift);
      fmpz_sparse_get_fmpz_poly(hdense, h);
      fmpz_poly_set_coeff_ui(shift, 1, 2);
      printf("shift = "); fmpz_poly_print_pretty(shift, "x"); printf("\n");
      fmpz_poly_compose(hdense, hdense, shift);
      fmpz_sparse_set_fmpz_poly(h, hdense);
      fmpz_sparse_rem_cyc_nmod(hmodq, h, e, q);
      printf("check: "); nmod_poly_print_pretty(hmodq, "x");
      printf("\n\n");
      fmpz_poly_clear(hdense);
      fmpz_poly_clear(shift);
    }

    fmpz_sparse_clear(f);
    fmpz_sparse_clear(g);
    fmpz_sparse_clear(h);
    fmpz_poly_clear(hmod);
    nmod_poly_clear(hmodq);
    nmod_poly_clear(hdiverse);
    fmpz_clear(D);
  }


  { /* test 2 */
    slong i;
    slong m = 5;
    slong t = 4;
    slong d = 100;
    slong logH = 10;
    ulong e = 11;
    ulong q = 101;
    
    fmpz_t dz;
    fmpz_sparse_t prod0, g;
    fmpz_poly_t prod1;
    nmod_poly_t prod2, temp;

    printf("---------------------test 2--------------------------\n");

    printf("Computing the product of %ld polynomials with degree %ld, %ld terms,\n and %ld-bit coeffs\n",
      m, d, t, logH);

    fmpz_init(dz);
    fmpz_sparse_init(prod0);
    fmpz_sparse_init(g);
    fmpz_poly_init(prod1);
    nmod_poly_init(prod2, q);
    nmod_poly_init(temp, q);

    fmpz_set_si(dz, d);
    fmpz_sparse_one(prod0);
    nmod_poly_one(prod2);

    for (i=0; i<m; ++i) {
      /* g = random polynomial */
      fmpz_sparse_randtest(g, rs, t, dz, logH);
      fmpz_sparse_mul(prod0, prod0, g);
      fmpz_sparse_rem_cyc_nmod(temp, g, e, q);
      nmod_poly_mul(temp, prod2, temp);
      nmod_poly_rem_cyc(prod2, temp, e);
    }

    printf("The full product is: ");
    fmpz_sparse_print_pretty(prod0, "x");
    printf("\n\n");

    fmpz_sparse_rem_cyc_dense(prod1, prod0, e);
    printf("Mod (x^%lu - 1) it equals: ", e);
    fmpz_poly_print_pretty(prod1, "x");
    printf("\n\n");

    printf("Also coeffs mod %lu it equals: ", q);
    nmod_poly_print_pretty(prod2, "x");
    printf("\n\n");

    fmpz_sparse_rem_cyc_nmod(temp, prod0, e, q);
    if (nmod_poly_equal(temp, prod2)) {
      printf("Mod check passed!\n\n");
    }
    else {
      printf("CHECK FAILED. The last one should be: ");
      nmod_poly_print_pretty(temp, "x");
      printf("\n\n");
      retval = 1;
    }

    fmpz_clear(dz);
    fmpz_sparse_clear(prod0);
    fmpz_sparse_clear(g);
    fmpz_poly_clear(prod1);
    nmod_poly_clear(prod2);
    nmod_poly_clear(temp);
  }

  { /* test 3 */
    /* computes ab - c for 3 random sparse polys and tests conversion */
    /* and I/O routines */

    slong i;
#define NPOLY (4)
    fmpz_t deg;
    fmpz_t tempz;
    fmpz_sparse_struct sparse[NPOLY];
    fmpz_poly_struct dense[NPOLY];
    const char* var = "variableName";
    char* ts1;
    char* ts2;
    int resp;
    slong dtest;

    printf("\n--------------------------test 3----------------------\n");

    fmpz_init_set_ui(deg, 50);
    fmpz_init(tempz);

    /* the first polys hold the result */
    fmpz_sparse_init(sparse+0);
    fmpz_poly_init(dense+0);

    for (i=1; i<3; ++i) {
      /* middle polys are random sparse polys */
      fmpz_sparse_init(sparse+i);
      fmpz_sparse_randtest(sparse+i, rs, 10, deg, 10);
      fmpz_poly_init(dense+i);
      fmpz_sparse_get_fmpz_poly(dense+i, sparse+i);
    }

    /* last polys are random dense */
    fmpz_sparse_init(sparse+3);
    fmpz_poly_init(dense+3);
    fmpz_poly_randtest_not_zero(dense+3, rs, 30, 20); /* unlucky */
    fmpz_poly_randtest_not_zero(dense+3, rs, 30, 20); /* unlucky */
    fmpz_poly_randtest_not_zero(dense+3, rs, 30, 20);

    fmpz_sparse_set_fmpz_poly(sparse+3, dense+3);

    /* p[0] = p[1] * p[2] */
    fmpz_sparse_mul(sparse+0, sparse+1, sparse+2);
    fmpz_poly_mul(dense+0, dense+1, dense+2);

    /* p[0] -= p[3] */
    fmpz_sparse_sub(sparse+0, sparse+0, sparse+3);
    fmpz_poly_sub(dense+0, dense+0, dense+3);

    fmpz_sparse_set_fmpz_poly(sparse+1, dense+0);
    fmpz_sparse_get_fmpz_poly(dense+1, sparse+0);

    if (fmpz_sparse_equal(sparse+0, sparse+1) 
        && fmpz_poly_equal(dense+0, dense+1)) {
      printf("equality check passed\n");
    }
    else {
      printf("ERROR: equality check failed. add some debugging here.\n");
      retval = 1;
    }

    ts1 = fmpz_sparse_get_str_pretty(sparse+0, var);
    ts2 = fmpz_poly_get_str_pretty(dense+1, var);
    if (strcmp(ts1, ts2) == 0) {
      printf("pretty string check passed\n");
    } else {
      printf("ERROR: pretty string check\n");
      retval = 1;
    }
    free(ts1);
    free(ts2);

    ts1 = fmpz_sparse_get_str(sparse+0);
    resp = fmpz_sparse_set_str(sparse+2, ts1);
    if (resp == 0 && fmpz_sparse_equal(sparse+0, sparse+2)) {
      printf("get/set string check passed\n");
    }
    else if (resp != 0) {
      printf("ERROR on set string, code %d\n", resp);
      retval = 1;
    }
    else {
      printf("ERROR: get/set string check\n");
      retval = 1;
    }
    free(ts1);

    /* A few more little getting/setting tests */
    fmpz_sparse_set(sparse+1, sparse+0);
    fmpz_sparse_sub(sparse+0, sparse+2, sparse+0);
    fmpz_sparse_add(sparse+0, sparse+1, sparse+0);
    fmpz_sparse_sub(sparse+0, sparse+1, sparse+2);
    if (fmpz_sparse_is_zero(sparse+0)) {
      printf("sub to zero check passed\n");
    } else {
      printf("ERROR in sub to zero check\n");
      retval = 1;
    }

    dtest = fmpz_poly_degree(dense+3) / 2;
    fmpz_sparse_set(sparse+0, sparse+3);
    fmpz_set_si(deg, dtest);
    fmpz_sparse_get_coeff(tempz, sparse+0, deg);
    if (fmpz_equal(tempz, fmpz_poly_get_coeff_ptr(dense+3, dtest))) {
      printf("get coeff check passed\n");
    } else {
      printf("ERROR in get coeff\n");
      retval = 1;
    }
    fmpz_set_si(tempz, 1985);
    fmpz_sparse_set_coeff(sparse+3, deg, tempz);
    fmpz_poly_set_coeff_fmpz(dense+3, dtest, tempz);
    fmpz_sparse_set_fmpz_poly(sparse+2, dense+3);

    fmpz_sparse_set(sparse+1, sparse+3);
    fmpz_poly_set(dense+1, dense+3);
    fmpz_zero(tempz);
    fmpz_sparse_set_coeff(sparse+1, deg, tempz);
    fmpz_poly_set_coeff_ui(dense+1, dtest, 0);
    fmpz_sparse_set_fmpz_poly(sparse+0, dense+1);

    if (fmpz_sparse_equal(sparse+3,sparse+2) && fmpz_sparse_equal(sparse+1,sparse+0)) {
      printf("set coeff checks passed\n");
    } else {
      printf("ERROR in set coeff\n");
      retval = 1;
    }


    fmpz_clear(deg);
    fmpz_clear(tempz);
    for (i=0; i<NPOLY; ++i) {
      fmpz_sparse_clear(sparse+i);
      fmpz_poly_clear(dense+i);
    }
  }

  flint_randclear(rs);
  return retval;
}
