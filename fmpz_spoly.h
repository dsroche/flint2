/*=============================================================================

    fmpz_spoly.h: Sparse univariate Laurent polynomials with 
    fmpz coefficients and fmpz exponents.

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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#ifndef FMPZ_SPOLY_H
#define FMPZ_SPOLY_H

#ifdef FMPZ_SPOLY_INLINES_C
#define FMPZ_SPOLY_INLINE FLINT_DLL
#else
#define FMPZ_SPOLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#include <string.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_vec.h"
#include "nmod_poly.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FMPZ_SPOLY_HEAP_XOVER (32)

/*  Type definitions *********************************************************/

typedef struct
{
    fmpz *coeffs;  /* Nonzero coefficients, in same order as the exponents. */
    fmpz *expons;  /* Term exponents, in reverse sorted order. */
    slong length;  /* Number of nonzero terms. */
    slong alloc;   /* Size that has been allocated for coeffs and expons. */
} fmpz_spoly_struct;

typedef fmpz_spoly_struct fmpz_spoly_t[1];

typedef struct
{
    fmpz_t q;
    ulong log2_order;
    fmpz* points;
    slong length;
} fmpz_spoly_bp_interp_basis_struct;

typedef fmpz_spoly_bp_interp_basis_struct fmpz_spoly_bp_interp_basis_t[1];

typedef struct
{
    fmpz* evals;
    const fmpz_spoly_bp_interp_basis_struct* basis;
} fmpz_spoly_bp_interp_eval_struct;

typedef fmpz_spoly_bp_interp_eval_struct fmpz_spoly_bp_interp_eval_t[1];

typedef struct
{
    nmod_t* cmods;
    ulong* shifts;
    nmod_t* emods;
    slong length;
    slong cimg_per_round, cimg_needed;
    slong eimg_per_round, eimg_needed;
} fmpz_spoly_sp_interp_basis_struct;

typedef fmpz_spoly_sp_interp_basis_struct fmpz_spoly_sp_interp_basis_t[1];

typedef struct
{
    nmod_poly_struct* evals;
    const fmpz_spoly_sp_interp_basis_struct* basis;
} fmpz_spoly_sp_interp_eval_struct;

typedef fmpz_spoly_sp_interp_eval_struct fmpz_spoly_sp_interp_eval_t[1];

extern const fmpz_t FMPZ_SPOLY_NEGATIVE_ONE;

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_spoly_init(fmpz_spoly_t poly);

FLINT_DLL void fmpz_spoly_init2(fmpz_spoly_t poly, slong alloc);

FLINT_DLL void fmpz_spoly_clear(fmpz_spoly_t poly);

FLINT_DLL void _fmpz_spoly_normalise(fmpz_spoly_t poly);

FLINT_DLL void _fmpz_spoly_reserve(fmpz_spoly_t poly, slong terms);

FMPZ_SPOLY_INLINE
void _fmpz_spoly_set_length(fmpz_spoly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; ++i)
        {
            _fmpz_demote(poly->coeffs + i);
            _fmpz_demote(poly->expons + i);
        }
    }
    poly->length = newlen;
}

/*  Polynomial parameters  ***************************************************/

FMPZ_SPOLY_INLINE 
slong fmpz_spoly_terms(const fmpz_spoly_t poly)
{
    return poly->length;
}

FMPZ_SPOLY_INLINE
const fmpz * fmpz_spoly_degree_ptr(const fmpz_spoly_t poly)
{
    if (poly->length == 0) return FMPZ_SPOLY_NEGATIVE_ONE;
    else return poly->expons + 0;
}

FMPZ_SPOLY_INLINE 
void fmpz_spoly_degree(fmpz_t res, const fmpz_spoly_t poly)
{
    if (poly->length > 0) fmpz_set(res, poly->expons + 0);
    else fmpz_set_si(res, -1);
}

FMPZ_SPOLY_INLINE 
slong fmpz_spoly_degree_si(const fmpz_spoly_t poly)
{
    if (poly->length > 0) return fmpz_get_si(poly->expons + 0);
    else return -1;
}

FMPZ_SPOLY_INLINE
const fmpz * fmpz_spoly_lowdeg_ptr(const fmpz_spoly_t poly)
{
    if (poly->length == 0) return FMPZ_SPOLY_NEGATIVE_ONE;
    else return poly->expons + (poly->length - 1);
}

FMPZ_SPOLY_INLINE 
void fmpz_spoly_lowdeg(fmpz_t res, const fmpz_spoly_t poly)
{
    if (poly->length > 0) fmpz_set(res, poly->expons + (poly->length - 1));
    else fmpz_set_si(res, 1);
}

FMPZ_SPOLY_INLINE 
slong fmpz_spoly_lowdeg_si(const fmpz_spoly_t poly)
{
    if (poly->length > 0) return fmpz_get_si(poly->expons + (poly->length - 1));
    else return 1;
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_poly(const fmpz_spoly_t poly)
{
    return poly->length == 0 || fmpz_sgn(poly->expons + (poly->length - 1)) >= 0;
}

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void fmpz_spoly_zero(fmpz_spoly_t poly);

FMPZ_SPOLY_INLINE
void fmpz_spoly_one(fmpz_spoly_t poly)
{
    fmpz_spoly_zero(poly);
    FLINT_ASSERT(poly->alloc >= 1);
    fmpz_init_set_si(poly->coeffs + 0, 1);
    fmpz_init_set_si(poly->expons + 0, 0);
    poly->length = 1;
}

FLINT_DLL void fmpz_spoly_set(fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2);

FMPZ_SPOLY_INLINE
void fmpz_spoly_set_fmpz_fmpz(fmpz_spoly_t poly, 
    const fmpz_t coeff, const fmpz_t expon)
{
    fmpz_spoly_zero(poly);
    if (!fmpz_is_zero(coeff)) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set(poly->coeffs + 0, coeff);
        fmpz_init_set(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_set_fmpz_si(fmpz_spoly_t poly, 
    const fmpz_t coeff, slong expon)
{
    fmpz_spoly_zero(poly);
    if (!fmpz_is_zero(coeff)) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set(poly->coeffs + 0, coeff);
        fmpz_init_set_si(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_set_si_fmpz(fmpz_spoly_t poly, 
    slong coeff, const fmpz_t expon)
{
    fmpz_spoly_zero(poly);
    if (coeff) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set_si(poly->coeffs + 0, coeff);
        fmpz_init_set(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_set_si_si(fmpz_spoly_t poly, 
    slong coeff, slong expon)
{
    fmpz_spoly_zero(poly);
    if (coeff) {
        FLINT_ASSERT(poly->alloc >= 1);
        fmpz_init_set_si(poly->coeffs + 0, coeff);
        fmpz_init_set_si(poly->expons + 0, expon);
        poly->length = 1;
    }
}

FLINT_DLL void fmpz_spoly_set_fmpz_poly(fmpz_spoly_t poly1, 
    const fmpz_poly_t poly2);

FLINT_DLL void fmpz_spoly_get_fmpz_poly(fmpz_poly_t out, 
    const fmpz_spoly_t in);

FLINT_DLL int fmpz_spoly_set_str(fmpz_spoly_t poly, const char * str);

FLINT_DLL char * fmpz_spoly_get_str(const fmpz_spoly_t poly);

FMPZ_SPOLY_INLINE
void fmpz_spoly_swap(fmpz_spoly_t poly1, fmpz_spoly_t poly2)
{
    FLINT_GENERIC_SWAP(fmpz *, poly1->coeffs, poly2->coeffs);
    FLINT_GENERIC_SWAP(fmpz *, poly1->expons, poly2->expons);
    FLINT_GENERIC_SWAP(slong, poly1->length, poly2->length);
    FLINT_GENERIC_SWAP(slong, poly1->alloc, poly2->alloc);
}

FLINT_DLL void fmpz_spoly_truncate(fmpz_spoly_t poly, const fmpz_t deg);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_set_trunc(fmpz_spoly_t res, 
    const fmpz_spoly_t poly, const fmpz_t deg);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_set_trunc_fmpz_poly(fmpz_poly_t res, 
    const fmpz_spoly_t poly, slong deg);

/*  Randomisation  ***********************************************************/

FLINT_DLL void fmpz_spoly_randtest(fmpz_spoly_t res, flint_rand_t state, 
     slong terms, const fmpz_t degree, mp_bitcnt_t bits);

FLINT_DLL void fmpz_spoly_randtest_kron(fmpz_spoly_t res, flint_rand_t state,
     slong terms, const fmpz_t degree, slong limit, mp_bitcnt_t bits, slong vars);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_randtest_unsigned(fmpz_spoly_t res, 
    flint_rand_t state, slong terms, const fmpz_t degree, mp_bitcnt_t bits);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_randtest_not_zero(fmpz_spoly_t res, 
    flint_rand_t state, slong terms, const fmpz_t degree, mp_bitcnt_t bits);

/*  Getting and setting coefficients  ****************************************/

FLINT_DLL slong _fmpz_spoly_index(const fmpz_spoly_t poly, const fmpz_t e);

FLINT_DLL slong _fmpz_spoly_index_si(const fmpz_spoly_t poly, slong e);

FLINT_DLL slong fmpz_spoly_get_coeff_si_si(const fmpz_spoly_t poly, 
    slong e);

FLINT_DLL slong fmpz_spoly_get_coeff_si_fmpz(const fmpz_spoly_t poly, 
    const fmpz_t e);

FLINT_DLL void fmpz_spoly_get_coeff_fmpz_si(fmpz_t res,
    const fmpz_spoly_t poly, slong e);

FLINT_DLL void fmpz_spoly_get_coeff(fmpz_t res,
    const fmpz_spoly_t poly, const fmpz_t e);

FLINT_DLL void fmpz_spoly_set_coeff_si_si(fmpz_spoly_t poly, slong c, slong e);

FLINT_DLL void fmpz_spoly_set_coeff_si_fmpz(fmpz_spoly_t poly, 
    slong c, const fmpz_t e);

FLINT_DLL void fmpz_spoly_set_coeff_fmpz_si(fmpz_spoly_t poly, 
    const fmpz_t c, slong e);

FLINT_DLL void fmpz_spoly_set_coeff(fmpz_spoly_t poly, 
    const fmpz_t c, const fmpz_t e);

FMPZ_SPOLY_INLINE
fmpz* fmpz_spoly_get_coeff_ptr(fmpz_spoly_t poly, const fmpz_t e)
{
    slong ind = _fmpz_spoly_index(poly, e);
    if (ind < 0) return NULL;
    else return poly->coeffs + ind;
}

FMPZ_SPOLY_INLINE
fmpz* fmpz_spoly_get_coeff_ptr_si(fmpz_spoly_t poly, slong e)
{
    slong ind = _fmpz_spoly_index_si(poly, e);
    if (ind < 0) return NULL;
    else return poly->coeffs + ind;
}

FMPZ_SPOLY_INLINE 
void fmpz_spoly_get_term(fmpz_t coeff, fmpz_t expon, 
    const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(coeff, poly->coeffs + i);
    fmpz_set(expon, poly->expons + i);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_get_term_coeff(fmpz_t res, const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(res, poly->coeffs + i);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_get_term_expon(fmpz_t res, const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    fmpz_set(res, poly->expons + i);
}

FMPZ_SPOLY_INLINE
fmpz* fmpz_spoly_get_term_coeff_ptr(const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return poly->coeffs + i;
}

FMPZ_SPOLY_INLINE
fmpz* fmpz_spoly_get_term_expon_ptr(const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return poly->expons + i;
}

FMPZ_SPOLY_INLINE
slong fmpz_spoly_get_term_coeff_si(const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return fmpz_get_si(poly->coeffs + i);
}

FMPZ_SPOLY_INLINE
slong fmpz_spoly_get_term_expon_si(const fmpz_spoly_t poly, slong i)
{
    FLINT_ASSERT(i >= 0 && i < poly->length);
    return fmpz_get_si(poly->expons + i);
}

/*  Comparison  **************************************************************/

FMPZ_SPOLY_INLINE 
int fmpz_spoly_equal(const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    return
        poly1->length == poly2->length &&
        _fmpz_vec_equal(poly1->expons, poly2->expons, poly1->length) &&
        _fmpz_vec_equal(poly1->coeffs, poly2->coeffs, poly1->length);
}

/* FIXME (not yet implemented) */
FLINT_DLL int fmpz_spoly_equal_trunc(const fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2, const fmpz_t n);

FMPZ_SPOLY_INLINE 
int fmpz_spoly_is_zero(const fmpz_spoly_t poly)
{
    return poly->length == 0;
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_unit(const fmpz_spoly_t poly)
{
    return poly->length == 1 && 
        fmpz_is_zero(poly->expons + 0) && 
        fmpz_is_pm1(poly->coeffs + 0);
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_term(const fmpz_spoly_t poly, 
    const fmpz_t c, const fmpz_t e)
{
    return (poly->length == 0 && fmpz_is_zero(c)) ||
        (poly->length == 1 && 
         fmpz_equal(poly->coeffs + 0, c) && 
         fmpz_equal(poly->expons + 0, e)
        );
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_term_fmpz_si(const fmpz_spoly_t poly, 
    const fmpz_t c, slong e)
{
    return (poly->length == 0 && fmpz_is_zero(c)) ||
        (poly->length == 1 && 
         fmpz_equal(poly->coeffs + 0, c) && 
         fmpz_equal_si(poly->expons + 0, e)
        );
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_term_si_fmpz(const fmpz_spoly_t poly, 
    slong c, const fmpz_t e)
{
    return (poly->length == 0 && (c == 0)) ||
        (poly->length == 1 && 
         fmpz_equal_si(poly->coeffs + 0, c) && 
         fmpz_equal(poly->expons + 0, e)
        );
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_is_term_si_si(const fmpz_spoly_t poly, slong c, slong e)
{
    return (poly->length == 0 && (c == 0)) ||
        (poly->length == 1 && 
         fmpz_equal_si(poly->coeffs + 0, c) && 
         fmpz_equal_si(poly->expons + 0, e)
        );
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_equal_fmpz(const fmpz_spoly_t poly, const fmpz_t c)
{
    return fmpz_spoly_is_term_fmpz_si(poly, c, 0);
}

/* FIXME (not yet implemented) */
FLINT_DLL int fmpz_spoly_equal_fmpz_poly(const fmpz_spoly_t spoly, 
    const fmpz_poly_t dpoly);

/*  Addition and subtraction  ************************************************/

FLINT_DLL void fmpz_spoly_add(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void fmpz_spoly_sub(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void fmpz_spoly_neg(fmpz_spoly_t res, const fmpz_spoly_t poly);

/*  Scalar multiplication and division  **************************************/

FLINT_DLL void fmpz_spoly_scalar_mul_ui(fmpz_spoly_t res,
        const fmpz_spoly_t poly, ulong c);

FLINT_DLL void fmpz_spoly_scalar_mul_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c);

FLINT_DLL void fmpz_spoly_scalar_mul(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_addmul(fmpz_spoly_t poly1,
    const fmpz_spoly_t poly2, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_submul(fmpz_spoly_t poly1,
    const fmpz_spoly_t poly2, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_fdiv_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_fdiv(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_tdiv_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_tdiv(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_mul_2exp(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong exp);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_fdiv_2exp(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong exp);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_tdiv_2exp(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong exp);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_mod(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_scalar_smod(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c);

/*  Bit packing  *************************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_bit_pack(fmpz_t res,
    const fmpz_spoly_t poly, mp_bitcnt_t bit_size);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_bit_unpack(fmpz_spoly_t res,
    const fmpz_t f, mp_bitcnt_t bit_size);

/*  Multiplication  **********************************************************/

FLINT_DLL void fmpz_spoly_new_mul_classical(fmpz_spoly_t res,
        const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void fmpz_spoly_mul_heaps(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void fmpz_spoly_mul_OS(fmpz_spoly_t res, flint_rand_t state, 
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void fmpz_spoly_mul_interp(fmpz_spoly_t res, flint_rand_t state, 
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2, slong terms);

FLINT_DLL void fmpz_spoly_mul_classical(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FMPZ_SPOLY_INLINE 
void fmpz_spoly_mul(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    if (poly1->length + poly2->length < FMPZ_SPOLY_HEAP_XOVER)
    {
        fmpz_spoly_mul_classical(res, poly1, poly2);
    }
    else fmpz_spoly_mul_heaps(res, poly1, poly2);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_sqr(fmpz_spoly_t res, const fmpz_spoly_t poly)
{
    fmpz_spoly_mul(res, poly, poly);
}

FLINT_DLL void _fmpz_spoly_mul_coeffs(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void _fmpz_spoly_mul_coeffs_slow(fmpz_spoly_t res, flint_rand_t state, 
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2, const fmpz * expons,
    slong len);

/*  Powering  ****************************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_pow_recurrence(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_pow_binomial(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_pow_binexp(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_pow(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_pow_trunc(fmpz_spoly_t res,
    const fmpz_spoly_t poly, ulong e, const fmpz_t n);

/*  Shifting  ****************************************************************/

FLINT_DLL void fmpz_spoly_shift_left(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t n);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_shift_left_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong n);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_shift_right(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t n);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_shift_right_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong n);

FLINT_DLL void _fmpz_spoly_vec_shift(fmpz_spoly_t poly, 
    slong start, slong end, slong dist);

FLINT_DLL void _fmpz_spoly_append_si(fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2, slong c, ulong e);

FLINT_DLL void _fmpz_spoly_append(fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2, const fmpz_t c, const fmpz_t e);

/*  Monomial multiplication and division *************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_mul_si_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c, slong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_mul_si_fmpz(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c, const fmpz_t e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_mul_fmpz_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c, slong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_mul_fmpz_fmpz(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c, const fmpz_t e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_fdiv_si_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c, slong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_fdiv_si_fmpz(fmpz_spoly_t res,
    const fmpz_spoly_t poly, slong c, const fmpz_t e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_fdiv_fmpz_si(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c, slong e);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_mon_fdiv_fmpz_fmpz(fmpz_spoly_t res,
    const fmpz_spoly_t poly, const fmpz_t c, const fmpz_t e);

/*  Norms  *******************************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_2norm(fmpz_t res, const fmpz_spoly_t poly);

FMPZ_SPOLY_INLINE
ulong fmpz_spoly_max_limbs(const fmpz_spoly_t poly)
{
    return _fmpz_vec_max_limbs(poly->coeffs, poly->length);
}

FMPZ_SPOLY_INLINE
slong fmpz_spoly_max_bits(const fmpz_spoly_t poly)
{
    return _fmpz_vec_max_bits(poly->coeffs, poly->length);
}

FMPZ_SPOLY_INLINE
slong fmpz_spoly_max_ebits(const fmpz_spoly_t poly)
{
    fmpz * lead = poly->expons + 0;
    fmpz * trail = poly->expons + (poly->length - 1);
    if (poly->length == 0) return 0;
    else if (poly->length == 1 || fmpz_cmpabs(lead, trail) >= 0)
    {
        return fmpz_bits(lead);
    }
    else return fmpz_bits(trail);
}

FMPZ_SPOLY_INLINE
ulong fmpz_spoly_max_elimbs(const fmpz_spoly_t poly)
{
    slong bc = fmpz_spoly_max_ebits(poly);
    if (bc >= 0) return (bc + FLINT_BITS - 1) / FLINT_BITS;
    else return (bc - FLINT_BITS + 1) / FLINT_BITS;
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_height(fmpz_t res, const fmpz_spoly_t poly)
{
    _fmpz_vec_height(res, poly->coeffs, poly->length);
}

FMPZ_SPOLY_INLINE
mp_bitcnt_t fmpz_spoly_height_bits(const fmpz_spoly_t poly)
{
    return FLINT_ABS(_fmpz_vec_max_bits(poly->coeffs, poly->length));
}

/*  Euclidean division  ******************************************************/

FMPZ_SPOLY_INLINE
void fmpz_spoly_divrem(fmpz_spoly_t Q, fmpz_spoly_t R,
    const fmpz_spoly_t A, const fmpz_spoly_t B)
{
    FLINT_ASSERT(0);
    /* FIXME this is just a placeholder! */
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_div(fmpz_spoly_t Q, 
    const fmpz_spoly_t A, const fmpz_spoly_t B)
{
    /* TODO better efficiency */
    fmpz_spoly_t temp;
    fmpz_spoly_init(temp);
    fmpz_spoly_divrem(Q, temp, A, B);
    fmpz_spoly_clear(temp);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_rem(fmpz_spoly_t R, fmpz_spoly_t A, fmpz_spoly_t B)
{
    /* TODO better efficiency */
    fmpz_spoly_t temp;
    fmpz_spoly_init(temp);
    fmpz_spoly_divrem(temp, R, A, B);
    fmpz_spoly_clear(temp);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_divrem_dense(fmpz_spoly_t Q, fmpz_poly_t R,
    const fmpz_spoly_t A, const fmpz_poly_t B)
{
    /* FIXME this is just a placeholder! */
    FLINT_ASSERT(0);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_div_dense(fmpz_spoly_t Q, fmpz_spoly_t A, fmpz_poly_t B)
{
    /* TODO better efficiency */
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_spoly_divrem_dense(Q, temp, A, B);
    fmpz_poly_clear(temp);
}

/* FIXME */
FLINT_DLL void fmpz_spoly_rem_dense(fmpz_poly_t R, fmpz_spoly_t A, fmpz_poly_t B);

FMPZ_SPOLY_INLINE 
void fmpz_spoly_rem_cyc(fmpz_spoly_t res, const fmpz_spoly_t poly, const fmpz_t e)
{
  if(poly == res)
  {
    fmpz_spoly_t temp;
    fmpz_spoly_init2(temp, poly->length);
    fmpz_spoly_rem_cyc(temp, poly, e);
    fmpz_spoly_swap(res, temp);
    fmpz_spoly_clear(temp);
  }
  else
  {
    fmpz_spoly_zero(res);
    _fmpz_spoly_reserve(res, poly->length);
    _fmpz_vec_scalar_mod_fmpz(res->expons, poly->expons, poly->length, e);
    _fmpz_vec_set(res->coeffs, poly->coeffs, poly->length);
    res->length = poly->length;
    _fmpz_spoly_normalise(res);
  }
}

FLINT_DLL void fmpz_spoly_rem_cyc_dense(fmpz_poly_t res,
    const fmpz_spoly_t poly, ulong e);

FLINT_DLL void fmpz_spoly_rem_cyc_nmod(nmod_poly_t res,
    const fmpz_spoly_t poly, ulong e, ulong q);

FLINT_DLL void fmpz_spoly_rem_cyc_mod_diverse(nmod_poly_t res,
    const fmpz_spoly_t poly, ulong a, ulong e, ulong q);

/*  Greatest common divisor  *************************************************/

FMPZ_SPOLY_INLINE
void fmpz_spoly_xgcd(fmpz_spoly_t r, 
    fmpz_spoly_t s, fmpz_spoly_t t,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    /* FIXME this is just a placeholder! */
    FLINT_ASSERT(0);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_gcd(fmpz_spoly_t res, 
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    /* TODO better efficiency */
    fmpz_spoly_t s, t;
    fmpz_spoly_init(s);
    fmpz_spoly_init(t);
    fmpz_spoly_xgcd(res, s, t, poly1, poly2);
    fmpz_spoly_clear(s);
    fmpz_spoly_clear(t);
}

FMPZ_SPOLY_INLINE 
void fmpz_spoly_lcm(fmpz_spoly_t res, 
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2)
{
    /* TODO better efficiency */
    fmpz_spoly_gcd(res, poly1, poly2);
    fmpz_spoly_div(res, poly1, res);
    fmpz_spoly_mul(res, res, poly2);
}

/*  Gaussian content  ********************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_content(fmpz_t res, const fmpz_spoly_t poly);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_primitive_part(fmpz_spoly_t res, const fmpz_spoly_t poly);

/*  Sparse interpolation ****************************************************/

FLINT_DLL void fmpz_spoly_bp_interp_basis_init(fmpz_spoly_bp_interp_basis_t res, 
    flint_rand_t state, slong terms, mp_bitcnt_t d, mp_bitcnt_t h);

FMPZ_SPOLY_INLINE
void fmpz_spoly_bp_interp_basis_clear(fmpz_spoly_bp_interp_basis_t res)
{
    fmpz_clear(res->q);
    if (res->length) _fmpz_vec_clear(res->points, res->length);
}


FMPZ_SPOLY_INLINE
void fmpz_spoly_bp_interp_eval_init(fmpz_spoly_bp_interp_eval_t res,
    const fmpz_spoly_bp_interp_basis_t basis)
{
    res->basis = basis;
    if (res->basis->length) res->evals = _fmpz_vec_init(res->basis->length);
    else res->evals = NULL;
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_bp_interp_eval_clear(fmpz_spoly_bp_interp_eval_t res)
{
    if (res->basis->length) _fmpz_vec_clear(res->evals, res->basis->length);
}

/* (redundant forward declaration) */
FLINT_DLL void fmpz_spoly_evaluate_powers(fmpz* res, slong len,
    const fmpz_spoly_t poly, const fmpz_t w, const fmpz_t p);

FMPZ_SPOLY_INLINE
void fmpz_spoly_bp_interp_eval(fmpz_spoly_bp_interp_eval_t res, const fmpz_spoly_t poly)
{
    fmpz_spoly_evaluate_powers(res->evals,
        res->basis->length, poly, res->basis->points + 1, res->basis->q);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_bp_interp_eval_set(fmpz_spoly_bp_interp_eval_t res,
    const fmpz_spoly_bp_interp_eval_t eval)
{
    FLINT_ASSERT(res->basis->length == eval->basis->length);
    _fmpz_vec_set(res->evals, eval->evals, res->basis->length);
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_bp_interp_eval_is_zero(const fmpz_spoly_bp_interp_eval_t eval)
{
    return _fmpz_vec_is_zero(eval->evals, eval->basis->length);
}

FLINT_DLL void fmpz_spoly_bp_interp_mul(fmpz_spoly_bp_interp_eval_t res,
    const fmpz_spoly_bp_interp_eval_t poly1,
    const fmpz_spoly_bp_interp_eval_t poly2);

FLINT_DLL void fmpz_spoly_bp_interp_addmul_si(fmpz_spoly_bp_interp_eval_t res,
    slong c,
    const fmpz_spoly_bp_interp_eval_t poly2);

FLINT_DLL void fmpz_spoly_bp_interp_pow(fmpz_spoly_bp_interp_eval_t res,
    const fmpz_spoly_bp_interp_eval_t poly, ulong pow);

FLINT_DLL int fmpz_spoly_bp_interp(fmpz_spoly_t res,
    const fmpz_spoly_bp_interp_eval_t eval);

FLINT_DLL void fmpz_spoly_sp_interp_basis_init(fmpz_spoly_sp_interp_basis_t res, 
    flint_rand_t state, slong terms, mp_bitcnt_t d, mp_bitcnt_t h);

FMPZ_SPOLY_INLINE
void fmpz_spoly_sp_interp_basis_clear(fmpz_spoly_sp_interp_basis_t res)
{
    if (res->length)
    {
        flint_free(res->cmods);
        flint_free(res->shifts);
        flint_free(res->emods);
    }
}

FLINT_DLL void fmpz_spoly_sp_interp_eval_init(fmpz_spoly_sp_interp_eval_t res,
    const fmpz_spoly_sp_interp_basis_t basis);

FLINT_DLL void fmpz_spoly_sp_interp_eval_clear(fmpz_spoly_sp_interp_eval_t res);

FLINT_DLL void fmpz_spoly_sp_interp_eval(fmpz_spoly_sp_interp_eval_t res, 
    const fmpz_spoly_t poly);

FMPZ_SPOLY_INLINE
void fmpz_spoly_sp_interp_eval_set(fmpz_spoly_sp_interp_eval_t res,
        const fmpz_spoly_sp_interp_eval_t eval)
{
    slong i;
    FLINT_ASSERT(res->basis == eval->basis);
    for (i = 0; i < res->basis->length; ++i)
    {
        nmod_poly_set(res->evals + i, eval->evals + i);
    }
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_sp_interp_eval_is_zero(const fmpz_spoly_sp_interp_eval_t eval)
{
    slong i;
    for (i = 0; i < eval->basis->length; ++i)
    {
        if (! nmod_poly_is_zero(eval->evals + i)) return 0;
    }
    return 1;
}

FLINT_DLL void fmpz_spoly_sp_interp_mul(fmpz_spoly_sp_interp_eval_t res,
    const fmpz_spoly_sp_interp_eval_t poly1,
    const fmpz_spoly_sp_interp_eval_t poly2);

FLINT_DLL void fmpz_spoly_sp_interp_addmul_si(fmpz_spoly_sp_interp_eval_t res,
    slong c,
    const fmpz_spoly_sp_interp_eval_t poly2);

FLINT_DLL void fmpz_spoly_sp_interp_pow(fmpz_spoly_sp_interp_eval_t res,
    const fmpz_spoly_sp_interp_eval_t poly, ulong pow);

FLINT_DLL int fmpz_spoly_sp_interp(fmpz_spoly_t res,
    const fmpz_spoly_sp_interp_eval_t eval);

/*  Divisibility testing  ***************************************************/

FMPZ_SPOLY_INLINE
int fmpz_spoly_divides(fmpz_spoly_t q, 
    const fmpz_spoly_t a, const fmpz_spoly_t b)
{
    /* TODO better efficiency */
    int res;
    fmpz_spoly_t temp;
    fmpz_spoly_init(temp);
    fmpz_spoly_divrem(q, temp, a, b);
    res = fmpz_spoly_is_zero(temp);
    fmpz_spoly_clear(temp);
    return res;
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_divides_dense(fmpz_spoly_t q, 
    const fmpz_spoly_t a, const fmpz_poly_t b)
{
    /* TODO better efficiency */
    int res;
    fmpz_poly_t temp;
    fmpz_poly_init(temp);
    fmpz_spoly_divrem_dense(q, temp, a, b);
    res = fmpz_poly_is_zero(temp);
    fmpz_poly_clear(temp);
    return res;
}

/*  Derivative  **************************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_derivative(fmpz_spoly_t res, const fmpz_spoly_t poly);

/*  Evaluation  **************************************************************/

FLINT_DLL void fmpz_spoly_evaluate(fmpz_t res,
    const fmpz_spoly_t f, const fmpz_t a);

FLINT_DLL void fmpz_spoly_evaluate_mod(fmpz_t res, const fmpz_spoly_t poly, 
    const fmpz_t a, const fmpz_t m);

FLINT_DLL ulong fmpz_spoly_evaluate_mod_ui(const fmpz_spoly_t poly, 
    ulong a, ulong m);

FLINT_DLL void fmpz_spoly_evaluate_powers(fmpz* res, slong len,
    const fmpz_spoly_t poly, const fmpz_t w, const fmpz_t p);

/*  Composition  *************************************************************/

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_compose(fmpz_spoly_t res,
    const fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

/* FIXME (not yet implemented) */
FLINT_DLL void fmpz_spoly_compose_dense(fmpz_spoly_t res,
    const fmpz_poly_t poly1, const fmpz_spoly_t poly2);

/*  Input and output  ********************************************************/

FLINT_DLL int fmpz_spoly_fprint(FILE * file, const fmpz_spoly_t poly);

FLINT_DLL int fmpz_spoly_fprint_pretty(FILE * file, 
    const fmpz_spoly_t poly, const char *x);

FMPZ_SPOLY_INLINE
int fmpz_spoly_print(const fmpz_spoly_t poly)
{
    return fmpz_spoly_fprint(stdout, poly);
}

FMPZ_SPOLY_INLINE
int fmpz_spoly_print_pretty(const fmpz_spoly_t poly, const char *x)
{
    return fmpz_spoly_fprint_pretty(stdout, poly, x);
}

FLINT_DLL int fmpz_spoly_fread(FILE * file, fmpz_spoly_t poly);

FMPZ_SPOLY_INLINE
int fmpz_spoly_read(fmpz_spoly_t poly)
{
    return fmpz_spoly_fread(stdin, poly);
}

FMPZ_SPOLY_INLINE
void fmpz_spoly_debug(const fmpz_spoly_t poly)
{
    flint_printf("(alloc = %wd, length = %wd,\n  coeffs = ", poly->alloc, poly->length);
    if (poly->coeffs) _fmpz_vec_print(poly->coeffs, poly->alloc);
    else flint_printf("NULL");
    flint_printf("\n  expons = ");
    if (poly->expons) _fmpz_vec_print(poly->expons, poly->alloc);
    else flint_printf("NULL");
    flint_printf("\n");
    fflush(stdout);
}

/*  OTHER  ******************************************************************/

FLINT_DLL void fmpz_van_prime(fmpz_t res, flint_rand_t state, slong support, 
    mp_bitcnt_t deg_bits, double gamma);

FLINT_DLL void fmpz_diff_prime(fmpz_t res, flint_rand_t state, slong support, 
    mp_bitcnt_t deg_bits, double gamma);

FLINT_DLL slong fmpz_spoly_sumcheck(fmpz ** res, const fmpz_spoly_t poly1, 
    const fmpz_spoly_t poly2);

FLINT_DLL slong fmpz_spoly_sumset(fmpz ** res, flint_rand_t state, const 
    fmpz_spoly_t poly1, const fmpz_spoly_t poly2);

FLINT_DLL void _fmpz_mod_poly_powmod_x_2exp(fmpz* res,
    const fmpz* poly, slong len, ulong k, const fmpz_t p);

FLINT_DLL slong _fmpz_mod_poly_binary_roots(fmpz* roots, fmpz* expons,
    const fmpz* poly, slong len, const fmpz_t theta, slong k, const fmpz_t p);

FLINT_DLL void _fmpz_spoly_transp_vandermonde_precomp(fmpz* bb, slong blen,
    const fmpz* vv_inv, fmpz_poly_struct * const * tree, const fmpz* tree_root,
    const fmpz* xx, slong len, const fmpz_t p);

FLINT_DLL void _fmpz_spoly_transp_vandermonde(fmpz* bb, slong blen,
    const fmpz* vv, const fmpz* xx, slong len, const fmpz_t p);

FLINT_DLL void _fmpz_spoly_transp_vandermonde_inv_precomp(fmpz* xx,
    fmpz_poly_struct * const * tree, const fmpz* tree_root,
    const fmpz* bb, slong len, const fmpz_t p);

FLINT_DLL void _fmpz_spoly_transp_vandermonde_inv(fmpz* xx,
    const fmpz* vv, const fmpz* bb, slong len, const fmpz_t p);

/* TODO: CRT would be nice, but we would need a fmpz_spoly_nmod type first. */

/* TODO: Finding all integer roots of an fmpz_spoly polynomial */

/* TODO: Computing low-degree and cyclotomic factors of an fmpz_spoly polynomial */

#ifdef __cplusplus
}
#endif

#endif /* FMPZ_SPOLY_H */
