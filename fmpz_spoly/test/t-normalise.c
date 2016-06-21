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

  Authored 2016 by Daniel S. Roche; US Government work in the public domain.

******************************************************************************/

#include <stdlib.h>

#include "flint.h"
#include "fmpz_spoly.h"
#include "fmpz_vec.h"
#include "fmpz.h"
#include "ulong_extras.h"

void init_copy(fmpz_spoly_t out, const fmpz_spoly_t in)
{
    fmpz_spoly_init2(out, in->alloc);
    out->length = in->length;
    _fmpz_vec_set(out->coeffs, in->coeffs, in->length);
    _fmpz_vec_set(out->expons, in->expons, in->length);
}

int check_valid(const fmpz_spoly_t a)
{
    slong i;

    if (a->length >= a->alloc)
        return 0;
    if (a->alloc < 1)
        return 0;
    
    for (i = 0; i < a->length; ++i)
    {
        if (fmpz_is_zero(a->coeffs + i))
            return 0;
        if (fmpz_cmp_si(a->expons + i, 0) < 0)
            return 0;
        if (i + 1 < a->length && 
            fmpz_cmp(a->expons + i, a->expons + (i + 1)) <= 0)
            return 0;
    }
    
    return 1;
}

typedef struct {
    slong ind;
    const fmpz* expon;
} eq_check_struct;

int eq_check_fn(const void* a, const void* b)
{
    return fmpz_cmp(((eq_check_struct*) b)->expon,
                    ((eq_check_struct*) a)->expon);
}

int check_equal(const fmpz_spoly_t a, const fmpz_spoly_t b)
{
    slong i, j;
    int pass = 1;
    eq_check_struct* sind;
    fmpz_t csum;

    fmpz_init(csum);
    sind = flint_malloc(b->length * sizeof *sind);

    for (i = 0; i < b->length; ++i)
    {
        sind[i].ind = i;
        sind[i].expon = b->expons + i;
    }

    qsort(sind, b->length, sizeof *sind, eq_check_fn);

    i = j = 0;
    while (j < b->length)
    {
        fmpz_set(csum, b->coeffs + sind[j].ind);

        for (++j; j < b->length && fmpz_equal(sind[j - 1].expon, sind[j].expon); ++j)
            fmpz_add(csum, csum, b->coeffs + sind[j].ind);

        if (!fmpz_is_zero(csum))
        {
            if (i < a->length && fmpz_equal(a->expons + i, sind[j - 1].expon)
                && fmpz_equal(a->coeffs + i, csum))
            {
                ++i;
            }
            else
            {
                pass = 0;
                break;
            }
        }
    }

    if (i < a->length || j < b->length)
        pass = 0;

    fmpz_clear(csum);
    flint_free(sind);

    return pass;
}

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("normalise....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); ++i)
    {
        fmpz_spoly_t a, b;
        slong cbits, terms, dbits;
        slong j;

        cbits = n_randbits(state, n_randint(state, 8)) + 1;
        terms = n_randbits(state, n_randint(state, 10)) + 1;
        dbits = n_randint(state, 8) + 1;

        fmpz_spoly_init(a);
        _fmpz_spoly_reserve(a, terms);

        for (j = 0; j < terms; ++j)
        {
            int choice = n_randint(state, 4);
            if (choice == 0 && j > 0)
            {
                /* duplicate previous term */
                slong prev = n_randint(state, j);
                fmpz_set(a->expons + j, a->expons + prev);
                fmpz_randtest(a->coeffs + j, state, cbits);
            }
            else if (choice == 1 && j > 0)
            {
                /* cancel previous term */
                slong prev = n_randint(state, j);
                fmpz_set(a->expons + j, a->expons + prev);
                fmpz_neg(a->coeffs + j, a->coeffs + prev);
            }
            else if (choice == 2)
            {
                /* make zero term */
                fmpz_randtest(a->expons + j, state, dbits);
                fmpz_abs(a->expons + j, a->expons + j);
                fmpz_zero(a->coeffs + j);
            }
            else
            {
                /* make random term */
                fmpz_randtest(a->expons + j, state, dbits);
                fmpz_abs(a->expons + j, a->expons + j);
                fmpz_randtest(a->coeffs + j, state, cbits);
            }
        }

        _fmpz_spoly_set_length(a, terms);
        
        init_copy(b, a);

        _fmpz_spoly_normalise(a);

        if (!check_valid(a))
        {
            flint_printf("FAIL (part 1, unsorted) %wd %wd %wd:\n", 
                    terms, cbits, dbits);
            fmpz_spoly_print(a), flint_printf("\n\n");
            abort();
        }

        if (!check_equal(a, b))
        {
            flint_printf("FAIL (part 1, unequal) %wd %wd %wd:\n", 
                    terms, cbits, dbits);
            fmpz_spoly_print(a), flint_printf("\n\n");
            fmpz_spoly_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_spoly_clear(a);
        fmpz_spoly_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
