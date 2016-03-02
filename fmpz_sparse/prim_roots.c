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

#include "fmpz_sparse.h"
#include "fmpz.h"

slong _fmpz_sparse_prim_roots(fmpz_t p, fmpz * qq, fmpz * ww, flint_rand_t state,
    slong len, mp_bitcnt_t p_bits, mp_bitcnt_t q_prod_bits)
{
    /* This function computes single random prime p, with the proper number
     * of bits, and corresponding primes qq, as well as an array
     * ww of elements in GF(q).
     * For each index i, it holds that p divides (qq[i]-1), and
     * that ww[i] is a p-th primitive root mod qq[i].
     * Furthermore, the product of primes in qq has at least q_prod_bits bits.
     * All arrays must be allocated by the caller, and the fmpz values
     * initialized beforehand. 
     * The len parameter is the length of qq and ww.
     * The returned value gives the actual number of primes q computed, 
     * which will be at most len.
     * In case len is too small to compute sufficiently many primes, -1 is
     * returned.
     */
    slong count = 0;
    fmpz_t qprod;
    ulong a = 0;
    mp_bitcnt_t qpbits;

    fmpz_randprime(p, state, p_bits, 0);
    fmpz_init_set_si(qprod, 1);
    fmpz_set_si(qq+0, 1);

    while (count < len && fmpz_bits(qprod) < q_prod_bits)
    {
        /* qq[count] = a*p + 1 and a is even */
        fmpz_addmul_ui(qq+count, p, 2); 
        a += 2;

        if (fmpz_is_prime(qq+count))
        {
            do
            {
                fmpz_randm(ww+count, state, qq+count);
                fmpz_powm_ui(ww+count, ww+count, a, qq+count);
            } while (fmpz_cmp_si(ww+count, 1) <= 0);
            fmpz_mul(qprod, qprod, qq+count);
            if (++count < len) {
                fmpz_set(qq+count, qq+(count-1));
            }
        }
    }

    qpbits = fmpz_bits(qprod);
    fmpz_clear(qprod);

    if (qpbits >= q_prod_bits) return count;
    else return -1;
}
