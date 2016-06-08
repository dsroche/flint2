#include "flint.h"
#include "math.h"
#include "ulong_extras.h"

ulong bitprimes1(slong numbits)
{
    ulong lower = UWORD(1) << (numbits - 1);
    ulong upper = 2 * (lower - 1) + 1;
    return n_prime_pi(upper) - n_prime_pi(lower);
}

ulong bitprimes2(slong numbits)
{
    if (numbits < 4) return 0;
    double lowerx = (UWORD(1) << (numbits - 1)) - 1;
    double upperx = (lowerx + 1.0) * 2.0;
    ulong lowerb = 1.25506 * lowerx / log(lowerx);
    ulong upperb = ((ulong) (upperx / log(upperx))) + 1;
    return upperb - lowerb;
}

ulong bitprimes3(slong numbits)
{
    if (numbits < 5) return 0;
    double x = UWORD(1) << (numbits - 1);
    return 3.0 * x / (5.0 * log(x));
}

ulong bitprimes4(slong numbits)
{
    ulong n = ((UWORD(1) << (numbits - 1)) - 1) * 2 + 1;
    ulong lo, hi;
    n_prime_pi_bounds(&lo, &hi, n);
    n = n / 2 + 1;
    n_prime_pi_bounds(&hi, &hi, n);
    if (hi > lo) return 0;
    else return lo - hi;
}

int main(int argc, char const * const * argv)
{
    slong bits;
    ulong total = 0;
    slong stopat = 0;
    ulong cur = 1;

    if (argc == 2) stopat = atoi(argv[1]);

    if (stopat < 2)
    {
        flint_printf("usage: %s max_bits\n", argv[0]);
        return 1;
    }
    else if (stopat > FLINT_BITS)
    {
        flint_printf("ERROR: can't go past %d bits\n", FLINT_BITS);
        return 2;
    }

    flint_printf("slong _n_prime_bits_table[%wd] = {0, 0", stopat + 1);

    for (bits = 2; bits <= stopat; ++bits)
    {
        ulong next_total;
        cur = cur * 2 + 1;
        next_total = n_prime_pi(cur);
        flint_printf(",");
        if (bits % 8 == 0) flint_printf("\n    ");
        else flint_printf(" ");
        flint_printf("%wu", next_total - total);
        total = next_total;
    }
    flint_printf("};\n");
    flint_printf("slong _n_prime_bits_tlen = %wd;\n", bits);

    return 0;
}
