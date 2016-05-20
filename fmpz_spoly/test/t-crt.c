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

    Copyright (C) 2010, 2012 William Hart
    Copyright (C) 2011, 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"

int
main(void)
{
  
    int i, result = 1;
    FLINT_TEST_INIT(state);
    
    
    flint_printf("crt from multiple mods....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
      fmpz_t res1, res2, a, b, c, total, temp;

      fmpz_init(a);
      fmpz_init(b);
      fmpz_init(total);
      fmpz_init(temp);
      fmpz_init(res1);
      fmpz_init(res2);

      fmpz_one(total);
      fmpz_set_ui(res2, 1000);

      fmpz_set_ui(a,1000);
      fmpz_mod_ui(c, a, 23);
      fmpz_mod_ui(b, a, 17);
      fmpz_mod_ui(a, a, 31);

      fmpz_set_ui(temp, 23);
      fmpz_CRT(res1, res1, total, c, temp, 0);
      fmpz_mul(total, total, temp);
      
      fmpz_set_ui(temp, 17);
      fmpz_CRT(res1, res1, total, b, temp, 0);
      fmpz_mul(total, total, temp);
      
      fmpz_set_ui(temp, 31);
      fmpz_CRT(res1, res1, total, a, temp, 0);
      fmpz_mul(total, total, temp);

      result = fmpz_equal(res1, res2);

      if (!result)
      {
        flint_printf("FAIL:\n");
        flint_printf("\nres1: "), fmpz_print(res1);
        flint_printf("\nres2: "), fmpz_print(res1);
        flint_printf("\n");
        abort();
      }
      
      fmpz_clear(a);
      fmpz_clear(b);
      fmpz_clear(c);
      fmpz_clear(temp);
      fmpz_clear(total);
      fmpz_clear(res1);
      fmpz_clear(res2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
