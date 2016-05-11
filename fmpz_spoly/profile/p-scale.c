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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   bits     Fixed bitsize of the coefficients
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Step size for the length
   deglo    Minimum degree
   deghi    Maximum degree
   degh     Step size for the degree
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define nalgs    2
#define img      1
#define imgname  "scale.ppm"

/*
   Write a binary 24-bit ppm image.
 */
int write_rgb_ppm(const char* file_name, unsigned char* pixels, 
                   unsigned int width, unsigned int height)
{
    FILE* file = fopen(file_name, "wb");
    if (file == NULL)
        return -1;
    flint_fprintf(file, "P6\n%d %d\n255\n", width, height);
    fwrite(pixels, sizeof(unsigned char), width * height * 3, file);
    fclose(file);
    return 0;
}

int
main(void)
{
    int i, j, len, deg, deghi;
    double T[200][1][nalgs];

    FLINT_TEST_INIT(state);
    
       
    deg = 0.0;
    deghi = 1000.0;
    for (j = 0; deg <= deghi; j++)
    {
      T[j][0][1] = deg;
      T[j][0][0] = (1000-deg);  
      deg += 5.0;
    }
    
    /*
       Print 2-D coloured image to file
     */
    if (img)
    {
        unsigned char * PIXELS;
        int k;
        
        PIXELS = (unsigned char *) flint_malloc(3 * 200 * 1 * sizeof(unsigned char));
        k = 0;
        for (i = 0; i < 200; i++)
        {
                double max = DBL_MIN, v[nalgs];
                int m;
                
                for (m = 0; m < nalgs; m++)
                {
                  v[m] = T[i][0][m];
                  if (v[m] > max)
                    max = v[m];
                }
                for (m = 0; m < nalgs; m++)
                {
                  v[m] = v[m] / max;
                  PIXELS[k++] = (unsigned char) (v[m] * 255);
                }
                for (; m < 3; m++)
                  PIXELS[k++] = (unsigned char) 0;    
        }

        k = write_rgb_ppm(imgname, PIXELS, 2, 200);
        flint_free(PIXELS);
        if (k)
        {
            flint_printf("Exception:  writing ppm image failed\n");
        }
    }

    flint_randclear(state);
    return 0;
}
