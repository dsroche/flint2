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

    Authored 2016 by A. Whitman Groves; US Government work in the public domain

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_spoly.h"
#include "ulong_extras.h"
#include "profiler.h"

/*
   Definitions for the parameters of the timing process.
   
   lenlo    Minimum length
   lenhi    Maximum length
   lenh     Step size for the length
   bitslo   Minimum bit size
   bitshi   Maximum bit size
   bitsh    Step size for the bit size
   cols     Number of different lengths
   rows     Number of different bit sizes
   cpumin   Minimum number of ms spent on each test case
   ncases   Number of test cases per point (length, bit size)
   nalgs    Number of algorithms
   img      Whether an RGB coloured image should be produced
   imgname  File name for image
 */

#define bits     512
#define lenlo    60
#define lenhi    98
#define lenh     1
#define deglo    500
#define deghi    14500
#define degh     500
#define rows     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define cols     ((deghi + 1 - deglo + (degh - 1)) / degh)
#define cpumin   10
#define ncases   1
#define nalgs    3
#define img      1
#define imgname  "poly_heaps_OS.ppm"

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
    int i, j, len, deg;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_spoly_t f, g, h;
    fmpz_poly_t x, y, z;
    fmpz_t degree;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(degree);
       
    fmpz_spoly_init(f);
    fmpz_spoly_init(g);
    fmpz_spoly_init(h);
    
    fmpz_poly_init(x);
    fmpz_poly_init(y);
    fmpz_poly_init(z);

    for (len = lenlo, i = 0; len <= lenhi; len += lenh, i++)
    {
        slong s[nalgs];
        for (deg = deglo, j = 0; deg <= deghi; deg += degh, j++)
        {
            int c, n, reps = 0;
            
            for (c = 0; c < nalgs; c++)
                s[c] = WORD(0);
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t t[nalgs];
                int l, loops = 1;
                
                /*
                   Construct random polynomials f and g
                 */
                {
                  fmpz_init_set_ui(degree, deg);
                  fmpz_spoly_randtest(f, state, deg - (len*deg*5)/deglo, degree, bits);
                  fmpz_spoly_randtest(g, state, deg - (len*deg*5)/deglo, degree, bits);
                }

                /*
                 * Construct dense polynomials from f and g
                */
                {
                  fmpz_spoly_get_fmpz_poly(x, f);
                  fmpz_spoly_get_fmpz_poly(y, g);
                }

              loop:

                timeit_start(t[0]);
                for (l = 0; l < loops; l++)
                    fmpz_poly_mul(z, x, y);
                timeit_stop(t[0]);
                
                timeit_start(t[1]);
                for (l = 0; l < loops; l++)
                    fmpz_spoly_mul_heaps(h, f, g);
                timeit_stop(t[1]);

                timeit_start(t[2]);
                for (l = 0; l < loops; l++)
                    fmpz_spoly_mul_OS(h, state, f, g);
                timeit_stop(t[2]);
                
                for (c = 0; c < nalgs; c++)
                    if (t[c]->cpu <= cpumin)
                    {
                        loops *= 10;
                        goto loop;
                    }
                
                for (c = 0; c < nalgs; c++)
                    s[c] += t[c]->cpu;
                reps += loops;
            }
            
            for (c = 0; c < nalgs; c++)
                T[i][j][c] = s[c] / (double) reps;
            
            if (s[0] <= s[1] && s[0] <= s[2])
                X[i][j] = 0;
            else if (s[1] <= s[2])
                X[i][j] = 1;
            else
                X[i][j] = 2;
            flint_printf("len = %d, deg = %d, winner = %d\n", deg - (len*deg*5)/deglo, deg, X[i][j]), fflush(stdout);
        }
        {
           slong sum = 0, c;
           for (c = 0; c < nalgs; c++)
              sum += s[c];
           flint_printf("len = %d, deg = %d, time = %wdms\n", deg - (len*deg*5)/deglo, deg, sum), fflush(stdout);
        }
    }
    
    /* 
       Print 2-D ASCII image of the winning algorithms
     */
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
            flint_printf("%d", X[i][j]);
        flint_printf("\n");
    }
    
    fmpz_spoly_clear(f);
    fmpz_spoly_clear(g);
    fmpz_spoly_clear(h);
    
    fmpz_poly_clear(x);
    fmpz_poly_clear(y);
    fmpz_poly_clear(z);
    
    fmpz_clear(degree);

    /*
       Print 2-D coloured image to file
     */
    if (img)
    {
        unsigned char * PIXELS;
        int k;
        
        PIXELS = (unsigned char *) flint_malloc(3 * rows * cols * sizeof(unsigned char));
        k = 0;
        for (i = 0; i < rows; i++)
        {
            for (j = 0; j < cols; j++)
            {
                double max = DBL_MIN, v[nalgs];
                int m;
                
                for (m = 0; m < nalgs; m++)
                {
                    v[m] = T[i][j][m] - T[i][j][X[i][j]];
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
        }

        k = write_rgb_ppm(imgname, PIXELS, cols, rows);
        flint_free(PIXELS);
        
        if (k)
        {
            flint_printf("Exception:  writing ppm image failed\n");
        }
    }

    flint_randclear(state);
    return 0;
}
