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

#define bits     100
#define lenlo    100
#define lenhi    10000
#define lenh     4950
#define deg      100
#define shiftlo  7
#define shifthi  20
#define shifth   6
#define nvars    2
#define rows     ((lenhi + 1 - lenlo + (lenh - 1)) / lenh)
#define cols     ((shifthi + 1 - shiftlo + (shifth - 1)) / shifth)
#define cpumin   1000
#define ncases   1
#define nalgs    3
#define img      1
#define imgname  "heaps_OS_poly_shift.ppm"

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
    slong i, j, len;
    ulong shift;
    int X[rows][cols];
    double T[rows][cols][nalgs];
    fmpz_spoly_t f, g, h;
    fmpz_poly_t fdense, gdense, hdense;
    fmpz_t degree;

    FLINT_TEST_INIT(state);
       
    fmpz_spoly_init(f);
    fmpz_spoly_init(g);
    fmpz_spoly_init(h);

    fmpz_init_set_ui(degree, deg);
    
    fmpz_poly_init(fdense);
    fmpz_poly_init(gdense);
    fmpz_poly_init(hdense);

    for (len = lenhi, i = 0; len >= lenlo; len -= lenh, i++)
    {
        slong times[nalgs];
        int c;

        j = 0;

        for (c = 0; c < nalgs; c++)
            times[c] = WORD(0);
        
        for (shift = shiftlo, j = 0; shift <= shifthi; shift += shifth, j++)
        {
            int n;
            int reps[nalgs];
            
            for (c = 0; c < nalgs; c++)
            {
                times[c] = WORD(0);
                reps[c] = 0;
            }
            
            for (n = 0; n < ncases; n++)
            {
                timeit_t timer;
                int l, loops;

                /*
                   Construct random sparse polynomials f and g
                 */
                {
                    fmpz_spoly_randtest_kron(f, state, len, degree, bits, shift, nvars);
                    fmpz_spoly_randtest_kron(g, state, len, degree, bits, shift, nvars);
                    fmpz_spoly_get_fmpz_poly(fdense, f);
                    fmpz_spoly_get_fmpz_poly(gdense, g);
                    if (fmpz_spoly_terms(f) != len)
                    {
                        flint_printf("ERROR: wrong number of terms\n");
                        abort();
                    }
                }

                /* mul_heaps */
                loops = 1;
                c = 0;
                do
                {
                    timeit_start(timer);
                    for (l = 0; l < loops; ++l)
                        fmpz_spoly_mul_heaps(h, g, f);
                    timeit_stop(timer);
                    times[c] += timer->cpu;
                    reps[c] += loops;
                    loops *= 2;
                }
                while (timer->cpu < cpumin);

                /* mul_OS */
                loops = 1;
                c = 1;
                do
                {
                    timeit_start(timer);
                    for (l = 0; l < loops; ++l)
                        fmpz_spoly_mul_OS(h, state, g, f);
                    timeit_stop(timer);
                    times[c] += timer->cpu;
                    reps[c] += loops;
                    loops *= 2;
                }
                while (timer->cpu < cpumin);

                /* fmpz_poly_mul */
                c = 2;
                if (shift >= 18)
                {
                    reps[c] = 1;
                    times[c] = WORD_MAX;
                }
                else
                {
                    loops = 1;
                    do
                    {
                        timeit_start(timer);
                        for (l = 0; l < loops; ++l)
                            fmpz_poly_mul(hdense, gdense, fdense);
                        timeit_stop(timer);
                        times[c] += timer->cpu;
                        reps[c] += loops;
                        loops *= 2;
                    }
                    while (timer->cpu < cpumin);
                }
            }
            
            for (c = 0; c < nalgs; c++)
                T[i][j][c] = times[c] / (double) reps[c];
            
            if (times[1] <= FLINT_MIN(times[0], times[2]))
                X[i][j] = 1;
            else if (times[0] <= times[2])
                X[i][j] = 0;
            else
                X[i][j] = 2;
           flint_printf("len = %wd, shift = %wu, winner = %d\n", len, shift, X[i][j]), fflush(stdout);
        }
        {
           slong sum = 0, c;
           for (c = 0; c < nalgs; c++)
              sum += times[c];
           flint_printf("len = %wd, shift = %wu, time = %wdms\n", len, shift, sum), fflush(stdout);
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
    fmpz_poly_clear(fdense);
    fmpz_poly_clear(gdense);
    fmpz_poly_clear(hdense);
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
                double v[nalgs];
                int m;
                
                for (m = 0; m < FLINT_MIN(3, nalgs); m++)
                {
                    v[m] = T[i][j][X[i][j]] / T[i][j][m];
                    PIXELS[k++] = (unsigned char) (v[m] * 255);
                }
                for(; m < 3; m++)
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
