/* wavelet/dwt.c
 *
 * Copyright (C) 2004 Ivo Alxneit
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* function dwt_step is based on the public domain function pwt.c
 * available from http://www.numerical-recipes.com
 */

#include "gsl_errno.h"
#include "gsl_wavelet.h"

#define ELEMENT(a,stride,i) ((a)[(stride)*(i)])

static int binary_logn (const size_t n);
static void dwt_step (const gsl_wavelet * w, float *a, size_t stride,
                      size_t n, int isign, gsl_wavelet_workspace * work);

static int
binary_logn (const size_t n)
{
  size_t ntest;
  size_t logn = 0;
  size_t k = 1;

  while (k < n)
    {
      k *= 2;
      logn++;
    }

  ntest = (1 << logn);

  if (n != ntest)
    {
      return -1;                /* n is not a power of 2 */
    }

  return logn;
}

static void
dwt_step (const gsl_wavelet * w, float *a, size_t stride, size_t n,
          gsl_wavelet_direction dir, gsl_wavelet_workspace * work)
{
  double ai, ai1;
  size_t i, ii;
  size_t jf;
  size_t k;
  size_t n1, ni, nh, nmod;

  for (i = 0; i < work->n; i++)
    {
      work->scratch[i] = 0.0;
    }

  nmod = w->nc * n;
  nmod -= w->offset;            /* center support */

  n1 = n - 1;
  nh = n >> 1;

  if (dir == forward)
    {
      for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
          ni = i + nmod;
          for (k = 0; k < w->nc; k++)
            {
              jf = n1 & (ni + k);
              work->scratch[ii] += w->h1[k] * ELEMENT (a, stride, jf);
              work->scratch[ii + nh] += w->g1[k] * ELEMENT (a, stride, jf);
            }
        }
    }
  else
    {
      for (ii = 0, i = 0; i < n; i += 2, ii++)
        {
          ai = ELEMENT (a, stride, ii);
          ai1 = ELEMENT (a, stride, ii + nh);
          ni = i + nmod;
          for (k = 0; k < w->nc; k++)
            {
              jf = (n1 & (ni + k));
              work->scratch[jf] += (w->h2[k] * ai + w->g2[k] * ai1);
            }
        }
    }

  for (i = 0; i < n; i++)
    {
      ELEMENT (a, stride, i) = work->scratch[i];
    }
}

int
gsl_wavelet_transform (const gsl_wavelet * w, 
                       float *data, size_t stride, size_t n,
                       gsl_wavelet_direction dir, 
                       gsl_wavelet_workspace * work)
{
  size_t i;

  if (work->n < n)
    {
      GSL_ERROR ("not enough workspace provided", GSL_EINVAL);
    }

  if (binary_logn (n) == -1)
    {
      GSL_ERROR ("n is not a power of 2", GSL_EINVAL);
    }

  if (n < 2)
    {
      return GSL_SUCCESS;
    }

  if (dir == forward)
    {
      for (i = n; i >= 2; i >>= 1)
        {
          dwt_step (w, data, stride, i, dir, work);
        }
    }
  else
    {
      for (i = 2; i <= n; i <<= 1)
        {
          dwt_step (w, data, stride, i, dir, work);
        }
    }

  return GSL_SUCCESS;
}



int
gsl_wavelet_transform_forward (const gsl_wavelet * w, 
                               float *data, size_t stride, size_t n,
                               gsl_wavelet_workspace * work)
{
  return gsl_wavelet_transform (w, data, stride, n, forward, work);
}

int
gsl_wavelet_transform_inverse (const gsl_wavelet * w, 
                               float *data, size_t stride, size_t n,
                               gsl_wavelet_workspace * work)
{
  return gsl_wavelet_transform (w, data, stride, n, backward, work);
}



