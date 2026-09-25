/* eigen/jacobi.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman
 */
/* Simple linear algebra operations, operating directly
 * on the gsl_vector and gsl_matrix objects. These are
 * meant for "generic" and "small" systems. Anyone
 * interested in large systems will want to use more
 * sophisticated methods, presumably involving native
 * BLAS operations, specialized data representations,
 * or other optimizations.
 */

/*
 * Modified by: J.-M. Nuzillard, April, 2002.
 * All calls to gsl_vector and gsl_matrix utility routines 
 * are replaced by regular C constructs.
 * The code is now self-contained but no index range checking
 * is achieved.
 * The dimension of the problem is given as first argument
 * to gsl_eigen_jacobi.
 * Returned values are 0 for success, 1 if max_rot is reached
 * and 2 if workspace memory allocation fails.
 *
 * The coordinates of the first eigenvector are
 * evec[0][0], evec[1][0], .... 
 * In other words, they are stored columnwise.
 *
 * file: OutlsdSrc/defs.h
 * Copyright(C)2000 CNRS - UMR 7312 - Jean-Marc Nuzillard
 *
 * This file is part of LSD.
 *
 * LSD is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * LSD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LSD; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdlib.h>
#include <math.h>

#define REAL double

#include "jacobi.h"
static void
jac_rotate(REAL **a,
           unsigned int i, unsigned int j, unsigned int k, unsigned int l,
           double *g, double *h,
           double s, double tau) ;


static void
jac_rotate(REAL **a,
           unsigned int i, unsigned int j, unsigned int k, unsigned int l,
           double *g, double *h,
           double s, double tau)
{
  *g = a[i][j] ;
  *h = a[k][l] ;
  a[i][j] = (*g) - s*((*h) + (*g)*tau);
  a[k][l] = (*h) + s*((*g) - (*h)*tau);
}

int
gsl_eigen_jacobi(unsigned int n,
                 REAL **a,
                 REAL *eval,
                 REAL **evec,
                 unsigned int max_rot, 
                 unsigned int *nrot)
{
    unsigned int i, j, iq, ip;
    double t, s;

    REAL *b = (REAL *) malloc(n * sizeof(REAL));
    REAL *z = (REAL *) malloc(n * sizeof(REAL));
    if(b == 0 || z == 0) {
      if(b != 0) free(b);
      if(z != 0) free(z);
	return 2 ;
    }

    /* Set eigenvectors to coordinate basis. */
    for(ip=0; ip<n; ip++) {
      for(iq=0; iq<n; iq++) {
        evec[ip][iq] = 0.0;
      }
      evec[ip][ip] = 1.0;
    }

    /* Initialize eigenvalues and workspace. */
    for(ip=0; ip<n; ip++) {
      REAL a_ipip = a[ip][ip];
      z[ip] = 0.0;
      b[ip] = a_ipip;
      eval[ip] = a_ipip;
    }

    *nrot = 0;

    for(i=1; i<=max_rot; i++) {
      REAL thresh;
      REAL tau;
      REAL g, h, c;
      REAL sm = 0.0;
      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          sm += fabs(a[ip][iq]);
        }
      }
      if(sm == 0.0) {
        free(z);
        free(b);
        return 0;
      }

      if(i < 4)
        thresh = 0.2*sm/(n*n);
      else
        thresh = 0.0;

      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          const REAL d_ip = eval[ip];
          const REAL d_iq = eval[iq];
	    const REAL a_ipiq = a[ip][iq];
          g = 100.0 * fabs(a_ipiq);
          if(   i > 4
             && fabs(d_ip)+g == fabs(d_ip)
             && fabs(d_iq)+g == fabs(d_iq)
	     ) {
            a[ip][iq] = 0.0;
          }
          else if(fabs(a_ipiq) > thresh) {
            h = d_iq - d_ip;
            if(fabs(h) + g == fabs(h)) {
              t = a_ipiq/h;
            }
            else {
              REAL theta = 0.5*h/a_ipiq;
              t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
              if(theta < 0.0) t = -t;
            }

            c   = 1.0/sqrt(1.0+t*t);
            s   = t*c;
            tau = s/(1.0+c);
            h   = t * a_ipiq;
            z[ip] -= h;
            z[iq] += h;
		eval[ip] = d_ip - h;
		eval[iq] = d_iq + h;
		a[ip][iq] = 0.0;

            for(j=0; j<ip; j++){
              jac_rotate(a, j, ip, j, iq, &g, &h, s, tau);
            }
            for(j=ip+1; j<iq; j++){
              jac_rotate(a, ip, j, j, iq, &g, &h, s, tau);
            }
            for(j=iq+1; j<n; j++){
              jac_rotate(a, ip, j, iq, j, &g, &h, s, tau);
            }
            for (j=0; j<n; j++){
              jac_rotate(evec, j, ip, j, iq, &g, &h, s, tau);
            }
            ++(*nrot);
          }
        }
      }
      for (ip=0; ip<n; ip++) {
        b[ip] += z[ip];
        z[ip]  = 0.0;
	eval[ip] = b[ip];
      }

      /* continue iteration */
    }

    return 1;
}

