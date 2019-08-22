/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * libBeyn.h -- header file for libBeyn, a simple implementation of
 *           -- Beyn's algorithm for nonlinear eigenproblems
 *           -- 
 *           -- This is packaged together with SCUFF-EM, but
 *           -- it is really a standalone independent entity
 *           -- for general-purpose use in solving nonlinear
 *           -- eigenproblems.
 */


#ifndef BEYN_H
#define BEYN_H

#include <complex.h>
#include <gsl/gsl_matrix.h>

/// User-supplied function that provides the operator M(z) whose "roots" are to be found.
typedef int (*beyn_function_M_t)(gsl_matrix_complex *target_M, complex double z, void *params);

/// (optional) User-supplied function that, given \f$ \hat V \f$, calculates \f$ M(z)^{-1} \hat V \f$.
typedef int (*beyn_function_M_inv_Vhat_t)(gsl_matrix_complex *target_M_inv_Vhat,
	       	const gsl_matrix_complex *Vhat, complex double z, void *params); 

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct BeynSolver
{
   int M;   // dimension of matrices
   int L;   // number of columns of VHat matrix

   gsl_vector_complex *Eigenvalues, *EVErrors;
   gsl_matrix_complex *Eigenvectors;
   gsl_matrix_complex *A0, *A1, *A0Coarse, *A1Coarse, *MInvVHat;
   gsl_matrix_complex *VHat;
   gsl_vector *Sigma, *Residuals;
   double complex *Workspace;

 } BeynSolver;

// constructor, destructor
BeynSolver *CreateBeynSolver(int M, int L);
void DestroyBeynSolver(BeynSolver *Solver);

// reset the random matrix VHat used in the Beyn algorithm
// 
void ReRandomize(BeynSolver *Solver, unsigned int RandSeed);

// for both of the following routines,
// the return value is the number of eigenvalues found,
// and the eigenvalues and eigenvectors are stored in the
// Lambda and Eigenvectors fields of the BeynSolver structure

// Beyn method for circular contour of radius R,
// centered at z0, using N quadrature points
//int BeynSolve(BeynSolver *Solver,
//              BeynFunction UserFunction, void *Params,
//              double complex z0, double R, int N);

// Beyn method for elliptical contour of horizontal, vertical
// radii Rx, Ry, centered at z0, using N quadrature points
int BeynSolve(BeynSolver *Solver,
              beyn_function_M_t M_function, beyn_function_M_inv_Vhat_t M_inv_Vhat_function, void *params,
              double complex z0, double Rx, double Ry, int N);

#endif // BEYN_H
