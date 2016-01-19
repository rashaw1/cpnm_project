/************************************************************************************************
 *
 * PURPOSE: To define routines for orthogonalising a subset of the basis functions given an overlap
 *          matrix for that subset, and interface this to the system of Gaussians (although this
 *          will work for any given overlap matrix).
 *
 * DATE           AUTHOR              CHANGES
 * ==================================================================================
 * 18/12/15       Robert Shaw         Original code.
 *
 *************************************************************************************************/

#ifndef ORTHOGONALISEHEADERDEF
#define ORTHOGONALISEHEADERDEF

#include <Eigen/Dense>

// Declare forward dependencies
class System;

// Declare constants
const int GRAM_SCHMIDT = 1;
const int CANONICAL = 2;
const int SYM_LOWDIN = 3;

// Declare routines

// Interface routine to orthogonalise the first n basis functions
// in a System, using whichever method specified (default is canonical).
Eigen::MatrixXd orthogonalise(System& sys, int n_, const int method = CANONICAL);

// Gram-Schmidt orthogonalisation - returns the matrix f, where
// the orthogonal functions are given by the rows of  f = P L^-1
// where the rows of P are the basis functions, and L is the
// lower-triangular matrix from the Cholesky decomposition
// of the overlap matrix, S
Eigen::MatrixXd gramSchmidt(Eigen::MatrixXd& S, Eigen::MatrixXd& P);

// Canonical orthogonalisation - returns f, as above, but where
// f = P W D^-1/2, with W the eigenvectors of S, and D the diagonal
// matrix of eigenvalues of S.
Eigen::MatrixXd canonical(Eigen::MatrixXd& S, Eigen::MatrixXd& P);

// Symmetric Lowdin orthogonalisations - returns f, as above, but with
// f = P S^-1/2, where the inverse square root of a matrix is calculated
// in the usual way as S^-1/2 = W D^-1/2 W
Eigen::MatrixXd symLowdin(Eigen::MatrixXd& S, Eigen::MatrixXd& P);

#endif
