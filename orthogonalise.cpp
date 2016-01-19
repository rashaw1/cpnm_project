/**********************************************************************************************
 *
 * PURPOSE: Implements orthogonalisation routines.
 *
 * DATE             AUTHOR             CHANGES          
 * ==============================================================================
 * 18/12/15         Robert Shaw        Original code.
 * 19/12/15         Robert Shaw        Sparse matrix unpacking added. 
 *
 **********************************************************************************************/

#include "orthogonalise.hpp"
#include "system.hpp"
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cmath>

// Interface to orthogonalise a subset of a System's basis functions
Eigen::MatrixXd orthogonalise(System& sys, int n_, const int method)
{

	// Check n_ is at least size 1, and not bigger than the number of
	// basis functions present in sys
	if (n_ <= 0){
		// Throw an error
		std::cerr << "Invalid number of functions to orthogonalise.\n"
				  << "Doing up to five instead.\n";

		// Set n_ to be the smaller of 5 and N
		n_ = ( sys.getN() >= 5 ? 5 : sys.getN() );
	} else if (n_ > sys.getN()) {
		// Throw an error
		std::cerr << "Too many basis functions to orthogonalise.\n"
				  << "\nDoing all available functions instead.\n";
		n_ = sys.getN();
	}
		
	// Form the overlap matrix of the first n_ basis functions in
	// system, by starting with a matrix full of zeroes, and then
	// looping through the vector of non-zero overlap integrals until
	// all relevant integrals are found
	
	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_, n_); // n_ x n_ matrix of zeroes
	
	// If all the non-diagonal overlap integrals are zero (it could happen!)
	// then nothing much needs to be done (they're already orthogonal)
	if ( sys.sInts.size() == sys.getN() ) {

		// S should just be the identity matrix
		for (int i = 0; i < n_; i++) S(i, i) = 1.0;

	} else {

		// Loop through the sInts, putting the non-zero
		// integrals into the overlap matrix
		int count = 1;
		int position = 0;
		for (int i = 0; i < n_; i++){
			for (int j = 0; j < i+1; j++){

				// Check we have not reached the end of sInts
				if (position < sys.sInts.size()) {

					// If at the correct index, insert integral
					if (sys.sIndices[position] == count) {
						S(j, i) = sys.sInts[position];
						S(i, j) = S(j, i); // Overlap matrix is symmetric
						// Move to next non-zero integral
						position++;
					}

				}

				// Move to next position in matrix
				count++;
			}
		}
	}

	// Now that the overlap matrix has been formed, call the correct
	// orthogonalisation routine. The coefficient matrix for these
	// systems will always be the identity matrix, as the basis functions
	// are the centres.
	Eigen::MatrixXd P = Eigen::MatrixXd::Identity(n_, n_);
	
	switch(method){

	case GRAM_SCHMIDT: {
		S = gramSchmidt(S, P);
		break;
	}
	case CANONICAL: {
		S = canonical(S, P);
		break;
	}
	case SYM_LOWDIN: {
		S = symLowdin(S, P);
		break;
	}
	default: {
		// Throw error
		std::cerr << "Unknown method requested.\n"
				  << "Defaulting to canonical.\n";
        S = canonical(S, P);
	}

	}

	return S;
}

// Gram-Schmidt orthogonalisation
Eigen::MatrixXd gramSchmidt(Eigen::MatrixXd& S, Eigen::MatrixXd& P)
{
	// Compute Cholesky decomposition, S = LL^T
	Eigen::LLT<Eigen::MatrixXd> lltOfS(S);
	// Retrieve lower triangular matrix, L, and compute its inverse
	Eigen::MatrixXd L = lltOfS.matrixL();
	L = L.inverse();

	// Return PL^-1
	return P*L;
}

// Canonical orthogonalisation
Eigen::MatrixXd canonical(Eigen::MatrixXd& S, Eigen::MatrixXd& P)
{
	// Calculate eigenvalues/eigenvectors of S (a self adjoint matrix)
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
	// Retrieve eigenvectors
	Eigen::MatrixXd W = solver.eigenvectors();

	// Retrieve eigenvalues as a diagonal matrix
	Eigen::MatrixXd D = solver.eigenvalues().asDiagonal();
	// Compute D^-1/2, simply by taking the inverse square
	// root of each element (as it is diagonal)
	for (int i = 0; i < D.rows(); i++){
		// All the eigenvalues should be positive and non-zero,
		// but it is worth checking at this stage
		if (D(i, i) > 0){
			D(i, i) = 1.0/sqrt(D(i, i));
		} else {
			// Throw an error
			std::cerr << "Non-positive eigenvalue "
					  << D(i, i)
					  << " found in canonical procedure.\n"
					  << "Setting to zero.\n";
			D(i, i) = 0.0;
		}
	}

	// Return PWD^-1/2
	return P*W*D;
}

// Symmetric Lowdin orthogonalisation
Eigen::MatrixXd symLowdin(Eigen::MatrixXd& S, Eigen::MatrixXd& P)
{
    // Calculate eigenvalues/eigenvectors of S (a self adjoint matrix)
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
	// Retrieve eigenvectors
	Eigen::MatrixXd W = solver.eigenvectors();
	
	// Retrieve eigenvalues as a diagonal matrix
	Eigen::MatrixXd D = solver.eigenvalues().asDiagonal();
	// Compute D^1/2
	for (int i = 0; i < D.rows(); i++){
		// Eigenvalues should be positive 
		// but should check here just in case
		if (D(i, i) > 0){
			D(i, i) = 1.0/sqrt(D(i, i));
		} else {
			// Throw an error
			std::cerr << "Non-positive eigenvalue "
					  << D(i, i)
					  << " found in sym. Lowdin procedure.\n"
					  << "Setting to zero.\n";
			D(i, i) = 0.0;
		}
	}

	// Return P S^-1/2 = P W D^-1/2 W^T
	return P*W*D*W.transpose();
}
