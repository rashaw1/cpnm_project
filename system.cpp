/***************************************************************************************
 *
 * PURPOSE: To implement class System
 *
 * DATE           AUTHOR             CHANGES
 * ======================================================================
 * 18/12/15       Robert Shaw        Original code.
 *
 ***************************************************************************************/

#include "system.hpp"
#include <iostream>

// Constructor
System::System(double THRESHOLD_) : N(0), zeroes(0), THRESHOLD(THRESHOLD_)
{
}

// Add a Gaussian
void System::addGaussian(Gaussian g_)
{
	gaussians.push_back(g_);
	N += 1;
}

// Calculate the overlap integrals
void System::calcOverlap()
{
	double currentIntegral;
	int currentIndex;
	
	// Loop over all unique pairs of Gaussians
	// (the overlap matrix is necessarily real, symmetric, positive definite)
	// i=j should give an overlap of 1 
	for (int i = 0; i < N; i++){
		for (int j = 0; j < i+1; j++){

			// Calculate integral
			currentIntegral = gaussians[i].overlap(gaussians[j]);

			// Check if lower than threshold
			if (currentIntegral < THRESHOLD) zeroes+=1;
			else {
				// Push the non-zero integral into the vector of overlap integrals
				sInts.push_back(currentIntegral);

				// Determine the correct matrix index given our packing system
				currentIndex = 1 + j + ( i*(i+1) )/2;
				sIndices.push_back(currentIndex);
			}
		}
	}
}

// Calculate the sparsity
double System::sparsity() const
{
	// Only the upper half of the matrix, and the diagonal, have
	// been calculated.
	double numEntries = 0.5*N*(N+1);

	// Return sparsity as a percentage
	return 100.0*( zeroes / numEntries );
}

// Overloaded equals operator
System& System::operator=(const System& other)
{
	// Copy N, threshold, zeroes
	N = other.N;
	THRESHOLD = other.THRESHOLD;
	zeroes = other.zeroes;

	// Deep copy gaussians
	gaussians.clear();
	for (int i = 0; i < other.gaussians.size(); i++)
		gaussians.push_back(other.gaussians[i]);

	return *this;
}
