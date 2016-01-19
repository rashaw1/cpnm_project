/*************************************************************************************
 * 
 * PURPOSE: To define a system of 'atoms', represented by gaussian function(s) placed
 *          at their locations, and the routines to calculate the overlap matrix,
 *          and determine its sparsity
 *
 * CONTAINS:
 *          data:
 *              N - the number of Gaussians
 *              gaussians - an array of Gaussians
 *              sInts - a vector of the non-zero overlap integrals - as the overlap matrix
 *                      is expected to be sparse, it is better to store only the non-zero
 *                      values, assuming that the 2N integers extra (see below) take
 *                      relatively little memory.  
 *              sIndices - vector of the matrix index of the corresponding sInt
 *              zeroes - a counter for the number of zero overlap elements
 *              THRESHOLD - the threshold under which the overlap integral is considered
 *                          to be zero.
 *          routines:
 *              calcOverlap() - calculates the overlap integrals, and at the same time,
 *                              the number of zeroes in the overlap matrix.
 *              sparsity() - determines the sparsity (percentage of zeroes) of the overlap
 *                           matrix
 *
 *
 * DATE         AUTHOR           CHANGES
 * ===========================================================================
 * 17/12/15     Robert Shaw      Original code. 
 * 
 ************************************************************************************/

#ifndef SYSTEMHEADERDEF
#define SYSTEMHEADERDEF

#include <vector>
#include "gaussian.hpp"

class System
{
private:
	int N, zeroes; // Number of gaussians, and zeroes in overlap matrix
	std::vector<Gaussian> gaussians; // Gaussian functions
	double THRESHOLD; // Threshold under which integrals are considered zero
public:
	std::vector<double> sInts; // All non-zero overlap integrals
	std::vector<int> sIndices; // The indices of the non-zero overlap integrals

	System(double THRESHOLD_); // Constructor

	// Accessors
	int getN() const { return N; }
	int getZeroes() const { return zeroes; }
	double getThreshold() const { return THRESHOLD; }
	const Gaussian& getGaussian(int i) const { return gaussians[i]; }
	
	void addGaussian(Gaussian g_); // Adds a Gaussian function to the System
	void calcOverlap(); // Calculates the overlap matrix, determines no. of zeroes
	double sparsity() const; // Calculates the sparsity of the overlap matrix

	// Overload the equals operator
	System& operator=(const System& other);
};

#endif
