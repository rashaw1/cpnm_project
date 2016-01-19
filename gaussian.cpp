/**************************************************************************************
 *
 * PURPOSE: Implements class Gaussian. 
 *
 * DATE         AUTHOR           CHANGES
 * ================================================================
 * 17/12/15     Robert Shaw      Original code. 
 *
 *************************************************************************************/

#include "gaussian.hpp"
#include <cmath> // For M_PI; pow, and exp functions
#include <numeric> // for inner_product
#include <iostream>

// Constructors

// Default constructor initialises zeta to 1 and positions to origin
Gaussian::Gaussian() : zeta(1.0)
{
	// Resize the pos vector of coordinates and set to zero
	pos.resize(3, 0.0);

	// Calculate the normalisation constant
	normalise();
}

// Explicit constructor
Gaussian::Gaussian(double zeta_, double x_, double y_, double z_) : zeta(zeta_)
{
	// Push the x-, y-, and z-coordinates into the pos vector
	pos.push_back(x_); // x = pos[0]
	pos.push_back(y_); // y = pos[1]
	pos.push_back(z_); // z = pos[2]

	// Make sure that zeta isn't zero
	if (zeta == 0) zeta = 1.0;
	
	// Calculate the normalisation constant
	normalise();
}

// Routines

// Calculate the square distance between two Gaussians - i
double Gaussian::dist2(const Gaussian& other) const
{
	// The square distance between vectors x, y
	// = |x - y|^2 = (x - y).(x - y)
	std::vector<double> dist = pos;
	for (int i = 0; i < 3; i++) dist[i] -= other.pos[i];
	
	return std::inner_product(dist.begin(), dist.end(), dist.begin(), 0.0);
}

// Calculate the normalisation constant (overridable)
void Gaussian::normalise()
{
	// The normalisation constant is (2*zeta/PI)^(3/4)
	// - see Helgaker et al, Molecular Electronic Structure Theory, 9.2.1
	norm = 2.0*zeta/M_PI;
	norm = pow(norm, 0.75);
}

// Calculate the overlap between two Gaussians
double Gaussian::overlap(const Gaussian& other) const
{
	// This is calculated using the Gaussian Product Rule
	// See Helgaker et al, Molecular Electronic Structure Theory, 9.2.3

	double p = zeta + other.zeta; // Total exponent
	double mu = zeta*other.zeta/p; // Reduced exponent
	double K = exp( -mu * dist2(other) ); // Pre-exponential factor
	p = M_PI /p; // Reuse the variable
	
	return K*norm*other.norm*pow(p, 1.5);
}



