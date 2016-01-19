/***************************************************************************
 *
 * PURPOSE: To define a gaussian function (or gaussian-type orbital) class
 *          containing the properties of the gaussian basis functions to be
 *          placed on each atom in the System, and their associated procedures.
 *
 * CONTAINS: 
 *         class Gaussian:
 *             data:
 *                pos - a vector of the coordinates of the position of the gaussian
 *                zeta - the exponent of the gaussian
 *                norm - the normalisation constant
 *             routines:
 *                normalise() - calculates the normalisation constant
 *                dist2(other_gaussian) - calculates the square distance
 *                              from this gaussian to the other_gaussian
 *                overlap(other_gaussian) - calculates the overlap of this
 *                              and other_gaussian
 *
 *
 * DATE        AUTHOR           CHANGES
 * ========================================================
 * 17/12/15    Robert Shaw      Original code. 
 *
 ***************************************************************************/

#ifndef GAUSSIANHEADERDEF
#define GAUSSIANHEADERDEF

#include <vector> // Vector class needed for position coordinates

class Gaussian
{
private:
	// Users should not be able to modify these unwittingly
	// - public accessors are provided instead
	std::vector<double> pos; // Coordinates
	double zeta; // Exponent
	double norm; // Normalisation constant
public:
	Gaussian(); // Default constructor - will set x, y, z, zeta to zero
	Gaussian(double zeta_, double x_, double y_, double z_); // Explicit constructor

	// Accessors
	// Get-accessors are const so that variables cannot be altered
	std::vector<double> getCoords() const { return pos; }
	double getZeta() const { return zeta; }
	double getNorm() const { return norm; }

	// Calculate square distance between two gaussians
	double dist2(const Gaussian& other) const;

	// The following are virtual for extensibility - other types of GTO
	// (e.g. p, d, ...) could then be derived from this class 

	// Calculate normalisation constant
	virtual void normalise();

	// Calculate overlap integral between two gaussians
	virtual double overlap(const Gaussian& other) const;
};
	
#endif
