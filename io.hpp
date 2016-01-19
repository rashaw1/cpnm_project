/*********************************************************************************************
 *
 * PURPOSE: To handle input and output, reading the input file, and formatting the output
 *          from the program (i.e. data and graphs)
 *
 * CONTAINS:
 *
 * DATE            AUTHOR            CHANGES 
 * =====================================================================================
 * 19/12/15        Robert Shaw       Original code.
 *
 **********************************************************************************************/

#ifndef IOHEADERDEF
#define IOHEADERDEF

#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <vector>

class System; // Forward declaration
class Gaussian;

// Returns the next command
std::vector<int> getNextCmd(std::ifstream& in, int lastCmd_);

// Identify a token
int findToken(std::string t_);

// Make the system from the geometry/basis specification
System makeSystem(std::ifstream& in);

// Add all basis functions of a particular atom type
void addAtomType(std::ifstream& in, System& sys, std::string atomtype_,
				 double exponent_, int geomstart_, int geomend_);

// Print details of the system to file
void printSystem(const System& sys, std::ofstream& out, bool printBasis = false);

// Print the non-zero overlap integrals, and their indices, to file
void printIntegrals(const System& sys, std::ofstream& out);

// Print the data for a sparse-graph
// fineness controls the size of `block' that the matrix is split into
// for example, a fineness of 100 will give block sizes of [N/100]
void printSparseGraph(const System& sys, std::ofstream& out, int fineness);

// Print the orthogonalisation results
void printOrthog(const System& sys, std::ofstream& out, const Eigen::MatrixXd& f, int orthogType);

// Print the details of a gaussian basis function
void printGaussian(const Gaussian& g, std::ofstream& out);
#endif 
