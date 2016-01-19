/**********************************************************************************************
 *
 * PURPOSE: To implement the IO procedures declared in io.hpp.
 *
 * DATE           AUTHOR             CHANGES 
 * ==================================================================================
 * 19/12/15       Robert Shaw        Original code.
 *
 **********************************************************************************************/

#include "io.hpp"
#include "system.hpp"
#include "gaussian.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>

// Match token to a particular command
int findToken(std::string t)
{
	// Get rid of extraneous spaces
	t.erase(std::remove(t.begin(), t.end(), ' '), t.end());

	// Make lowercase
	std::transform(t.begin(), t.end(), t.begin(), ::tolower);

	// Find the correct command
	int rval = 0;
	if (t == "basis") { rval = 1; }
	else if (t == "geom") { rval = 2; }
	else if (t == "threshold") { rval = 3; }
	else if (t == "print") { rval = 4; }
	else if (t == "orthog") { rval = 5; }
	else if (t == "integrals") { rval = 6; }
	else if (t == "sparsegraph") { rval = 7; }
	else if (t == "canonical") { rval = 8; }
	else if (t == "gramschmidt") { rval = 9; }
	else if (t == "symlowdin") { rval = 10; }

	return rval;
}

// Read the basis, geom, and threshold, to make the system
System makeSystem(std::ifstream& in)
{
	double threshold = 1e-4; // Default threshold value
	int geomstart = 0; int geomend = 0;
	int basisstart = 0; int basisend = 0;

	// Read line by line and parse
	// to find the starting and ending points of the
	// geometry and basis specifications, and the
	// threshold value (if given).
	std::string line, token;
	int pos; int linecount = 0;
	while (std::getline(in, line)) {
		// Erase any comments, indicated by !
		pos = line.find('!');
		if (pos != std::string::npos){
			line.erase(pos, line.length());
		}

		// Tokenise
		pos = line.find(',');
		if (pos != std::string::npos){
			token = line.substr(0, pos);

			// Find token
			switch(findToken(token)){
			case 1: { // Basis
				// Find the start and end points of
				// the basis specfication
				basisstart = linecount + 1;
				basisend = basisstart -1;
				while (line != "basisend"){
					if (std::getline(in, line)){
						basisend++;
						linecount++;
					}
				}
				break;
			}
			case 2: { // Geometry
				// Find start and end of geometry 
				geomstart = linecount + 1;
				geomend = geomstart - 1;
				while (line != "geomend"){
					if (std::getline(in, line)) {
						geomend++;
						linecount++;
					}
				}
				break;
			}
			case 3: { // Threshold
				threshold = std::stod(line.substr(pos+1, line.length()));
				break;
			}
   			}
		}
		linecount ++;
	}

	// Make the system
	System sys(threshold);

	// Read the basis and geometry, if they've been specified correctly
	if ( (basisend - basisstart) > 0 && (geomend - geomstart) > 0) {
		// Rewind to beginning of file
		in.clear();
		in.seekg(0, std::ios::beg);

		std::vector<double> exponents;
		std::vector<std::string> atomtypes;

		for (int l = 0; l < basisend; l++) {
			std::getline(in, line);
			if (l > basisstart - 1) {
				pos = line.find(','); // Find first token, the atom type
				token = line.substr(0, pos); // Tokenise

				// Get rid of whitespace
				token.erase(std::remove(token.begin(), token.end(), ' '), token.end());

				// Assign atom type and exponent
				atomtypes.push_back(token);
				exponents.push_back(std::stod(line.substr(pos+1, line.length())));
			}
		}

		// Add all the atoms in atomtypes
		for (int i = 0; i < atomtypes.size(); i++)
			addAtomType(in, sys, atomtypes[i], exponents[i], geomstart, geomend);

	} else {
		// Throw error
		if (basisend - basisstart <= 0) 
			std::cerr << "No basis specification was given!\n";
		
		if (geomend - geomstart <= 0)
			std::cerr << "No geometry specification was given!\n";
	}

	return sys;
}

// Loop through geometry specification
// and add all atoms of a particular type to a system
void addAtomType(std::ifstream& in, System& sys, std::string atomtype_,
				 double exponent_, int geomstart_, int geomend_)
{
	// Rewind to beginning of file
	in.clear();
	in.seekg(0, std::ios::beg);

	// Read the geometry line by line
	double x, y, z;
	std::string curratom, line;
	int pos;

	for (int l = 0; l < geomend_; l++) {
		std::getline(in, line);
		if (l > geomstart_ - 1) {
			pos = line.find(',');
			curratom = line.substr(0, pos);
			if (curratom == atomtype_) {
				// Read the geometry and add a basis function to sys
				line.erase(0, pos+1);
				pos = line.find(',');
				x = std::stod(line.substr(0, pos-1));
				line.erase(0, pos+1);
				pos = line.find(',');
				y = std::stod(line.substr(0, pos));
				line.erase(0, pos+1);
				z = std::stod(line);

				Gaussian g(exponent_, x, y, z);
				sys.addGaussian(g);
			}
		}
	}
}

// Return the next command
// -1 = Error
// 0 = no more commands
// 1 = print integrals
// 2, fineness = print sparsegraph
// 3, n = canonical orthog. first n funcs
// 4, n = gram-schmidt orthog. first n funcs
// 5, n = sym. lowdin orthog. first n funcs
std::vector<int> getNextCmd(std::ifstream& in, int lastCmd_)
{
    // Rewind to beginning of file
	in.clear();
	in.seekg(0, std::ios::beg);
	
	std::vector<int> cmd;
	
	// Read line by line until find the correct command
	std::string line, token;
	int pos; int cmdcount = 0;
	bool found = false;
	
	// Get to end of geometry first
	while(!found && std::getline(in, line)){
		if (line == "geomend")
			found = true;
	}
	
	while (cmdcount < lastCmd_+1 && std::getline(in, line)) {
		// Erase any comments
		pos = line.find('!');
		if (pos != std::string::npos){
			line.erase(pos, line.length());
		}

		// Tokenise
		pos = line.find(',');
		if (pos != std::string::npos){
			token = line.substr(0, pos);
			int id = findToken(token);
			// Find token
			switch(id){
			case 4: { // Print command
				cmdcount++;

				// Is this the next command?
				if (cmdcount == lastCmd_+1){
					// Get the rest of the command
				    token = line.substr(pos+1, line.length());

					int newpos = token.find(',');
					if( newpos != std::string::npos) token = token.substr(0, newpos);
					
					switch(findToken(token)){
					case 6: { // Print overlap integrals
						cmd.push_back(1);
						break;
					}
					case 7: { // Print sparse graph data
						cmd.push_back(2);
						token = line.substr(pos+newpos+2, line.length());
						cmd.push_back(std::stoi(token));
						break;
					}
					default: {
						std::cerr << "Invalid print command\n";
						cmd.push_back(-1);
					}
					}
				}
				break;
			}
			case 5: { // Orthog command
				cmdcount++;

				// Is this the next command?
				if (cmdcount == lastCmd_+1){
					// Get the next part of the command
					line.erase(0, pos+1);
					pos = line.find(',');
					if (pos != std::string::npos){
						token = line.substr(0, pos);
						int nfuncs = std::stoi(line.substr(pos+1, line.length()));
											   
						switch(findToken(token)){
						case 8: { // Canonical
							cmd.push_back(3);
							break;
						}
						case 9: { // Gram-Schmidt
							cmd.push_back(4);
							break;
						}
						case 10: { // Symmetric Lowdin
							cmd.push_back(5);
							break;
						}
						default: {
							std::cerr << "Invalid orthog command.\n";
							cmd.push_back(-1);
						}
						}

						cmd.push_back(nfuncs);
					} else {
						std::cerr << "No orthogonalisation method specified!\n";
						cmd.push_back(-1);
					}
				}
				break;
			}
			default: {
				if (id > 5 || id < 1){
					std::cerr << "Command not found. " << id << "\n";
					cmd.push_back(-1);
				}
			}
			}
		} 

	}

	// Check if no more commands were found
	if (cmdcount <= lastCmd_) cmd.push_back(0);
	
	return cmd;
}

// Print details of the system to file
void printSystem(const System& sys, std::ofstream& out, bool printBasis)
{
	int N = sys.getN();
	out << "\nThis system has " << N << " basis functions\n"
		<< "It's sparsity is: " << sys.sparsity() << " percent\n"
		<< "with a threshold of " << sys.getThreshold() << "\n\n"
		<< "That is equivalent to " << sys.getZeroes()
		<< " zeroes out of " << (N*(N+1))/2 << " possible unique integrals.\n\n";

	// Print details of the basis functions, if wanted
	if (printBasis) {
		out << "LIST OF BASIS FUNCTIONS\n\n"
			<< std::setw(12) << "Zeta"
			<< std::setw(12) << "Norm"
			<< std::setw(12) << "x"
			<< std::setw(12) << "y"
			<< std::setw(12) << "z\n"
			<< std::string(60, '.') << "\n";
		out << std::setprecision(6);
		for (int i = 0; i < N; i++)
			printGaussian(sys.getGaussian(i), out);
	}
}

// Print details of a Gaussian basis function
void printGaussian(const Gaussian& g, std::ofstream& out)
{
	std::vector<double> coords = g.getCoords();
	out << std::setw(12) << g.getZeta()
		<< std::setw(12) << g.getNorm()
		<< std::setw(12) << coords[0]
		<< std::setw(12) << coords[1]
		<< std::setw(12) << coords[2] << "\n";
}
											   
// Print the non-zero integrals to file
void printIntegrals(const System& sys, std::ofstream& out)
{
	int N = sys.getN();
	out << "NON-ZERO INTEGRALS: " << (N*(N+1))/2 - sys.getZeroes() << "\n\n";
	
	out << std::setw(8) << "Row"
		<< std::setw(8) << "Column"
		<< std::setw(20) << "Integral\n"
		<< std::string(36, '.') << "\n";

	out << std::setprecision(8);

	// Loop through rows and columns in order
	int currIndex = 1; 
	int position = 0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < i+1; j++){
			if (sys.sIndices[position] == currIndex){
				out << std::setw(8) << j+1
					<< std::setw(8) << i+1
					<< std::setw(20) << sys.sInts[position] << "\n";

				position++;
				if (position >= sys.sInts.size()) {
					j = i+1;
					i = N;
				}
			}
			currIndex++;
		}
	}
}			

// Print out the sparse graph data
// This data is used to make a greyscale graph in which
// the density of the matrix is represented by darkness
// fineness controls the size of the blocks, where for example
// 100 would give blocks of size [N/100]
void printSparseGraph(const System& sys, std::ofstream& out, int fineness)
{
	int N = sys.getN();
	
	// Fineness only makes sense if positive
	if (fineness < 1) {
		std::cerr << "Invalid choice of fineness - must be > 0.\n";
		fineness = 1;
	}
    int blocksize = (N/fineness > 0 ? N/fineness : 1);
	int matrixSize = N/blocksize + (N-N/blocksize > 0 ? 1 : 0);

	// Make a matrix of densities, all entries initialised to zero 
	std::vector<std::vector<double> > densities(matrixSize, std::vector<double>(matrixSize, 0.0));
	
	// Loop through the overlap matrix, summing the blocks
    int currIndex = 1;
	int position = 0;
	int currRow = 0;
	int currCol = 0;
	for (int i = 0; i < N; i++){
		currCol = i/blocksize; // Will give nearest integer below
		for (int j = 0; j < i+1; j++) {
			currRow = (j)/blocksize;
			if (sys.sIndices[position] == currIndex){
				densities[currRow][currCol] += sys.sInts[position];

				position++;
				if (position >= sys.sInts.size()) {
					j = i+1;
					i = N;
				}
			}
			currIndex++;
		}
	}

	// Print this to file in the format
	// row    col    density
	out << std::setprecision(6);
	for(int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			// Symmetrise the matrix (below diagonal will be zero otherwise)
			if ( i > j ) densities[i][j] = densities[j][i]; 

			out << std::setw(15) << i
				<< std::setw(15) << j
				<< std::setw(15) << densities[i][j] << "\n";
		}
	}		
}

// Print the results of the orthogonalisation procedure to file
void printOrthog(const System& sys, std::ofstream& out, const Eigen::MatrixXd& f, int orthogType)
{
	int nfuncs = f.rows();

	std::string oType;
	switch(orthogType){
	case 1: { oType = "CANONICAL"; break; }
	case 2: { oType = "GRAM SCHMIDT"; break; }
	case 3: { oType = "SYMMETRIC LOWDIN"; break; }
	default: oType = "UNKNOWN";
	}
	out << oType << " ORTHOGONALISATION RESULTS\n\n";
	
	// First print out details of the relevant Gaussians
	out << "BASIS FUNCTIONS\n"
		<< std::setw(12) << "Zeta"
		<< std::setw(12) << "Norm"
		<< std::setw(12) << "x"
		<< std::setw(12) << "y"
		<< std::setw(12) << "z\n"
		<< std::string(60, '.') << "\n";
	out << std::setprecision(4);
	for (int i = 0; i < nfuncs; i++)
		printGaussian(sys.getGaussian(i), out);

	// Now print the coefficients of the orthonormalised functions
	// in same order as above
	out << "\n\nFUNCTION SPECIFICATION";

	for (int i = 0; i < nfuncs; i++){
		out << "\nFUNCTION " << i+1 << " COEFFICIENTS:\n";
		for (int j = 0; j < nfuncs; j++)
			out << f(j, i) << "\n";
	}
	
}




											   
							
