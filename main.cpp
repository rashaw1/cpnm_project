/****************************************************************************************
 *
 * PURPOSE: Main program, forms a system from input, and does requested commands
 *
 * DATE            AUTHOR             CHANGES
 * ==============================================================================
 * 19/12/15        Robert Shaw        Original code.
 *
 ****************************************************************************************/

#include "system.hpp"
#include "io.hpp"
#include "orthogonalise.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
	int program = 0;

	// Get input file from first command line argument
	if(argc < 2){ // No input file given
		std::cerr << "Usage: ./main [input_file]\n";
		program = -1;
	} else {
		std::string ifname = argv[1]; // Input filename
		std::string ofname = ifname; // Output file prefix
		std::size_t pos = ofname.find('.');
		if (pos != std::string::npos) { // Cut off extension
			ofname.erase(pos, ofname.length());
		}

		// Open the input file
		std::ifstream input(ifname);
		// Check it opened successfully
		if (!input.is_open()){
			std::cerr << "Failed to open input file.\n";
			program = -1;
		} else {
			
			// Make the system
			System sys = makeSystem(input);

			// Calculate the overlap integrals
			sys.calcOverlap();

			// Open main output file and print system details
			std::ofstream output(ofname + ".out");
			printSystem(sys, output, true);

			// Do all the optional commands
			int lastcmd = 0;
			int flag = 1;
			int orthog = 0;
			std::vector<int> currcmd;
			Eigen::MatrixXd f;
			while(flag > 0){
				currcmd = getNextCmd(input, lastcmd);
				switch(currcmd[0]){
				case 1: { // Print the overlap integrals
					std::ofstream intout(ofname + ".ints");
					printIntegrals(sys, intout);
					intout.close();
					break;
				}
				case 2: { // Print the sparse graph data
					std::ofstream sparseout(ofname + ".sparse");
					printSparseGraph(sys, sparseout, currcmd[1]);
					sparseout.close();
					break;
				}
				case 3: { // Canonical orthogonalisation
					f = orthogonalise(sys, currcmd[1], CANONICAL);
					orthog = 1;
					break;
				}
				case 4: { // Gram-Schmidt orthogonalisation
					f = orthogonalise(sys, currcmd[1], GRAM_SCHMIDT);
					orthog = 2;
					break;
				}
				case 5: { // Symmetric Lowdin orthogonalisation
					f = orthogonalise(sys, currcmd[1], SYM_LOWDIN); 
					orthog = 3;
					break;
				}
				case -1: { // Error
					output << "\nErroneous command given.\n";
					flag = 0;
					program = -1;
					break;
				}
				default: { // No more commands
					output << "\nProgram finished.\n";
					flag = 0;
				}
				}
				lastcmd++;
			}

			// Print orthogonalisation data if needed
			if (orthog > 0) {
				std::ofstream orthogout(ofname + ".orthog");
			    printOrthog(sys, orthogout, f, orthog);
				orthogout.close();
			}
			
			output.close();
		}
	}
	
	return program;
}
	
