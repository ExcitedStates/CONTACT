/*
    CONTACT: COntact Networks Through Alternate Conformation Transitions
    Henry van den Bedem, Gira Bhabha, Kun Yang, Peter Wright, James S. Fraser
    e-mail: vdbedem@stanford.edu, james.fraser@ucsf.edu

        Copyright (C) 2013 Stanford University

	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to do
	so, subject to the following conditions: 

	This entire text, including the above copyright notice and this permission notice
	shall be included in all copies or substantial portions of the Software. 

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
	IN THE SOFTWARE.

    
*/

#include "CreateGraph_stats.hpp"
#include <cstdlib>

int main(int argc, char** argv) {
    Dynamic_Graph myGraph;
    myGraph.SetParameter(std::string(argv[2]), std::atof(argv[3]), std::atoi(argv[4]), std::string(argv[5]));
    myGraph.ReadFile(std::string(argv[1]));
    /*
    std::string ID1, ID2;
    int modNo;
    std::cout << "Model No: " << std::endl;
    std::cin >> modNo;
    std::cout << "Residue ID: " << std::endl;
    std::cin >> ID1;
    std::cout << "Residue ID: " << std::endl;
    std::cin >> ID
    */
    myGraph.RunChains(1);
    return 0;
}
