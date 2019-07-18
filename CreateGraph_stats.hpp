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

#ifndef _CREATE_GRAPH_
#define _CREATE_GRAPH_
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstdio>
#include <fstream>
#include <cfloat>
#include <stack>
#include <list>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <mmdb/mmdb_manager.h>

struct Node {
    int ResInd;
    PPCAtom atoms;
    int atomNo;
    AltLoc altLoc;
};

class Dynamic_Graph {
public:
    Dynamic_Graph() {}
    bool ReadFile(std::string);
    bool SetParameter(std::string, double, int, std::string);
    bool RunChains(int MolNo);
    bool CreateGraph(int MolNo, ChainID);
    void FindPath(std::multimap<int, int>&, std::vector<Node>&, std::vector<std::vector<double> >&, std::vector<std::vector<bool> >&, const std::vector<std::vector<double> >&, const std::vector<std::vector<std::pair < std::string, std::string > > >& );
    void BuildPath(std::vector<Node>&, std::multimap<int, int>&, std::vector<std::vector<int> >&, std::vector<std::vector<bool> >&, std::vector<std::vector<double> >&, std::vector<std::pair<int, int> >, const std::vector<std::vector<double> >&, std::vector<double>&, const std::vector<std::vector<std::pair < std::string, std::string > > >& );
    void ProProcess();
private:
    PCMMDBManager myPCMMDBManager;
    int uddHandle_relief, uddHandle_norelief;
    bool isMChain(PCAtom);
    int ModelNo;
    std::vector<int> graph;
    std::vector<PCResidue> pcresidues;
    std::string file;
    double ratio;
    int length;
    bool isSideChain;
};

#endif
