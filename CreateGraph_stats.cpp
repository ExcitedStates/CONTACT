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
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <mmdb/mmdb_tables.h>

#define UNDEFINED -1
#define OFFSET 3

double gPRCNTL;
int gNpaths, gNfailedpaths;
using namespace std;

double norm(double x1, double y1, double z1,
            double x2, double y2, double z2) {
    double a, b, c;
    a = x1 - x2;
    b = y1 - y2;
    c = z1 - z2;
    return sqrt(a * a + b * b + c * c);
}

double getPercentile(vector<double>& arr, const double prcntl )
{
  double value = 0;
  if ( arr.size() > 0 ) {
    sort ( arr.begin(), arr.end() );
    int idx = int(arr.size()*prcntl);
    value = arr[idx];
  }
  else
   value = 0.;
  cout << prcntl << " Percentile of array is: " << value << endl;
  return value;
}

bool Dynamic_Graph::isMChain(PCAtom atom) {
    if(atom->CheckID(" CA ", " C", "*") || atom->CheckID(" N  ", " N", "*") || atom->CheckID(" H  ", " H", "*") || atom->CheckID(" HA ", " H", "*") || atom->CheckID(" O  ", " O", "*") || atom->CheckID(" C  ", " C", "*")) 
        return true;
    return false;
}

void Dynamic_Graph::BuildPath(vector<Node>& nodes, multimap<int, int>& ResInd2NodeID, vector<vector<int> >& CollisionPair, vector<vector<bool> >& graph_indicator, vector<vector<double> >& graph, vector<pair<int, int> > path, const vector<vector<double> >& ovrlp, vector<double>& relief, const vector<vector<pair<string,string> > >& atmnames ) {
        pair<int, int> EndNode = path[path.size() - 1];
        bool flag = true;
        if(path.size() < length) {
            int head = EndNode.first;
            int tail = EndNode.second;
            set<int> pathSet;
            set<int> resSet;
            for(int i = 0; i < path.size(); ++i) {
                    pathSet.insert(path[i].first);
                    pathSet.insert(path[i].second);
                    resSet.insert(nodes[path[i].first].ResInd);
                    resSet.insert(nodes[path[i].second].ResInd);
            }
            pair<multimap<int, int>::iterator, multimap<int, int>::iterator> iter;
            iter = ResInd2NodeID.equal_range(nodes[tail].ResInd);
            for(multimap<int, int>::iterator it = iter.first; it != iter.second; ++it) {
                    if((*it).second == tail) continue;
                    
		    bool allsigns = false;
                    bool sign = false;
                    for(int ind = 0; ind < path.size(); ++ind) {
                            sign = sign || /*(graph_indicator[path[ind].first][(*it).second]*/ (ovrlp[path[ind].first][(*it).second] <= gPRCNTL);
			    //cout << ovrlp[path[ind].first][(*it).second] << " ";
                    }
		    allsigns = sign || allsigns;
                    if(!sign /*&& pcresidues[nodes[(*it).second].ResInd]->GetSeqNum() == 361*/) {
			    relief.push_back ( ovrlp[path[path.size() - 1].first][(*it).second] );
                            //relief_altLoc.push_back ( nodes[(*it).second].altLoc );
                            if(accumulate(CollisionPair[(*it).second].begin(), CollisionPair[(*it).second].end(), 0)) {
                                    for(int k = 0; k < CollisionPair[(*it).second].size(); ++k) {
                                            if(CollisionPair[(*it).second][k] == 1) {
                                                    if(!resSet.count(nodes[k].ResInd)) {
    							flag = false;
                                                        vector<pair<int, int> > tmpPath = path;
                                                        tmpPath.push_back(make_pair((*it).second, k));
                                                        BuildPath(nodes, ResInd2NodeID, CollisionPair, graph_indicator, graph, tmpPath, ovrlp, relief, atmnames );
                                                    }
                                            }
                                    }
                            }
                    }
	    if (allsigns ) {
		    gNfailedpaths++;
		    int N_norelief;
		    pcresidues[nodes[(*it).second].ResInd]->GetUDData (uddHandle_norelief, N_norelief );
		    N_norelief++;
		    pcresidues[nodes[(*it).second].ResInd]->PutUDData(uddHandle_norelief, N_norelief );
	            }
            }
        }
        if(flag/*&& pcresidues[nodes[path[0].first].ResInd]->GetSeqNum() == 409*/) {
                FILE* pFile = fopen(file.c_str(), "a");
                pair<multimap<int, int>::iterator, multimap<int, int>::iterator> iter;
                iter = ResInd2NodeID.equal_range(nodes[path[path.size() - 1].second].ResInd);
                for(multimap<int, int>::iterator it = iter.first; it != iter.second; ++it) {
                        if((*it).second != path[path.size() - 1].second && /*!graph_indicator[path[path.size() - 1].first][(*it).second]*/ 
                           ovrlp[path[path.size() - 1].first][(*it).second] > gPRCNTL) {// && !resSet.count(nodes[(*it).second].ResInd)) {
                                fprintf(pFile, "PATH(ModelNo, ChainID, SeqNum(InsertionCode), AltLoc (clash) --> AltLoc (relief), clash with prvs, relief, clash atom (prvs rsd), clash atom (current rsd) \n");
				if ( gNpaths++ > 1000000 ) {std::cout << "Paths exceeded 1,000,000\n"; exit(1);}
                                string atmi = "-";
                                string atmj = "-";
		    	        int N_relief;
                                for(int i = 0; i < path.size(); ++i) {
                                        if ( i>0 ) {
                                          if ( path[i-1].first < path[i-1].second ) {
                                            atmi = atmnames[path[i-1].first][path[i-1].second].first;
                                            atmj = atmnames[path[i-1].first][path[i-1].second].second;
                                          }
                                          else {
                                            atmi = atmnames[path[i-1].first][path[i-1].second].second;
                                            atmj = atmnames[path[i-1].first][path[i-1].second].first;
                                          }
                                        }
                                        fprintf(pFile, "%d, %s, %d(%s), %s --> %s, %f, %f, %s, %s\n", 
                                            pcresidues[nodes[path[i].first].ResInd]->GetModelNum(),  // Model Number
                                            pcresidues[nodes[path[i].first].ResInd]->GetChainID(), // Chain
                                            pcresidues[nodes[path[i].first].ResInd]->GetSeqNum(), // Residue Number
                                            pcresidues[nodes[path[i].first].ResInd]->insCode, // Insertion code
                                            i>0 ? nodes[path[i-1].second].altLoc : nodes[path[i].first].altLoc, // Second of prvs pair in path is clashing aternate conformer of this residue
                                            i>0 ? nodes[path[i].first].altLoc : "-", // First of current pair relieves the previous clash with this residue
                                            i>0 ? ovrlp[path[i-1].first][path[i-1].second] : 0,  // Clash value between previous and this residue
                                            i>0 ? relief[i-1]: 0, // Relief value by changing alternate values of this residue
                                            atmi.c_str(), // Atom of previous residue 
                                            atmj.c_str() ) ; // Atom of current residue clashing
		    	                    pcresidues[nodes[path[i].first].ResInd]->GetUDData (uddHandle_relief, N_relief );
		    	                    N_relief++;
		    	                    pcresidues[nodes[path[i].first].ResInd]->PutUDData(uddHandle_relief, N_relief );
                                }
                                if ( path[path.size() - 1].first < path[path.size() - 1].second ) {
                                  atmi = atmnames[path[path.size() - 1].first][path[path.size() - 1].second].first;
                                  atmj = atmnames[path[path.size() - 1].first][path[path.size() - 1].second].second;
                                }
                                else {
                                  atmi = atmnames[path[path.size() - 1].first][path[path.size() - 1].second].second;
                                  atmj = atmnames[path[path.size() - 1].first][path[path.size() - 1].second].first;
                                }
                                fprintf(pFile, "%d, %s, %d(%s), %s --> %s, %f, %f, %s, %s\n\n", 
                                    pcresidues[nodes[(*it).second].ResInd]->GetModelNum(), 
                                    pcresidues[nodes[(*it).second].ResInd]->GetChainID(),
                                    pcresidues[nodes[(*it).second].ResInd]->GetSeqNum(),
                                    pcresidues[nodes[(*it).second].ResInd]->insCode,
                                    nodes[path[path.size() - 1].second].altLoc,
                                    nodes[(*it).second].altLoc, 
                                    ovrlp[path[path.size() - 1].first][path[path.size() - 1].second],
                                    ovrlp[path[path.size() - 1].first][(*it).second],
                                    atmi.c_str(),
                                    atmj.c_str() );//, nodes[(*it).second].altLoc);                              
		    	            pcresidues[nodes[(*it).second].ResInd]->GetUDData (uddHandle_relief, N_relief );
		    	            N_relief++;
		    	            pcresidues[nodes[(*it).second].ResInd]->PutUDData(uddHandle_relief, N_relief );
                        }
                }
                //printf("%d, %s, %d(%s), %s\n", pcresidues[nodes[path[path.size() - 1].second].ResInd]->GetModelNum(), pcresidues[nodes[path[path.size() - 1].second].ResInd]->GetChainID(), pcresidues[nodes[path[path.size() - 1].second].ResInd]->GetSeqNum(), pcresidues[nodes[path[path.size() - 1].second].ResInd]->insCode, nodes[path[path.size() - 1].second].altLoc);
                fclose(pFile);
        }
}

void Dynamic_Graph::FindPath(multimap<int, int>& ResInd2NodeID, vector<Node>& nodes, vector<vector<double> >& graph, vector<vector<bool> >& graph_indicator, const vector<vector<double> >& ovrlp, const vector < vector <pair<string,string> > >& atmnames ) {
        vector<vector<int> > CollisionPair;
        CollisionPair.resize(nodes.size());
        for(int i = 0; i < CollisionPair.size(); ++i) {
                CollisionPair[i].resize(nodes.size());
                for(int j = 0; j < CollisionPair[i].size(); ++j) {
                        CollisionPair[i][j] = 0;
                }
        }
        for(int i = 0; i < graph.size(); ++i) {
                for(int j = i + 1; j < graph.size(); ++j) {
                        if(graph_indicator[i][j] && strcmp(nodes[i].altLoc, "") && strcmp(nodes[j].altLoc, "")) {
                                CollisionPair[i][j] = 1;
                                CollisionPair[j][i] = 1;
                        }
                }
        }
        printf("Begin Path Finder for %s\n", file.c_str());
        for(int i = 0; i < graph.size(); ++i) {
                for(int j = 0; j < graph.size(); ++j) {
                        if(graph_indicator[i][j] && strcmp(nodes[i].altLoc, "") && strcmp(nodes[j].altLoc, "")) {
                                pair<multimap<int, int>::iterator, multimap<int, int>::iterator> ret;
                                ret = ResInd2NodeID.equal_range(nodes[i].ResInd);
                                //bool flag = false;
                                bool flag = true; //Bug fix? Next few lines are unnecessary, and discard valid alternates.
                                for(multimap<int, int>::iterator it = ret.first; it != ret.second; ++it) {
                                        if(/*!graph_indicator[(*it).second][j]*/ ovrlp[(*it).second][j] > gPRCNTL) {
                                                flag = true;
                                                break;
                                        }
                                }  
                                if(flag) {                              
                                        vector<pair<int, int> > path;
					vector<double> relief;
                                        path.push_back(make_pair(i, j));
                                        BuildPath(nodes, ResInd2NodeID, CollisionPair, graph_indicator, graph, path, ovrlp, relief, atmnames );
                                }        
                        }
                }
        }
}

bool Dynamic_Graph::SetParameter(string outfile, double rat, int len, string isSC) {
    file = outfile;
    ratio = rat;
    length = len;
    if(isSC == "T") {
        isSideChain = true;
    }
    else {
        isSideChain = false;
    }
    return true;
}

bool isSameRes(PCResidue res1, PCResidue res2) {
    //return res1->GetModelNum() == res2->GetModelNum() && string(res1->GetChainID()) == string(res2->GetChainID()) && res1->GetSeqNum() == res2->GetSeqNum() && string(res1->GetInsCode()) == string(res2->GetInsCode());
    return res1 == res2;
}

bool Dynamic_Graph::ReadFile(string filename) {
	int RC,lcount;
	char S[500];
	InitMatType();
	gNpaths = 0;
	gNfailedpaths = 0;
	PCMMDBManager myMMDBManager = NULL;
	myMMDBManager = new CMMDBManager();
	myMMDBManager->SetFlag(MMDBF_PrintCIFWarnings);
  	RC = myMMDBManager->ReadPDBASCII(filename.c_str());
	if(RC) {
		printf(" ***** ERROR #%i READ:\n\n %s\n\n",RC,GetErrorDescription(RC));
		myMMDBManager->GetInputBuffer(S,lcount);
    		if (lcount >= 0) 
      			printf("       LINE #%i:\n%s\n\n", lcount, S);
    		else if (lcount == -1)
      			printf("       CIF ITEM: %s\n\n", S);
    		delete myMMDBManager;	
    		return false;
	}
	myPCMMDBManager = myMMDBManager;
	return true;
}

bool Dynamic_Graph::RunChains(int MolNo) {
    if(MolNo < 1 || MolNo > myPCMMDBManager->GetNumberOfModels()) {
        printf("MolNo should be between 1 and %d\n", myPCMMDBManager->GetNumberOfModels());
        return false;
    }
    uddHandle_relief = myPCMMDBManager->RegisterUDInteger ( UDR_RESIDUE, "relief"   );
    uddHandle_norelief = myPCMMDBManager->RegisterUDInteger ( UDR_RESIDUE, "norelief"   );
    PCModel Mod = myPCMMDBManager->GetModel(MolNo);
    int NumOfChains = Mod->GetNumberOfChains();
    remove(file.c_str());
    printf("NUMBER OF CHAINS = %d\n", NumOfChains);
    for(int i = 0; i < NumOfChains; ++i) {
        printf("CHAINNO = %d\n", i);
        PCChain Chain = Mod->GetChain(i);
        ChainID ChID;
        strcpy(ChID, string(Chain->GetChainID()).c_str());
        CreateGraph(MolNo, ChID);
    }
    ProProcess();
}

bool Dynamic_Graph::CreateGraph(int MolNo, ChainID ChID) {
    PCModel pCModel = myPCMMDBManager->GetModel(MolNo);
    int SelHnd = myPCMMDBManager->NewSelection();
    myPCMMDBManager->Select(
        SelHnd,
        STYPE_RESIDUE,
        MolNo,
        ChID,
        ANY_RES, "*",
        ANY_RES, "*",
        "!HOH",
        "*",
        "*",
        "*",
        SKEY_NEW
                    );
    PPCResidue pPCRes;

    int nPCRes;
    myPCMMDBManager->GetSelIndex(SelHnd, pPCRes, nPCRes);
    if ( nPCRes < 2 )
      return false;
    pcresidues.resize(nPCRes);
    for(int i = 0; i < nPCRes; ++i) {
        pcresidues[i] = pPCRes[i];
    }
    printf("Number of Residues = %d\n", nPCRes);
    
    std::vector<Node> nodes;
    string notMCstr = "!N,CA,C,O,H,HA";
    for(int i = 0; i < nPCRes; ++i) {
	pcresidues[i]->PutUDData ( uddHandle_norelief, 0 );
	pcresidues[i]->PutUDData ( uddHandle_relief, 0 );
        PPCAtom AtomTable;
        int NumOfAtoms;
        pcresidues[i]->GetAtomTable(AtomTable, NumOfAtoms);
        set<string> locs;
        for(int j = 0; j < NumOfAtoms; ++j) {
            if ( isSideChain && AtomTable[j]->CheckID(" H  ", " H", "*") ) continue; // Main-chain H can be sole atom in multiple confs. Creates multiple nodes for this residue. Problem for side-chain only.
            string loc = string(AtomTable[j]->altLoc);
            locs.insert(loc);
        }
        if(locs.size() > 1) {
            locs.erase("");
        }
        for(set<string>::iterator ite = locs.begin(); ite != locs.end(); ++ite) {
            Node TmpNode;
            TmpNode.ResInd = i;
            strcpy(TmpNode.altLoc, (*ite).c_str());  
            nodes.push_back(TmpNode);          
        }
    }
    for(int i = 0; i < nodes.size(); ++i) {
        int SelHnd = myPCMMDBManager->NewSelection();
        myPCMMDBManager->Select(
            SelHnd,
            STYPE_ATOM,
            pcresidues[nodes[i].ResInd]->GetModelNum(),
            pcresidues[nodes[i].ResInd]->GetChainID(),
            pcresidues[nodes[i].ResInd]->GetSeqNum(), "*",
            pcresidues[nodes[i].ResInd]->GetSeqNum(), "*",
            "*",
            isSideChain ? notMCstr.c_str() : "*",
            "*",
            (string(",") + string(nodes[i].altLoc)).c_str(),
            SKEY_NEW
                    );
        PPCAtom ppCAtom = NULL;
        int atomNo;
        myPCMMDBManager->GetSelIndex(SelHnd, ppCAtom, atomNo);

        nodes[i].atoms = ppCAtom;
        nodes[i].atomNo = atomNo;
    }
    multimap<int, int> ResInd2NodeID;
    for(int i = 0; i < nodes.size(); ++i) {
        ResInd2NodeID.insert(make_pair(nodes[i].ResInd, i));
    }

    map<string, double> myVdWaalsRadius;
    for(int i = 0; i < nElementNames; ++i) {
        myVdWaalsRadius.insert(make_pair(string(ElementName[i]), VdWaalsRadius[i]));
    }
    myVdWaalsRadius[" H"] = 1.17;
    myVdWaalsRadius[" O"] = 1.40;
    myVdWaalsRadius[" C"] = 1.75;
    std::vector<vector<double> > graph;
    std::vector<vector<double> > overlaps (nodes.size());
    std::vector<vector<bool> > graph_indicator;
    vector<vector<pair<string,string> > > atomnames;
    graph.resize(nodes.size());
    graph_indicator.resize(nodes.size());
    atomnames.resize ( nodes.size());
    for(int i = 0; i < graph.size(); ++i) {
        graph[i].resize(nodes.size());
        graph_indicator[i].resize(nodes.size());
	overlaps[i].resize(nodes.size());
        atomnames[i].resize ( nodes.size());
    }
    for(int i = 0; i < nodes.size(); ++i) {
        for(int j = 0; j < nodes.size(); ++j) {
            graph[i][j] = DBL_MAX;
            graph_indicator[i][j] = false;
        }
    }
    char S1[128], S2[128];
    for(int i = 0; i < nodes.size(); ++i) {
        for(int j = i + 1; j < nodes.size(); ++j) {
            if(nodes[i].ResInd != nodes[j].ResInd) {
                vector<double> atomDist;
		vector<double> alloverlaps;
                vector<pair<string,string> > atomName;
		        bool flag = false;
                for(int ii = 0; ii < nodes[i].atomNo; ++ii) {
                    for(int jj = 0; jj < nodes[j].atomNo; ++jj) {
                        bool CS;
                        if(isSideChain) {
                            CS = !isMChain(nodes[i].atoms[ii]) && !isMChain(nodes[j].atoms[jj]);
                        }
                        else {
                            //CS = (abs(double(pcresidues[nodes[i].ResInd]->GetSeqNum()) - double(pcresidues[nodes[j].ResInd]->GetSeqNum())) < 1.5 && (!isMChain(nodes[i].atoms[ii]) || !isMChain(nodes[j].atoms[jj]))) || (abs(double(pcresidues[nodes[i].ResInd]->GetSeqNum()) - double(pcresidues[nodes[j].ResInd]->GetSeqNum())) > 1.5 && pcresidues[nodes[i].ResInd]->isMainchainHBond(pcresidues[nodes[j].ResInd]));
                            if(abs(double(pcresidues[nodes[i].ResInd]->GetSeqNum()) - double(pcresidues[nodes[j].ResInd]->GetSeqNum())) < 1.5) {
                                if(!isMChain(nodes[i].atoms[ii]) || !isMChain(nodes[j].atoms[jj])) {
                                    CS = true;
                                }
                                else CS = false;
                            }
                            else {
                                 /*if(pcresidues[nodes[i].ResInd]->isMainchainHBond(pcresidues[nodes[j].ResInd])) {
                                    if(!isMChain(nodes[i].atoms[ii]) || !isMChain(nodes[j].atoms[jj])) {
                                        CS = true;
                                    }
                                    else CS = false;
                                }
                                else */ CS = true;
                            }
                        }
			if ( string(pcresidues[nodes[i].ResInd]->GetResName()) == "CYS" && string(pcresidues[nodes[j].ResInd]->GetResName()) == "CYS" && string(nodes[i].atoms[ii]->name) == " SG " && string(nodes[j].atoms[jj]->name) == " SG " )
			    {
				//CS = false;
				std::cout << "Found SG" << std::endl;
			    }
                        if(CS) {
                        //if(!isMChain(nodes[i].atoms[ii]) && !isMChain(nodes[j].atoms[jj])) {
                        	double rii, rjj;
                        	string resNameii, resNamejj;
                        	string atomNameii, atomNamejj;
                        	resNameii = string(nodes[i].atoms[ii]->GetResName());
                        	resNamejj = string(nodes[j].atoms[jj]->GetResName());
                        	atomNameii = string(nodes[i].atoms[ii]->name);
                        	atomNamejj = string(nodes[j].atoms[jj]->name);
                        	if(resNameii == "PHE" && (atomNameii == " HD1" || atomNameii == " HE1" || atomNameii == " HZ " || atomNameii ==  " HE2" || atomNameii == " HD2")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "TRP" && (atomNameii == " HD1" || atomNameii == " HE1" || atomNameii == " HZ2" || atomNameii ==  " HH2" || atomNameii == " HZ3" || atomNameii == " HE3")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "TYR" && (atomNameii == " HD1" || atomNameii == " HE1" || atomNameii ==  " HE2" || atomNameii == " HD2")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "ASN" && (atomNameii == "HD21" || atomNameii == "HD22")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "GLN" && (atomNameii == "HE21" || atomNameii == "HE22")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "LYS" && (atomNameii == " HZ1" || atomNameii == " HZ1" || atomNameii == " HZ3")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "ARG" && (atomNameii == " HE " || atomNameii == "HH11" || atomNameii == "HH12" || atomNameii ==  "HH21" || atomNameii == "HH22")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "SER" && (atomNameii == " HG ")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "CYS" && (atomNameii == " HG ")) {
                        		rii = 1.0;
                        	}
                        	else if(resNameii == "THR" && (atomNameii == " HG1")) {
                        		rii = 1.0;
                        	}
                        	else {
                        		if(atomNameii == " C  ") {
                        			rii = 1.65;
                        		}
                        		else if(atomNameii == " H  ") {
                        			rii = 1.00;
                        		}
                        		else {
                        			rii = myVdWaalsRadius[string(nodes[i].atoms[ii]->element)]; 
                        		}
                        	}
                        	
                        	if(resNamejj == "PHE" && (atomNamejj == " HD1" || atomNamejj == " HE1" || atomNamejj == " HZ " || atomNamejj ==  " HE2" || atomNamejj == " HD2")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "TRP" && (atomNamejj == " HD1" || atomNamejj == " HE1" || atomNamejj == " HZ2" || atomNamejj ==  " HH2" || atomNamejj == " HZ3" || atomNamejj == " HE3")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "TYR" && (atomNamejj == " HD1" || atomNamejj == " HE1" || atomNamejj ==  " HE2" || atomNamejj == " HD2")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "ASN" && (atomNamejj == "HD21" || atomNamejj == "HD22")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "GLN" && (atomNamejj == "HE21" || atomNamejj == "HE22")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "LYS" && (atomNamejj == " HZ1" || atomNamejj == " HZ1" || atomNamejj == " HZ3")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "ARG" && (atomNamejj == " HE " || atomNamejj == "HH11" || atomNamejj == "HH12" || atomNamejj ==  "HH21" || atomNamejj == "HH22")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "SER" && (atomNamejj == " HG ")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "CYS" && (atomNamejj == " HG ")) {
                        		rjj = 1.0;
                        	}
                        	else if(resNamejj == "THR" && (atomNamejj == " HG1")) {
                        		rjj = 1.0;
                        	}
                        	else {
                        		if(atomNamejj == " C  ") {
                        			rjj = 1.65;
                        		}
                        		else if(atomNamejj == " H  ") {
                        			rjj = 1.00;
                        		}
                        		else {
                        			rjj = myVdWaalsRadius[string(nodes[j].atoms[jj]->element)]; 
                        		}
                        	}
                            double leastDist = rii + rjj;
                            double tmpatomDist = norm(nodes[i].atoms[ii]->x, nodes[i].atoms[ii]->y, nodes[i].atoms[ii]->z, 
                                                nodes[j].atoms[jj]->x, nodes[j].atoms[jj]->y, nodes[j].atoms[jj]->z);  
                            atomDist.push_back(tmpatomDist);
                            atomName.push_back ( make_pair ( atomNameii, atomNamejj ) );
			    alloverlaps.push_back ( tmpatomDist / leastDist );
                            //std::printf("%lf\n", ratio);
                            if(tmpatomDist < 1.0001 * leastDist) {
    			                flag = true;			
    			            } 
			            }
	                }
                }
                if(atomDist.size() != 0) graph[i][j] = graph[j][i] = *min_element(atomDist.begin(), atomDist.end());
		if(atomDist.size() != 0) {
		  overlaps[i][j] = overlaps[j][i] = *min_element(alloverlaps.begin(), alloverlaps.end());
                  vector<double>::iterator dit = min_element(alloverlaps.begin(), alloverlaps.end());
                  atomnames[i][j] = atomnames[j][i] = atomName[dit-alloverlaps.begin()];
                }
		graph_indicator[i][j] = graph_indicator[j][i] = flag;
            }
        }
    }

    std::vector<double> clash_ratios;
    for ( int i=0; i<overlaps.size(); ++i )
      for ( int j=i+1; j<overlaps.size(); ++j )
      {
        if ( graph_indicator[i][j] && overlaps[i][j] < 1.0001 )
        {
          //std::cout << overlaps[i][j] << std::endl;
          clash_ratios.push_back ( overlaps[i][j] );
        }
      }
    
    //cout << "Size of clash_ratios: " << clash_ratios.size() << endl;
    double clash_prcntl = getPercentile ( clash_ratios, ratio );
    gPRCNTL = clash_prcntl;

    for ( int i=0; i<overlaps.size(); ++i )
      for ( int j=i+1; j<overlaps.size(); ++j )
        if ( overlaps[i][j] >  clash_prcntl )
          graph_indicator[i][j] = graph_indicator[j][i] = false;

FindPath(ResInd2NodeID, nodes, graph, graph_indicator, overlaps, atomnames );
    //ProProcess();

  
   ofstream WriteNorelief("relief_by_residue.txt");
   int totalpaths = 0;
   int totalnopaths = 0;
   WriteNorelief << "residue,no_relief,relief" << endl;
   for(int i = 0; i < nPCRes; ++i) {
     int N_norelief, N_relief;
     pcresidues[i]->GetUDData (uddHandle_norelief, N_norelief );
     pcresidues[i]->GetUDData (uddHandle_relief, N_relief );
     WriteNorelief << pcresidues[i]->GetSeqNum() << "," <<  double(N_norelief)/double(gNfailedpaths) << "," << double(N_relief)/double(gNpaths) << endl;
     totalnopaths += N_norelief;
     totalpaths += N_relief;
   }
   WriteNorelief.close();
   cout << "Total Paths: " << totalnopaths << " " << totalpaths << endl;
   return true;
}

int isContained(vector<string>& str1, vector<string>& str2) {
	bool result = true;

	if(str1.size() >= str2.size()) {
		for(int i = 0; i < str1.size() - str2.size() + 1; ++i) {
			if(str1[i] == str2[0]) {
				bool flag = true;
				for(int j = 0; j < str2.size(); ++j) {
					if(str1[i + j] != str2[j]) {
						flag = false;
						break;
					}
				}
				if(flag) return 1;
			}
		}
		return 0;
	}
	else {
		for(int i = 0; i < str2.size() - str1.size() + 1; ++i) {
			if(str2[i] == str1[0]) {
				bool flag = true;
				for(int j = 0; j < str1.size(); ++j) {
					if(str2[i + j] != str1[j]) {
						flag = false;
						break;
					}
				}
				if(flag) return -1;
			}
		}
		return 0;
	}	
}
#include <sstream>
void Dynamic_Graph::ProProcess() {
    cout << "Begin Post Process:\n";
	list<vector<string> > PathTable;
	ifstream PathFile(file.c_str());
	if(!PathFile.is_open()) {
		cout << "Cannot open file!\n";
		return;
	}
	string tmpstr2;
	vector<string> OnePath;
	while(getline(PathFile, tmpstr2)) {
		string tmpstr, tokens;
                std::istringstream tstr(tmpstr2);
		for (int i = 0; i < 3; ++i) {
			  tstr >> tokens;
                          tmpstr += tokens;
			  tmpstr += " ";
		}
		if(tmpstr2 == string("")) {
			if(PathTable.size() == 0) {
				PathTable.push_back(OnePath);
			}
			else {
				list<vector<string> >::iterator it = PathTable.begin();
				bool flag = true;
				while(it != PathTable.end()) {
					int sign = isContained(*it, OnePath);
					if(sign == -1) {
					    list<vector<string> >::iterator tmp = it;
					    ++it;
						PathTable.erase(tmp);
					}
					else if(sign == 1) {
						flag = false;
						break;
					}
					else {
						++it;
					}
				}
				if(flag) {
					PathTable.push_back(OnePath);
				}
			}
			OnePath.clear();
		}
		else if(tmpstr != "PATH(ModelNo, ChainID, SeqNum(InsertionCode), ") {
			OnePath.push_back(tmpstr);
		}
	}
	PathFile.close();
	cout << gNpaths << " " << gNfailedpaths << endl;
	ofstream WritePath((file+"RSDS_ONLY").c_str());
	cout << "Number of Paths = " << PathTable.size() << "\n\n\n";
	for(list<vector<string> >::iterator it = PathTable.begin(); it != PathTable.end(); ++it) {
		WritePath << "PATH(ModelNo, ChainID, SeqNum(InsertionCode), AltLoc):" << "\n";
		for(int j = 0; j < (*it).size(); ++j) {
			WritePath << (*it)[j] << "\n";
		}
		WritePath << "\n";
	}
	WritePath.close();
}
