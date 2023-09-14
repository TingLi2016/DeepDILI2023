#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "Elements.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Calculator1.h"

bitset<2000> Calculator1::used;
Calculator1::Calculator1(Molecule & _mol)
:mol(_mol){	
	RBF = SCBOHD = NMB = NRB = NDB = NAB = 0;
	numAtoms = mol.getAtomsNum(false);
}

Calculator1::~Calculator1(void){}

void Calculator1::AromaticBondConverter(){
	try{
		if(mol.isFeaturesReady(AromaticBondConvert)) return;

		vector<double> fingerPrint(18, -999); 
		vector<short> nr(13, 0); // nr[i] is the number of rings with i (i > 2) atoms 
		nr[0] = mol.getBondsNum(false) - numAtoms + 1;	

		mol.setRingStatus(true);
		// process rings: find all base cycles
		if(nr[0] > 0) ringSearch(); 

		// Convert to Aromatic Bonds
		convertAromaticBonds(nr);

		/* number of 3-membered ~ 12-membered rings printout to output file */
		for (short i = 2; i <= 12; i++) 
			fingerPrint[i - 2] = nr.at(i);

		// Change Bond Topology
		changeBondTopology();

		// Find NCRS
		fingerPrint[12] = (double)findNCRS();

		// Find NMB  NRB  RBF  NDB		
		findBondProperties();
		fingerPrint[11] = (double)NMB;		
		fingerPrint[13] = (double)NRB;
		fingerPrint[14] = (double)RBF;
		fingerPrint[15] = (double)NDB;
		fingerPrint[16] = (double)NAB;
		fingerPrint[17] = (double)SCBOHD;
		mol.setFingerprint(fingerPrint, AromaticBondConvert);
	}catch (const exception & e){
		cout << e.what() << " thrown from Calculator1" << endl;
	}
}

// Fuse all base cycles together, and get total number of rings
// This version can not give right result for rings that share more than three atoms
int Calculator1::findNCRS(){
	vector< vector<short> > & rings = mol.getRings();
	int NCRS = rings.size();	
	if(NCRS > 1){		
		vector< set<short> > currentLevel;
		vector< vector<char> > currentP;
		vector< vector<char> > nextP;
		vector<char> twoP(NCRS, 0);

		for(short i = rings.size() - 1; i >= 0; i--){
			currentLevel.push_back(set<short>(rings[i].begin(), rings[i].end()));
			currentP.push_back(vector<char>(NCRS, 0));
			currentP.back()[i] = 1;
		}

		vector< set<short> > nextLevel;
		short level = 0;

		while(currentLevel.size() > 0){
			short nr = currentLevel.size();	
			level++;
					
			nextLevel.clear();
			nextP.clear();

			for(short i = 0; i < nr; i++){
				for(short j = i + 1; j < nr; j++){
					set<short> two(currentLevel[i].begin(), currentLevel[i].end());
					two.insert(currentLevel[j].begin(), currentLevel[j].end());
					int cL = 0;
					for(short m = twoP.size() - 1; m >= 0; m--){
						twoP[m] = currentP[i][m] > 0 || currentP[j][m] > 0 ? 1:0;
						cL += twoP[m];
					}
					if(cL > level + 1) continue; // not a next level, skip

					if(currentLevel[i].size() + currentLevel[j].size() - two.size() > 1){ //they share at least two atoms
						// check existent of two
						bool find = false;
						for(short n = nextLevel.size() - 1; n >= 0; n--){
							if(two.size() == nextLevel[n].size()){
								bool isSame = true;
								for(short m = twoP.size() - 1; m >= 0; m--){
									if(twoP[m] != nextP[n][m]){
										isSame = false;
										break;
									}
								}
								if(isSame){
									find = true;
									break;
								}
							}
						}	
						if(!find){
							nextLevel.push_back(two);	
							nextP.push_back(twoP);							
						}
					}
				}
			}

			NCRS += nextLevel.size();
			swap(currentLevel, nextLevel);
			swap(currentP, nextP);			
			// Number of circus. If there are too many circus, MAX_RING_FUSE - 1 is assigned.
			if(NCRS >= MAX_RING_FUSE - 1){
				NCRS = MAX_RING_FUSE - 1;
				break;
			}
		}		
	}
	return NCRS;
}
/*

pair<short, short> getRight(vector<short> & ring, short idx){
	short idx1 = idx + 1;
	if(idx1 >= (short)ring.size()){
		return pair<short, short>(ring[0], 0);
	}else{
		return pair<short, short>(ring[idx1], idx1);
	}
}

pair<short, short> getLeft(vector<short> & ring, short idx){
	short idx1 = idx - 1;
	if(idx1 < 0){
		return pair<short, short>(ring.back(), ring.size() - 1);
	}else{
		return pair<short, short>(ring[idx1], idx1);
	}
}

void fuseTwo(vector<short> & r1, vector<short> & r2, vector<short> & fused){
	// find a shared atom
	for(short m = r1.size() - 1; m >= 0; m--){
		cout << r1[m] << " ";
	}
	cout  <<  endl;
	
	for(short m = r2.size() - 1; m >= 0; m--){
		cout << r2[m] << " ";
	}
	cout << endl;

	short i, j;	
	for(i = r1.size() - 1; i >= 0; i--){
		bool find = false;
		for(j = r2.size() - 1; j >= 0; j--){			
			if(r1[i] == r2[j]){ // find one
				find = true;
				break;	
			}
		}
		if(find) break;
	}
	
	// extend i and j, find the border atoms
	pair<short, short> r1_L = getLeft(r1, i);
	pair<short, short> r2_L = getLeft(r2, j);

	pair<short, short> r2_R = getRight(r2, j);
	pair<short, short> r1_R = getRight(r1, i);

	//case 1:
	if(r1_L.first == r2_L.first || r1_R.first == r2_R.first){
		short step1 = 0;
		while(r1_L.first == r2_L.first){
			step1++;
			r1_L = getLeft(r1, r1_L.second);
			r2_L = getLeft(r2, r2_L.second);
		}
		r1_L = getRight(r1, r1_L.second);
		short step2 = 0;
		while(r1_R.first == r2_R.first){
			step2++;
			r1_R = getRight(r1, r1_R.second);
			r2_R = getRight(r2, r2_R.second);
		}
	
		short nr1 = r1.size() - step1 - step2 + 1;
		for(short m = 0; m < nr1; m++){
			fused.push_back(r1_L.first);
			r1_L = getLeft(r1, r1_L.second);
		}

		short nr2 = r2.size() - step1 - step2 - 1;
		for(short m = 0; m < nr2; m++){
			fused.push_back(r2_R.first);
			r2_R = getRight(r2, r2_R.second);
		}
		set<short> s1(fused.begin(), fused.end());
		if(s1.size() < fused.size()){
			short xx = 1;
		}
		for(short m = fused.size() - 1; m >= 0; m--){
			cout << fused[m] << " ";
		}
		cout << endl;
		cout << endl;
		return;
	}
	//case 2:
	if(r1_L.first == r2_R.first || r1_R.first == r2_L.first){
		short step1 = 0;
		while(r1_L.first == r2_R.first){
			step1++;
			r1_L = getLeft(r1, r1_L.second);
			r2_R = getRight(r2, r2_R.second);
		}
		r1_L = getRight(r1, r1_L.second);		
		short step2 = 0;
		while(r1_R.first == r2_L.first){
			step2++;
			r1_R = getRight(r1, r1_R.second);
			r2_L = getLeft(r2, r2_L.second);
		}

		short nr1 = r1.size() - step1 - step2 + 1;
		for(short m = 0; m < nr1; m++){
			fused.push_back(r1_L.first);
			r1_L = getLeft(r1, r1_L.second);
		}

		short nr2 = r2.size() - step1 - step2 - 1;
		for(short m = 0; m < nr2; m++){
			fused.push_back(r2_L.first);
			r2_L = getLeft(r2, r2_L.second);
		}

		set<short> s1(fused.begin(), fused.end());
		if(s1.size() < fused.size()){
			short xx = 1;
		}
		for(short m = fused.size() - 1; m >= 0; m--){
			cout << fused[m] << " ";
		}
		cout << endl;
		return;
	}		
}

int Calculator1::findNCRS(){
	vector< vector<short> > & rings = mol.getRings();
	int NCRS = rings.size();	
	if(NCRS > 1){		
		vector< vector<short> > currentLevel;
		vector< vector<char> > currentP;
		vector< vector<char> > nextP;
		vector<char> twoP(NCRS, 0);

		for(short i = rings.size() - 1; i >= 0; i--){
			currentLevel.push_back(vector<short>(rings[i].begin(), rings[i].end()));
			currentP.push_back(vector<char>(NCRS, 0));
			currentP.back()[i] = 1;
		}

		vector< vector<short> > nextLevel;
		short level = 0;

		while(currentLevel.size() > 0){
			short nr = currentLevel.size();	
			level++;
					
			nextLevel.clear();
			nextP.clear();

			for(short i = 0; i < nr; i++){
				for(short j = i + 1; j < nr; j++){
					set<short> two(currentLevel[i].begin(), currentLevel[i].end());
					two.insert(currentLevel[j].begin(), currentLevel[j].end());
					int cL = 0;
					for(short m = twoP.size() - 1; m >= 0; m--){
						twoP[m] = currentP[i][m] > 0 || currentP[j][m] > 0 ? 1:0;
						cL += twoP[m];
					}
					if(cL > level + 1) continue; // not a next level, skip

					if(currentLevel[i].size() + currentLevel[j].size() - two.size() > 1){ //they share at least two atoms
						// check existent of two
						bool find = false;
						for(short n = nextLevel.size() - 1; n >= 0; n--){
							if(two.size() == nextLevel[n].size()){
								bool isSame = true;
								for(short m = twoP.size() - 1; m >= 0; m--){
									if(twoP[m] != nextP[n][m]){
										isSame = false;
										break;
									}
								}
								if(isSame){
									find = true;
									break;
								}
							}
						}	
						if(!find){
							vector<short> fused;							
							fuseTwo(currentLevel[i], currentLevel[j], fused);
							short idx = -1;
							if(!isAlreadyFound(currentLevel, fused, idx)){
								nextLevel.push_back(fused);	
								nextP.push_back(twoP);							
							}
						}
					}
				}
			}

			NCRS += nextLevel.size();
			swap(currentLevel, nextLevel);
			swap(currentP, nextP);			
			// Number of circus. If there are too many circus, MAX_RING_FUSE - 1 is assigned.
			if(NCRS >= MAX_RING_FUSE - 1){
				NCRS = MAX_RING_FUSE - 1;
				break;
			}
		}		
	}
	return NCRS;
}
*/
void Calculator1::findBondProperties(){	
	short bondNum = mol.getBondsNum(false);
	for(short i = 0; i < bondNum; i++){
		char bondType = mol.getBond(i).getBondType();
		if(mol.getBond(i).getBondTopology() == Bond::CHAIN_BOND && bondType == Bond::SINGLE_BOND) NRB++;
		if(bondType == Bond::DOUBLE_BOND){
			NDB++;
			SCBOHD += 2.0;
		}
		if(bondType == Bond::SINGLE_BOND) SCBOHD += 1.0;
		if(bondType == Bond::TRIPLE_BOND) SCBOHD += 3.0;
		if(bondType == Bond::AROMATIC_BOND){
			SCBOHD += 1.5;
			NAB++;
		}
		if(bondType != Bond::SINGLE_BOND) NMB++;		
	}

	if(mol.getBondsNum(false) != 0)
		RBF = double(NRB) / mol.getBondsNum(false);
}

// Change Bond Topology
void Calculator1::changeBondTopology(){
	short bondNum = mol.getBondsNum(false);
	vector< vector<short> > & rings = mol.getRings();
	short nr = rings.size();
	for(short i = 0; i < bondNum; i++) 
		mol.getBond(i).setBondTopology(Bond::CHAIN_BOND);

	bool btir1, btir2;
	for(short i = 0; i < bondNum; i++){
		bool inr = false;		
		short a1 = mol.getBond(i).getFirstAtomNumber();
		short a2 = mol.getBond(i).getSecondAtomNumber();
		for (short j = 0; j < nr; j++){
			btir1 = btir2 = false;
			for(short m = rings[j].size() - 1; m >= 0; m--){
				if (rings[j][m] == a1) btir1 = true;
				else if (rings[j][m] == a2) btir2 = true;
				if(btir1 && btir2){
					inr = true;
					break;
				}
			}
			if (inr) break;
		}
		if (inr) mol.getBond(i).setBondTopology(1);
	}
}

// Convert to Aromatic Bonds
void Calculator1::convertAromaticBonds(vector<short> & nr){
	vector< vector<short> > & rings = mol.getRings();
	for (short i = rings.size() - 1; i >= 0; i--){
		short nri = rings[i].size();
		if(nri >= 2 && nri <= 12)
			nr[nri]++;
	}

	// Convert to Aromatic Bonds
	vector< vector<char> > & connBonds = mol.getConnectedBondTypes(false); // original bond types
	bool bondtypeModified = false;
	for (short i = rings.size() - 1; i >= 0; i--){
		short nri = rings[i].size();
		if(nri == 6){
			bool r = true;
			for(short j = 0; j < nri; j++){
				bool btir = false;
				for(short k = connBonds[rings[i][j]].size() - 1; k >= 0; k--){
					if(connBonds[rings[i][j]][k] == Bond::DOUBLE_BOND){
						btir = true;
						break;
					}
				}
				if (!btir){
					r = false;
					break;
				}
			}

			if(r){	
				for(short j = 0; j < nri; j++){
					if(mol.getAtom(rings[i][j]).getAtomSymbol() != Elements::CARBON){
						r = false;
						break;
					}
				}
			}

			if(r){
				short bondNum = mol.getBondsNum(false);
				bool trace_1, trace_2;
				for(short j = 0; j < bondNum; j++){
					trace_1 = trace_2 = false;
					short a1 = mol.getBond(j).getFirstAtomNumber();
					short a2 = mol.getBond(j).getSecondAtomNumber();
					for(short m = 0; m < nri; m++){
						if(a1 == rings[i][m]) trace_1 = true;
						else if(a2 == rings[i][m]) trace_2 = true;
						if(trace_1 && trace_2) break;
					}
					if(trace_1 && trace_2){
						mol.getBond(j).setBondType(Bond::AROMATIC_BOND);
						bondtypeModified = true;
					}
				}
				nr[2]++;
			}
		}
	}

	if(bondtypeModified) // if bond types has been modified here, connectedBondTypes must be reset;
		connBonds.clear(); 
}

void Calculator1::removeUnrelatedAtoms(){
	bool changed = true;
	while(changed){
		changed = false;
		for(short i = connAtoms.size() - 1; i >= 0; i--){
			if(connAtoms[i].size() == 1){
				for(short j = connAtoms[connAtoms[i][0]].size() - 1; j >= 0; j--){
					if(connAtoms[connAtoms[i][0]][j] == i){
						connAtoms[connAtoms[i][0]].erase(connAtoms[connAtoms[i][0]].begin() + j);
						break;
					}
				}
				connAtoms[i].clear();
				changed = true;
			}
		}		
	}
}

// find all back edges using DFS algorithm
void Calculator1::dfsBackEdges( void ){
	mark.assign(numAtoms, WHITE);
	short current = -1;	
	stack<short> myQueue;

	for(short i = 0; i < numAtoms; i++){
		if(connAtoms[i].size() > 1){
			myQueue.push(i);
			break;
		}
	}

	while(!myQueue.empty()){
		current = myQueue.top();
		mark[current] = GRAY;		
		myQueue.pop();
		// check end node in the path
		for(short j = connAtoms[current].size() - 1; j >= 0; j--){
			switch(mark[connAtoms[current][j]]){
			case WHITE:
				myQueue.push(connAtoms[current][j]);					
				mark[connAtoms[current][j]] = GRAY;
				break;
			case GRAY:
				edges.push_back(pair<short, short>(current, connAtoms[current][j]));
				break;
			}
		}
		mark[current] = BLACK;
	}
}

// Find all shortest paths from back edge
void Calculator1::bfsShortestRing(pair<short, short> & edge, vector< vector<short> > &paths){
	paths.clear();
	deque< short > Queue;	
	p.assign(numAtoms, -1);
	mark.assign(numAtoms, WHITE);
	p[edge.first] = edge.first;
	mark[edge.first] = GRAY;
	len.assign(numAtoms, -1);
	short time = 0;

	for(short j = connAtoms[edge.first].size() - 1; j >= 0; j--){
		if(connAtoms[edge.first][j] != edge.second){
			Queue.push_back(connAtoms[edge.first][j]);
			mark[connAtoms[edge.first][j]] = GRAY;
			p[connAtoms[edge.first][j]] = edge.first;
			len[connAtoms[edge.first][j]] = time;
		}
	}
	
	short current = -1;	
	short findTime = -1;
	while(!Queue.empty()){		
		current = Queue.front();
		
		if(findTime > -1 && len[current] >= findTime) break;
		mark[current] = GRAY;
		Queue.pop_front();	
		time = len[current] + 1;
		// check end node in the path			
		for(short j = connAtoms[current].size() - 1; j >= 0; j--){
			if(mark[connAtoms[current][j]] == WHITE){
				if(connAtoms[current][j] == edge.second) {// find the shorest path
					paths.push_back(vector<short> () );
					vector<short> &path = paths.back();
					path.push_back(edge.second);				
					while(current != edge.first){
						path.push_back(current);					
						current = p[current];			
					}
					path.push_back(edge.first);
									
					if(findTime < 0) findTime = time;
					break;
				}else{
					Queue.push_back(connAtoms[current][j]);
					p[connAtoms[current][j]] = current;
					mark[connAtoms[current][j]] = GRAY;
					len[connAtoms[current][j]] = time;
				}
			}
		}
	}
}

void Calculator1::sharedBackEdge(short eID, vector<short> &path, vector<short> & sharedEdges){
	sharedEdges.clear();
	for(short i = edges.size() - 1; i >= 0; i--){
		if(eID == i) continue;
		set<short> two(path.begin(), path.end());
		two.insert(edges[i].first);
		two.insert(edges[i].second);
		if(two.size() == path.size()){
			sharedEdges.push_back(i);
		}
	}
}

bool Calculator1::isAlreadyFound(vector<vector<short> > & rings, vector<short> & path, short & idx){
	for(short i = rings.size() - 1; i >= 0; i--){
		if(path.size() == rings[i].size()){
			set<short> two(path.begin(), path.end());
			two.insert(rings[i].begin(), rings[i].end());
			if(two.size() == path.size()){
				idx = i;
				return true;
			}
		}
	}
	return false;
}

inline void Calculator1::cleanUpeSet(vector< set<pair<short, short> > > & eSet){
	set<pair<short, short> > single;
	bool changed = true;
	while(changed){
		changed = false;
		for(short i = eSet.size() - 1; i >= 0; i--){
			if(eSet[i].size() == 1){
				single.insert(eSet[i].begin(), eSet[i].end());
				eSet[i].clear();
				changed = true;
			}else{
				vector< set<pair<short, short> >::iterator > remove;
				for(set<pair<short, short> >::iterator it = eSet[i].begin(); it != eSet[i].end(); it++){
					if(single.find(*it) != single.end()) remove.push_back(it);
				}
				for(short j = remove.size() - 1; j >= 0; j--) eSet[i].erase(remove[j]);
				changed = remove.size() > 0;
			}
		}
	}
}

inline bool Calculator1::brokenSearch(vector<vector<short> > & rings, pair<short, short> & e1, pair<short, short> & e2){
	// break e1 then search e2
	short bid = e1.first;
	if(e1.first == e2.first || e1.first == e2.second) bid = e1.second;
	vector<short> cbid = connAtoms[bid];
	connAtoms[bid].clear();

	vector< vector<short> > paths;
	bfsShortestRing(e2, paths);
	connAtoms[bid] = cbid;
	short idx;
	for(short j = paths.size() - 1; j >= 0; j--){
		if(!isAlreadyFound(rings, paths[j], idx)){
			rings.push_back(paths[j]);
			return true;
		}					
	}
	return false;
}

void Calculator1::ringSearch(){	
	try{
		vector< vector<short> > & rings = mol.getRings();
		size_t nr = mol.getBondsNum(false) - numAtoms + 1;
		connAtoms = mol.getConnectedAtomIDs(false);	

		//iterated delete all single connections
		removeUnrelatedAtoms();

		// find all back edges (each one means a cycle)
		dfsBackEdges();

		vector< set<pair<short, short> > > eSet;
		vector< vector<short> > paths;
		vector<short> s;
		for(short i = edges.size() - 1; i >= 0; i--){
			if(rings.size() == nr)
				break;
			bfsShortestRing(edges[i], paths);
			for(short j = paths.size() - 1; j >= 0; j--){			
				sharedBackEdge(i, paths[j], s);
				short idx = -1;
				if(!isAlreadyFound(rings, paths[j], idx)){
					idx = rings.size();
					if(idx == nr) break;
					rings.push_back(paths[j]);
					eSet.push_back(set<pair<short, short> >());
				}
				eSet[idx].insert(edges[i]);
				for(short m = s.size() - 1; m >= 0; m--)
					eSet[idx].insert(edges[s[m]]);
			}
		}

		if(rings.size() != nr){ // ido second step
			cleanUpeSet(eSet);
			for(short i = eSet.size() - 1; i >= 0; i--){
				if(nr == rings.size()) break;
				short es = eSet[i].size();
				if(es >= 2){
					set<pair<short, short> >::iterator it = eSet[i].begin();
					while(true){
						pair<short, short> e1 = *it;
						it++;
						if(it != eSet[i].end()){
							pair<short, short> e2 = *(it);
							brokenSearch(rings, e1, e2);
							if(nr > rings.size())
								brokenSearch(rings, e2, e1);							
						}else{
							break;
						}
						if(nr == rings.size())	break;
					}
				}
			}
			if(nr != rings.size())
				throw mold2Exception("Can not find all base cycles.");
		}
	}catch(exception & ex){
		cout << ex.what() << endl;
	}
}