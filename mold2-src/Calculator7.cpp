#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Constants.h"
#include "Calculator7.h"


Calculator7::Calculator7(Molecule & _mol):mol(_mol){
	numAtoms = _mol.getAtomsNum(false);	
	doAll();
}

Calculator7::~Calculator7(void){
}

void Calculator7::PathWalk(){
	try{
		if(mol.isFeaturesReady(::PathWalk)) return;
		if(mol.isAllPathFailed()){
			vector<double> tmp(4, ::numeric_limits<float>::signaling_NaN());
			mol.setFingerprint(tmp, ::PathWalk);
			return;
		}

		/* EXP of Path-distance / Walk-distance over all atoms: EXP */
		vector<double> fingerprint(4, 0);

		/* -------------------------------- paths 2~5 ------------------------------------- */
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector< vector<int> > & path = mol.getPath();
		for(short i = 0; i < numAtoms; i++){
			/* ---------------------------- EXP2 ----------------------------------------- */
			vector<short> aw1, aw2;
			aw1.assign(connAtoms[i].begin(), connAtoms[i].end());
			short na = connAtoms[i].size();
			if(path[0][i] + na > 0)
				fingerprint[0] += double(path[0][i]) / (path[0][i] + na);

			gotoNextLevel(aw1, aw2);
			/* ---------------------------- EXP3 ----------------------------------------- */
			gotoNextLevel(aw2, aw1);		
			short nai = aw1.size();
			if(nai > 0)
				fingerprint[1] += double(path[1][i]) / nai;

			/* ---------------------------- EXP4 ----------------------------------------- */
			gotoNextLevel(aw1, aw2);
			nai = aw2.size();
			if(nai > 0)
				fingerprint[2] += double(path[2][i]) / nai;

			/* ---------------------------- EXP5 ----------------------------------------- */
			gotoNextLevel(aw2, aw1);
			nai = aw1.size();		
			if(nai > 0)
				fingerprint[3] += double(path[3][i]) / nai;
		}

		/* ------------------------- calculate EXP2 ~ EXP5 -------------------------------- */
		fingerprint[0] /= numAtoms;
		fingerprint[1] /= numAtoms;
		fingerprint[2] /= numAtoms;
		fingerprint[3] /= numAtoms;		
		mol.setFingerprint(fingerprint, ::PathWalk);
		path.clear();
	}catch(exception &ex){
		cout << ex.what() << " thrown in PathWalk" << endl;
	}
}


void Calculator7::KierPathIndex(){
	try{
		if(mol.isFeaturesReady(::KierPathIndex)) return;
		if(mol.isAllPathFailed()) {
			vector<double> tmp(4, ::numeric_limits<float>::signaling_NaN());
			mol.setFingerprint(tmp, ::KierPathIndex);
			return;
		}

		/* atom2bonds[i][j]--number of bonds type j connected to the atom i */
		vector< vector<short> > atom2bonds(numAtoms, vector<short>(mol.m_maxBondType, 0));

		/* atomic Covalent Radius: 	i is order of atom
		cr[i]: sp-N -- cr[i] = 1: sp; cr[i] = 2: sp2; cr[i] = 3: sp3; */
		vector<short> cr(numAtoms, 0);
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		/* get type of bond that on the shortest path */
		DistancePathMatrix &DPM = mol.getDistancePathMatrix();
		for(short i = 0; i < numAtoms; i++){
			for(short j = i; j < numAtoms; j++){
				vector<short> &pij = DPM.getPath(i, j);
				for(short k = pij.size() - 1; k >= 0; k--){
					if(atom2bonds[pij[k]][0] != 0)	/* bonds of the atom count already -- does not count again */
						continue;
					else{	/* bonds of the atoms does not count before */
						short nAtoms = connAtoms[pij[k]].size();
						for(short m = 0; m < nAtoms; m++){
							atom2bonds[pij[k]][0]++;	/* make a flag */
							atom2bonds[pij[k]][connBond[pij[k]][m]]++;	/* add type of the bond */
						}
					}
				}
			}
		}

		/* ---------------- get Covalent Radius of atoms on shortest path that saving location in Constants::CovalentRadius[] --------------- */

		for(short i = 0; i < numAtoms; i++){
			char symbol = mol.getAtom(i).getAtomSymbol();
			if(symbol == Elements::CARBON){	/* symbol C */
				if(atom2bonds[i][3] == 1 || atom2bonds[i][2] == 2){
					cr[i] = 105;	/* value of sp saving location in Constants::CovalentRadius[] */
				}else if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){
					cr[i] = 104; /* value of sp2 saving location in Constants::CovalentRadius[] */
				}else{
					cr[i] = 103;	/* value of sp3 saving location in Constants::CovalentRadius[] */
				}
			}else if(symbol == Elements::NITROGEN){	/* symbol N */
				if(atom2bonds[i][3] == 1 || atom2bonds[i][2] == 2){
					cr[i] = 108;	/* value of sp saving location in Constants::CovalentRadius[] */
				}else if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){
					cr[i] = 107;	/* value of sp2 saving location in Constants::CovalentRadius[] */
				}else{
					cr[i] = 106;	/* value of sp3 saving location in Constants::polarizability_value[] */
				}
			}else if(symbol == Elements::OXYGEN){ /* symbol O */
				if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){
					cr[i] = 110; /* value of sp2 saving location in Constants::CovalentRadius[] */
				}else{
					cr[i] = 109; /* value of sp3 saving location in Constants::CovalentRadius[] */
				}
			}else if(symbol == Elements::PHOSPHORUS){ /* symbol P */
				if(atom2bonds[i][4] > 0 || atom2bonds[i][2] > 0){
					cr[i] = 112; /* value of sp2 saving location in Constants::CovalentRadius[] */
				}else{
					cr[i] = 111; /* value of sp3 saving location in Constants::CovalentRadius[] */
				}
			}else if(symbol == Elements::SULFUR){	/* symbol S */
				if(atom2bonds[i][4] > 0 || atom2bonds[i][2] > 0){
					cr[i] = 114; /* value of sp2 saving location in Constants::CovalentRadius[] */
				}else{
					cr[i] = 113; /* value of sp3 saving location in Constants::CovalentRadius[] */
				}
			}else{
				cr[i] = symbol; /* others atoms' covalent radius saving location in Constants::CovalentRadius[] */
			}		
		}

		/* -------------------------------- alpha -------------------------------------- */
		double a = 0.0;
		for(short i = 0; i < numAtoms; i++){
			a += Constants::CovalentRadius[cr[i]] / 0.77 - 1;
		}

		/* --------- calculate paths equal to 1, 2, and 3 that between any atom --------- */

		/* --------------------------------- k2pi ---------------------------------------- */
		short numBonds = mol.getBondsNum(false);
		vector<double> kpi(3, 0);
		if(pow((numBonds + a), 2) == 0)
			kpi[0] = Constants::ERROR_SIGNAL;
		else
			kpi[0] = pow((numAtoms + a - 1), 2) * (numAtoms + a) / pow((numBonds + a), 2);

		/* --------------------------------- k2pi ---------------------------------------- */
		vector<int> & p = mol.getP();
		if(pow((p[1] + a), 2) == 0)
			kpi[1] = Constants::ERROR_SIGNAL;
		else if(numAtoms < 3)
			kpi[1] = 0;
		else 
			kpi[1] = pow((numAtoms + a - 2), 2) * (numAtoms + a -1) / pow((p[1] + a), 2);


		/* --------------------------------- k3pi ---------------------------------------- */

		if(pow((p[2] + a), 2) == 0)
			kpi[2] = Constants::ERROR_SIGNAL;
		else{ /* for odd number of atoms and number of atoms > 3 */
			if(numAtoms > 3 && numAtoms % 2 != 0){
				kpi[2] = pow((numAtoms + a - 3), 2) * (numAtoms + a - 1) / pow((p[2] + a), 2);
			}/* for even number of atoms and number of atoms > 3 */
			else if(numAtoms > 3 && numAtoms % 2 == 0){
				kpi[2] = pow((numAtoms + a - 2), 2) * (numAtoms + a -3) / pow((p[2] + a), 2);
			}
		}

		/* --------------------------------- mfi ---------------------------------------- */
		double mfi = kpi[0] * kpi[1] / numAtoms;
		vector<double> fingerprint;

		fingerprint.insert(fingerprint.begin(), kpi.begin(), kpi.end());
		fingerprint.push_back(mfi);
		mol.setFingerprint(fingerprint, ::KierPathIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in KierPathIndex" << endl;
	}
}

void Calculator7::AllPathConnectivity(){
	try{
		if(mol.isFeaturesReady(::AllPathConnectivity)) return;
		if(mol.isAllPathFailed()) {
			vector<double> tmp(31, ::numeric_limits<float>::signaling_NaN());
			mol.setFingerprint(tmp, ::AllPathConnectivity);
			return;
		}

		vector<float> & valence = mol.getValence();
		vector<short> & pqn = mol.getPqn();
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector<double> & vc = mol.getVC();
		vector<double> & vvc = mol.getVVC();
		vector<double> & pqvc = mol.getPQVC();

		for(short i = 0; i < numAtoms; i++){
			short natoms = connAtoms[i].size();	
			if(natoms > 0){
				double tmpVal = 1.0 / sqrt(double(natoms));
				vc[0] += tmpVal; /* connectivity order-0 Index */					
				pqvc[0] += tmpVal / pqn[i];	/* first order atoms */
			}
			if(valence[i] > 0)
				vvc[0] += 1.0 / sqrt(valence[i]);	/* valence vertex connectivity order-0 index */
		}
		
		pqvc[0] *= 2;	/* principal quantum vertex connectivity order-0 index */

		/* -------------- calculate the connectivity order-1 Index ----------------- */
		short numBonds = mol.getBondsNum(false);	
		for(short i = 0; i < numBonds; i++){
			Bond &bond = mol.getBond(i);	
			short ai = bond.getFirstAtomNumber();
			short b = bond.getSecondAtomNumber();

			/* calculate the vertex connectivity order-1 Index */
			short nai = connAtoms[ai].size();
			short nb = connAtoms[b].size();
			double tmpVal = 1.0 /sqrt(double(nai) * nb);
			vc[1] += tmpVal;
			vvc[1] += 1.0 / sqrt(valence[ai] * valence[b]);	/* valence vertex connectivity order-0 index */
			/*  principal quantum vertex connectivity order-1 index */
			pqvc[1] += tmpVal / (pqn[ai] * pqn[b]);
		}
		pqvc[1] *= 4;	/* principal quantum vertex connectivity order-1 index */

		/* --- calculate the average vertex connectivity index with order-0 ~ order-5 and average valence connectivity index with order-0 ~ order-5 ----- */
		vector<double> avc(6, 0);
		vector<double> avvc(6, 0); 
		avc[0] = vc[0] / numAtoms;		/* the average vertex connectivity order-0 index */
		avvc[0] = vvc[0] / numAtoms;	/* the average valence connectivity Index chi-0 */

		if(numBonds != 0){
			avc[1] = vc[1] / numBonds;		/* the average vertex connectivity order-1 index */
			avvc[1] = vvc[1] / numBonds;	/* the average valence vertex connectivity Index chi-1 */
		}

		vector<int> & a = mol.getA();

		for(int i = 2; i < 6; i++){
			if(a[i] > 0){
				avc[i] = vc[i] / a[i];	/* the average vertex connectivity order-2 index */
				avvc[i] = vvc[i] / a[i];	/* the average valence vertex connectivity order-2 Index */
			}
		}
		/* ---------------------- principal quantum vertex connectivity index with order-1 ~ order-5 ------------------------- */

		pqvc[2] *= 8;	/* principal quantum vertex connectivity order-2 index */
		pqvc[3] *= 16;	/* principal quantum vertex connectivity order-3 index */
		pqvc[4] *= 32;	/* principal quantum vertex connectivity order-4 index */
		pqvc[5] *= 64;	/* principal quantum vertex connectivity order-5 index */

		/* ---------------------- aromaticity valence vertex connectivity order-1 index: AV1 ------------------------- */
		double av1 = numBonds > 0 ? 3 * vvc[1] / numBonds : 3 * vvc[1];	/* aromaticity valence vertex connectivity order-1 index: AV1 */

		vector<double> fingerprint;
		fingerprint.insert(fingerprint.end(), vc.begin(), vc.end());
		fingerprint.insert(fingerprint.end(), avc.begin(), avc.end());
		fingerprint.insert(fingerprint.end(), vvc.begin(), vvc.end());
		fingerprint.insert(fingerprint.end(), avvc.begin(), avvc.end());
		fingerprint.insert(fingerprint.end(), pqvc.begin(), pqvc.end());
		fingerprint.push_back(av1);
		mol.setFingerprint(fingerprint, ::AllPathConnectivity);
		vc.clear();
		pqvc.clear();
		vvc.clear();
	}catch(exception &ex){
		cout << ex.what() << " thrown in AllPathConnectivity" << endl;
	}
}


void Calculator7::gotoNextLevel(vector<short>& currentLevel, vector<short>& nextLevel){
	short nai = currentLevel.size();
	vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
	nextLevel.clear();
	for(short k = 0; k < nai; k++){
		nextLevel.insert(nextLevel.end(), connAtoms[currentLevel[k]].begin(), connAtoms[currentLevel[k]].end());	/* atoms on the level 4 */
	}
}

void Calculator7::MaximumPathIndex(){
	try{
		if(mol.isFeaturesReady(::MaximumPathIndex)) return;
		if(mol.isAllPathFailed()){
			vector<double> tmp(3, ::numeric_limits<float>::signaling_NaN());
			mol.setFingerprint(tmp, ::MaximumPathIndex);
			return;
		}

		int k = 0;
		vector<short> & detourMatrixData = mol.getDetourMatrix().getData();

		for(short i = detourMatrixData.size() - 1; i >= 0; i--){
			k += detourMatrixData[i];	/* calculate the Maximum Path Index */
		}

		/* ------------------------ calculate Wiener Type Maximum Path Index: wmpi ------------ */
		vector<short> & detourMatrix = mol.getDetourMatrix().getData();
		double h = 0.0;
		int t = 0;/* Wiener Type Maximum Path Index */	
		for(int m = detourMatrix.size() - 1; m >= 0; m--){
			if(detourMatrix[m] != 0){
				int dmp = detourMatrix[m] * (detourMatrix[m] + 1); /* get from Maximum-Path Matrix */
				t += dmp;	/* calculate Wiener Type Maximum Path Index */
				h += 1.0f / dmp;	/* reciprocal Maximum-Path Matrix */
			}
		}	

		/*------------------------- calculate reciprocal Wiener Type Maximum Path Index ------------ */
		/*k: Maximum Path Index; t: Wiener Type Maximum Path Index; h: reciprocal Wiener Type Maximum Path Index */
		vector<double> fingerprint(3, 0);
		fingerprint[0] = k;
		fingerprint[1] = t / 2;
		fingerprint[2] = h * 2;

		mol.setFingerprint(fingerprint, ::MaximumPathIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in MaximumPathIndex" << endl;
	}
}

void Calculator7::MinMaxPathIndex(){
	try{
		if(mol.isFeaturesReady(::MinMaxPathIndex)) return;
		if(mol.isAllPathFailed()){
			mol.setFingerprint(::numeric_limits<float>::signaling_NaN(), ::MinMaxPathIndex);
			return;
		}

		double k = 0.0;
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		vector<short> & detourMatrix = mol.getDetourMatrix().getData();
		for(int m = detourMatrix.size() - 1; m >= 0; m--){
			if(detourMatrix[m] > 0)
				k += double(distanceMatrix[m]) / detourMatrix[m];
		}
		mol.setFingerprint(k, ::MinMaxPathIndex);	
	}catch(exception &ex){
		cout << ex.what() << " thrown in MinMaxPathIndex" << endl;
	}
}

void Calculator7::AllPathWiener(){
	try{
		if(mol.isFeaturesReady(::AllPathWiener)) return;
		int k = mol.isAllPathFailed() ? ::numeric_limits<int>::signaling_NaN() : sum_element(mol.getAllPath().getData());
		mol.setFingerprint(k, ::AllPathWiener);
	}catch(exception &ex){
		cout << ex.what() << " thrown in AllPathWiener" << endl;
	}
}

void Calculator7::D_D_RingAndCircuitsIndex(){
	try{
		if(mol.isFeaturesReady(::D_D_RingIndex)) return;
		if(mol.isAllPathFailed()){
			vector<double> tmp(20, ::numeric_limits<float>::signaling_NaN());
			mol.setFingerprint(tmp, ::D_D_RingIndex);
			return;
		}

		vector<float> dd_sums(numAtoms, 0);
		vector< vector<float> > dd(numAtoms, vector<float> (numAtoms, 0));
		short sizeDDr = numAtoms > 13? numAtoms:13;
		vector<float> ddr(sizeDDr, 0);
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		vector<short> & detourMatrix = mol.getDetourMatrix().getData();
		/* ----------- Distance+Detour path Quotient Matrix and rows sums of the matrix ----------- */
		size_t m = 0;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j < numAtoms; j++){
				if(detourMatrix[m] != 0){
					dd[j][i] = dd[i][j] = float(distanceMatrix[m]) / detourMatrix[m];
				}
				m++;
			}
		}

		for(short i = 0; i < numAtoms; i++) dd_sums[i] = sum_element(dd[i]);

		/* ------- calculated the distance/detour ring index of order 3 ~ 12 ------- */
		// get rings copy not refer to original one
		vector< vector<short> > rings = mol.getRings();
		short nr = rings.size();
		for(short i = 0; i < nr; i++){
			short nri = rings[i].size() - 1;
			for(short j = nri - 1; j >= 0; j--){
				ddr[nri] += dd_sums[rings[i][j]];
			}
		}

		/* printout results of Distance+Detour path on ring Index of order 3 ~ 12 to output file */
		vector<double> fingerprint;
		for(short i = 3; i < 13; i++)
			fingerprint.push_back(ddr[i]);

		// calculateD_D_RingAndCircuitsIndex
		/***************************** search new circuits ********************************************/
		for(short i = 0; i < nr; i++){ /* first ring */
			for(short m = i + 1; m < nr; m++){/* the rings are behind of the first ring */
				/* ---------------------- building new circuit ------------------- */
				set<short> new_circuit(rings[m].begin(), rings[m].end());
				// set clear will clear the reference ringi. so I use following code
				new_circuit.insert(rings[i].begin(), rings[i].end());			
				/* check the new_circuit saved in my_rings or no */
				if(rings[i].size() + rings[m].size() - new_circuit.size() > 1){ /* found a new circuit */
					rings.push_back(vector<short> (new_circuit.begin(), new_circuit.end()));
				} 					
			} 
		} 

		/***************************** end of search new circuits ********************************************/
		/* ------- calculated the distance/detour ring index of order 3 ~ 12 ------- */
		ddr.assign(ddr.size(), 0);
		nr = rings.size();
		for(short i = 0; i < nr; i++){
			short nri = rings[i].size() - 1;
			for(short j = 0; j < nri; j++){
				ddr[nri] += dd_sums[rings[i][j]];
			}
		}

		/* printout results of Distance+Detour path on ring Index of order 3 ~ 12 to output file */
		for(short i = 3; i < 13; i++)
			fingerprint.push_back(ddr[i]);
		mol.setFingerprint(fingerprint, ::D_D_RingIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in D_D_RingIndex" << endl;
	}
}

void Calculator7::doAll(void ){
	PathWalk();
	KierPathIndex();
	AllPathConnectivity();
	MinMaxPathIndex();
	MaximumPathIndex();
	AllPathWiener();
	D_D_RingAndCircuitsIndex();
}