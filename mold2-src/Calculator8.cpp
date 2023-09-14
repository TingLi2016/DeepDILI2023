#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Constants.h"
#include "Calculator8.h"

Calculator8::Calculator8(Molecule& _mol)
:mol(_mol){
	numAtoms = mol.getAtomsNum(false);
	dmew.resize(3, SymmetricMatrix<float> (numAtoms));
	zmvpw.resize(4, SymmetricMatrix<float> (numAtoms));
	pv.resize(numAtoms, vector<short>(2, 0));
	laplacianMatrix.resize(numAtoms, vector<char>(numAtoms, 0));
	vertexDistanceCounts.resize(numAtoms, 0);
	initializeVertexDistanceCounts();

	initializeBW_MW_VW();
	initializePV();
	initializePW();
	initializeLaplacianMatrix();
	doAll();
}

Calculator8::~Calculator8(void){
}

void Calculator8::initializeLaplacianMatrix(){
	vector< vector<short> > & connectedAtomIDs =mol.getConnectedAtomIDs(false);		
	for(short i = 0; i < numAtoms; i++){
		for(short j = 0; j < numAtoms; j++){
			if(i == j)
				laplacianMatrix[i][j] = connectedAtomIDs[i].size();
			else
				for(short n = connectedAtomIDs[i].size() - 1; n >= 0; n--)
					laplacianMatrix[i][connectedAtomIDs[i][n]] = -1;
		}
	}
}

void Calculator8::initializeBW_MW_VW(){
	vector<short> & atomic_num = mol.getAtomicNum(false);
	vector<float> atomic_weight(numAtoms, 0);
	vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
	vector<float> vdw(numAtoms, 0);	/* Van Der Waals */
	/* Distance Matrix of Electronegativity Weighted */
	vector< vector<float> > en(numAtoms, vector<float>(3, 0));

	//TODO: the initialization of en is very weird in the original code
	for(short i = 0; i < numAtoms; i++) { // get atomic number of the atom		
		int symbol = mol.getAtom(i).getAtomSymbol();
		atomic_weight[i] = Constants::atom_weight[symbol];
		vdw[i] = Constants::vdw_radius[symbol];
		en[i][0] = Constants::elect_negativities[symbol][0];
		en[i][1] = Constants::elect_negativities[symbol][1];
		en[i][2] = Constants::elect_negativities[symbol][2];
	}

	vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
	for(short i = 0; i < numAtoms; i++){
		for(short j = i; j < numAtoms; j++){
			if(i == j){ /* calculte the diagonal entries first */
				zmvpw[ZW].setValueToCell(i, i, 1.0f - (6.0f / atomic_num[i])); /* the diagonal entries: C-atomic number / atomic number of atom-i */
				zmvpw[MW].setValueToCell(i, i, 1.0f - (12.0112f / atomic_weight[i])); /* the diagonal entries: C-weight / atomic weight of atom-i */
				zmvpw[VW].setValueToCell(i, i, 1.0f - (1.70f / vdw[i])); /* the diagonal entries: C-weight / van der Waals weight of the atom-i */
				dmew[0].setValueToCell(i, i, 1.0f - (Constants::elect_negativities[5][0] / en[i][0])); /* the diagonal entries: C-weight / electronegativity weight of atom-i */
				dmew[1].setValueToCell(i, i, 1.0f - (Constants::elect_negativities[5][1] / en[i][1])); /* the diagonal entries: C-weight / electronegativity weight of atom-i */
				dmew[2].setValueToCell(i, i, 1.0f - (Constants::elect_negativities[5][2] / en[i][2])); /* the diagonal entries: C-weight / electronegativity weight of atom-i */
			}else{	/* calculte others data of the Bonds Weighted distance matrix */
				vector<short> & DPM = mol.getDistancePathMatrix().getPath(i, j);
				for(short k = DPM.size() - 1; k > 0; k--){ /* running over all atoms that on the shortest path */
					/* two vertices incident to the considered bond */
					float nzw = (float)atomic_num[DPM[k]] * atomic_num[DPM[k - 1]];
					float nmw = atomic_weight[DPM[k]] * atomic_weight[DPM[k - 1]];
					float nvw = vdw[DPM[k]] * vdw[DPM[k - 1]];
					float ndmew0 = en[DPM[k]][0] * en[DPM[k - 1]][0];
					float ndmew1 = en[DPM[k]][1] * en[DPM[k - 1]][1];
					float ndmew2 = en[DPM[k]][2] * en[DPM[k - 1]][2];
					
					for(short m = connAtoms[DPM[k]].size() - 1; m >= 0; m--){
						/* look for the next atom of the path */
						if(connAtoms[DPM[k]][m] == DPM[k - 1]){
							/* if the bond is an Aromatic bond */
							float p = connBond[DPM[k]][m];
							if(connBond[DPM[k]][m] == Bond::AROMATIC_BOND)
								p = 1.5;	/* value of the aromatic bond equal to 1.5 */
							zmvpw[ZW].cellValuePlus(i, j, (36.0f / nzw) * (1.0f / p));	/* calculate the Bonds Weighted distance matrix */
							zmvpw[MW].cellValuePlus(i, j, (144.2689f / nmw) * (1.0f / p));	/* calculate the Mass Weighted Matrix */
							zmvpw[VW].cellValuePlus(i, j, (2.89f / nvw) * (1.0f / p));
							dmew[0].cellValuePlus(i, j, (Constants::elect_negativities[5][0] * Constants::elect_negativities[5][0] / ndmew0) * (1.0f / p));
							dmew[1].cellValuePlus(i, j, (Constants::elect_negativities[5][1] * Constants::elect_negativities[5][1] / ndmew1) * (1.0f / p));
							dmew[2].cellValuePlus(i, j, (Constants::elect_negativities[5][2] * Constants::elect_negativities[5][2] / ndmew2) * (1.0f / p));
						}
					}					
				} /* end of FOR-k */
			}
		}
	} /* end of FOR-i */
}

void Calculator8::initializePW(){
	/* -------------------------- calculate the Polarizability Weighted Distance Matrix ------------------------------- */
	vector< vector<short> > & connAtoms =mol.getConnectedAtomIDs(false);
	vector< vector<char> > & connBond =mol.getConnectedBondTypes(false);
	for(short i = 0; i < numAtoms; i++){
		for(short j = i; j < numAtoms; j++){
			if(pv[i][1] == 0) /* we have not the polarizability-value of the atom-i */
				continue;	  /* ignore the atom */
			else{
				if(i == j){ /* calculte the diagonal entries first */
					/* the diagonal entries: C-weight / atomic weight of atom-i */
					zmvpw[PW].setValueToCell(i, i, 1 - (Constants::polarizability_value[pv[i][0]] / Constants::polarizability_value[pv[i][1]])); /* the polarizability-value of atom-C using same sp-N of the atoms that on the shortest path */
				}else{	/* calculte others data of the matrix */
					vector<short> & DPM = mol.getDistancePathMatrix().getPath(i, j);
					for(short k = DPM.size() - 1; k > 0; k--){ /* running over all atoms that on the shortest path */
						short atomk1 = DPM[k - 1];
						short atomk = DPM[k];
						if(Constants::polarizability_value[pv[atomk1][1]] == 0 ||	/* we have not the polarizability-value of the atom[i][j][k] */
							Constants::polarizability_value[pv[atomk][1]] == 0)		/* we have not the polarizability-value of the atom[i][j][k-1] */
							continue;   /* ignore the atom */
						else{
							/* two vertices incident to the considered bond */
							float n = Constants::polarizability_value[pv[atomk][1]] * Constants::polarizability_value[pv[atomk1][1]];
							short na = connAtoms[atomk].size();
							for(short m = 0; m < na; m++){
								/* look for the next atom of path */
								if(connAtoms[atomk][m] == atomk1){
									/* if the bond is an Aromatic bond */
									float p = connBond[atomk][m];
									if(connBond[atomk][m] == Bond::AROMATIC_BOND)
										p = 1.5;	/* value of the aromatic bond equal to 1.5 */

									/* the Constants::polarizability_value of atom-C using same sp-N of atoms that on the shortest path */
									float q = Constants::polarizability_value[pv[atomk][0]] * Constants::polarizability_value[pv[atomk1][0]];
									zmvpw[PW].cellValuePlus(i, j, (q / n) * (1.0f / p));	/* calculate the Polarizability weighted distance matrix */
								}
							}
						}						
					} /* end of FOR-k */
				} // end of esle -- i != j
			}
		}
	} /* end of FOR-i */
}

void Calculator8::initializePV(){
	/* atomic polarizability: i is order of atom
	pv[i][0]: sp-N -- 1=sp3; 2=sp2; 3=sp;
	pv[i][1]: polarizability-value of atom-i saving location in Constants::polarizability_value[] */

	vector< vector<char> > & connBond =mol.getConnectedBondTypes(false);
	/* atom2bonds[i][j]--number of bonds type j connected to the atom i */
	vector< vector<char> > atom2bonds(numAtoms, vector<char>(mol.m_maxBondType, 0));

	/* get type of bond that on the shortest path */
	DistancePathMatrix &DPM =mol.getDistancePathMatrix();
	for(short i = 0; i < numAtoms; i++){
		for(short j = i; j < numAtoms; j++){
			vector<short> &pij = DPM.getPath(i, j);
			for(short k = pij.size() - 1; k > 0; k--){
				if(atom2bonds[pij[k]][0] != 0)	/* bonds of the atom count already -- does not count again */
					continue;
				else{	/* bonds of the atoms does not count before */
					short connAtoms = connBond[pij[k]].size();
					for(short m = 0; m < connAtoms; m++){
						atom2bonds[pij[k]][0]++;	/* make a flag */
						atom2bonds[pij[k]][connBond[pij[k]][m]]++;	/* add type of the bond */
					}
				}
			}
		}
	}

	/* ---------------- get polarizability-value of shortest path atoms saving location in Constants::polarizability_value[] --------------- */
	for(short i = 0; i < numAtoms; i++){
		char symbol = mol.getAtom(i).getAtomSymbol();
		if(symbol == Elements::CARBON){	/* symbol C */
			if(atom2bonds[i][3] == 1 || atom2bonds[i][2] == 2){
				pv[i][0] = pv[i][1] = 3;	/* value of sp saving location in Constants::polarizability_value[] */
			}
			if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){
				pv[i][0] = pv[i][1] = 2; /* value of sp2 saving location in Constants::polarizability_value[] */
			}else{
				pv[i][0] = pv[i][1] = 1;	/* value of sp3 saving location in Constants::polarizability_value[] */
			}
		}else if(symbol == Elements::NITROGEN){	/* symbol N */
			if(atom2bonds[i][3] == 1 || atom2bonds[i][2] == 2){
				pv[i][1] = 6;	/* value of sp saving location in Constants::polarizability_value[] */
				pv[i][0] = 3;	/* sp */
			}
			if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){
				pv[i][1] = 5;	/* value of sp2 saving location in Constants::polarizability_value[] */
				pv[i][0] = 2;	/* sp2 */
			}else{
				pv[i][1] = 4;	/* value of sp3 saving location in Constants::polarizability_value[] */
				pv[i][0] = 1;	/* sp3 */
			}
		}else if(symbol == Elements::OXYGEN){	/* symbol O */
			if (atom2bonds[i][2] > 0){
				pv[i][1] = 8; /* value of sp2 saving location in Constants::polarizability_value[] */
				pv[i][0] = 3;
			}else{
				pv[i][1] = 7; /* value of sp3 saving location in Constants::polarizability_value[] */
				pv[i][0] = 2;
			}
		}else if(symbol == Elements::FLOURINE){	/* symbol F */
			pv[i][1] = 12; /* value of sp saving location in Constants::polarizability_value[] */
			pv[i][0] = 1;
		}else if(symbol == Elements::PHOSPHORUS){	/* symbol P */
			pv[i][1] = 11; /* value of sp3 saving location in Constants::polarizability_value[] */
			pv[i][0] = 1;
		}else if(symbol == Elements::SULFUR){	/* symbol S */
			if (atom2bonds[i][2] > 0){
				pv[i][1] = 10; /* value of sp2 saving location in Constants::polarizability_value[] */
				pv[i][0] = 2;
			}else{
				pv[i][1] = 9; /* value of sp3 saving location in Constants::polarizability_value[] */
				pv[i][0] = 1;
			}
		}else if(symbol == Elements::CLORINE){	/* symbol Cl */
			pv[i][1] = 13; /* value of sp saving location in Constants::polarizability_value[] */
			pv[i][0] = 1;
		}else if(symbol == Elements::BROMINE){	/* symbol Br */
			pv[i][1] = 14; /* value of sp saving location in Constants::polarizability_value[] */
			pv[i][0] = 1;
		}else if(symbol == Elements::IODINE){	/* symbol I */
			pv[i][1] = 15; /* value of sp saving location in Constants::polarizability_value[] */
			pv[i][0] = 1;
		}		
	}
}
void Calculator8::calculateWeight(SymmetricMatrix<float> & xw, double &le_ev, double &sum_ev, double &abs_sum_ev){
	Array2D<double> a(numAtoms, numAtoms);		

	for(short i = 0; i < numAtoms; i++){
		a[i][i] = xw.getValueInCell(i, i);
		for(short j = i + 1; j < numAtoms; j++){
			a[j][i] = a[i][j] = xw.getValueInCell(i, j); /* Heteroatoms and Multiple bonds weighted Distance Matrix */				
		}
	}		
	vector<double> ev = getEigenvalue(a);
	for(short i = 0; i < numAtoms; i++){
		sum_ev += ev[i]; /* SEigHMB, SEigMD, SEigVD, SEigPD */
		abs_sum_ev += fabs(ev[i]); /* AEigHMB, AEigMD, AEigVD, AEigPD */
		if(ev[i] > 1.0e-9 && ev[i] > le_ev)
			le_ev = ev[i];
	}
}


void Calculator8::SumWeighted(){
	try{
		if(mol.isFeaturesReady(::SumWeighted)) return;

		vector<double> le_ev(7, 0), sum_ev(7, 0), abs_sum_ev(7, 0);
		for(short i = 0; i < 3; i++){
			calculateWeight(zmvpw[i], le_ev[i], sum_ev[i], abs_sum_ev[i]);/* Heteroatoms and Multiple bonds weighted Distance Matrix */
			calculateWeight(dmew[i], le_ev[i + 4], sum_ev[i + 4], abs_sum_ev[i + 4]);/* Distance Matrix of Electronegativity Weighted */
		}
		calculateWeight(zmvpw[3], le_ev[3], sum_ev[3], abs_sum_ev[3]);/* Heteroatoms and Multiple bonds weighted Distance Matrix */

		vector<double> fingerprint;
		/* printout the results of MEigHMB, MEigMD, MEigVD, MEigPD,MEigePsD, MEigeSD, MEigeARD to output file */
		fingerprint.insert(fingerprint.end(), le_ev.begin(), le_ev.end());

		/* printout the results of SEigHMB, SEigMD, SEigVD, SEigPD, SEigePsD, SEigeSD, SEigeAeRD to output file */
		fingerprint.insert(fingerprint.end(), sum_ev.begin(), sum_ev.end());	

		/* printout the results of AEigHMB, AEigMD, AEigVD, AEigPD, AEigePsD, AEigeSD, AEigeARD to output file */
		fingerprint.insert(fingerprint.end(), abs_sum_ev.begin(), abs_sum_ev.end());	
		mol.setFingerprint(fingerprint, ::SumWeighted);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SumWeighted"<< endl;
	}
}

void Calculator8::BalabanAVDCI_BalabanTypeWeighted(){
	try{
		if(mol.isFeaturesReady(::BalabanAVDCI_BalabanTypeWeighted)) return;

		/* sigma of atom-i: the vertex distance degree of the ith atom */
		vector<int> & sigma = mol.getDistanceMatrix(false).getSigma();
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);

		/*------------------------- calculate Average vertex distance connectivity index: AVDCI ------------------- */
		double k = 0.0;
		for(short i = 0; i < numAtoms; i++){	/* running over all atoms */
			for(short j = connAtoms[i].size() - 1; j >= 0; j--){	/* first order atoms */
				if(connAtoms[i][j] != i){	/* does not recount atoms */
					/* calculate the (sigma of atom-i) * (sigma of the first order atoms) */
					k += 1.0 / sqrt(double(sigma[i]) * sigma[connAtoms[i][j]]);
				}
			}
		}
		short numBonds = mol.getBondsNum(false);
		short nr = numBonds - numAtoms + 1;	/* number of rings */
		k *= numBonds / (nr + 1.0);	/* Average vertex distance connectivity index: AVDCI */

		vector<double> fingerprint(8, 0);
		fingerprint[0] = k / 2;	/* Average vertex distance connectivity index: AVDCI */

		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);

		/* ---------------------- f: the conventional bond order -------------------- */
		vector<float> bond(numAtoms, 0);
		for(short i = 0; i < numAtoms; i++){
			for(short j = connBond[i].size() - 1; j >= 0; j--){
				if(connBond[i][j] == Bond::AROMATIC_BOND)	/* aromatic bond */
					bond[i] += 1.5f;
				else
					bond[i] += connBond[i][j];
			}
			bond[i] += 1.0f;
		}

		for(short i = 0; i < numAtoms; i++){	/* running over all atoms */
			for(short j = connAtoms[i].size() - 1; j >= 0; j--){	/* the first order atoms */
				short k = connAtoms[i][j];
				double sisk = sigma[i] * sigma[k];
				for(short wid = 0; wid < 3; wid++){	
					/* [1]: Balaban heteroatoms bonds weighted index: BHBWI */
					/* [2]: Balaban mass weighted index: BMWI */
					/* [3] Balaban van der Waals weighted index: BVDWWI */
					fingerprint[wid + 1] += 1.0 / sqrt(sisk / fabs((zmvpw[wid].getValueInCell(i, i) + 1.0) * bond[i] * (zmvpw[wid].getValueInCell(k, k) + 1.0) * bond[k]));

					/* [4]: Atomic electronegativities weighted with Pauling-Scale index: EWPSI */
					/* [5]: Atomic electronegativities weighted with Sanderson-Scale index: EWSSI */
					/* [6]: Atomic electronegativities weighted with Allred-Rochow Scale index: EWARSI */
					fingerprint[wid + 4] += 1.0 / sqrt(sisk / fabs((1.0 + dmew[wid].getValueInCell(i, i)) * bond[i] * (1 + dmew[wid].getValueInCell(k, k)) * bond[k]));
				}
				/* Balaban polarizability weighted index: BPWI */
				fingerprint[7] += 1.0 / sqrt(sisk / ((fabs(zmvpw[3].getValueInCell(i, i)) + 1.0) * bond[i] * (fabs(zmvpw[3].getValueInCell(k, k)) + 1.0) * bond[k]));
			}
		}
		mol.setFingerprint(fingerprint, ::BalabanAVDCI_BalabanTypeWeighted);
	}catch(exception &ex){
		cout << ex.what() << " thrown in BalabanAVDCI_BalabanTypeWeighted"<< endl;
	}
}

// Polarizability weighted distance matrix
void Calculator8::PolarizabilityWeighted(){
	try{
		if(mol.isFeaturesReady(::PolarizabilityWeighted)) return;

		float pw = zmvpw[PW].cumulativeSum() / 2;
		mol.setFingerprint(pw, ::PolarizabilityWeighted);
	}catch(exception &ex){
		cout << ex.what() << " thrown in PolarizabilityWeighted"<< endl;
	}
}

//Bonds Weighted distance matrix
void Calculator8::BW_MW_VW_EWeighted(){
	try{
		if(mol.isFeaturesReady(::BW_MW_VW_EWeighted)) return;

		vector<double> fingerprint;
		fingerprint.push_back(zmvpw[ZW].cumulativeSum() / 2);
		fingerprint.push_back(zmvpw[MW].cumulativeSum() / 2);
		fingerprint.push_back(zmvpw[VW].cumulativeSum() / 2);
		fingerprint.push_back(dmew[0].cumulativeSum() / 2);
		fingerprint.push_back(dmew[1].cumulativeSum() / 2);
		fingerprint.push_back(dmew[2].cumulativeSum() / 2);
		mol.setFingerprint(fingerprint, ::BW_MW_VW_EWeighted);	
	}catch(exception &ex){
		cout << ex.what() << " thrown in BW_MW_VW_EWeighted"<< endl;
	}
}

void Calculator8::BalabanShortPathIndex(){
	try{
		if(mol.isFeaturesReady(::BalabanShortPathIndex)) return; 

		vector<int> sigma = mol.getDistanceMatrix(false).getSigma();	
		double bspi = 0.0;
		DistancePathMatrix & DPM = mol.getDistancePathMatrix();
		for(short i = 0; i < numAtoms; i++){
			double eweight = 0;
			for(short j = i; j < numAtoms; j++){
				double r = 1.0;
				vector<short> & pij = DPM.getPath(i, j);
				for(short k = pij.size() - 1; k > 0; k--){
					// product of edge weight
					r *= 1.0 / sqrt(double(sigma[pij[k]]) * sigma[pij[k-1]]);
				}
				eweight += r; /* sum the product of edge weight on each shortest path */
			}
			bspi += eweight;	/* Balaban Short-Path index: BSPI */
		}
		mol.setFingerprint(bspi, ::BalabanShortPathIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in BalabanShortPathIndex"<< endl;
	}
}

void Calculator8::SecondMoharIndex(){
	try{
		if(mol.isFeaturesReady(::SecondMoharIndex)) return;

		Array2D<double> a(numAtoms, numAtoms, 0.0);
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				a[i][j] = laplacianMatrix[i][j];/* get Lapacian matrix and the data start from 1 to N */
			}
		}

		vector<double> ev = getEigenvalue(a);	
		double min_ev = numeric_limits<double>::max();
		for(short i = ev.size() - 1; i >= 0; i--){
			if(ev[i] > 0.0001 && ev[i] < min_ev)
				min_ev = ev[i];
		}
		double ti2 = 0.0;
		if(numAtoms * min_ev != 0)
			ti2 = 4.0 / (numAtoms * min_ev); /* Second Mohar Index: TI2 */
		mol.setFingerprint(ti2, ::SecondMoharIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SecondMoharIndex"<< endl;
	}
}

void Calculator8::SpanningTree(){
	try{
		if(mol.isFeaturesReady(::SpanningTree)) return;

		Array2D<double> a(numAtoms, numAtoms);
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j< numAtoms; j++){
				/* get Lapacian matrix and the data start from 1 to N */
				a[i][j] = laplacianMatrix[i][j];
			}
		}

		vector<double> ev = getEigenvalue(a); /* function call: calculate the eigenvalues */
		double stlv = 1;
		for(short i = 0; i < numAtoms; i++){
			if(ev[i] > 1e-9){ /* take the positive eigenvalue */
				stlv *= ev[i];
			}
		}
		stlv /= numAtoms;
		stlv = fabs(log(stlv));
		mol.setFingerprint(stlv, ::SpanningTree);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SpanningTree"<< endl;
	}
}


void Calculator8::LaplacianMatrix(){
	try{
		if(mol.isFeaturesReady(::LaplacianMatrix)) return;

		short numBonds =mol.getBondsNum(false);
		vector< vector<short> > rings = mol.getRings();
		short nr = rings.size();		

		DistanceMatrix<float> resistance_matrix;
		resistance_matrix.resize(numAtoms);
		vector< vector<short> > ring_path_both(nr);

		DistancePathMatrix &DPM = mol.getDistancePathMatrix();
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j < numAtoms; j++){
				int ring_bond = 0;
				bool doit = false;
				vector<char> used(numAtoms, 0);
				vector<short> ring_path;
				vector<short> & pij = DPM.getPath(i, j);
				short path_length = pij.size();

				for(short k = path_length - 1; k > 0; k--){
					bool find = false;					
					for(short m = 0; m < numBonds; m++){
						Bond & bond = mol.getBond(m);
						short first = bond.getFirstAtomNumber();
						short second = bond.getSecondAtomNumber();
						if(((first == pij[k] && second == pij[k - 1]) ||
							(first == pij[k - 1] && second == pij[k])) &&
							bond.getBondTopology() == 1){
								find = true; /* found ring path */
								break;
						}
					}
					if (!find){					
						resistance_matrix.cellValuePlus(i, j, 1.0f);
						if(ring_bond == 0)
							continue;
						else doit = true;
					}else{
						if(k == 1) doit = true;
						ring_bond++;
						used[pij[k]] = used[pij[k - 1]] = 1;
						if(ring_path.empty()){
							ring_path.push_back(pij[k]);
						}
						ring_path.push_back(pij[k - 1]);
					}

					if(doit){
						double dum = 1.0 / ring_bond;
						for(short n = 0; n < nr; n++){
							ring_path_both[n].clear();
							for(short p = rings[n].size() - 1; p >= 0; p--){
								for(short q = ring_path.size() - 1; q >= 0; q--){
									if(rings[n][p] == ring_path[q]){
										ring_path_both[n].push_back(ring_path[q]);
										break;
									}
								}
							}
						}
						find = false;
						short len = ring_path.size();
						for(short n = 0; n < nr; n++){
							if(len == ring_path_both[n].size()){
								dum += 1.0/(rings[n].size() - ring_bond);
								find = true;
							}
						}
						if(!find){
							int v = 0;
							for(short n = 0; n < nr; n++){
								if(ring_path_both[n].size() > 1){
									for(short p = rings[n].size() - 1; p >=0; p--){
										if(used[rings[n][p]] == 0){
											v++;
											used[rings[n][p]] = 1;
										}
									}
								}
							}
							dum += 1.0 /(v+1);
						}
						resistance_matrix.cellValuePlus(i, j,  (float)(1.0 / dum));					
						ring_bond = 0;
						doit = false;
						used.assign(used.size(), 0);
						ring_path.clear();
					}
				}/* end FOR-k */
			}/* end FOR-j */
		}/* end FOR-i */

		double ilm =  sum_element(resistance_matrix.getData());
		
		// First No-Zeor eigenvalue of Laplacian Matrix: elm1
		double t = ((double)numBonds) / numAtoms;
		double elm1 = numBonds == 0 ? 0 : 2.0 * ilm * log10(t);
		vector<double> fingerprint(2, 0);
		fingerprint[0] = ilm;
		fingerprint[1] = elm1;
		mol.setFingerprint(fingerprint, ::LaplacianMatrix);
	}catch(exception &ex){
		cout << ex.what() << " thrown in LaplacianMatrix"<< endl;
	}
}

void Calculator8::initializeVertexDistanceCounts(){	
	/* initial vertex_distance_counts[i] to 0 */
	vector<short> & distanceMatrix =mol.getDistanceMatrix(false).getData();
	/* calculate sum of vertex distance counts */
	for(int i = distanceMatrix.size() - 1; i >= 0; i--){ 
		if(distanceMatrix[i] != 0)
			vertexDistanceCounts[distanceMatrix[i]] += 2;
	}
	/* count sigma_i and saving in vertex_distance_counts[0] */
	vertexDistanceCounts[0] = numAtoms * numAtoms;
}

void Calculator8::SchultzMTI_VVD(){
	try{
		if(mol.isFeaturesReady(::SchultzMTI_VVD)) return;

		//Schultz type Molecular Topolgical Index of valence vertex degrees: SMTIV
		double smtiv = 0.0;
		/* printout elements of distance_matrix[][].path_count */
		short numAtoms =mol.getAtomsNum(false);
		vector<float> & cbo =mol.getConventionalBondOrder();
		vector<int> & sigma =mol.getDistanceMatrix(false).getSigma();
		for(short i = 0; i < numAtoms; i++){
			/* ---------------------- calculate the conventional bond order ---------------- */
			double smtiv_valence = cbo[i];

			/* ------------------------ calculate valence electrons ------------------------ */
			Atom & atom =mol.getAtom(i);
			double delta_v = valenceElections(atom.getAtomSymbol(), atom.getCharge(), i, smtiv_valence);

			/* ------------ calculate the Schultz type Molecular Topolgical Index of valence vertex degrees: SMTIV --------- */
			smtiv += (delta_v + sigma[i]) * delta_v;
		}
		mol.setFingerprint(smtiv, ::SchultzMTI_VVD);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SchultzMTI_VVD"<< endl;
	}
}

double Calculator8::valenceElections(char symbol, char charge, short i, double Zagreb_Mx_valence){
	int lone_pair_eletrons = 0;
	vector< vector<char> > & connBond =mol.getConnectedBondTypes(false);
	double ev = Zagreb_Mx_valence + 1;
	if(symbol == Elements::NITROGEN){/* symbol N */
		if(charge == 0) lone_pair_eletrons = 2;
		ev = (Zagreb_Mx_valence + lone_pair_eletrons) + 1;
	}else if(symbol == Elements::OXYGEN) {/* symbol O */
		lone_pair_eletrons = 4;
		ev = (Zagreb_Mx_valence + lone_pair_eletrons) + 1;
	}else if(symbol == Elements::FLOURINE){/* symbol F */
		if(charge == 0) lone_pair_eletrons = 6;
		ev = (Zagreb_Mx_valence + lone_pair_eletrons) + 1;
	}else if(symbol == Elements::SILICON){/* symbol Si */
		ev = 13.0;
	}else if(symbol == Elements::PHOSPHORUS){/* symbol P */		
		if(Zagreb_Mx_valence <= 3) lone_pair_eletrons = 2;
		ev = 14.0;
	}else if(symbol == Elements::SULFUR){/* symbol S */
		if(sum_element(connBond[i]) == 4) lone_pair_eletrons = 2;
		ev = 15.0;		
	}else if(symbol == Elements::CLORINE){/* symbol Cl */
		ev = 16.0;
		lone_pair_eletrons = 6;
	}else if(symbol == Elements::BROMINE){/* symbol Br */		
		if(charge == 0) lone_pair_eletrons = 6;
		ev = 34.0;
	}else if(symbol == Elements::IODINE){/* symbol I */		
		if(charge == 0) lone_pair_eletrons = 6;
		ev = 52.0;
	}

	double delta_v = (Zagreb_Mx_valence + lone_pair_eletrons) / (ev - (Zagreb_Mx_valence + lone_pair_eletrons));
	
	return delta_v;
}

void Calculator8::VertexDegree(){
	try{
		if(mol.isFeaturesReady(::VertexDegree)) return;

		/* calculate eta_i, sigma, and sum of vertex distance counts of first order */
		double lvdpc = 0;
		double s = 1.0;
		double h = 0.0;
		double nti = 1.0;

		vector< vector<short> > & connAtoms =mol.getConnectedAtomIDs(false);
		vector<int> & sigma =mol.getDistanceMatrix(false).getSigma();
		int sigmaI =mol.getDistanceMatrix(false).sumAllCells();
		for(short i = 0; i < numAtoms; i++){
			short nX = connAtoms[i].size();
			if(nX > 0){
				s *= nX; /* Narumi-type topological index */
				h += 1.0 / nX; /* harmonic topological index */
			}else{ // one atom case
				h = 1.0;
			}

			/* Log of vertex distance path count: LVDPC */
			if(sigma[i] == 0)          /* the item just has one atom -- atom's distance=0 */
				lvdpc = 0;  /* printout a spical signal */
			else
				lvdpc += (double) log((double)sigma[i]); /* calculate the Log of vertex distance path count: LVDPC */
		}

		nti = log(s * 1.0);	/* Narumi topological index: NTI */

		double hti = numAtoms / h; /* Harmonic topological index: HTI */

		double gti = pow(s, 1.0 / numAtoms); /* calculate Geometric topological index: GTI */

		/* average vertex distance path count: AVDPC */
		double avdpc = ((double)sigmaI) / numAtoms;

		/* calculate Polarity number */
		int tdc3 = numAtoms < 4 ? 0:vertexDistanceCounts[3] / 2;    /* get TDC3 */

		/* calculate the Balaban-type of mean square vertex distance index: BMSVDI */
		int m = 0;
		vector<short> & distanceMatrix =mol.getDistanceMatrix(false).getData();
		for(int i = distanceMatrix.size() - 1; i >= 0; i--){
			m += distanceMatrix[i] * distanceMatrix[i];		
		}
		m *= 2;

		double bmsvdi = numAtoms == 1 ? sqrt(double(m)) : sqrt(double(m)) / (numAtoms * (numAtoms - 1)); /* calculate BMSVDI */

		vector<double> fingerprint(7, 0.0);
		fingerprint[0] = nti;
		fingerprint[1] = hti;
		fingerprint[2] = gti;
		fingerprint[3] = tdc3;
		fingerprint[4] = lvdpc;
		fingerprint[5] = avdpc;
		fingerprint[6] = bmsvdi;
		mol.setFingerprint(fingerprint, ::VertexDegree);
	}catch(exception &ex){
		cout << ex.what() << " thrown in VertexDegree"<< endl;
	}
}

void Calculator8::MeanTotalDistance(){
	try{
		if(mol.isFeaturesReady(::MeanTotalDistance)) return;

		/* --- calculate the number of distances with value g in the triangular Destance submatrix --- */
		int g = 0;
		for(short i = 1; i < numAtoms; i++)
			if(vertexDistanceCounts[i] != 0)	/* number of different distance between atoms */
				g++;	/* count number of vertex distance degree over all atoms */

		/* calculate the Wiener Index */
		double w =mol.getDistanceMatrix(false).sumAllCells() / 2; /* calculate w as half of total topological distance */

		/* ---------------- claculate VDCEI, VDCMI, TVDCEI, TVDCMI -------------- */
		double vdcei = 0.0;
		double vdcmi = 0.0;
		double tvdcei = 0.0;
		double tvdcmi = 0.0;
		for(short i = 1; i <= g; i++){
			/* ---------------- vertex distance count equality index: VDCEI -------------- */
			double r = double(vertexDistanceCounts[i]) / (numAtoms * (numAtoms - 1));
			vdcei += -r * log10(r) / Constants::LOG20;

			/* ---------------- vertex distance count magnitude index: VDCMI -------------- */
			double s = i / w;
			double tmpV = -vertexDistanceCounts[i] / 2.0 / Constants::LOG20;
			vdcmi += tmpV * s * log10(s);

			/* ---------------- total vertex distance count equality index: TVDCEI -------------- */
			tvdcei += tmpV * log10(vertexDistanceCounts[i]/2.0);

			/* ---------------- total vertex distance count magnitude index: TVDCMI -------------- */
			tvdcmi += tmpV * i * log10(double(i));
		}
		double r = numAtoms * (numAtoms - 1.0) / 2.0;
		if(r > 0){
			tvdcei += r * log10(r) / Constants::LOG20;	/* TVDCEI */
			tvdcmi += w * log10(w) / Constants::LOG20;	/* TVDCMI */
		}
		vector<double> fingerprint(4, 0);
		fingerprint[0] = vdcei;
		fingerprint[1] = vdcmi;
		fingerprint[2] = tvdcei;
		fingerprint[3] = tvdcmi;
		mol.setFingerprint(fingerprint, ::MeanTotalDistance);
	}catch(exception &ex){
		cout << ex.what() << " thrown in MeanTotalDistance"<< endl;
	}
}

void Calculator8::Zagreb_M2_v(){
	try{
		if(mol.isFeaturesReady(::Zagreb_M2_v)) return;

		double Zagreb_M2_v = 0.0; // first Zagreb index by valence vertex degrees: Zagreb_M2_v
		vector< vector<char> > & connBond =mol.getConnectedBondTypes(false);
		vector< vector<short> >& connAtoms =mol.getConnectedAtomIDs(false);
		// valence electrons of Zagreb_M2_v without H
		for(short i = 0; i < numAtoms; i++){
			double zagreb_M2_valence = 0;
			int aromatic_bond_count1 = 0; // count aromatic bonds
			for(short j = connBond[i].size() - 1; j >= 0; j--){ /* second order atoms */
				/* ---------------------- calculate the conventional bond order ---------------- */
				int aromatic_bond_count = 0;
				/* compare the second order atoms' connection-bonds */
				short atomID = connAtoms[i][j];
				for(short k = connBond[atomID].size() - 1; k >= 0; k--){
					if(connBond[atomID][k] == Bond::AROMATIC_BOND){ /* aromatic bond */
						if(aromatic_bond_count < 2){ /* the aromatic bond does not shared with other */
							zagreb_M2_valence += 1.5; /* index of aromatic bond equal to 1.5 */
							aromatic_bond_count++; /* count aromatic bonds */
						}else
							zagreb_M2_valence += 1; /* if the aromatic bond shared with other */
					}
					else /* others bond type */
						zagreb_M2_valence += connBond[atomID][k];
				}

				/* take out atom-i's bonds that connection the second order atoms */
				if(connBond[i][j] == Bond::AROMATIC_BOND){ /* aromatic bond */
					if(aromatic_bond_count1 < 2){
						zagreb_M2_valence -= 1.5;  /* index of aromatic bond equal to 1.5 */
						aromatic_bond_count1++; /* count aromatic bonds */
					}else
						zagreb_M2_valence -= 1; /* the aromatic bond shared with other */
				}
				else /* others bond type */
					zagreb_M2_valence -= connBond[i][j];
			}

			/* calculate valence electrons */
			Atom & atom =mol.getAtom(i);
			double delta_v = valenceElections(atom.getAtomSymbol(), atom.getCharge(), i, zagreb_M2_valence);

			Zagreb_M2_v += delta_v * delta_v; /* calculate and collection Zagreb_M2_v */
		}
		mol.setFingerprint(Zagreb_M2_v,::Zagreb_M2_v);
	}catch(exception &ex){
		cout << ex.what() << " thrown in Zagreb_M2_v"<< endl;
	}
}

void Calculator8::Zagreb(){
	try{
		if(mol.isFeaturesReady(::Zagreb)) return;

		/* calculate the first Zagreb Index M1: Zagreb_M1 */
		short numBonds =mol.getBondsNum(false);
		int Zagreb_M1 = 2 * numBonds; // First Zagreb Index M1: Zagreb_M1
		if(numAtoms > 2)
			Zagreb_M1 += vertexDistanceCounts[2];

		/* calculate the first Zagreb Index M2: Zagreb_M2 */
		int Zagreb_M2 = 0; // second Zagreb index: Zagreb_M2
		DistanceMatrix<short> & distanceMatrix =mol.getDistanceMatrix(false);
		for(short i = 0; i < numAtoms; i++){
			int n = 0;
			for(short j = 0; j< numAtoms; j++){
				/* compare adjacent atoms to the ith atom as second order path count */
				/* follow the definition of Molecular Descriptors using the formula of
				path-count=2 */
				if(distanceMatrix.getValueInCell(i, j) == 2)
					n++; /* collection the second order path count */
			}
			Zagreb_M2 += n * n; /* calculate Zagreb_M2 */
		}

		/* calculate Zagreb_M1_v */
		double Zagreb_M1_v = 0.0; //first Zagreb index by valence vertex degrees: Zagreb_M1_v

		vector<float> & cbo =mol.getConventionalBondOrder();
		for(short i = 0; i < numAtoms; i++){
			/* ---------------------- calculate the conventional bond order ---------------- */
			double Zagreb_M1_valence = cbo[i];

			/* ------------------------- calculate valence electrons ------------------------ */
			Atom & atom =mol.getAtom(i);
			// valence vertex degree
			double delta_v = valenceElections(atom.getAtomSymbol(), atom.getCharge(), i, Zagreb_M1_valence);
			Zagreb_M1_v += delta_v * delta_v; /* calculate and collection Zagreb_M1_v */
		}

		/* Quadratic index: Q_index */
		/* Vertex degree topological index: VDTI */
		int	vdti = 3 - 2 * numAtoms + Zagreb_M1 / 2;

		vector<double> fingerprint(4, 0);
		fingerprint[0] = Zagreb_M1;
		fingerprint[1] = Zagreb_M1_v;
		fingerprint[2] = Zagreb_M2;
		fingerprint[3] = vdti;
		mol.setFingerprint(fingerprint,::Zagreb);
	}catch(exception &ex){
		cout << ex.what() << " thrown in Zagreb"<< endl;
	}
}

void Calculator8::doAll(){
	BW_MW_VW_EWeighted();
	PolarizabilityWeighted();
	BalabanAVDCI_BalabanTypeWeighted();
	SumWeighted();
	BalabanShortPathIndex();
	SecondMoharIndex();
	SpanningTree();
	LaplacianMatrix();

	Zagreb();
	MeanTotalDistance();
	VertexDegree();
	Zagreb_M2_v();
	SchultzMTI_VVD();
}
