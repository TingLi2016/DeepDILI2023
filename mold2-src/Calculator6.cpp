#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Constants.h"
#include "Calculator6.h"

Calculator6::Calculator6(Molecule & _mol)
:mol(_mol){
	doAll();
}

Calculator6::~Calculator6(void){
}

void Calculator6::VanderwaalsVol(void){
	try{
		if(mol.isFeaturesReady(::vanderwaalsVol)) return;

		double SVDW_C;	/* sum of atomic Van Der Waals Carbon-Scele: SVDW_C */
		double MVDW_C;	/* mean of atoic van der Waals Carbon-scale: MVDW_C */

		double Pe;		/* Pauling electronegativities */
		double SEP_C;	/* sum of atomic electronegativites Pauling-Scale on Carbon: SEP_C */
		double MEP_C;	/* mean of atomic electronegativities Pauling-scaled on Carbon: MEP_C */

		double Se;		/* Sanderson electronegativities */
		double sum_Se;	/* sum of atomic electrongativities Sanderson-scaled on Carbon: SES_C */
		double mean_Se;	/* mean of atoic electronegativity Sanderson-scaled on Carbon: MES_C */

		double Ae;		/* Allred-Rochow electronegativities */
		double SEAR_C;	/* sum of atomic  electrongativities Allred-Rochow-scaled on Carbon: SEAR_C */
		double MEAR_C;	/* mean of atomic  electrongativity Allred-Rochow-scaled on Carbon: MEAR_C */

		/* Normal descriptor calculation. */
		SVDW_C =0;
		SEP_C = 0;
		sum_Se = 0;
		SEAR_C = 0;

		short numAtoms = mol.getAtomsNum(true);
		vector<short> & atoms = mol.countElements();
		for(short j = 0; j < Elements::TOTALELEMENTS; j++){
			if(atoms[j] > 0){
				//calculate sum of atomic Van Der Walls volume 
				double radius = Constants::vdw_radius[j] / 1.70; // scaled on carbon atom, vdw_radius_carbon = 1.70 
				double volume = radius * radius * radius * atoms[j]; // atomic van der waals volume
				SVDW_C += volume;

				// calculate sum of atomic electronegativities
				// calculate sum of atomic Pauling electronegativities 
				Pe = (Constants::elect_negativities[j][1] / 2.55) * atoms[j]; // scaled on carbon atom 
				SEP_C += Pe;
				// calculate sum of atomic Sanderson electronegativities 
				Se = (Constants::elect_negativities[j][2] / 2.47) * atoms[j]; // scaled on carbon atom 
				sum_Se += Se;
				// calculate sum of atomic Allred-Rochow electronegativities 
				Ae = (Constants::elect_negativities[j][3] / 2.5) * atoms[j]; // scaled on carbon atom 
				SEAR_C += Ae;
			}
		}

		// calculate and print out mean atomic van der waals volume 
		MVDW_C = SVDW_C / numAtoms;
		MEP_C = SEP_C / numAtoms;
		mean_Se = sum_Se / numAtoms;
		MEAR_C = SEAR_C / numAtoms;

		//print out sum of atomic van der waals volume 
		vector<double> fingerprint(8, 0.0);
		fingerprint[0] = SVDW_C;
		fingerprint[1] = MVDW_C;
		fingerprint[2] = SEP_C;
		fingerprint[3] = MEP_C;
		fingerprint[4] = sum_Se;
		fingerprint[5] = mean_Se;
		fingerprint[6] = SEAR_C;
		fingerprint[7] = MEAR_C;
		mol.setFingerprint(fingerprint, ::vanderwaalsVol);
	}catch(exception & ex){
		cout << ex.what() << " thrown in vanderwaalsVol"<< endl;
	}
}


void Calculator6::AtomsAndBonds(void){
	try{
		if(mol.isFeaturesReady(::AtomsAndBonds))return;

		vector<short> & atoms = mol.countElements();
		vector<double> fingerPrint(12, 0.0f);
		short numBonds = mol.getBondsNum(true);
		short numAtoms = mol.getAtomsNum(true);
		fingerPrint[0] = numAtoms;  // total number of atoms 
		fingerPrint[1] = numAtoms - atoms[Elements::symbolToId("H")]; // number of non-H atoms
		fingerPrint[2] = numBonds;  // number of bonds
		fingerPrint[3] = numBonds - atoms[Elements::symbolToId("H")]; // number of non-H bonds
		fingerPrint[4] = numBonds - numAtoms + 1; // number of rings

		/* count number of double, triple, aromatic, bonds in each molecular */
		for(short i = 0; i < numBonds; i++){
			if(mol.getBond(i).getBondType() == Bond::TRIPLE_BOND){
				fingerPrint[5] += 1.0; // number of triple-bonds
			}
		}

		/* computing and printout Halogen atomsch */

		fingerPrint[6] = atoms[Elements::symbolToId("F")] + 
			atoms[Elements::symbolToId("Cl")] + 
			atoms[Elements::symbolToId("Br")] + 
			atoms[Elements::symbolToId("I")] + 
			atoms[Elements::symbolToId("At")];

		/* Molecular size index: MSI */
		fingerPrint[7] = numAtoms * log10(double(numAtoms)) / Constants::LOG20;

		/*------------ Atomic Composition Index: ACI and print it to output file ------------------*/
		double tempaci = 0.0; /* hold part of aci */
		for(short i = atoms.size() - 1; i >= 0; i--){		
			if(atoms[i] > 0){ /* found same type of atoms storing in the cell */
				tempaci += atoms[i] * (log10(double(atoms[i])));
				/* computing and print: Mean Value of Atomic Composition Index: MACI */
				fingerPrint[9] += -double(atoms[i]) / numAtoms * (log10(double(atoms[i])) - log10(double(numAtoms)));
			}
		}
		tempaci /= Constants::LOG20;
		fingerPrint[9] /= Constants::LOG20;

		fingerPrint[8] = fingerPrint[7] - tempaci; // aci
		/* computing and printout Branch index: BI */
		vector<short> temp(numAtoms, 0);
		for(short i = 0; i < numBonds; i++){
			Bond & bond = mol.getBond(i);
			temp[bond.getFirstAtomNumber()]++;  /* search, collection and sort first_atoms_number */
			temp[bond.getSecondAtomNumber()]++; /* search, collection and sort second_atom_number */
		}

		/* computing the SIGMA-BONDs(the sigma bonds as vertex degree), BI and MSCI */
		double number = 1;
		for(short i = 0; i < numAtoms; i++){
			short n = temp[i]-(mol.getAtom(i).getHydrogenCountPlus1() - 1); // computing SIGMA BONDs
			if( n > 0) number *= n; // factorial of sigma-bonds; calculating for MSCI
			if(n > 2) fingerPrint[10] += n - 2;   // all the vertex degrees greater than 2, computing BI 
		}

		/* Molecular structure connectivity index: MSCI */
		number = fabs(number);
		fingerPrint[11] = 1.0 / sqrt(number);  /* calculating MSCI */
		mol.setFingerprint(fingerPrint, ::AtomsAndBonds);
	}catch(exception & ex){
		cout << ex.what() << " thrown in AtomsAndBonds"<< endl;
	}
}

void Calculator6::polarizability(void){
	try{
		if(mol.isFeaturesReady(::polarizability))return;

		/* Normal descriptor calculation. */
		short numAtoms = mol.getAtomsNum(true);
		//i is the atom order number, j is the bond type
		vector< vector<short> > & atom2bonds = mol.getAtom2bonds(true);

		double pol=0; // polarizability: this function does not contain a NITROGEN. It's different from one in Molecule. But I don't know why
		for(short i = 0; i < numAtoms; i++){
			char symbol = mol.getAtom(i).getAtomSymbol();
			if(symbol == Elements::CARBON){
				if(atom2bonds[i][3] == 1 || atom2bonds[i][2] == 2){ // use SP
					pol += Constants::polarizability_value[3];
				}
				if(atom2bonds[i][4] > 0 || atom2bonds[i][2] == 1){ // use SP2
					pol += Constants::polarizability_value[2];
				}else{ // use SP3
					pol += Constants::polarizability_value[1];
				}
			}else if(symbol == Elements::HYDROGEN){
				pol += Constants::polarizability_value[0];
			}else if(symbol == Elements::OXYGEN){
				if (atom2bonds[i][2]>0){ //sp2
					pol += Constants::polarizability_value[8];
				}else{ //sp3
					pol += Constants::polarizability_value[7];
				}
			}else if(symbol == Elements::SULFUR){
				if (atom2bonds[i][2]>0){ //sp2
					pol += Constants::polarizability_value[10];
				}else{ //sp3
					pol += Constants::polarizability_value[9];
				}
			}else if(symbol == Elements::FLOURINE){
				pol += Constants::polarizability_value[12];
			}else if(symbol == Elements::CLORINE){
				pol += Constants::polarizability_value[13];
			}else if(symbol == Elements::BROMINE){
				pol += Constants::polarizability_value[14];
			}else if(symbol == Elements::IODINE){
				pol += Constants::polarizability_value[15];
			}
		}

		/* scaled on Carbon atom: C-sp3 */
		pol /= Constants::polarizability_value[1];

		/* calculation mean atomic polarizability */
		double pol_mean; //average polarizability
		if(numAtoms != 0)
			pol_mean= pol / numAtoms;
		else
			pol_mean= Constants::ERROR_SIGNAL;
		vector<double> fingerprint(2, 0);
		fingerprint[0] = pol;
		fingerprint[1] = pol_mean;
		mol.setFingerprint(fingerprint,::polarizability);
	}catch(exception & ex){
		cout << ex.what() << " thrown in polarizability"<< endl;
	}
}

void Calculator6::ElectrotopologicalVariation(void){
	try{
		if(mol.isFeaturesReady(::ElectrotopologicalVariation))return;

		short numAtoms = mol.getAtomsNum(false);
		vector<short> & pqn = mol.getPqn();
		vector<float> & valence = mol.getValence();
		/* -------------------------- intrinsic state --------------------------------------- */
		vector<double> is(numAtoms);
		vector< vector<short> >& connAtoms = mol.getConnectedAtomIDs(false);
		for(short i = 0; i < numAtoms; i++){
			is[i] = ((4.0/pqn[i]/pqn[i]) * valence[i] + 1);
			if(numAtoms > 1) is[i] /= connAtoms[i].size();
		}

		/* ------------------------- electrotopological state index: SES ---------------------- */
		double ses = 0.0;
		vector<double> s(is);
		DistanceMatrix<short> & distanceMatrix = mol.getDistanceMatrix(false);
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				short nP = distanceMatrix.getValueInCell(i, j) + 1;
				s[i] += (is[i] - is[j]) / (nP * nP);
			}

			/* --------------- sum electrotopological states: SES ------------------- */
			ses += s[i]; 
		}

		double mvvenv = -1.0e-9;	/* initial to a very small negative value */
		double mvvepv = -1.0e-9;	/* initial to a very small negative value */
		double maenv = 0.0;
		vector<double> minus(numAtoms, 0);
		for(short i = 0; i < numAtoms; i++){
			minus[i] = s[i] - is[i];
			maenv += fabs(minus[i]);						/* electrotopological variation */
		}
		mvvepv = *max_element(minus.begin(), minus.end());
		mvvenv = *min_element(minus.begin(), minus.end())*(-1.0);

		/*------------------------- calculate ETI ------------------------ */
		double eti = 0.0;
		short numBonds = mol.getBondsNum(false);
		for(short i = 0; i < numBonds; i++){	/* running over all bonds */
			Bond &bond = mol.getBond(i);	
			short k = bond.getFirstAtomNumber();
			short j = bond.getSecondAtomNumber();
			double skj = s[k] * s[j];
			if(skj > 0.0000001)					/* take the positive value */
				eti += 1.0 / sqrt(skj);	/* calculate the first part of ETI */
		}
		int nr = numBonds - numAtoms + 1;	/* number of rings */
		eti *= numBonds / (nr + 1.0);	/* calculate ETI */

		/* ------------------------- mean electrotopological sitate index: MESI ---------------------- */
		double mesi = ses / numAtoms;
		vector<double> fingerprint(6, Constants::ERROR_SIGNAL);
		fingerprint[0] = mvvenv;
		fingerprint[1] = mvvepv;
		fingerprint[2] = maenv;
		fingerprint[3] = eti;
		fingerprint[4] = ses;
		fingerprint[5] = mesi;
		mol.setFingerprint(fingerprint, ::ElectrotopologicalVariation);
	}catch(exception & ex){
		cout << ex.what() << " thrown in ElectrotopologicalVariation"<< endl;
	}
}

void Calculator6::EigenvalueLowest(void){
	try{
		if(mol.isFeaturesReady(::eigenvalueLowest))return;
		short numAtoms = mol.getAtomsNum(false);
		//burdenmatrix[][]: the off-diagonal elements of burdenmatrix[][] is bond-type dependent
		Array2D<double> burdenmatrix(numAtoms, numAtoms, 0.001);
		initializeBurdenMatrix_Lowest(burdenmatrix);

		vector<short> & atomicNum = mol.getAtomicNum(false);
		for(short i = 0; i < numAtoms; i++)
			burdenmatrix[i][i] = atomicNum[i]; //atomic order numbers

		//jacobi starts at 1 not 0. copy burdenmatrix[i][j] to burdVDW[i][j] and burdSand[][]
		Array2D<double> burdVDW(numAtoms, numAtoms);
		Array2D<double> burdSand(numAtoms, numAtoms);		
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				burdSand[i][j] = burdVDW[i][j] = burdenmatrix[i][j];//copy to burdVDW[][]
			}
		}

		for(short i = 0; i < numAtoms; i++){
			int symbol = mol.getAtom(i).getAtomSymbol();
			burdVDW[i][i] = Constants::vdw_volume[symbol]; //diagonal elements put vdw volume
			burdSand[i][i] = Constants::elect_negativities[symbol][1];
			burdenmatrix[i][i] = symbol + 1; //atomic order numbers
		}

		//then call jacobi to get the eigenvalue
		vector<double> eigenV = getEigenvalue(burdenmatrix);
		vector<double> eigenVvdw = getEigenvalue(burdVDW);
		vector<double> eigenVsand = getEigenvalue(burdSand);

		//sort eigenvalue, in order from min to max
		sort(eigenV.begin(), eigenV.end());
		sort(eigenVvdw.begin(), eigenVvdw.end());
		sort(eigenVsand.begin(), eigenVsand.end());

		///////////////////////////////////////////////////
		vector<double> fingerprint(24, 0);
		short k = 0;
		for(short i = eigenV.size() - 1; i >= 0 && k < 8; i--)
			fingerprint[k++] = eigenV[i];

		k = 8;
		for(short i = eigenVvdw.size() - 1; i >= 0 && k < 16; i--)
			fingerprint[k++] = eigenVvdw[i];

		k = 16;
		for(short i = eigenVsand.size() - 1; i >= 0 && k < 24; i--)
			fingerprint[k++] = eigenVsand[i];
		mol.setFingerprint(fingerprint, ::eigenvalueLowest);
	}catch(exception & ex){
		cout << ex.what() << " thrown in eigenvalueLowest"<< endl;
	}
}

void Calculator6::EigenvaluePolLowest(void){
	try{
		if(mol.isFeaturesReady(::eigenvaluePolLowest))return;
		short numAtoms = mol.getAtomsNum(false);
		//create burdenmatrix[][]
		Array2D<double> burdenmatrix(numAtoms, numAtoms, 0.001);
		initializeBurdenMatrix_Lowest(burdenmatrix);
		vector< vector<short> > & a2b = mol.getAtom2bonds(false);
		//the off-diagonal elements of burdenmatrix[][] is bond-type dependent
		for(short i = 0; i < numAtoms; i++){
			short vid = getPolarizabilityValueID(a2b, mol.getAtom(i).getAtomSymbol(), i);
			if(vid > 0)
				burdenmatrix[i][i] = Constants::polarizability_value[vid];		
		}

		/* then call jacobi to get the eigenvalue: eigen[] */
		vector<double> eigenV = getEigenvalue(burdenmatrix);

		//sort eigenvalue, in order from min to max
		sort(eigenV.begin(), eigenV.end());
		vector<double> fingerprint(8, 0);
		short k = 0;
		short size = eigenV.size() - 1;
		for(short i = size > 7 ? 7 : size; i >= 0; i--)
			fingerprint[k++] = eigenV[i];
		mol.setFingerprint(fingerprint, ::eigenvaluePolLowest);
	}catch(exception & ex){
		cout << ex.what() << " thrown in eigenvaluePolLowest"<< endl;
	}
}

void Calculator6::EigenvalueHighest(void){
	try{
		if(mol.isFeaturesReady(::eigenvalueHighest)) return;

		short numAtoms = mol.getAtomsNum(false);
		//create burdenmatrix[][]
		Array2D<double> burdenmatrix(numAtoms, numAtoms, 0.0);
		//the off-diagonal elements of burdenmatrix[][] is bond-type dependent
		initializeBurdenMatrix_Highest(burdenmatrix);

		Array2D<double> burdVDW(numAtoms, numAtoms, 0.0);
		Array2D<double> burdSand(numAtoms, numAtoms, 0.0);		
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				burdSand[i][j] = burdVDW[i][j] = burdenmatrix[i][j];//copy to burdVDW[][]
			}
		}

		for(short i = 0; i < numAtoms; i++){
			int symbol = mol.getAtom(i).getAtomSymbol();
			burdVDW[i][i] = Constants::vdw_volume[symbol]; //diagonal elements put vdw volume
			burdSand[i][i] = Constants::elect_negativities[symbol][1];
			burdenmatrix[i][i] = symbol + 1; //atomic order numbers
		}

		//then call jacobi to get the eigenvalue
		vector<double> eigenV = getEigenvalue(burdenmatrix);
		vector<double> eigenVvdw = getEigenvalue(burdVDW);
		vector<double> eigenVsand = getEigenvalue(burdSand);

		//absolute eigenvalue
		for(short i = 0; i< numAtoms; i++){
			eigenV[i] = fabs(eigenV[i]);
			eigenVvdw[i] = fabs(eigenVvdw[i]);
			eigenVsand[i] = fabs(eigenVsand[i]);
		}

		//sort eigenvalue, in order from min to max
		sort(eigenV.begin(), eigenV.end());
		sort(eigenVvdw.begin(), eigenVvdw.end());
		sort(eigenVsand.begin(), eigenVsand.end());

		///////////////////////////////////////////////////
		vector<double> fingerprint(24, 0);
		short k = 0;
		for(short i = eigenV.size() - 1; i >= 0 && k < 8; i--)
			fingerprint[k++] = eigenV[i];

		k = 8;
		for(short i = eigenVvdw.size() - 1; i >= 0 && k < 16; i--)
			fingerprint[k++] = eigenVvdw[i];

		k = 16;
		for(short i = eigenVsand.size() - 1; i >= 0 && k < 24; i--)
			fingerprint[k++] = eigenVsand[i];
		mol.setFingerprint(fingerprint, ::eigenvalueHighest);
	}catch(exception & ex){
		cout << ex.what() << " thrown in eigenvalueHighest"<< endl;
	}
}

void Calculator6::EigenvaluePolHighest(void){
	try{
		if(mol.isFeaturesReady(::eigenvaluePolHighest))return;

		short numAtoms = mol.getAtomsNum(false);
		//create burdenmatrix[][]
		Array2D<double> burdenmatrix(numAtoms, numAtoms, 0.0);
		initializeBurdenMatrix_Highest(burdenmatrix);
		vector< vector<short> > & a2b = mol.getAtom2bonds(false);
		for(short i = 0; i < numAtoms; i++){
			short vid = getPolarizabilityValueID(a2b, mol.getAtom(i).getAtomSymbol(), i);
			if(vid >= 0)
				burdenmatrix[i][i] = Constants::polarizability_value[vid];		
		}

		/* then call jacobi to get the eigenvalue: eigen[] */
		vector<double> eigenV = getEigenvalue(burdenmatrix);

		//absolute eigenvalues 
		for(short i = 0; i< numAtoms; i++)
			eigenV[i] = fabs(eigenV[i]); //absolute value

		//sort eigenvalue, in order from min to max
		sort(eigenV.begin(), eigenV.end());
		vector<double> fingerprint(8, 0);
		short k = 0;
		for(short i = eigenV.size() - 1; i >= 0 && k < 8; i--)
			fingerprint[k++] = eigenV[i];
		mol.setFingerprint(fingerprint, ::eigenvaluePolHighest);
	}catch(exception & ex){
		cout << ex.what() << " thrown in eigenvaluePolHighest"<< endl;
	}
}

void Calculator6::initializeBurdenMatrix_Lowest(Array2D<double>& burdenmatrix){	
	short numBonds = mol.getBondsNum(false);
	vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);

	for(short i = 0; i < numBonds; i++){
		Bond & bond = mol.getBond(i);
		short m = bond.getFirstAtomNumber();
		short n = bond.getSecondAtomNumber();
		char bondType = bond.getBondType();

		double value = 0.15;
		if(bondType != Bond::AROMATIC_BOND) value = 0.1 * bondType;

		burdenmatrix[n][m] = burdenmatrix[m][n] = value;

		if(connAtoms[n].size() == 1 || connAtoms[m].size() == 1){
			burdenmatrix[m][n] += 0.01;
			burdenmatrix[n][m] += 0.01;
		}
	}
}

void Calculator6::initializeBurdenMatrix_Highest(Array2D<double>& burdenmatrix){
	//the off-diagonal elements of burdenmatrix[][] is bond-type dependent
	short numBonds = mol.getBondsNum(false);
	for(short i = 0; i < numBonds; i++){
		Bond & bond = mol.getBond(i);
		short m = bond.getFirstAtomNumber();
		short n = bond.getSecondAtomNumber();
		char bondType = bond.getBondType();

		double value = sqrt(1.5);
		if(bondType != Bond::AROMATIC_BOND)
			value = sqrt(double(bondType));

		burdenmatrix[n][m] = burdenmatrix[m][n] = value;
	}
}

short Calculator6::getPolarizabilityValueID(vector< vector<short> > & a2b, char symbol , short i){
	short vid = -1;	
	if(symbol == Elements::CARBON){ /* the atom is a C */
		if(a2b[i][3] == 1 || a2b[i][2] == 2){
			vid = 3; // use SP				
		}
		if(a2b[i][4] > 0 || a2b[i][2] == 1){
			vid = 2; // use SP2
		}else{
			vid = 1; // use SP3
		}
	}else if(symbol == Elements::HYDROGEN){ /* the atom is a H */
		vid = 0;
	}else if(symbol == Elements::OXYGEN){ /* the atom is an O */
		if (a2b[i][2]>0){
			vid = 8; //sp2
		}else{
			vid = 7; //sp3
		}
	}else if(symbol == Elements::SULFUR){ /* the atom is a S */
		if (a2b[i][2]>0){
			vid = 10; //sp2
		}else{
			vid = 9; //sp3
		}
	}else if (symbol == Elements::PHOSPHORUS){ /* the atom is a P */
		vid = 11;
	}else if (symbol == Elements::FLOURINE){ /* the atom is a F */
		vid = 12;
	}else if (symbol == Elements::CLORINE){ /* the atom is a Cl */
		vid = 13;
	}else if (symbol == Elements::BROMINE){ /* the atom is a Br */
		vid = 14;
	}else if (symbol == Elements::IODINE){ /* the atom is an I */
		vid = 15;
	}
	return vid;
}

vector<float> Calculator6::calculatePol(void){
	short numAtoms = mol.getAtomsNum(true);
	vector<float> pol(numAtoms, 0);
	vector< vector<short> > & a2b = mol.getAtom2bonds(true);
	for(short i = 0; i < numAtoms; i++){
		short vid = getPolarizabilityValueID(a2b, mol.getAtom(i).getAtomSymbol(), i);
		if(vid >= 0)
			pol[i] = Constants::polarizability_value[vid];
	}
	return pol;
}

void Calculator6::Geary_MoranAutoCorrelation(){
	try{
		if(mol.isFeaturesReady(::Geary_MoranAutoCorrelation)) return;

		short numAtoms = mol.getAtomsNum(true);

		//count mean molecular weight, vdw volume and Sanderson electronegativities including H
		vector<short> atomics = mol.countElements();
		double molecular_weight = 0.0;
		double volume = 0.0;
		double Se_elect = 0.0;
		for(vector<short>::size_type i = 0; i < atomics.size(); i++){
			if(atomics[i] > 0){
				molecular_weight += (Constants::atom_weight[i] * atomics[i]);
				volume += (Constants::vdw_volume[i] * atomics[i]);
				Se_elect += (Constants::elect_negativities[i][1] * atomics[i]);
			}
		}

		vector<float> pol = calculatePol();// mol.getPol();  //TODO:It's weird:Geary_pol no N
		double polar = sum_element(pol);

		//calculate the average value including H
		double mean_m = molecular_weight / numAtoms;
		double mean_v = volume / numAtoms;
		double mean_se = Se_elect / numAtoms;
		double mean_p = polar / numAtoms;

		// check the value in distance_matrix[][] and calculate Geary autocorrelation

		double gtsa_m_down = 0;
		double gtsa_vdw_down = 0;
		double gtsa_se_down = 0;
		double gtsa_p_down = 0;
		//count how many H connected to the atom i
		// I don't understand this part
		vector<double> wx(numAtoms, 0), wy(numAtoms, 0), volx(numAtoms, 0), voly(numAtoms, 0), se(numAtoms, 0);
		short noHnumAtoms = mol.getAtomsNum(false);
		vector<char> & num_of_H = mol.getCon_to_H();

		for(short m = 0; m < numAtoms ; m++){
			for(short n=0; n < numAtoms; n++){
				if(n < noHnumAtoms){
					int x = mol.getAtom(n).getAtomSymbol();
					wx[n] = Constants::atom_weight[x] + num_of_H[n]; // atomic weight of atom i + num of H connected to it
					gtsa_m_down +=  (wx[n] - mean_m)*(wx[n] - mean_m); //calculate autocorelation length-1 down part

					volx[n] = Constants::vdw_volume[x] + Constants::vdw_volume[0] * num_of_H[n]; // van der waals volume of atom i + atom H connected to i
					gtsa_vdw_down += (volx[n] - mean_v) * (volx[n] - mean_v); //calculate autocorelation down part

					se[n] = Constants::elect_negativities[x][1] + Constants::elect_negativities[0][1] * num_of_H[n];//Sanderson electronegativities of atom i + atom H connected to i
					gtsa_se_down += (se[n] - mean_se) * (se[n] - mean_se); //calculate autocorelation down part
					gtsa_p_down += (pol[n] - mean_p) * (pol[n] - mean_p);
				} // end for_n statement
			}
			int x = mol.getAtom(m).getAtomSymbol();
			wy[m] = Constants::atom_weight[x] + num_of_H[m];
			voly[m] = Constants::vdw_volume[x] + Constants::vdw_volume[0] * num_of_H[m];
			se[m] = Constants::elect_negativities[x][1] + Constants::elect_negativities[0][1] * num_of_H[m];
		}

		vector<double> gtsa_m_up(Constants::MAX_CONNECTIONS, 0), 
			gtsa_vdw_up(Constants::MAX_CONNECTIONS, 0), 
			gtsa_se_up(Constants::MAX_CONNECTIONS, 0),
			gtsa_p_up(Constants::MAX_CONNECTIONS, 0);
		vector<short> pair(Constants::MAX_CONNECTIONS, 0);

		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		size_t k = 0;
		for(short i = 0; i < noHnumAtoms ; i++){	
			for(short j = i + 1; j < noHnumAtoms; j++){
				short dis = distanceMatrix[k++];
				if(dis > 0 && dis < Constants::MAX_CONNECTIONS){
					gtsa_m_up[dis] += (wx[i] - wy[j]) * (wx[i] - wy[j]) * 2; //calculate autocorrelation length-1 upper part
					gtsa_vdw_up[dis] += (volx[i] - voly[j]) * (volx[i] - voly[j]) * 2;
					gtsa_se_up[dis] += (se[i] - se[j]) * (se[i] - se[j]) * 2;
					gtsa_p_up[dis] += (pol[i] - pol[j]) * (pol[i] - pol[j]) * 2;
					pair[dis] += 2;
				}
			}
		}

		double mtsa_m_down = 0;
		double mtsa_vdw_down = 0;
		double mtsa_se_down = 0;
		double mtsa_p_down = 0;

		for(short m = 0; m < noHnumAtoms ; m++){
			wy[m] = 0;
			voly[m] = 0;
			se[m] = 0;
			for(short n=0; n < noHnumAtoms; n++){
				wx[n] = 0;
				volx[n] = 0;
				se[n] = 0;
				int x = mol.getAtom(n).getAtomSymbol();
				wx[n] = Constants::atom_weight[x] + num_of_H[n]; // atomic weight of atom i + num of H connected to it
				mtsa_m_down +=  (wx[n] - mean_m)*(wx[n] - mean_m); //calculate autocorelation length-1 down part

				volx[n] = Constants::vdw_volume[x] + Constants::vdw_volume[0] * num_of_H[n]; // van der waals volume of atom i + atom H connected to i
				mtsa_vdw_down += (volx[n] - mean_v) * (volx[n] - mean_v); //calculate autocorelation down part

				se[n] = Constants::elect_negativities[x][1] + Constants::elect_negativities[0][1] * num_of_H[n];//Sanderson electronegativities of atom i + atom H connected to i
				mtsa_se_down += (se[n] - mean_se) * (se[n] - mean_se); //calculate autocorelation down part
				mtsa_p_down += (pol[n] - mean_p) * (pol[n] - mean_p);

			}
			int x = mol.getAtom(m).getAtomSymbol();
			wy[m] = Constants::atom_weight[x] + num_of_H[m];
			voly[m] = Constants::vdw_volume[x] + Constants::vdw_volume[0] * num_of_H[m];
			se[m] = Constants::elect_negativities[x][1] + Constants::elect_negativities[0][1] * num_of_H[m];
		}

		vector<short> pair1(Constants::MAX_CONNECTIONS, 0);
		vector<double> mtsa_m_up(Constants::MAX_CONNECTIONS, 0), 
			mtsa_vdw_up(Constants::MAX_CONNECTIONS, 0), 
			mtsa_se_up(Constants::MAX_CONNECTIONS, 0),
			mtsa_p_up(Constants::MAX_CONNECTIONS, 0);

		k = 0;
		for(short i = 0; i < noHnumAtoms ; i++){	
			for(short j = i + 1; j < noHnumAtoms; j++){
				short dis = distanceMatrix[k++];
				if(dis > 0 && dis < Constants::MAX_CONNECTIONS){
					mtsa_m_up[dis] += (wx[i] - mean_m) * (wy[j] - mean_m) * 2; //calculate autocorrelation length-1 upper part
					mtsa_vdw_up[dis] += (volx[i] - mean_v) * (voly[j] - mean_v) * 2;
					mtsa_se_up[dis] += (se[i] - mean_se) * (se[j] - mean_se) * 2;
					mtsa_p_up[dis] += (pol[i] - mean_p)* (pol[j] - mean_p) * 2;
				}
			}
		}

		//calculate GATS sum
		vector<double> GTSA_m(Constants::MAX_CONNECTIONS, 0), 
			GTSA_vdw(Constants::MAX_CONNECTIONS, 0), 
			GTSA_se(Constants::MAX_CONNECTIONS, 0),
			GTSA_p(Constants::MAX_CONNECTIONS, 0);

		//calculate MATS sum
		vector<double> MTSA_m(Constants::MAX_CONNECTIONS, 0), 
			MTSA_vdw(Constants::MAX_CONNECTIONS, 0), 
			MTSA_se(Constants::MAX_CONNECTIONS, 0),
			MTSA_p(Constants::MAX_CONNECTIONS, 0);

		for (short i=1; i< Constants::MAX_CONNECTIONS; i++){
			if(pair[i] !=0){
				if(gtsa_m_down != 0) GTSA_m[i] = (gtsa_m_up[i] * (numAtoms - 1) * numAtoms)/(2*gtsa_m_down * pair[i]);
				if(mtsa_m_down != 0) MTSA_m[i] = (mtsa_m_up[i] * noHnumAtoms * noHnumAtoms)/(mtsa_m_down * pair[i]);

				if(gtsa_vdw_down != 0) GTSA_vdw[i] = (gtsa_vdw_up[i] * (numAtoms - 1)* numAtoms)/(2*gtsa_vdw_down * pair[i]);
				if(mtsa_vdw_down != 0) MTSA_vdw[i] = (mtsa_vdw_up[i] * noHnumAtoms * noHnumAtoms)/(mtsa_vdw_down * pair[i]);

				if(gtsa_se_down != 0) GTSA_se[i] = (gtsa_se_up[i] * (numAtoms - 1) * numAtoms)/(2*gtsa_se_down * pair[i]);
				if(mtsa_se_down != 0) MTSA_se[i] = (mtsa_se_up[i] * noHnumAtoms * noHnumAtoms)/(mtsa_se_down * pair[i]);

				if(gtsa_p_down != 0) GTSA_p[i] = (gtsa_p_up[i] * (numAtoms - 1) *  numAtoms)/(2*gtsa_p_down * pair[i]);
				if(mtsa_p_down != 0) MTSA_p[i] = (mtsa_p_up[i] * noHnumAtoms * noHnumAtoms)/(mtsa_p_down * pair[i]);		
			}
		}

		vector<double> fingerprint;
		fingerprint.insert(fingerprint.end(), GTSA_m.begin() + 1, GTSA_m.end());
		fingerprint.insert(fingerprint.end(), GTSA_vdw.begin() + 1, GTSA_vdw.end());
		fingerprint.insert(fingerprint.end(), GTSA_se.begin() + 1, GTSA_se.end());
		fingerprint.insert(fingerprint.end(), GTSA_p.begin() + 1, GTSA_p.end());

		fingerprint.insert(fingerprint.begin(), MTSA_m.begin() + 1, MTSA_m.end());
		fingerprint.insert(fingerprint.begin(), MTSA_vdw.begin() + 1, MTSA_vdw.end());
		fingerprint.insert(fingerprint.begin(), MTSA_se.begin() + 1, MTSA_se.end());
		fingerprint.insert(fingerprint.begin(), MTSA_p.begin() + 1, MTSA_p.end());
		mol.setFingerprint(fingerprint, ::Geary_MoranAutoCorrelation);
	}catch(exception & ex){
		cout << ex.what() << " thrown in Geary_MoranAutoCorrelation"<< endl;
	}
}

void Calculator6::InformationConent0_5(){
	try{
		if(mol.isFeaturesReady(::InformationContent0_5)) return;

		short numAtoms = mol.getAtomsNum(true);
		/* ----------------------------- calculate tatol bonds ---------------------------- */
		short numBonds = mol.getBondsNum(true);
		vector<float> bonds(numBonds, 0);
		for(short i = 0; i < numBonds; i++){
			short aromatic_bond_count = 0;
			if(mol.getBond(i).getBondType() == Bond::AROMATIC_BOND){
				if(aromatic_bond_count < 2){ /* has no one share the aromatic bond */
					bonds[i] += 1.5;/* add aromatic bond's value */
					aromatic_bond_count++;  /* count aromatic bond */
				}else
					bonds[i] += 1;  /* if the aromatic bond shared with other, value=1 */
			}else
				/* others bonds' type collection */
				bonds[i] += mol.getBond(i).getBondType();
		}
		vector<char> H_con(numAtoms, 0);
		for(int i = 0; i < numAtoms; i++)
			if(mol.getAtom(i).getAtomSymbol() == Elements::HYDROGEN)
				H_con[i] = 1;

		vector< vector<short> > &connH = mol.getConnectedAtomIDs(true);
		vector< vector<char> > & connBondsH = mol.getConnectedBondTypes(true);
		vector< vector<short> > &conn = mol.getConnectedAtomIDs(false);	
		vector< vector<int> > acm(6, vector<int>(numAtoms, 0));
		vector<short> & atomic_num = mol.getAtomicNum(true);

		for(short i = 0; i < 6; i++){		
			for(short j = 0; j < numAtoms; j++){
				short naj = connH[j].size();
				int m = naj * 1000000;	/* vertex degree, made a bigger number for easy calculate */
				int n = 1;
				int p = 0;
				switch(i){
				case 0:	/* 0-order */
					acm[i][j] = atomic_num[j];
					break;
				case 1:	/* 1-order */
					if(j < (short)conn.size())
						n += conn[j].size() * 10000;
					for(short k = 0; k < naj; k++){
						if(H_con[j] > 0)
							p += atomic_num[connH[j][k]];
					}
					acm[i][j] = acm[0][j] + m + n + p;
					break;
				case 2: /* 2-order */
					for(short k = 0; k < naj; k++){
						for(short q = connH[connH[j][k]].size() - 1; q >= 0; q--){
							n += connBondsH[connH[j][k]][q];
							p += atomic_num[connH[connH[j][k]][q]];
						}
					}
					acm[i][j] = acm[1][j] + m + n + p;
					break;
				case 3: /* 3-order */
					for(short k = 0; k < naj; k++){
						for(short q = connH[connH[j][k]].size() - 1; q >= 0; q--){
							for(short s = connH[connH[connH[j][k]][q]].size() - 1; s >= 0; s--){
								p += atomic_num[connH[connH[connH[j][k]][q]][s]];
							}
						}
					}
					acm[i][j] = acm[2][j] + p;
					break;
				case 4: /* 4-order */
					for(short k = 0; k < naj; k++){
						for(short q = connH[connH[j][k]].size() - 1; q >= 0; q--){
							for(short s = connH[connH[connH[j][k]][q]].size() - 1; s >= 0; s--){
								for(short t = connH[connH[connH[connH[j][k]][q]][s]].size() - 1; t >= 0; t--){
									p += atomic_num[connH[connH[connH[connH[j][k]][q]][s]][t]];
								}								
								n += connH[connH[connH[connH[j][k]][q]][s]].size();
							}
						}
					}
					acm[i][j] = acm[3][j] + n + p;
					break;
				case 5: /* 5-order */
					for(short k = 0; k < naj; k++){
						for(short q = connH[connH[j][k]].size() - 1; q >= 0; q--){
							for(short s = connH[connH[connH[j][k]][q]].size() - 1; s >= 0; s--){
								for(short t = connH[connH[connH[connH[j][k]][q]][s]].size() - 1; t >= 0; t--){
									for(short u = connH[connH[connH[connH[connH[j][k]][q]][s]][t]].size() - 1; u >= 0; u--){
										p += atomic_num[connH[connH[connH[connH[connH[j][k]][q]][s]][t]][u]];
									}
								}
							}
						}
					}
					acm[i][j] += acm[3][j] + p;
				}	
			}
		} /* end of FOR-i */

		/* -------------------------- calculate Information Content 0 ~ 5 ----------------------- */
		vector< vector<double> > ico(5, vector<double>(6, 0));
		double tmpV = sum_element(bonds);
		double v = tmpV > 0 ? log10(tmpV) / log10(2.0) : 0;
		for(short k = 0; k < 6; k++){
			/* ------------------- calculate the equivalence calsses ----------------- */
			int m = 0;
			vector< vector<short> > ec;
			for(short i = 0; i < numAtoms; i++){
				if(acm[k][i] != 0){	/* does not empty */
					/* working the first element */
					ec.push_back(vector<short>(1, i + 1));

					for(short j = i+1; j < numAtoms; j++){	/* working other elements */
						if(acm[k][i] == acm[k][j]){	/* found same element */
							ec[m].push_back(j + 1);	/* elements of the class */
							acm[k][j] = 0;	/* atom-j count already, set it to 0 for does not recount */
						}
					}
					acm[k][i] = 0;	/* atom-i count already, set it to 0 for does not recount */
					m++;	/* number of equivalence class */
				}
			}

			/* ---------------------- calculate ICO_0 ~ ICO_5 ------------------------------ */
			int np = ec.size();
			for(short i = 0; i < np; i++){
				double r = double(ec[i].size()) / numAtoms;
				ico[0][k] += -r * log10(r) / Constants::LOG20;
			}
			/* ---------------------- calculate TICO_0 ~ TICO_5 ---------------------------- */
			ico[1][k] = ico[0][k] * numAtoms;
			/* ---------------------- calculate SICO_0 ~ SICO_5 ---------------------------- */
			ico[2][k] = ico[0][k] * Constants::LOG20;
			if(numAtoms > 1)
				ico[2][k] /= log10(double(numAtoms));
			/* ---------------------- calculate CICO_0 ~ CICO_5 ---------------------------- */
			ico[3][k] = log10(double(numAtoms)) / Constants::LOG20 - ico[0][k];

			/* ---------------------- calculate BICO_0 ~ BICO_5 ---------------------------- */
			if(v != 0.0)
				ico[4][k] = ico[0][k] / v;
		}
		vector<double> fingerprint;
		/* printout all of the results to output file */
		for(short i = 0; i < 6; i++){
			for(short j = 0; j < 5; j++){
				fingerprint.push_back(ico[j][i]);
			}
		}
		mol.setFingerprint(fingerprint, ::InformationContent0_5);
	}catch(exception & ex){
		cout << ex.what() << " thrown in InformationContent0_5"<< endl;
	}
}

void Calculator6::SP_number(){
	try{
		if(mol.isFeaturesReady(::SP_Number)) return;
		short numAtoms = mol.getAtomsNum(false);
		vector<char> cType(numAtoms, 0);/* atomic conventional bond order of atomic symbol C
										cType[1]: number of sp
										cType[2]: number of sp2
										cType[3]: number of sp3
										cType[4]: number of aromatic bonds
										cType[5]: number of sp =c= bonds */

		/* ------------ count number of conventional bond order of atomic C ------------- */
		vector< vector<short> > & atom2bonds = mol.getAtom2bonds(false);
		for(short i = 0; i < numAtoms; i++){
			if(mol.getAtom(i).getAtomSymbol() == Elements::CARBON){ /* found atomic C */
				if(atom2bonds[i][4] >= 1){ /* one or more aromatic bond(s) */
					cType[i] = 4; /* aromatic bond */
				}else if(atom2bonds[i][3] == 1){ /* one triple bond */
					cType[i] = 1; /* SP */
				}else if(atom2bonds[i][2] == 2){ /* two double bonds */
					cType[i] = 5; /* sp2 =C= bond */
				}else if(atom2bonds[i][2] == 1){ /* sp2: the atom only has one double bond */
					cType[i] = 2; /* SP2 */
				}else if(atom2bonds[i][1] >= 1){
					cType[i] = 3; /* SP3 */
				}
			}
		}

		short num_allenes, num_p_sp3, num_unArom_sp2, num_Arom_sp2, num_t_sp, num_nt_sp, num_q_sp3, num_t_sp3, num_s_sp3, num_qr_sp3,num_tr_sp3, num_sr_sp3, num_p_sp2, num_s_sp2, num_t_sp2;
		num_allenes = num_p_sp3 = num_unArom_sp2 = num_Arom_sp2 = num_t_sp = num_nt_sp = num_q_sp3 = num_t_sp3 = num_s_sp3 = num_qr_sp3 = num_tr_sp3 = num_sr_sp3 = num_p_sp2 = num_s_sp2 = num_t_sp2 = 0;
		vector<char> & num_of_H = mol.getCon_to_H();
		vector<char> & ring = mol.getInRing();
		for(short i = 0; i < numAtoms; i++){
			switch(cType[i]){
		case 5: /* ----------- allenes groups ------------- */
			num_allenes++; /* number of allenes groups: n=C= */
			break;
		case 4: /* -------------- aromatic C --------------- */
			if(num_of_H[i] == 1)	/* there has one atomic H connected it */
				num_unArom_sp2++;	/* number of unsubstituted aromatic C(sp2) */
			else
				num_Arom_sp2++;		/* number of substituted aromatic C(sp2) */
			break;
		case 3:/* ------------------ sp3 ------------------- */
			if(ring[i] == 0){ /* atom i not in ring */
				if(num_of_H[i] == 0)		/* atom-i does not connected to atomic H */
					num_q_sp3++;			/* number of total quaternary C(sp3): nCq */
				else if(num_of_H[i] == 1)	/* atom-i connected one H */
					num_t_sp3++;			/* number of total tertiary C(sp3): nCt */
				else if(num_of_H[i] == 2)	/* atom-i connected two H */
					num_s_sp3++;			/* number of total secondary C(sp3): nCs */
				else if(num_of_H[i] == 3)	/* atom-i connected three H */
					num_p_sp3++;			/* number of total primary C(sp3): nCp */
			}else if(ring[i] == 1){ /* atom i on ring */
				if(num_of_H[i] == 0)	/* atom i does not connected to atomic H */
					num_qr_sp3++;		/* number of ring quaternary C(sp3): nCrR2 */
				else if(num_of_H[i] == 1) /* atom i has one H connected to it */
					num_tr_sp3++;		  /* number of ring tertiary C(sp3): nCrHR */
				else if(num_of_H[i] == 2) /* atom i has two H connected to it */
					num_sr_sp3++;		  /* number of ring secondary C(sp3): nCrH2 */
			}		
			break;
		case 2: /* ---------------- sp2 ------------------ */
			if(ring[i] == 0){ /* atom i not in ring */
				short na = mol.getConnectedAtomIDs(false)[i].size();
				if(na == 1)	/* the atom only has one bond */
					num_p_sp2++;	/* number of primary C(sp2): n=CH2 */
				else if(na == 2)	/* the atom has two bonds */
					num_s_sp2++;	/* number of secondary C(sp2): n=CHR */
				else if(na == 3)	/* the atom has three bonds */
					num_t_sp2++;	/* number of tertiary C(sp2): n=CR2 */
			}			
			break;
		case 1: /* ----------------- sp ------------------ */
			if(num_of_H[i] == 1) /* atom i has one H connected to it */
				num_t_sp++;		/* number of terminal C(sp): n#CH */
			else
				num_nt_sp++;	/* number of non-terminal C(sp): n#CR */
			break;
			}
		}

		vector<double> fingerprint;
		fingerprint.push_back(num_p_sp3);
		fingerprint.push_back(num_s_sp3);
		fingerprint.push_back(num_t_sp3);
		fingerprint.push_back(num_q_sp3);
		fingerprint.push_back(num_sr_sp3);
		fingerprint.push_back(num_tr_sp3);
		fingerprint.push_back(num_qr_sp3);
		fingerprint.push_back(num_unArom_sp2);
		fingerprint.push_back(num_Arom_sp2);
		fingerprint.push_back(num_p_sp2);
		fingerprint.push_back(num_s_sp2);
		fingerprint.push_back(num_t_sp2);
		fingerprint.push_back(num_allenes);
		fingerprint.push_back(num_t_sp);
		fingerprint.push_back(num_nt_sp);
		mol.setFingerprint(fingerprint, ::SP_Number);
	}catch(exception & ex){
		cout << ex.what() << " thrown in SP_Number"<< endl;
	}
}

void Calculator6::Empirical(void){
	try{
		if(mol.isFeaturesReady(::Empirical)) return;

		short numAtoms = mol.getAtomsNum(true);
		vector<char> con_to_H(numAtoms, 0);	/* number of atom-i connected H */

		/* ------------------ count how many H connected to the atom-i ---------------- */
		vector< vector<short> > & connAtomsH = mol.getConnectedAtomIDs(true);
		for(short i = 0; i < numAtoms; i++){
			for(short j = connAtomsH[i].size() - 1; j >= 0; j--){
				if(mol.getAtom(connAtomsH[i][j]).getAtomSymbol() == Elements::HYDROGEN){ /* atom-j is H */
					con_to_H[i]++; /* number of connecting to H plus one */
				}
			}
		}

		/* ------------------ distinguish atomic symbol of atom-i -------------- */
		short carbons = 0;/* atomic symbol C */
		short oh = 0;/* -OH */
		short sh = 0;/* -SH */
		short nh = 0;/* -NH */

		short noHnumAtoms = mol.getAtomsNum(false);
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		for(short i = 0; i < noHnumAtoms; i++){
			char symbol = mol.getAtom(i).getAtomSymbol();
			if(symbol == Elements::CARBON){	/* found atomic C */
				carbons++;
			}else if(symbol == Elements::NITROGEN){	/* the atom-i is N */
				if(con_to_H[i] >= 1 && connBond[i].size() == 1) /* atom-i has only bond and one H */
					nh++;	/* -NH */
			}else if(symbol == Elements::OXYGEN){	/* the atom-i is O */
				if(con_to_H[i] >= 1 && connBond[i].size() == 1 && /* atom-i has only bond and one or more H */
					mol.getAtom(connAtoms[i][0]).getAtomSymbol() != Elements::SULFUR) /* the O does not connecting to S */
					oh++;	/* -OH */
			}else if(symbol == Elements::SULFUR){	/* the atom-i is S */
				if(con_to_H[i] >= 1 && connBond[i].size() == 1) /* atom-i has only bond and one H */
					sh++;	/* -SH */
			}
		}

		/* calculate the number of independent cycles */
		/* ----------- 1. check the compound is saturated or unsaturated; 2. count aromatic bonds ----------- */
		double a = 0;
		short aromatic = 0; /* aromatic */
		for(short i = 0; i < mol.getBondsNum(false); i++){
			Bond & bond = mol.getBond(i);
			if(bond.getBondType() == Bond::AROMATIC_BOND){ /* the bond-i is an aromatic bond */
				a += 1.5;	/* value of aromatic is equal to 1.5 */
				aromatic++; /* count aromatic bonds */
			}else
				a += bond.getBondType();
		}

		double uicb = 0; /* saturated compound */
		if(a != mol.getBondsNum(false)){/* unsaturated compound */	
			if(a > 0)
				uicb = log10(a) / log10(2.0);
		}


		/* ------------------------- calculate hydrophilic factor index: HFI ------------------------------ */
		double hfi = oh + sh + nh;	/* hydrophilic groups */
		if(noHnumAtoms != 0)
			a = sqrt((2.0 * hfi) / (noHnumAtoms * noHnumAtoms));

		hfi = 1 + hfi;
		hfi *= log10(hfi) / log10(2.0);
		if(noHnumAtoms > 0)
			hfi += carbons * ((1.0/noHnumAtoms) * (log10(1.0/noHnumAtoms)/log10(2.0))) + a;
		hfi /= log10(1.0 + noHnumAtoms) / log10(2.0);

		/* ------------------------- calculate aromatic bonds ratio: ABR ------------------------------ */
		double abr = 0;
		if(mol.getBondsNum(false) != 0)
			abr = double(aromatic) / mol.getBondsNum(false);

		/* --------------------- printout results to output file ----------------------- */
		vector<double> fingerprint;
		fingerprint.push_back(uicb);
		fingerprint.push_back(hfi);
		fingerprint.push_back(abr);
		mol.setFingerprint(fingerprint, ::Empirical);
	}catch(exception & ex){
		cout << ex.what() << " thrown in Empirical"<< endl;
	}
}

void Calculator6::AutoCorrelation(void){
	try{
		if(mol.isFeaturesReady(::autocorrelation)) return;

		short numAtoms = mol.getAtomsNum(true);
		/*----- atomic weight, vdw volume, atomic elect-negativities, polarizability -----*/
		vector<double> TSA_m(Constants::MAX_CONNECTIONS, 0), TSA_vdw(Constants::MAX_CONNECTIONS, 0), TSA_se(Constants::MAX_CONNECTIONS, 0), TSA_p(Constants::MAX_CONNECTIONS, 0);

		/* check the value that is distance between two atoms and calculate autocorrelation */
		vector<double> w(numAtoms, 0), vol(numAtoms, 0), se(numAtoms, 0), p(numAtoms, 0);
		vector<float> & pol = mol.getPol();
		for(short i = 0; i< numAtoms; i++){
			int ai = mol.getAtom(i).getAtomSymbol();
			w[i] = Constants::atom_weight[ai] / Constants::atom_weight[5]; /* atomic weight scaled on carbon atom */ 
			vol[i] = Constants::vdw_volume[ai] / Constants::vdw_volume[5]; /* vdw volume scaled on carbon atom */ 
			se[i] = Constants::elect_negativities[ai][1] / Constants::elect_negativities[5][1]; /* atomic elect-negativities scaled on carbon atom */ 
			p[i] = pol[i] / Constants::polarizability_value[1]; /* atomic polarizability scaled on carbon atom_SP3 */
		}

		vector<short> & distanceMatrix = mol.getDistanceMatrix(true).getData();
		size_t k = 0;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j < numAtoms; j++){
				/* check the distance that between atom-i and atom-j */
				int dis = distanceMatrix[k++];
				if(dis > 0 && dis < Constants::MAX_CONNECTIONS){
					TSA_m[dis] += w[j] * w[i];
					TSA_vdw[dis] += vol[j] * vol[i]; 
					TSA_se[dis] += se[j] * se[i];
					TSA_p[dis] += p[j] * p[i];
				}
			} /* end for_j statement */
		}/* end for_i statement */

		vector<double> fingerprint;
		for (short i=1; i<Constants::MAX_CONNECTIONS; i++)
			fingerprint.push_back( TSA_m[i]);

		for (short i=1; i<Constants::MAX_CONNECTIONS; i++)
			fingerprint.push_back( TSA_vdw[i]);

		for (short i=1; i<Constants::MAX_CONNECTIONS; i++) 
			fingerprint.push_back( TSA_se[i]);

		for (short i=1; i<Constants::MAX_CONNECTIONS; i++) 
			fingerprint.push_back( TSA_p[i]);
		mol.setFingerprint(fingerprint, ::autocorrelation);
	}catch(exception & ex){
		cout << ex.what() << " thrown in autocorrelation"<< endl;
	}
}
void Calculator6::doAll(void){
	polarizability();	
	AtomsAndBonds();
	VanderwaalsVol();

	ElectrotopologicalVariation();
	EigenvalueLowest();
	EigenvaluePolLowest();
	EigenvalueHighest();
	EigenvaluePolHighest();	
	Geary_MoranAutoCorrelation();
	SP_number();
	Empirical();
	InformationConent0_5();	
	AutoCorrelation();
}
