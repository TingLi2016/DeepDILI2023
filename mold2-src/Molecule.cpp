#include "stdafx.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Atom.h"
#include "Bond.h"
#include "Elements.h"
#include "Molecule.h"
#include "Constants.h"
#include "Calculator1.h"
#include "Features.h"

Molecule::Molecule(void):
m_nChiralFlag(0), m_nStextEntries(0), m_ringIsReady(false), m_useSearchDoAllPath(false), m_numAtoms(0), m_numBonds(0){
	m_allPathFailed = false;
	m_maxBondType = 5;
}

Molecule::~Molecule(void){}

void Molecule::clear(void){
	m_nChiralFlag = 0;
	m_nStextEntries = 0;
	m_ringIsReady = false;
	m_useSearchDoAllPath = false;
	m_allPathFailed = false;
	m_deletePart = 0;

	m_FieldName.clear();
	m_ID = -1;
	path.clear();
	p.clear();
	a.clear();
	vc.clear();
	vvc.clear();
	pqvc.clear();

	clearPathsSet();
	pathLen.clear();

	allPath.clear();
	detourMatrix.clear();

	rings.clear();
	m_Header.clear();
  m_headerMoleculeId.clear();

	m_fingerprint.assign(TOTALFEATURES, ::numeric_limits<float>::signaling_NaN());
	m_FeatureStatus.assign(TOTALGROUPS, false);

	atomic_num.clear();
	atom2bonds.clear();
	pol.clear();
	inRing.clear(); // used to check if atom i is on ring
	valence.clear();
	pqn.clear();
	conventionalBondOrder.clear();
	con_to_H.clear();
	aromatic.clear();
	isA.clear();

	m_atom_list.clear();
	m_bond_list.clear();
	connectedAtomIDs.clear();
	connectedBondTypes.clear();
	distanceMatrix.clear();

	distancePath.clear();
	distancePathMatrix.clear();
	m_elementCounts.clear();
	m_maxBondType = 5;

	// only this five data structure are generated for both H-added and H-deleted moleculess
	connectedAtomIDs_H.clear();
	connectedBondTypes_H.clear();
	distanceMatrix_H.clear();
	atom2bonds_H.clear();
	atomic_num_H.clear();

	m_numAtoms = 0; // without H
	m_numBonds = 0; // without H
}

void Molecule::setNumAtoms(short nAtoms){
	m_numAtoms = nAtoms;
}

void Molecule::setNumBonds(short nBonds){
	m_numBonds = nBonds;
}

vector< vector<char> > & Molecule::getIsA(){
	if(isA.empty()){
		short numAtoms = getAtomsNum(false);
		isA.resize(numAtoms, vector<char>(10, 0));
		for(short i = 0; i < numAtoms; i++){
			char symbol = getAtom(i).getAtomSymbol();
			/* ------------------ distinguish atomic symbol of atom-i -------------- */
			if(symbol == Elements::CARBON){	/* found atomic a C */
				isA[i][C] = 1;	/* atom-i is C */
			}else if(symbol == Elements::NITROGEN){	/* the atom-i is a N */
				isA[i][N] = 1;	/* atom-i is N */
			}else if(symbol == Elements::OXYGEN){	/* the atom-i is an O */
				isA[i][O] = 1;	/* atom-i is O */
			}else if(symbol == Elements::PHOSPHORUS){	/* the atom-i is a P */
				isA[i][P] = 1; /* atom-i is P */
			}else if(symbol == Elements::SULFUR){	/* the atom-i is a S */
				isA[i][S] = 1;/* atom-i is S */
			}else if(symbol == Elements::FLOURINE || symbol == Elements::CLORINE || symbol == Elements::BROMINE || symbol == Elements::IODINE){
				isA[i][X] = 1;
				if(symbol == Elements::FLOURINE){
					isA[i][F] = 1;
				}else if(symbol == Elements::CLORINE){
					isA[i][Cl] = 1;
				}else if(symbol == Elements::BROMINE){
					isA[i][Br] = 1;
				}else if(symbol == Elements::IODINE){
					isA[i][I] = 1;
				}
			}
		}
	}
	return isA;
}

vector<char> & Molecule::getCon_to_H(void){
	if(con_to_H.empty()){
		/* ------------------ count how many H connected to the atom-i ---------------- */
		short na = getAtomsNum(true);
		con_to_H.resize(na, 0);
		vector< vector<short> > & connAtoms = getConnectedAtomIDs(true);
		for(short i = 0; i < na; i++){
			for(short j = connAtoms[i].size() - 1; j >= 0; j--){
				if(getAtom(connAtoms[i][j]).getAtomSymbol() == Elements::HYDROGEN){
					con_to_H[i]++; /* number of connecting to H plus one */
				}
			}
		}
	}
	return con_to_H;
}

vector<char> & Molecule::getIsAromatic(void){ // H-deleted
	if(aromatic.empty()){
		short numAtoms = getAtomsNum(false);
		aromatic.resize(numAtoms, 0);
		vector< vector<char> > & connBond = getConnectedBondTypes(false);
		for(short i = 0; i < numAtoms; i++){
			/* ----------------- check atom-i is aromatic atom or no ---------------- */
			for(short j = connBond[i].size() - 1; j >= 0; j--){ /* running over all bonds of atom-i */
				if(connBond[i][j] == Bond::AROMATIC_BOND){ /* aromatic bond */
					aromatic[i] = 1; /* the atom-i is aromatic atom and made a flag as YES */
					break;
				}
			}
		}
	}
	return aromatic;
}

vector<double> Molecule::getFeatures(short startLocation){
	vector<double> features;
	short len = FeatureGroupLength[startLocation];
	features.assign(m_fingerprint.begin() + startLocation, m_fingerprint.begin() + startLocation + len);
	return features;
}

bool Molecule::isFeaturesReady(short startLocation){
	return m_FeatureStatus[startLocation];
}

template<class T> void Molecule::setFingerprint(vector<T> & subFeatures, short startLocation){
	short size = subFeatures.size();
	for(short i = 0; i < size; i++){
		m_fingerprint[FeatureLocation[startLocation] + i] = (float)subFeatures[i];
	}
	m_FeatureStatus[startLocation] = true;
}

void Molecule::setFingerprint(double subFeature, short location){
	m_fingerprint[FeatureLocation[location]] = (float)subFeature;
}

void Molecule::printFeatures(ostream & out){
	for(short i = 0; i < TOTALFEATURES; i++){
		if(abs(m_fingerprint[i]) < numeric_limits<float>::infinity())
			out << "\t" << m_fingerprint[i];
		else
			out << "\t" << numeric_limits<float>::max();
	}
	out << endl;

	out.flush();
}

void Molecule::showFingerprint(ostream & out, short id){
	short length = id == 72? 1:FeatureLocation[id + 1] - FeatureLocation[id];
	for(int i = 0; i < length; i++)
		out << m_fingerprint[FeatureLocation[id] + i] << "\t";
	out << endl;
}

short Molecule::getPolarizabilityValueID(char symbol , short i){ //use H-added
	short vid = -1;
	vector< vector<short> > & a2b = getAtom2bonds(true);
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
	}else if(symbol == Elements::NITROGEN){ /* the atom is an N */
		if(a2b[i][3] == 1 || a2b[i][2] == 2){
			vid = 6; // use SP
		}
		if(a2b[i][4] > 0 || a2b[i][2] == 1){
			vid = 5; // use SP2
		}else{
			vid = 4; // use SP3
		}
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

void Molecule::initPqn_valence(void){
	short numAtoms = getAtomsNum(false);
	pqn.assign(numAtoms, 2);
	vector<float> & cbo = getConventionalBondOrder();
	valence.assign(cbo.begin(), cbo.end());
	vector< vector<char> > & connBond = getConnectedBondTypes(false);
	for(short i = 0; i < numAtoms; i++){
		/* ---------------------- calculate the conventional bond order ---------------- */
		char symbol = getAtom(i).getAtomSymbol();
		/* ------------------ calculate valence electrons --------------------------- */
		/* valence vertex degree = Sigma-electrons + Pi-electrons + lone pair electrons */

		if(symbol == Elements::NITROGEN){/* symbol N */
			valence[i] += 2;
		}else if(symbol == Elements::OXYGEN){/* symbol O */
			valence[i] += 4;
		}else if(symbol == Elements::FLOURINE){/* symbol F */
			valence[i] += 6;
		}
		/*valence vertex degree(Si) = valence vertex degree / (its atomic number - valence vertex degree - 1)*/
		else if(symbol == Elements::SILICON){/* symbol Si */
			valence[i] = valence[i] / (13 - valence[i]);
			pqn[i] = 3;
		}
		/* P, S, Cl, Br, I: valence vertex degree = (Sigma-electrons + Pi-electrons + lone pair electrons) /
		   (its atomic number - (Sigma-electrons + Pi-electrons + lone pair electrons) - 1) */
		else if(symbol == Elements::PHOSPHORUS){/* symbol P */
			if(valence[i] > 3)
				valence[i] = valence[i] / (14 - valence[i]);
			else
				valence[i] = (valence[i] + 2) / (12 - valence[i]);

			pqn[i] = 3;
		}else if(symbol == Elements::SULFUR){/* symbol S */
			short n = 0;
			short nBonds = connBond[i].size();
			for(short j = 0; j < nBonds; j++) /* compare bonds' type */
				n += connBond[i][j];
			if(n <= 2)
				valence[i] = (valence[i] + 4) / (11 - valence[i]);
			else if(n == 4)
				valence[i] = (valence[i] + 2) / (13 - valence[i]);
			else
				valence[i] = valence[i] / (15 - valence[i]);

			pqn[i] = 3;
		}else if(symbol == Elements::CLORINE){/* symbol Cl */
			valence[i] = (valence[i] + 6) / (10 - valence[i]);
			pqn[i] = 3;
		}else if(symbol == Elements::BROMINE){/* symbol Br */
			valence[i] = (valence[i] + 6) / (28 - valence[i]);
			pqn[i] = 4;
		}else if(symbol == Elements::IODINE){ /* symbol I */
			valence[i] = (valence[i] + 6) / (46 - valence[i]);
			pqn[i] = 5;
		}
	}
}

vector<float> & Molecule::getConventionalBondOrder(void){
	if(conventionalBondOrder.empty()){
		short numAtoms = getAtomsNum(false);
		vector< vector<char> > & connBond = getConnectedBondTypes(false);
		conventionalBondOrder.resize(numAtoms, 0);
		for(short i = 0; i < numAtoms; i++){
			int aromatic_bond_count = 0; //count aromatic bonds
			for(short j = connBond[i].size() - 1; j >= 0; j--){ /* compare bonds' type */
				if( connBond[i][j] == Bond::AROMATIC_BOND){ /* aromatic bond */
					if(aromatic_bond_count < 2){/* the aromatic bond does not shared with other */
						conventionalBondOrder[i] += 1.5;/* index of aromatic bond equal to 1.5 */
						aromatic_bond_count++;/* count aromatic bonds */
					}else
						conventionalBondOrder[i] += 1;/* if the aromatic bond shared with other */
				}else/* others bond type */
					conventionalBondOrder[i] += connBond[i][j];
			}
		}
	}
	return conventionalBondOrder;
}

vector<short> & Molecule::getPqn(void){
	if(pqn.empty()) initPqn_valence();
	return pqn;
}

vector<float> & Molecule::getValence(){
	if(valence.empty()) initPqn_valence();
	return valence;
}

vector<char> & Molecule::getInRing(){
	if(inRing.empty()){
		short numAtoms = getAtomsNum(false);
		inRing.resize(numAtoms, 0);
		vector< vector<short> > & rings = getRings();
		for(short j = rings.size() - 1; j >= 0; j--){
			for(short k = rings[j].size() - 1; k >= 0; k--){
				inRing[rings[j][k]] = 1; /* made a flag as 1 on the ring atom */
			}
		}
	}
	return inRing;
}

void Molecule::setRingStatus(bool status){
	m_ringIsReady = status;
}

vector< vector<short> > & Molecule::getRings(void){
	if(!m_ringIsReady){
		Calculator1 cal(*this);
		cal.AromaticBondConverter();
	}
	return rings;
}

vector< vector<short> > & Molecule::getAtom2bonds(bool withH){
	vector< vector<short> > & a2b = withH ? atom2bonds_H : atom2bonds;

	if(a2b.empty()){
		short numAtoms = getAtomsNum(withH);
		a2b.resize(numAtoms, vector<short> (m_maxBondType, 0));
		for(short i = getBondsNum(withH) - 1; i >= 0; i--){
			Bond &bond = getBond(i);
			a2b[bond.getFirstAtomNumber()][bond.getBondType()]++;
			a2b[bond.getSecondAtomNumber()][bond.getBondType()]++;
		}
	}
	return a2b;
}

vector<float> & Molecule::getPol(){
	if(pol.empty()){
		short numAtoms = getAtomsNum(true);
		pol.resize(numAtoms, 0);
		for(short i = 0; i < numAtoms; i++){
			short vid = getPolarizabilityValueID(getAtom(i).getAtomSymbol(), i);
			if(vid >= 0)
				pol[i] = Constants::polarizability_value[vid];
		}
	}
	return pol;
}

vector<short> & Molecule::getAtomicNum(bool withH){
	vector<short> & an = withH ? atomic_num_H : atomic_num;
	if(an.empty()) {
		short numAtoms = getAtomsNum(withH);
		an.resize(numAtoms);
		for(short i = 0; i < numAtoms; i++) { // get atomic number of the atom
			an[i] = getAtom(i).getAtomSymbol() + 1;	/* atomic number of the atom-i */
		}
	}
	return an;
}

void Molecule::setChiralFlag(short nChiralFlag){
	m_nChiralFlag = nChiralFlag;
}

void Molecule::setStextEntriesNum(short nStextEntries){
	m_nStextEntries = nStextEntries;
}

DistancePathMatrix & Molecule::getDistancePathMatrix(){
	if(distancePathMatrix.empty()) calculateDistanceMatrix(false);
	return distancePathMatrix;
}

DistanceMatrix<short> & Molecule::getDistanceMatrix(bool withH){
	DistanceMatrix<short> & DM = withH ? distanceMatrix_H : distanceMatrix;
	if(DM.empty()) calculateDistanceMatrix(withH);
	return DM;
}

vector< vector<short> >& Molecule::getConnectedAtomIDs(bool withH)
{
	vector< vector<short> > & connAtoms = withH ? connectedAtomIDs_H : connectedAtomIDs;
	if(connAtoms.empty()) initializeConnectedAtomIDs(withH);
	return connAtoms;
}

vector< vector<char> >& Molecule::getConnectedBondTypes(bool withH){
	vector< vector<char> > & connBonds = withH ? connectedBondTypes_H : connectedBondTypes;
	if(connBonds.empty()) initializeConnectedBondTypes(withH);
	return connBonds;
}

DistanceMatrix<int> & Molecule::getDistancePath(void)
{
	if(distancePath.empty()) calculateDistanceMatrix(false);
	return distancePath;
}

void Molecule::pushHeader(string & line)
{
  m_Header.append(line);
}

void Molecule::pushAtom(Atom& atom)
{
	m_atom_list.push_back(atom);
}

void Molecule::pushBond(Bond& bond)
{
	m_bond_list.push_back(bond);
}

short Molecule::getAtomsNum(bool withH)
{
	if(withH) return m_atom_list.size();
	else return m_numAtoms;
}

short Molecule::getBondsNum(bool withH){
	if(withH) return m_bond_list.size();
	else return m_numBonds;
}

Atom & Molecule::getAtom(short index) // 0-based
{
	return m_atom_list[index];
}

Bond & Molecule::getBond(short index){
	return m_bond_list[index];
}

void Molecule::preprocess(void){
	m_deletePart = 0;
	m_fingerprint.assign(TOTALFEATURES, Constants::ERROR_SIGNAL);
	m_FeatureStatus.assign(TOTALGROUPS, false);
	vector< set<short> > jointSetArray;
	set<short>::iterator it1, it2;
	for(short i = this->getAtomsNum(false) - 1; i >= 0; i--){ // each atom is a independent set
		set<short> newSet;
		newSet.insert(i);
		jointSetArray.push_back(newSet);
	}

	for(short i = getBondsNum(false) - 1; i >= 0; i--){
		short a1 = getBond(i).getFirstAtomNumber();
		short a2 = getBond(i).getSecondAtomNumber();

		short idx1 = -1;
		short idx2 = -1;
		for(short j = jointSetArray.size() - 1; j >= 0 ; j--){
			it1 = jointSetArray.at(j).find(a1);
			if(it1 != jointSetArray.at(j).end()){
				idx1 = j;
				break;
			}
		}

		for(short j = jointSetArray.size() - 1; j >= 0 ; j--){
			it2 = jointSetArray.at(j).find(a2);
			if(it2 != jointSetArray.at(j).end()){
				idx2 = j;
				break;
			}
		}

		if(idx1 == -1 && idx2 == -1){ // not exist, add a new set
			set<short> newSet;
			newSet.insert(a1);
			newSet.insert(a2);
			jointSetArray.push_back(newSet);
		}else if(idx1 > -1 && idx2 > -1){ // both exist, combine two sets
			if(idx1 != idx2){
				jointSetArray.at(idx1).insert(jointSetArray.at(idx2).begin(), jointSetArray.at(idx2).end());
				jointSetArray.erase(jointSetArray.begin() + idx2);
			}
		}else{
			if(idx1 > -1){
				jointSetArray.at(idx1).insert(a2);
			}else{
				jointSetArray.at(idx2).insert(a1);
			}
		}
	}

	if(jointSetArray.size() > 1 || getAtomsNum(false) != jointSetArray.size()){ // has at least two parts
		// find the largest part (H is not counted)
		vector<short> size(jointSetArray.size(), 0);
		for(short i = jointSetArray.size() - 1; i >= 0 ; i--){
			set<short> s1 = jointSetArray.at(i);
			for(set<short>::iterator it = s1.begin(); it != s1.end(); it++){
				if(getAtom(*it).getAtomSymbol() != Elements::HYDROGEN){
					size[i]++;
				}
			}
		}

		short maxN = size[0];
		int maxSetIdx = 0;
		for(short i = size.size() - 1; i > 0 ; i--){
			if(maxN < size[i]){
				maxN = size[i];
				maxSetIdx = i;
			}
		}
		set<short> keepIDs = jointSetArray.at(maxSetIdx);
		set<short> remove_atomIDs;
		for(short i = getAtomsNum(false) - 1; i >= 0; i--){
			if(keepIDs.find(i) == keepIDs.end()){ // remove this atom
				remove_atomIDs.insert(i);
			}
		}
		if(remove_atomIDs.size() > 0){
			removeAtomsAndBonds(remove_atomIDs);
			m_deletePart = 1;
		}
	}
	removeH();
}

void Molecule::removeH(){
	set<short> remove_atoms;
	for(short i = getAtomsNum(false) - 1; i >= 0; i--){
		int x = m_atom_list[i].getCharge();
		if(m_atom_list[i].getAtomSymbol() == Elements::HYDROGEN){ // remove an H atom with one connection only
			short cN = 0;
			for(short j = m_bond_list.size() - 1; j >= 0; j--){
				if(m_bond_list[j].getFirstAtomNumber() == i || m_bond_list[j].getSecondAtomNumber() == i){
					cN++;
					if(cN >= 2) break;
				}
			}
			if(cN == 1) remove_atoms.insert(i);
		}
	}
	if(getAtomsNum(false) != 2 || remove_atoms.size() != 2) // H-H case
		removeAtomsAndBonds(remove_atoms);
	getRings();	// ring and aromatic check before do any analysis
	addHydrogen(); // add Hydrogen
}

void Molecule::removeAtom(short index){
	m_atom_list.erase(m_atom_list.begin() + index);
	m_numAtoms--;
}

void Molecule::removeBond(short index){
	m_bond_list.erase(m_bond_list.begin() + index);
	m_numBonds--;
}

void Molecule::removeAtomsAndBonds(set<short> & atomIDs){
	set<short>::iterator it1, it2;
	short nAtoms = getAtomsNum(false);
	vector<short> atomIDMap(nAtoms);
	for(short i = nAtoms - 1; i >= 0; i--){
		atomIDMap[i] = i;
	}

	// re-assign atom id
	stack<short> removeIDs;
	for(short i = 0; i < nAtoms; i++){
		if(atomIDs.find(i) != atomIDs.end()){ // remove this atom
			for(short j = i + 1; j < nAtoms; j++){
				atomIDMap[j]--;
			}
			removeIDs.push(i);
		}
	}

	//remove atoms
	while (!removeIDs.empty()){
		removeAtom(removeIDs.top());
		removeIDs.pop();
	}

	// remove bonds
	for(short i = 0; i < getBondsNum(false); i++){
		short first = m_bond_list[i].getFirstAtomNumber();
		short second = m_bond_list[i].getSecondAtomNumber();
		it1 = atomIDs.find(first);
		it2 = atomIDs.find(second);
		if( it1 != atomIDs.end() || it2 != atomIDs.end()){
			removeBond(i);
			i--;
		}else{
			m_bond_list[i].setFirstAtomNumber(atomIDMap[first]);
			m_bond_list[i].setSecondAtomNumber(atomIDMap[second]);
		}
	}
}

// If the molecule have Hydrogen and does not display, add Hydrogen atoms
void Molecule::addHydrogen(){
	for(short i = getAtomsNum(false) - 1; i >= 0; i--){
		Atom & atom = getAtom(i);
		if(atom.getHydrogenCountPlus1() != 0) continue;

		if(atom.getValence() == 15) atom.setHydrogenCountPlus1(1);

		double sum = 0.0;
		for(short j = getBondsNum(false) - 1; j >= 0 ; j--){
			Bond & bond = getBond(j);
			if(bond.getFirstAtomNumber() == i || bond.getSecondAtomNumber() == i){
				if(bond.getBondType() >= Bond::AROMATIC_BOND){ // For aromatic bonds.
					sum += 1.5;
				}else{
					sum += bond.getBondType();
				}
			}
		}

		if(atom.getValence() == 0) atom.setValence((char) (sum + 0.1));

		/* One positive or negative charge use 1 valence. */
		sum += fabs(double(Constants::defined_charge[atom.getCharge()]));

		/* Hydrogen count + 1. */
		int j = atom.getAtomSymbol();
		short k = 0;

		/* N, P, An, Sb, Bi have a positive charge make their total valences to be 5. */
		if(Constants::defined_charge[atom.getCharge()] == 1 && (j == 6 || j == 14 || j == 32 || j == 50 || j == 82)) k = 1;

		while(k < Constants::atom_valences[j][0]){
			k++;
			atom.setHydrogenCountPlus1(Constants::atom_valences[j][k] - (char) (sum + 0.1) + 1);
			if(atom.getHydrogenCountPlus1() >= 1) break;
		}

		/* Some metals can give out bonds more than their valences. */
		if(atom.getHydrogenCountPlus1() < 1) atom.setHydrogenCountPlus1(1);

		/* Single atom molecules. Only add Hydrogen atoms to B, C, N, O, F, Si, P, S, Cl, Br, I*/
		if(getAtomsNum(false) == 1 && atom.getHydrogenCountPlus1() > 1){
			if(j < 4 || (j > 8 && j < 13) || j > 16 || j == 34 || j == 52) atom.setHydrogenCountPlus1(1);
		}
	}

	/* Add Hydrogen atoms to Mol file. */
	for(short i = getAtomsNum(false) - 1; i >= 0; i--){
		if(getAtom(i).getValence() != 15){
			char nH = getAtom(i).getHydrogenCountPlus1() - 1;
			for(short j = 0; j < nH; j++){
				Atom atom;
				atom.setHydrogenCountPlus1(1);
				atom.setValence(1);
				pushAtom(atom);

				Bond bond;
				bond.setFirstAtomNumber(i);
				bond.setSecondAtomNumber(getAtomsNum(true) - 1);
				pushBond(bond);
			}
		}
	}

	/* Recheck the number of H atoms for each atom, in case there are pre-added H atoms (such as, explicit Hs). */
	checkHydrogenNums();

	//Call following functions after adding H
	countElements();
	getMoleculeWeight();
}

void Molecule::print(ostream & out, bool hasH, bool isSDF)
{
	//output head lines
	out << m_Header<< endl;
	short atomNum = getAtomsNum(hasH);
	for(short j = 0; j < atomNum; j++){
		Atom & atom = getAtom(j);
		if(hasH){
			atom.setHydrogenCountPlus1(0);
			atom.setValence(0);
		}
		out.precision(4);
		out << fixed;
		out << setw(10) << atom.getX();
		out << setw(10) << atom.getY();
		out << setw(10) << atom.getZ();
		out << " ";
		out << setw(3) << Elements::idToSymbol(atom.getAtomSymbol());
		out << setw(2) << (short)atom.getMassDifference();
		out << setw(3) << (short)atom.getCharge();
		out << setw(3) << (short)atom.getAtomStereoParity();
		out << setw(3) << (short)atom.getHydrogenCountPlus1();
		out << setw(3) << (short)atom.getStereoCareBox();
		out << setw(3) << (short)atom.getValence();
		out << setw(3) << (short)atom.getH0_designator();
		out << setw(6) << ""; // two n
		out << setw(3) << (short)atom.getAtomMappingNumber();
		out << setw(3) << (short)atom.getInversionRetentionFlag();
		out << setw(3) << (short)atom.getExactChangeFlag() << endl;
	}

	short bondNum = getBondsNum(hasH);
	for(short j = 0; j < bondNum; j++){
		Bond & bond = getBond(j);
		out << setw(3) << bond.getFirstAtomNumber() + 1;
		out << setw(3) << bond.getSecondAtomNumber() + 1;
		out << setw(3) << (short)bond.getBondType();
		out << setw(3) << (short)bond.getBondStereo();
		out << setw(3) << "";
		out << setw(3) << (short)bond.getBondTopology();
		out << setw(3) << (short)bond.getReactingCenterStatus() << endl;
	}
	out << "M  END" << endl;
	if(isSDF)
		out << "$$$$" << endl;
}

void Molecule::checkHydrogenNums(void)
{
	for(short i = getAtomsNum(true) - 1; i >= 0; i--) getAtom(i).setHydrogenCountPlus1(1);

	for(short i = getBondsNum(true) - 1; i >= 0; i--){
		Bond &bond = getBond(i);
		if(getAtom(bond.getFirstAtomNumber()).getAtomSymbol() == Elements::HYDROGEN){
			Atom & atom = getAtom(bond.getSecondAtomNumber());
			atom.setHydrogenCountPlus1(atom.getHydrogenCountPlus1() + 1);
		}else if(getAtom(bond.getSecondAtomNumber()).getAtomSymbol() == Elements::HYDROGEN){
			Atom & atom = getAtom(bond.getFirstAtomNumber());
			atom.setHydrogenCountPlus1(atom.getHydrogenCountPlus1() + 1);
		}
	}
}

vector<short>& Molecule::countElements(void){
	if(m_elementCounts.empty()){
		m_elementCounts.resize(Elements::TOTALELEMENTS, 0);

		for(short i = getAtomsNum(true) - 1; i >= 0; i--){
			short atomID = getAtom(i).getAtomSymbol();
			m_elementCounts[atomID] += 1;
		}
		setFingerprint(m_elementCounts, CountElements);
	}
	return m_elementCounts;
}

vector<double> Molecule::getMoleculeWeight(void){
	vector<double> molecular_weight(2, 0);
	short numAtoms = getAtomsNum(true);
	for(short i = 0; i < numAtoms; i++){
		molecular_weight[0] += Constants::atom_weight[getAtom(i).getAtomSymbol()];
	}

	/* computing average_molecular_weight */
	molecular_weight[1] = molecular_weight[0] / numAtoms;
	setFingerprint(molecular_weight, MolecularWeight);
	return molecular_weight;
}

void Molecule::setHeaderMoleculeId(string hmolId)
{
  m_headerMoleculeId = hmolId;
}

string Molecule::getHeaderMoleculeId() { return m_headerMoleculeId; }


inline void Molecule::ConnectivityOrder2Index(PathNode * p_path, double & vci, double & vvci, double & pqvci, int & a){
	double mms, mmv, mm;
	mms = mmv = mm = 1.0;
	vector<float> & valence = getValence();
	vector<short> & pqn = getPqn();
	vector< vector<short> > & connAtoms = getConnectedAtomIDs(false);
	PathNode * p = p_path;
	while(p != NULL){
		mm *= connAtoms[p->id].size(); // product of all path atoms
		mmv *= valence[p->id];	// product of valence vertex degree for valence connectivity index
		mms *= pqn[p->id];		// product of principal quantum numbers
		p = p->parent;
	}
	vci += 1.0 /sqrt(mm);		// connectivity order-3 Index
	vvci += 1.0 / sqrt(mmv);		// valence vertex connectivity order-3 Index
	pqvci += 1.0 /sqrt(mm) / mms;	// principal quantum vertex connectivity order-3 index
	a++;
}

DistanceMatrix<int> & Molecule::getAllPath(void){
	if(allPath.empty()) calculateDistanceMatrix(false);
	return allPath;
}

DistanceMatrix<short> & Molecule::getDetourMatrix(void){
	if(detourMatrix.empty()) calculateDistanceMatrix(false);
	return detourMatrix;
}
vector< vector<int> > & Molecule::getPath(void){
	if(path.empty()) calculateDistanceMatrix(false);
	return path;
}

vector<int> & Molecule::getP(){
	if(p.empty()) calculateDistanceMatrix(false);
	return p;
}
vector<int> & Molecule::getA(){
	if(a.empty()) calculateDistanceMatrix(false);
	return a;
}
vector<double> & Molecule::getVC() {
	if(vc.empty()) calculateDistanceMatrix(false);
	return vc;
}
vector<double> & Molecule::getVVC() {
	if(vvc.empty()) calculateDistanceMatrix(false);
	return vvc;
}
vector<double> & Molecule::getPQVC() {
	if(pqvc.empty()) calculateDistanceMatrix(false);
	return pqvc;
}

inline short findPathLength(PathNode * path){
	PathNode * pt = path->parent;
	short len = 0;
	while(pt != NULL){
		len++;
		pt = pt->parent;
	}
	return len;
}

inline void Molecule::savePaths(vector< PathNode * > &pathList, vector<short> & len, FILE * pFile){
	int nPath = pathList.size();
	fwrite(&nPath, sizeof(int), 1, pFile);
	for(int i = pathList.size() - 1; i >= 0; i--){
		fwrite(&(len[i]), sizeof(short), 1, pFile);
		PathNode * pt = pathList[i];
		while(pt != NULL){
			fwrite(&(pt->id), sizeof(short), 1, pFile);
			pt = pt->parent;
		}
	}
	if(ftell(pFile) > ::numeric_limits<long>::max())
		throw mold2Exception("cache file oversize");
}

void Molecule::calculatDistanceMatrixOnly(){
	short numAtoms = getAtomsNum(false);
	vector<short> shortestPath;
	for(short i = 0; i < numAtoms; i++){
		for(short j = i+1; j < numAtoms; j++){
			findShortestPath(false, i, j, shortestPath);
			short len = shortestPath.size() - 1;
			distanceMatrix.setValueToCell(i, j, len);
			int val = (len * len +len) >> 1;
			distancePath.setValueToCell(i, j, val);
			distancePathMatrix.setPath(i, j, shortestPath);
		}
	}
}

void Molecule::calculateDistanceMatrix(bool withH){
	short numAtoms = getAtomsNum(withH);
	if(withH){
		distanceMatrix_H.resize(numAtoms);
		vector<short> shortestPath;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i+1; j < numAtoms; j++){
				distanceMatrix_H.setValueToCell( i, j, findShortestPath(true,i, j, shortestPath));
			}
		}
	}else{
		distanceMatrix.resize(numAtoms);
		distancePathMatrix.resize(numAtoms);
		distancePath.resize(numAtoms);
		allPath.resize(numAtoms);
		detourMatrix.resize(numAtoms);

		path.resize(4, vector<int>(numAtoms, 0));
		p.resize(3, 0);

		a.resize(6, 0);
		vc.resize(6, 0);
		vvc.resize(6, 0);
		pqvc.resize(6, 0);

		vector< PathNode * > pathList;
		pathSets.clear();
		int nTimes = (numAtoms * numAtoms - numAtoms) / 2;
		pathLen.resize(nTimes);

		int minLen, pathID, maxLen, sumLen, nPaths, m, val;
		int totalPaths = 0;
		short i, j;
		int nPathCount = 0;
		vector<short> currentPathLen;
		int timesCount = 0;


		for(i = 0; i < numAtoms; i++){
			for(j = i+1; j < numAtoms; j++){
				if(numAtoms > 50){
					vector<short> shortestPath;
					findShortestPath(false, i, j, shortestPath); // find the shortest path only
					pathList.clear();
					PathNode * p_path = new PathNode(i, NULL);
					for(int xm = 1; xm < shortestPath.size(); xm++){
						p_path = new PathNode(shortestPath[xm], p_path);
					}
					pathList.push_back(p_path);

					minLen = shortestPath.size() - 1;

					distanceMatrix.setValueToCell(i, j, minLen);
					val = (minLen * minLen + minLen) >> 1;
					distancePath.setValueToCell(i, j, val);
					distancePathMatrix.setPath(i, j, shortestPath);

					// detour matrix
					detourMatrix.setValueToCell(i, j, minLen);
					// all path matrix
					allPath.cellValuePlus(i, j, minLen);

					vector<short> & len = m_useSearchDoAllPath ? currentPathLen:pathLen[nPathCount++];
					len.resize(1);
					len[0] = minLen;

					// initializePathWalk
					if(minLen > 1 && minLen < 6){
						path[minLen - 2][i]++;	// number of paths
						path[minLen - 2][j]++;	// number of retrograde the paths
					}

					// initializeKierPathIndex
					switch(minLen){
					case 1: // path equal to one between any atom
						p[0]++;
						break;
					case 2: // path equal to two between any atom
						p[1]++;
						break;
					case 3: // path equal to three between any atom
						p[2]++;
					}

					// initialzieAllPathConnectivity
					if(minLen > 1 && minLen < 6)
						ConnectivityOrder2Index(pathList[0], vc[minLen], vvc[minLen], pqvc[minLen], a[minLen]);

					pathSets.push_back(pathList);
					continue;
				}else{
					dfSearch(i, j, pathList); // find all paths between i and j;
				}
				timesCount++;
				minLen = numeric_limits<int>::max();
				pathID = -1;
				maxLen = sumLen = 0;
				nPaths = pathList.size() - 1;
				vector<short> & len = m_useSearchDoAllPath ? currentPathLen:pathLen[nPathCount++];
				len.resize(nPaths + 1);
				for(m = nPaths; m >= 0; m--){
					len[m] = findPathLength(pathList[m]);
					if(minLen > len[m]){
						minLen = len[m];
						pathID = m; // distance path
					}

					if(maxLen < len[m]) maxLen = len[m];
					sumLen += len[m];
				}

				totalPaths += sumLen;
				distanceMatrix.setValueToCell(i, j, minLen);
				val = (minLen * minLen + minLen) >> 1;
				distancePath.setValueToCell(i, j, val);
				PathNode * p_path = pathList[pathID];
				vector<short> shortest(len[pathID] + 1);
				while(p_path != NULL){
					shortest[minLen--] = p_path->id;
					p_path = p_path->parent;
				}
				distancePathMatrix.setPath(i, j, shortest);

				// detour matrix
				detourMatrix.setValueToCell(i, j, maxLen);
				// all path matrix
				allPath.cellValuePlus(i, j, sumLen);

				// initializePathWalk
				for(m = nPaths; m >= 0; m--){ // part is the same as KierPathIndex
					if(len[m] > 1 && len[m] < 6){
						path[len[m] - 2][i]++;	// number of paths
						path[len[m] - 2][j]++;	// number of retrograde the paths
					}
				}

				// initializeKierPathIndex
				for(m = nPaths; m >= 0; m--){
					switch(len[m]){
					case 1: // path equal to one between any atom
						p[0]++;
						break;
					case 2: // path equal to two between any atom
						p[1]++;
						break;
					case 3: // path equal to three between any atom
						p[2]++;
					}
				}

				// initialzieAllPathConnectivity
				for(m = nPaths; m >= 0; m--){
					// calculate the connectivity order-2 Index
					if(len[m] > 1 && len[m] < 6)
						ConnectivityOrder2Index(pathList[m], vc[len[m]], vvc[len[m]], pqvc[len[m]], a[len[m]]);
				}

				pathSets.push_back(pathList);
			}
		}

		setFingerprint(totalPaths, ::AllPathWiener);
		MolecularPathCount(NULL);
	}
}


void Molecule::initializeConnectedAtomIDs(bool withH){
	vector< vector<short> > & connAtoms = withH ? connectedAtomIDs_H : connectedAtomIDs;
	short numAtoms = getAtomsNum(withH);
	connAtoms.resize(numAtoms);
	for(short i = getBondsNum(withH) - 1; i >= 0; i--){
		Bond & bond = getBond(i);
		short atom1 = bond.getFirstAtomNumber();
		short atom2 = bond.getSecondAtomNumber();
		connAtoms[atom1].push_back(atom2);
		connAtoms[atom2].push_back(atom1);
	}
}

void Molecule::initializeConnectedBondTypes(bool withH){
	vector< vector<char> > & connBonds = withH ? connectedBondTypes_H : connectedBondTypes;

	short numAtoms = getAtomsNum(withH);
	connBonds.resize(numAtoms);

	for(short i = getBondsNum(withH) - 1; i >= 0; i--){
		Bond & bond = getBond(i);
		short atom1 = bond.getFirstAtomNumber();
		short atom2 = bond.getSecondAtomNumber();
		char bondType = bond.getBondType();
		connBonds[atom1].push_back(bondType);
		connBonds[atom2].push_back(bondType);
	}
}


vector< vector< PathNode * > > & Molecule::getPathSets(void){
	if(pathSets.empty() && !m_useSearchDoAllPath && getAtomsNum(false) > 1) calculateDistanceMatrix(false);
	return pathSets;
}

void Molecule::calculateMoleculePathCountFromArrays(){
	size_t nPath = 0;
	vector< vector<PathNode *> > & pathSets = getPathSets();
	short numAtoms = getAtomsNum(false);
	double r, ri, bapi;
	vector< vector<char> > & connBond = getConnectedBondTypes(false);
	vector< vector<short> > & connAtoms = getConnectedAtomIDs(false);
	vector<int> & apath = allPath.getSigma();

	for(short i = 0; i < numAtoms; i++){
		for(short j = i+1; j < numAtoms; j++){
			vector<PathNode *> & pathList = pathSets[nPath];
			for(int k = pathList.size() - 1; k >= 0 ; k--){	// order 1~numAtoms
				short len = pathLen[nPath][k];
				// the path start and end atoms: order 2~numAtoms
				path_length[i][len]++;

				// ---------- calculate molecular multiple path count and Randic indexr: RI -----------
				r = ri = bapi = 1.0;
				PathNode * xp = pathList[k];
				PathNode * xp1 = xp->parent;
				while(xp1){
					short akm = xp->id;
					short akm1 = xp1->id;
					xp = xp1;
					xp1 = xp1->parent;
					short nakm = connAtoms[akm].size();
					for(int n = nakm - 1; n >= 0; n--){
						// find the edge that between the pair atoms
						if(connAtoms[akm][n] == akm1){
							if(connBond[akm][n] == Bond::AROMATIC_BOND) // aromatic bond
								r *= 1.5; // aromatic bond product to 1.5
							else r *= connBond[akm][n]; // others bonds
						}
					}
					// ------------------------- Randic indexr: RI ---------------------------------------
					double txt = sqrt(double(nakm) * connAtoms[akm1].size());
					if(txt != 0) ri *= 1.0 / txt;

					// -------------- Balaban ID number calculated by all path: bapi -----------------------
					txt = sqrt(double(apath[akm]) * apath[akm1]);
					if(txt != 0) bapi *= 1.0 / txt;
				} // end of FOR-m
				mpath[len] += (float)r; // sum the product of all bonds' type on each path for multiple bond orders
				cweight[len] += (float)ri;	// sum the product of edge connectivity's weight of each path
				eweight[len] += (float)bapi;	// sum the product of edge weight
			}
			nPath++;
		}
	}
	clearPathsSet();
}

void Molecule::clearPathsSet(){
	for(int i = pathSets.size() - 1; i >= 0; i--){
		clearPathList(pathSets[i]);
	}
	pathSets.clear();
}

void Molecule::calculateMoleculePathCountFromFile(FILE * pFile){
	short len;
	int nPath;
	short numAtoms = getAtomsNum(false);
	short * pathX = new short[numAtoms];
	vector< vector<char> > & connBond = getConnectedBondTypes(false);
	vector< vector<short> > & connAtoms = getConnectedAtomIDs(false);
	vector<int> & apath = allPath.getSigma();
	double r, ri, bapi;
	fseek(pFile, 0, SEEK_SET);
	for(short i = 0; i < numAtoms; i++){
		for(short j = i+1; j < numAtoms; j++){
			fread(&nPath, sizeof(int), 1, pFile);
			for(int k = nPath - 1; k >= 0 ; k--){	// order 1~numAtoms
				fread(&len, sizeof(short), 1, pFile);
				fread(pathX, sizeof(short), len + 1, pFile);
				// the path start and end atoms: order 2~numAtoms
				path_length[i][len]++;

				// ---------- calculate molecular multiple path count and Randic indexr: RI -----------
				r = ri = bapi = 1.0;
				for(short m = len; m > 0; m--){
					short m1 = m - 1;
					short nakm = connAtoms[pathX[m]].size();
					for(int n = nakm - 1; n >= 0; n--){
						// find the edge that between the pair atoms
						if(connAtoms[pathX[m]][n] == pathX[m1]){
							if(connBond[pathX[m]][n] == Bond::AROMATIC_BOND) // aromatic bond
								r *= 1.5; // aromatic bond product to 1.5
							else r *= connBond[pathX[m]][n]; // others bonds
						}
					}
					// ------------------------- Randic indexr: RI ---------------------------------------
					double txt = sqrt(double(nakm) * connAtoms[pathX[m1]].size());
					if(txt != 0) ri *= 1.0 / txt;

					// -------------- Balaban ID number calculated by all path: bapi -----------------------
					txt = sqrt(double(apath[pathX[m]]) * apath[pathX[m1]]);
					if(txt != 0) bapi *= 1.0 / txt;
				} // end of FOR-m
				mpath[len] += (float)r; // sum the product of all bonds' type on each path for multiple bond orders
				cweight[len] += (float)ri;	// sum the product of edge connectivity's weight of each path
				eweight[len] += (float)bapi;	// sum the product of edge weight
			}
			nPath++;
		}
	}
	delete [] pathX;
	fclose(pFile);
}

void Molecule::MolecularPathCount(FILE * pFile){
	if(isFeaturesReady(::MolecularPathCount)) return;
	if(m_allPathFailed){
		vector<double> tmp(23, ::numeric_limits<float>::signaling_NaN());
		setFingerprint(tmp, ::MolecularPathCount);
		return;
	}

	// calculate the all-path: distance of each atom to anothers atoms
	short numAtoms = getAtomsNum(false);

	path_length.assign(numAtoms, vector<int>(numAtoms, 0));
	short mSize = numAtoms > 10 ? numAtoms:11;
	mpath.assign(mSize, 0);
	cweight.assign(numAtoms, 0);
	eweight.assign(numAtoms, 0);

	if(pFile) calculateMoleculePathCountFromFile(pFile);
	else calculateMoleculePathCountFromArrays();

	// ------------------ calculate path length of molecular, CBI, RI --------------------------
	double cbi = sum_element(mpath) + numAtoms;
	double ri = sum_element(cweight);
	double bapi = sum_element(eweight);
	short pSize = numAtoms > 11 ? numAtoms : 11;
	vector<int> px(pSize, 0);

	for(short i = 0; i < numAtoms; i++){
		for(short j = 2; j < numAtoms; j++)
			px[j] += path_length[i][j];
	}
	px[0] = numAtoms;
	px[1] = getBondsNum(false);

	ri += numAtoms;	// Randic indexr: RI
	bapi += numAtoms; // Balaban All-Path index: BAPI

	// ------------------ calculate molecular topological all path index: MTAPI --------------------------
	double mtapi = sum_element(px);

	// ---------------- ratio of convention bonds with total path counts: RCBTPC -----------------
	double rcbtpc = 0.0;
	if(mtapi != 0)
		rcbtpc = cbi / mtapi;

	// ---------- ratio of difference of conventional bonds and total path counts: RDCTPC -------------
	double rdctpc = (cbi - mtapi) / 88.6880604752573;

	// ================================ Printout results =================================

	vector<double> fingerprint;
	// printout results of molecular path length of atom-i in order 2~10 to output file
	for(short i = 2; i <= 10; i++)
		fingerprint.push_back(px[i]);

	//printout results of molecular topological  multiple path index of order 03~10 to output file
	for(short i = 3; i <= 10; i++)
		fingerprint.push_back(mpath[i]);

	// printout results of MTAPI, CBI, RCBTBC, RDCTPC, RI, BAPI to output file
	fingerprint.push_back(mtapi);
	fingerprint.push_back(cbi);
	fingerprint.push_back(rcbtpc);
	fingerprint.push_back(rdctpc);
	fingerprint.push_back(ri);
	fingerprint.push_back(bapi);
	setFingerprint(fingerprint, ::MolecularPathCount);
}


// Breath First Search
int Molecule::findShortestPath(bool withH, short from, short to, vector<short> &shortestPath){
	shortestPath.clear();
	vector< vector<short> > & connAtoms = getConnectedAtomIDs(withH);
	short numAtoms = getAtomsNum(withH);
	short current = -1;
	deque< short > myQueue;
	myQueue.push_back(from);
	vector<char> mark(numAtoms, 0);
	vector<short> parent(numAtoms, -1);
	parent[from] = from;

	while(!myQueue.empty()){
		current = myQueue.front();
		mark[current] = GRAY;
		myQueue.pop_front();

		// check end node in the path
		for(short j = connAtoms[current].size() - 1; j >= 0; j--){
			if( mark[connAtoms[current][j]] == WHITE){
				if(connAtoms[current][j] == to) {// find the shorest path
					shortestPath.push_back(to);
					while(current != from){
						shortestPath.push_back(current);
						current = parent[current];
					}
					shortestPath.push_back(from);
					return shortestPath.size() - 1;
				}else{
					myQueue.push_back(connAtoms[current][j]);
					parent[connAtoms[current][j]] = current;
					mark[connAtoms[current][j]] = GRAY;
				}
			}
		}
	}
	return -1;
}

inline void deleteBranch(PathNode * endP){
	while(endP != NULL && endP->children == 0){
		PathNode * tp = endP;
		endP = endP->parent;
		delete tp;
		if(endP != NULL) endP->children--;
	}
}

inline void Molecule::clearPathList(vector< PathNode * > & pathList){
	for(vector< PathNode *>::iterator it = pathList.begin(); it != pathList.end(); it++){
		deleteBranch(*it);
	}
	pathList.clear();
}

inline bool inPath(short id, PathNode *p_path){
	PathNode *p = p_path->parent;
	while(p != NULL){
		if(p->id == id) return true;
		else p = p->parent;
	}
	return false;
}

void Molecule::dfSearch(short from, short to, vector< PathNode * > & paths)
{
	if(m_useSearchDoAllPath) clearPathList(paths);
	else paths.clear();

	stack< PathNode * > myQueue;
	myQueue.push(new PathNode(from, NULL));

	vector< vector<short> > & connAtoms = getConnectedAtomIDs(false);
	PathNode *p_path;
	while(!myQueue.empty()){
		p_path = myQueue.top();
		myQueue.pop();
		if(paths.size() < 10){
			// check end node in the path
			for(short j = connAtoms[p_path->id].size() - 1; j >= 0; j--){
				if(connAtoms[p_path->id][j] == to) {// find a new path
					paths.push_back(new PathNode(to, p_path));
					continue;
				}

				if(!inPath(connAtoms[p_path->id][j], p_path)){ // not in the path, add a new path
					myQueue.push(new PathNode(connAtoms[p_path->id][j], p_path));
				}
			}
		}
		deleteBranch(p_path);
	}
}