#include "stdafx.h"

#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "SDfile.h"

#include "Elements.h"

// The constants listed below are defined in Symyx SDfile v2000 format
//  CTfile formats Symyx June 2010. [sch - min bond line length shortened to 12 to support sdf output from RDKit]
const int SDfile::m_CountsLineLength = 39;
const int SDfile::m_AtomBlockLength = 69;
const int SDfile::m_BondBlockLength = 12;

SDfile::SDfile(void)
{
	m_MoleculeCount = 0;
	m_in = NULL;
}

SDfile::~SDfile(void)
{
}


// Check whether reach the end of SDfile
bool SDfile::hasNext()
{
	if(m_in -> eof()){
		return false;
	}else{
		readHeader(*m_in);
		return m_in -> good();
	}
}

void cleanHeaderMoleculeName(string& s)
{
  s.erase(std::remove_if(s.begin(), s.end(), [](char c) { return c == '\r' || c == '\t'; } ), s.end());
}

void SDfile::readHeader(istream & in){
	string line;
	m_header.clear();
	getline(in, line);
  cleanHeaderMoleculeName(line);
  m_headerMoleculeId = line;
	m_header.append(line + "\n");
	getline(in, line);
	m_header.append(line + "\n");
	getline(in, line);
	m_header.append(line + "\n");
}

// Check whether reach the end of SDfile
void SDfile::getMoleculeFrom(istream & in, vector<Molecule> & mols){
	mols.clear();
	string line = "$$$$";
	while(line == "$$$$"){
		Molecule mol;
		getline(in, line);
    cout << "Got line: " << line << endl; // TODO: Remove
		trim(line);
		mol.m_ID = atoi(line.c_str());
    cout << "made id" << endl; // TODO: Remove

		readHeader(in);
    cout << "read header" << endl; // TODO: Remove
		readMol(in, mol);
    cout << "read mol" << endl; // TODO: Remove
		getline(in, line);
		getline(in, line);

		//preprecess molecule to remove smaller disjointed parts
		mols.push_back(mol);
	}
}

void SDfile::readMol(istream &in, Molecule & mol){
	mol.pushHeader(m_header);
  mol.setHeaderMoleculeId(m_headerMoleculeId);
	string line;
	getline(in, line);
	if(line.size() < m_CountsLineLength){
		jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
		throw mold2Exception("The counts line of the Mol/SD-file is invalid.", mol.getHeaderMoleculeId());
	}
	if(line.substr(36) == "v3000"){
		jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
		throw mold2Exception("Unsupported version v3000 encountered in Mol/SD-file.", mol.getHeaderMoleculeId());
	}

	int nAtoms = atoi( line.substr(0, 3).c_str() );
	int nBonds = atoi( line.substr(3, 3).c_str() );
	if(nAtoms < 1){
		jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
		throw mold2Exception("No identifiable structure found in Mol/SD-file.", mol.getHeaderMoleculeId());
	}
	mol.setNumAtoms(nAtoms);
	mol.setNumBonds(nBonds);
	mol.pushHeader(line);
	mol.setChiralFlag(atoi(line.substr(12, 3).c_str()));
	mol.setStextEntriesNum(atoi(line.substr(15, 3).c_str()));

	// read atoms block
	for(int i = 0; i < nAtoms; i++){
		getline(in, line);
		if(!in.good() || line.size() < m_AtomBlockLength){
			jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
			if(in.good()){
				throw mold2Exception("Atom block is invalid in Mol/SD-file.", mol.getHeaderMoleculeId());
			}
			throw mold2Exception("Could not read atoms block in Mol/SD-file.", mol.getHeaderMoleculeId());
		}
		Atom atom;
		atom.setX((float)atof(line.substr(0, 10).c_str()));
		atom.setY((float)atof(line.substr(10, 10).c_str()));
		atom.setZ((float)atof(line.substr(20, 10).c_str()));

		string symbol = line.substr(31, 3);
		trim(symbol);
		char atomSymbol = Elements::symbolToId(symbol);
		if(atomSymbol == Elements::NOTFOUND){ // Check unknown elements
			jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
			throw mold2Exception("Unknown element symbol encountered in Mol/SD-file.", mol.getHeaderMoleculeId());
		}
		atom.setAtomSymbol(atomSymbol);
		atom.setMassDifference(atoi(line.substr(34,2).c_str()));
		atom.setCharge(atoi(line.substr(36, 3).c_str()));
		atom.setStereoParity(atoi(line.substr(39, 3).c_str()));
		atom.setHydrogenCountPlus1(atoi(line.substr(42, 3).c_str()));
		atom.setStereoCareBox(atoi(line.substr(45, 3).c_str()));
		atom.setValence(atoi(line.substr(48, 3).c_str()));
		atom.setH0designator(atoi(line.substr(51, 3).c_str()));
		atom.setAtomMappingNumber(atoi(line.substr(60, 3).c_str()));
		atom.setInversionRetentionFlag(atoi(line.substr(63, 3).c_str()));
		atom.setExactChangeFlag(atoi(line.substr(66, 3).c_str()));
		mol.pushAtom(atom);
	}

	// read bonds block
	for(int i = 0; i < nBonds; i++){
		getline(in, line);
		if(!in.good() || line.size() < m_BondBlockLength){
			jumpToMolfileEnd(in, mol.m_FieldName, mol.m_ID);
			if(in.good()){
				throw mold2Exception("Invalid bond block in Mol/SD-file.", mol.getHeaderMoleculeId());
			}
			throw mold2Exception("Could not read bond block in Mol/SD-file.", mol.getHeaderMoleculeId());
		}
		Bond bond;
		bond.setFirstAtomNumber(atoi(line.substr(0, 3).c_str()) - 1);
		bond.setSecondAtomNumber(atoi(line.substr(3, 3).c_str()) - 1);
		int bt = atoi(line.substr(6, 3).c_str());
		if(bt >= mol.m_maxBondType) mol.m_maxBondType = bt + 1;
		bond.setBondType(bt);
		bond.setBondStereo(atoi(line.substr(9, 3).c_str()));
    // Let remaining fields default to 0 if not present (RDKit output support).
    size_t len = line.size();
    bond.setBondTopology(len > 15 ? atoi(line.substr(15, 3).c_str()) : 0);
		bond.setReactingCenterStatus(len > 18 ? atoi(line.substr(18,3).c_str()) : 0);
		mol.pushBond(bond);
	}
}

// get next molfile
void SDfile::getNextMolecule(Molecule & mol){
	mol.clear();
	readMol(*m_in, mol);
	jumpToMolfileEnd(*m_in, mol.m_FieldName, mol.m_ID);
	//preprecess molecule to remove smaller disjointed parts
	mol.preprocess();
}

bool SDfile::openSDfile(istream & is)
{
	m_in = &is;
	return true;
}

// close file
void SDfile::closeSDfile(void)
{
}

void SDfile::jumpToMolfileEnd(istream &in, string & fieldName, int &ID)
{
	string line;
	m_MoleculeCount++;
	while(in.good()){
		getline(in, line);
		if(fieldName == "" && line.size() > 0){
			bool findID = false;
			if(line.substr(0, 7) == ">  <ID>"){
				fieldName = "ID";
				findID = true;
			}else if(line.substr(0, 12) == ">  <AUTO_ID>"){
				fieldName = "AUTO_ID";
				findID = true;
			}

			// read ID
			if(findID){
				getline(in, line);
				trim(line);
				ID = atoi(line.c_str());
			}
		}
		if(line.substr(0, 4) == "$$$$")
			break;
	}
}

