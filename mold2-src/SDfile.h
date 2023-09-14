#pragma once

class Molecule;
class SDfile
{
public:
	SDfile(void);
	~SDfile(void);
	// Check whether reach the end of SDfile
	bool hasNext(void);
	// get next molfile
	void getNextMolecule(Molecule & mol);
	// open a SDfile
	bool openSDfile(istream & in);
	// close file
	void closeSDfile(void);
	void getMoleculeFrom(istream & in, vector<Molecule> & mol);

private:
	int m_MoleculeCount;
	string m_header;
  string m_headerMoleculeId;
	// SDfile stream
	istream* m_in;

	const static int m_CountsLineLength;
	const static int m_AtomBlockLength;
	const static int m_BondBlockLength;
	void jumpToMolfileEnd(istream &in, string & fieldName, int &ID);
	void readMol(istream &in, Molecule & mol);
	void readHeader(istream & in);
};
