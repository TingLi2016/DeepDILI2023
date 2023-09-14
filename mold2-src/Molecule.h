#pragma once

class PathNode{
public:
	PathNode(short _id, PathNode * _p){
		id = _id;
		parent = _p;
		if(_p != NULL) {
			parent->children++;
		}
		children = 0;
	}
	short id; // atom ID
	char children; // children's number
	PathNode * parent; //point to parent
};


class Molecule
{
private:
	vector< vector<int> > path;
	vector<int> p;
	vector<int> a;
	vector<double> vc;
	vector<double> vvc;
	vector<double> pqvc;
	vector< vector<PathNode *> > pathSets;
	vector< vector<short> > pathLen;
	DistanceMatrix<int> allPath;
	DistanceMatrix<short> detourMatrix;
	bool m_allPathFailed;

	vector< vector<int> > path_length;
	vector<float> mpath;
	vector<float> cweight;
	vector<float> eweight;

	short m_nChiralFlag;
	short m_nStextEntries;
	bool m_ringIsReady;
	bool m_useSearchDoAllPath;
	vector< vector<short> > rings;
	string m_Header;
  string m_headerMoleculeId;
	short m_numAtoms; // without H
	short m_numBonds; // without H

	vector<float> m_fingerprint;
	vector<bool> m_FeatureStatus;

	vector<short> atomic_num;
	vector< vector<short> > atom2bonds;
	vector<float> pol;
	vector<char> inRing; // used to check if atom i is on ring
	vector<float> valence;
	vector<short> pqn;
	vector<float> conventionalBondOrder;
	vector<char> con_to_H;
	vector<char> aromatic;
	vector< vector<char> > isA;

	vector<Atom> m_atom_list;
	vector<Bond> m_bond_list;
	vector< vector<short> > connectedAtomIDs; /* atom IDs connect to each atom */
	vector< vector< char > > connectedBondTypes; /* bond types connect to each atom */
	DistanceMatrix<short> distanceMatrix;

	vector< vector<short> > connectedAtomIDs_H; /* atom IDs connect to each atom */
	vector< vector< char > > connectedBondTypes_H; /* bond types connect to each atom */
	DistanceMatrix<short> distanceMatrix_H;
	vector< vector<short> > atom2bonds_H;
	vector<short> atomic_num_H;


	DistanceMatrix<int> distancePath;
	DistancePathMatrix distancePathMatrix;
	vector<short> m_elementCounts;

	inline void ConnectivityOrder2Index(PathNode *path, double & vci, double & vvci, double & pqvci, int & a);
	inline void savePaths(vector< PathNode * > &pathList, vector<short> & len, FILE * pFile);
	void calculateMoleculePathCountFromFile(FILE * pFile);
	void calculateMoleculePathCountFromArrays();
	void clearPathsSet();

	void calculateDistanceMatrix(bool withH);
	void calculatDistanceMatrixOnly();
	void MolecularPathCount(FILE * pFile);
	void initializeConnectedAtomIDs(bool withH);
	void initializeConnectedBondTypes(bool withH);
	void initPqn_valence(void);
	vector< vector< PathNode * > > & getPathSets(void);
	short getPolarizabilityValueID(char symbol , short i);

	void removeAtomsAndBonds(set<short> & atomIDs);
	int findShortestPath(bool withH, short from, short to, vector<short> & shortestPath);
	void dfSearch(short from, short to, vector< PathNode *> & pathList);
	inline void clearPathList(vector< PathNode * > & pathList);

	void checkHydrogenNums(void);
	void removeH(void);
	void removeBond(short index);
	void removeAtom(short index);
	void addHydrogen(void);

public:
	static const short H_ADDED = 2;
	static const short H_DELETED = 1;
	static const short H_ORIGINAL = 0;
	char m_deletePart;
	char m_maxBondType;

	static const char C = 0;
	static const char N = 1;
	static const char O = 2;
	static const char P = 3;
	static const char S = 4;
	static const char F = 5;
	static const char Cl = 6;
	static const char Br = 7;
	static const char I = 8;
	static const char X = 9; //(X = F, Cl, Br, I; A = Aromatics

	string m_FieldName;
	int m_ID;

	Molecule(void);
	~Molecule(void);
	void clear(void);

	void print(ostream & out, bool hasH, bool isSDF);
	void pushHeader(string & line);
	void printFeatures(ostream & out);
	void showFingerprint(ostream & out, short id);

	bool isFeaturesReady(short startLocation);
	vector<double> getFeatures(short startLocation);

	DistancePathMatrix & getDistancePathMatrix();
	bool isAllPathFailed(){return this->m_allPathFailed;}
	DistanceMatrix<short> & getDistanceMatrix(bool withH);
	vector< vector<short> >& getConnectedAtomIDs(bool withH);
	vector< vector<char> >& getConnectedBondTypes(bool withH);
	DistanceMatrix<int> & getDistancePath(void);
	vector<float> & getConventionalBondOrder(void);
	vector< vector<short> > & getRings(void);
	vector<char> & getCon_to_H(void);
	vector<char> & getIsAromatic(void);
	vector< vector<char> > & getIsA();
	DistanceMatrix<int> & getAllPath(void);
	DistanceMatrix<short> & getDetourMatrix(void);
	vector< vector<int> > & getPath(void);
	vector< vector<short> > & getAtom2bonds(bool withH);
	vector<char> & getInRing(void);
	vector<float> & getValence(void);
	vector<short> & getPqn(void);
	vector<float> & getPol(void);
	vector<int> & getP();
	vector<int> & getA();
	vector<double> & getVC();
	vector<double> & getVVC();
	vector<double> & getPQVC();

	void setNumAtoms(short nAtoms);
	void setNumBonds(short nBonds);

	void setChiralFlag(short nChiralFlag);
	void setStextEntriesNum(short nStextEntries);
	void setRingStatus(bool status);
	template<class T> void setFingerprint(vector<T> & subFeatures, short startLocation);
	void setFingerprint(double subFeature, short location);

	void pushAtom(Atom& atom);
	void pushBond(Bond& bond);
	void preprocess(void);
	short getAtomsNum(bool withH);
	short getBondsNum(bool withH);
	Atom & getAtom(short index);
	Bond & getBond(short index);

  void setHeaderMoleculeId(string hmolId);
  string getHeaderMoleculeId();

	// generate signature
	vector<short> & countElements(void);
	vector<double> getMoleculeWeight(void);
	vector<short> & getAtomicNum(bool withH);
};
