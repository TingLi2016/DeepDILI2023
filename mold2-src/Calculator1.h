#pragma once

class Calculator1
{
public:
	Calculator1(Molecule & _mol);
	~Calculator1(void);	
	void AromaticBondConverter(void);
	
private:
	Molecule & mol;
	static bitset<2000> used;
	static const short MAX_RING_FUSE = 200;

	int	NMB, NRB, NDB, NAB;
	double RBF, SCBOHD;
	vector< vector<short> > connAtoms;
	vector< pair<short, short> > edges;

	vector<short> p;
	vector<char> mark;
	vector<short> len;
	short numAtoms;
	
	int findNCRS();	
	void findBondProperties(void);
	void changeBondTopology(void);
	void convertAromaticBonds(vector<short> & nr);	
	void ringSearch(void);

	void removeUnrelatedAtoms(void);
	void dfsBackEdges();	
	inline bool brokenSearch(vector<vector<short> > & pathList, pair<short, short> & e1, pair<short, short> & e2);
	inline void cleanUpeSet(vector< set<pair<short, short> > > & eSet);
	bool isAlreadyFound(vector<vector<short> > & pathList, vector<short> & path, short & idx);
	void bfsShortestRing( pair<short, short> & edge, vector< vector<short> > &paths);
	void sharedBackEdge(short currentEdgeID, vector<short> &path, vector<short> & sharedEdges);
};
