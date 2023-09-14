#pragma once

class Calculator7
{
	Molecule & mol;
	short numAtoms;		
	void gotoNextLevel(vector<short>& currentLevel, vector<short>& nextLevel);

public:
	Calculator7(Molecule & _mol);
	~Calculator7(void);
	void PathWalk();
	void KierPathIndex();
	void AllPathConnectivity();
	void MinMaxPathIndex();
	void MaximumPathIndex();
	void AllPathWiener();
	void D_D_RingAndCircuitsIndex();
	void doAll(void );
};
