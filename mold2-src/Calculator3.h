#pragma once

class Calculator3
{
private:
	Molecule & mol;
	vector<char> & con_to_H;
	vector< vector<char> > & isA;
	vector< vector<char> > & connBond;
	vector< vector<short> > & connAtoms;
	vector<char> & aromatic;
		
	vector<char> bonds;
	short cc, c2, c3, co, c2o, c2c, cx, c2x, c3x, car, arc, arx, al, ar, al2x, sp3, oc, oal, oar, oh, o1, o2, o3, ro, nx, n2x, sr, s2r, s2o, ssr, x2, x2x, nc2o, ncc, h, m;
	vector<char> used;
	vector<short> results;

	void IisC(short j, short bij, short aij);
	void IisO(short i, short bij, short aij);
	void IisN(short i, short bij, short aij);
	void IisS(short i, short bij, short aij);
	void calculate_C_results(short i);
	void calculate_H_results(short i);
	void calculate_O_results(short i);
	void calculate_N_results(short i);
	void calculate_S_results(short i);	
public:
	Calculator3(Molecule &_mol);
	~Calculator3(void);
	void AtomCentredFragments(void);
};
