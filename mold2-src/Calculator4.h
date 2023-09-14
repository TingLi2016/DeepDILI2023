#pragma once

class Calculator4
{
public:
	Calculator4(Molecule & _mol);
	~Calculator4(void);
	void calculateMLogP(void);
private:
	Molecule & mol;
	vector< vector<char> > & connBond;
	vector< vector<short> > & connAtoms;
	vector<char> & inRing;
	vector<char> & aromatic;
	vector< vector<short> > & rings;
	vector< vector<char> > & isA;
	vector<char> & con_to_H;

	int countOfUnsaturatedBonds(void);
	vector<float> atomicCount(void);
	void stepOne(double & fprx, double & fqn, double & fncs, double & famp, short & nno2);
	short countIrng(void);
	void calculateIalk(short i, short & ialk);
	void calculateIbetal(short i, short & ibetal);
	void checkIhb5_6ring(short i,short &w, short & ihb);
	void calculateNpol(short i, short & w, short & npol);
	void ihbStepOne(short i, short & w, short & ihb);
};
