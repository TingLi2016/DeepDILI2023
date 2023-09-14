#pragma once

class Calculator2
{
private:
	vector<char> acharge;
	void transformCharge();
	vector<short> results;
	short oh, c3n, c2n, c2nh, cn, cnn, cn1, cn2, cn2h, co, c2o, coc, cs, c2s, cs1, x, ph, pho, phs, phn2, chain, chaino, chains, chainn2;
	short cc, nn, n2n, nc, no, n2o, so, s2o;	/* number of C-C */

public:
	Calculator2(Molecule & mol);
	~Calculator2(void);
	void AliphaticAromatic(void);
private:
	Molecule & mol;
	vector<char> & con_to_H;
	vector< vector<char> > & isA;
	vector<char> & aromatic;
	vector< vector<char> > & connBond;
	vector< vector<short> > & connAtoms;
	void C_nbi2(short i);
	void C_nbi3(short i);
	void C_ring(short i, short nai);
	void N_search(short i, short nai);
	void O_search(short i, short nai);
	void S_search(short i, short nais);
	void P_search(short i, short nai);
};
