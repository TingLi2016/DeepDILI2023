#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Calculator2.h"


Calculator2::Calculator2(Molecule & _mol)
:mol(_mol), 
con_to_H(mol.getCon_to_H()),
isA(mol.getIsA()),
aromatic(mol.getIsAromatic()),
connBond(mol.getConnectedBondTypes(false)),
connAtoms(mol.getConnectedAtomIDs(false))
{
	acharge.resize(_mol.getAtomsNum(false), 0);
	results.resize(Elements::TOTALELEMENTS + 1, 0);
	oh = c3n = c2n = c2nh = cn = cnn = cn1 = cn2 = 0;
	cn2h = co = c2o = coc = cs = c2s = cs1 = x = 0;
	ph = pho = phs = phn2 = chain = chaino = chains = 0;
	chainn2 = cc = nn = n2n = nc = no = n2o = so = s2o = 0;
	AliphaticAromatic();
}

Calculator2::~Calculator2(void){
}

void Calculator2::transformCharge(){
	short noHnumAtoms = mol.getAtomsNum(false);
	for(short i = 0; i < noHnumAtoms; i++){	
		Atom & atom = mol.getAtom(i);
		char ch = atom.getCharge();
		if(ch != 4 && ch != 0)
			acharge[i] = 4 - ch;
	}
}

void Calculator2::C_nbi2(short i){
	for(short j = 0; j < 2; j++){
		if(isA[connAtoms[i][j]][Molecule::N] == 1){// C*N
			if(connBond[i][j] == Bond::SINGLE_BOND) cn++;	// C-N
			else if(connBond[i][j] == Bond::DOUBLE_BOND){ //C=NH
				if(connAtoms[connAtoms[i][j]].size() == 1 && con_to_H[connAtoms[i][j]] == 1) c2nh++;
			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
				if(connAtoms[connAtoms[i][j]].size() == 2){
					if(aromatic[connAtoms[connAtoms[i][j]][k]] == 1) phn2++; // ph-n=
					else if(connAtoms[connAtoms[i][j]][k] != i) chainn2++; // chainn2++;		
				}
			}
			}else if(connBond[i][j] == Bond::TRIPLE_BOND){
				if(connAtoms[connAtoms[i][j]].size() == 1) c3n++; //C3N
			}

			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){ // C*N-OH
				if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::O] == 1 && con_to_H[connAtoms[connAtoms[i][j]][k]] == 1) oh++;
			}
		}else if(isA[connAtoms[i][j]][Molecule::O] == 1){ // C*O
			if(connBond[i][j] == Bond::DOUBLE_BOND && connAtoms[connAtoms[i][j]].size() == 1) c2o++;
			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){					
				if(connBond[i][j] == Bond::SINGLE_BOND && connAtoms[connAtoms[i][j]].size() == 2){
					if(aromatic[connAtoms[connAtoms[i][j]][k]] == 1) pho++; // ph-o-
					else if(connAtoms[connAtoms[i][j]][k] != i) chaino++; // chain-o-
				}
			}
		}else if(isA[connAtoms[i][j]][Molecule::S] == 1){ //C*S
			if(connBond[i][j] == Bond::DOUBLE_BOND && connAtoms[connAtoms[i][j]].size() == 1) c2s++;
			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){	
				if(connBond[i][j] == Bond::SINGLE_BOND && connAtoms[connAtoms[i][j]].size() == 2 ){
					if(aromatic[connAtoms[connAtoms[i][j]][k]] == 1) phs++; //ph-s- 
					else if(connAtoms[connAtoms[i][j]][k] != i) chains++; // chains++;
				}
			}
		}else if(isA[connAtoms[i][j]][Molecule::C] == 1){ // C*C
			if(connBond[i][j] == Bond::SINGLE_BOND){ // C-C
				if(aromatic[connAtoms[i][j]] == 1) ph++; // C-C=
				else chain++;
			}
		}			
	} /* end of FOR-j */

	/* --------------- count components and get results ----------------------- */
	if(ph == 1){
		if(c2o == 1 && con_to_H[i] == 1) results[32]++;	/* nCOHPh */
		if(c3n == 1) results[48]++;	/* nCNPh */
		if(c2nh == 1) results[50]++;	/* nC=NPh */
		if(cn == 1 && oh == 1) results[56]++;	/* nCNOHPh */
	}
	if(chaino == 1 && c3n == 1) results[1]++;	/* nOCN */
	if(pho == 1 && c3n == 1) results[2]++;	/* nOCNPh */
	if(chainn2 == 1){
		if(c2o == 1) results[3]++;	/* nNCO */
		if(c2s == 1) results[7]++;	/* nNCS */	
	}
	if(chains == 1 && c3n == 1) results[5]++;	/* nSCN */
	if(phs == 1 && c3n == 1) results[6]++;	/* nSCNPh */
	if(phn2 == 1){
		if(c2o == 1) results[4]++;	/* nNCOPh */
		if(c2s == 1) results[8]++;	/* nNCSPh */
	}
	if(chain == 1){
		if(c2o == 1 && con_to_H[i] == 1) results[31]++;	/* nCOH */
		if(c3n == 1) results[47]++;	/* nCN */
		if(c2nh == 1) results[49]++;	/* nC=N */
		if(cn == 1 && oh == 1) results[55]++;	/* nCNOH */
	}	
}


void Calculator2::C_nbi3(short i){
	for(short j = 0; j < 3; j++){
		if(isA[connAtoms[i][j]][Molecule::X] == 1) x++;	// C-X
		else if(isA[connAtoms[i][j]][Molecule::O] == 1){ //C*O
			short nbij = connAtoms[connAtoms[i][j]].size();
			if(nbij == 1){
				if(connBond[connAtoms[i][j]][0] == Bond::DOUBLE_BOND) c2o++;// C=O
				else if(connBond[connAtoms[i][j]][0] == Bond::SINGLE_BOND) co++;//C-O
			}else if(nbij == 2){
				for(short k = nbij - 1; k >= 0; k--){
					if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1 && connBond[connAtoms[i][j]][k] == Bond::SINGLE_BOND && connAtoms[connAtoms[i][j]][k] != i){
						coc++;	// C-O-C
						if(aromatic[connAtoms[connAtoms[i][j]][k]] == 1) pho++;	// C-O-Ph
					}
				}
			}
		}else if(isA[connAtoms[i][j]][Molecule::N] == 1){ //C*N
			if(connBond[i][j] == Bond::SINGLE_BOND){
				if(con_to_H[connAtoms[i][j]] == 2) cn2h++;	//C-N=2H
				else cn++;	// C-N
			}
			short nbij = connAtoms[connAtoms[i][j]].size();
			if(nbij == 2){
				if(connBond[connAtoms[i][j]][0] == Bond::SINGLE_BOND) cn1++; // C-N-
				if(connBond[i][j] == Bond::SINGLE_BOND){
					for(short k = nbij - 1; k >= 0; k--){ // C-N-N							
						if( isA[connAtoms[connAtoms[i][j]][k]][Molecule::N] == 1 && connAtoms[connAtoms[i][j]][k] != i && connBond[connAtoms[i][j]][k] == 1) cnn++;
					}
				}
			}else if(nbij == 3 && connBond[i][j] == Bond::SINGLE_BOND) cn2++; // C-N= 
		}else if(isA[connAtoms[i][j]][Molecule::S] == 1){ //C*S
			short nbij = connAtoms[connAtoms[i][j]].size();
			if(nbij == 1){
				if(connBond[i][j] == Bond::SINGLE_BOND) cs++; // C-S
				else if(connBond[i][j] == Bond::DOUBLE_BOND) c2s++; // C=S
			}else if(nbij == 2 && connBond[i][j] == Bond::SINGLE_BOND) cs1++; // C-S-			
		}else if(isA[connAtoms[i][j]][Molecule::C] == 1  && aromatic[connAtoms[i][j]] != 1) chain++; // chain				

		if(aromatic[connAtoms[i][j]] == 1) ph++; // aromatic 

	} /* end of FOR-j */

	/* --------------- count components and get results ----------------------- */
	if(ph == 1){
		if(c2o == 1){
			if(co == 1) results[10]++;	/* nCOOHPh */
			if(coc == 1) results[12]++;	/* nCOORPh */
			if(cn2h == 1) results[14]++;	/* nCONH2Ph */
			if(cn1 == 1) results[16]++;	/* nCONHRPh */
			if(cn2 == 1) results[18]++;	/* nCONR2Ph */
			if(x == 1) results[22]++;	/* nCOXPh */
			if(cs1 == 1) results[28]++;	/* nCOSRPh */
			if(chain == 1) results[34]++;	/* nCOPh */
			if(cnn == 1) results[36]++;	/* nCONNPh */
		}
		if(c2s == 1){
			if(co == 1) results[24]++;	/* nCSOHPh */
			if(cs == 1) results[26]++;	/* nCSSHPh */
			if(cs1 == 1) results[28]++;	/* nCSSRPh */
		}		
	}
	if(ph == 0){
		if(chain == 1){
			if(c2s == 1){
				if(cs1 == 1) results[27]++;	/* nCSSR */
				if(c2s == 1 && cs == 1) results[25]++;	/* nCSSH */
				if(c2s == 1 && co == 1) results[23]++;	/* nCSOH */
			}
			if(c2o == 1){
				if(coc == 1) results[11]++;	/* nCOOR */
				if(co == 1) results[9]++;	/* nCOOH */
				if(cnn == 1) results[35]++;	/* nCONN */			
				if(cs1 == 1) results[27]++;	/* nCOSR */			
				if(x == 1) results[21]++;	/* nCOX */
				if(cn2 == 1) results[17]++;	/* nCONR2 */	
				if(cn1 == 1) results[15]++;	/* nCONHR cn2==2: N has other two compositions to connecting it */			
				if(cn2h == 1) results[13]++;	/* nCONH2 */			
			}
		}else if(chain == 2 && c2o == 1) results[33]++;	/* nCO */
	}
	if(coc == 1 && pho == 0 && c2o == 1 && cn == 1) results[19]++;	/* nOCONPh */
	if(pho == 1 && c2o == 1 && cn == 1) results[20]++;	/* nOCON */
}


void Calculator2::C_ring(short i, short nai){
	vector<char> & inRing = mol.getInRing();
	if(inRing[i] == 1){ // i in a ring
		for(short k = 0; k < nai; k++){
			if(inRing[connAtoms[i][k]] == 1 && connBond[i][k] == Bond::DOUBLE_BOND){ 
				cn2++;
				break;
			}
		}
	}

	for(short j = 0; j < nai; j++){		
		if(isA[connAtoms[i][j]][Molecule::C] == 1){ // C
			if(connAtoms[i][j] == Bond::SINGLE_BOND) cc++; // C-C
			else if(connAtoms[i][j] == Bond::DOUBLE_BOND) c2n++;	// C=C
			else if(connAtoms[i][j] == Bond::TRIPLE_BOND) c3n++;	//C#C
		}else if(isA[connAtoms[i][j]][Molecule::X] == 1){ // X	halogens
			x++; 
			if(nai >= 2){
				short nn = 0;
				for(short n = 0; n < nai; n++){	/* look for halogens */
					if(isA[connAtoms[i][n]][Molecule::X] == 1) nn++;
				}

				if(nn == 1){ /* the atom-i only connected one halogens */
					nn = 0;
					for(short n = 0; n < nai; n++){	/* research all connection atoms of atom-i */
						if(isA[connAtoms[i][n]][Molecule::C] == 1 && connBond[i][n] == Bond::DOUBLE_BOND){
							for(short k = connAtoms[connAtoms[i][n]].size() - 1; k >= 0; k--){
								if(isA[connAtoms[connAtoms[i][n]][k]][Molecule::C] == 1 && connAtoms[connAtoms[i][n]][k] != i && connBond[connAtoms[i][n]][k] == Bond::SINGLE_BOND){/* atom-k connected atom-n with single bond */
									for(short m = connBond[connAtoms[connAtoms[i][n]][k]].size() - 1; m >= 0; m--){
										if(connBond[connAtoms[connAtoms[i][n]][k]][m] >= Bond::DOUBLE_BOND){
											results[101]++;	/* nCconjX */
											nn = -1;
											break;
										}
									}
								}
								if(nn == -1) break;
							}
						}
						if(nn == -1) break;
					}
				}
			}	
		} /* end of if(halogens[aij] == 1) */
	} /* end of FOR-j */

	/* --------------------------- count components and get results --------------------------- */
	if(cc == 1){
		if(con_to_H[i] == 2 && x == 1) results[88]++;	/* nRCH2X */
		if(x == 3) results[97]++;	/* nRCX3 */
		if(con_to_H[i] == 1 && x == 2) results[94]++;	/* nRCHX2 */
	}
	if(cc == 2 && con_to_H[i] == 1 && x == 1) results[89]++;	/* nR2CHX */
	if(cc == 3 && x == 1) results[90]++;	/* nR3CX */
	if(con_to_H[i] == 1 && c2n == 1 && x == 1) results[91]++;	/* nR=CHX */
	if(c2n == 1){
		if(cc == 1 && x == 1) results[92]++;	/* nR=CRX */
		if(x == 2) results[96]++;	/* nR=CX2 */	
	}
	if(c3n == 1 && x == 1) results[93]++;	/* nR#CX */	
	if(cc == 2 && x == 2) results[95]++;	/* nR2CX2 */	
	if(aromatic[i] == 1 && x == 1) results[98]++;	/* nPhX */
	if(inRing[i] == 1){
		if(x == 1){
			results[99]++;	/* nCXr */
			if(cn2 == 1) results[100]++;	/* nCXr= */
		}
	}
}

void Calculator2::N_search(short i, short nai){
	cn = c2o = cc = nn = n2n = nc = no = n2o = oh = ph = 0;
	/* the N has only one bond: single bond, it is also connecting to two atomic H */
	if(nai == 1 && connBond[i][0] == Bond::SINGLE_BOND && isA[connAtoms[i][0]][Molecule::C] == 1 &&	con_to_H[i] == 2){
		if(aromatic[connAtoms[i][0]] == 1) results[38]++;	// connecting to aromatic
		else results[37]++;	/* connecting to chain */
	}

	short m = 0;
	for(short j = 0; j < nai; j++){
		if(isA[connAtoms[i][j]][Molecule::C] == 1 && connBond[i][j] == Bond::SINGLE_BOND){
			if(aromatic[connAtoms[i][j]] == 1) ph++;	/* aromatic */
			else nc++;	/* chain */
		}else if(isA[connAtoms[i][j]][Molecule::N] == 1){
			if(con_to_H[connAtoms[i][j]] == 2) nn++;
			else if(con_to_H[connAtoms[i][j]] == 1 && connBond[connAtoms[i][j]][0] == Bond::SINGLE_BOND) n2n++;
		}else if(isA[connAtoms[i][j]][Molecule::O] == 1){ /* atom-j is O */
			if(connBond[i][j] == Bond::SINGLE_BOND){ /* atom-i and the O connecting with single bond */
				no++;	/* N-O */
				if(con_to_H[connAtoms[i][j]] == 1 && nai == 1) oh++; // OH
			}else if(connBond[i][j] == Bond::DOUBLE_BOND) n2o++; // N=O
		}else if(isA[connAtoms[i][j]][Molecule::C] == 1 && connAtoms[connAtoms[i][j]].size() == 3){
			cn++;	/* C-N */
			for(short k = 2; k >= 0; k--){
				short aijk = connAtoms[connAtoms[i][j]][k];
				char bijk = connBond[connAtoms[i][j]][k];
				if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1 && connBond[connAtoms[i][j]][k] == Bond::SINGLE_BOND) cc++;	// C-C
				else if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::O] == 1 && connBond[connAtoms[i][j]][k] == Bond::DOUBLE_BOND) c2o++; // C=O
			}
		}
		if(connBond[i][j] == Bond::SINGLE_BOND) m++;	/* single bond plus one */
	} /* end of FOR-j */

	if(m == nai) results[103]++;	/* nHAcc */
	/* --------------- count components and get results ----------------------- */
	if(ph == 0){
		if(nc > 1 && con_to_H[i] == 1) results[39]++;	/* nNHR */
		if(nc >= 3 && c2o == 0) results[41]++;			/* nNR2 */
		if(nc >= 1){
			if(nn >= 1) results[43]++;	/* nN-N */
			results[51] += acharge[i];	/* nN+ */
			if(con_to_H[i] == 1 && oh >= 1) results[53]++;	/* nNHOH */
			if(n2o >= 1) results[57]++;	/* nNNOx */
			if(n2o >= 1) results[59]++;	/* nNO */	
			if(n2o == 2) results[61]++;		/* nNO2 */
			if(n2o == 1 && no >= 1) results[61]++;	/* nNO2 */
		}		
	}else if(ph >= 1){
		if(nc >= 1 && con_to_H[i] == 1) results[40]++;	/* nNHRPh */
		if(nc == 2 && c2o == 0) results[42]++;		/* nNR2Ph */
		if(nn >= 1) results[44]++;	/* nN-NPh */
		results[52] += acharge[i];	/* nN+Ph */
		if(con_to_H[i] == 1 && oh >= 1) results[54]++;	/* nNHOHPh */
		if(n2o >= 1) results[58]++;	/* nNNOxPh */		
		if(no >= 1) results[60]++;	/* nNOPh */
		if(n2o == 2) results[62]++;		/* nNO2Ph */
		else if(n2o == 1 && no >= 1) results[62]++;	/* nNO2Ph */
	}
	if(cn == 2 && c2o == 2 && cc == 2) results[63]++;	/* nN(CO)2 */
	results[102] += con_to_H[i];	/* nHDon */
}

void Calculator2::O_search(short i, short nai){
	co = ph = x = 0;
	if(con_to_H[i] == 2) results[73]++;		/* nH2O */
	results[102] += con_to_H[i];	/* nHDon */

	short m = 0;
	for(short j = 0; j < nai; j++){ /* check all atoms that connecting to atom-i */
		short aij = connAtoms[i][j];
		if(nai == 1 && con_to_H[i] == 1 && isA[connAtoms[i][0]][Molecule::C] == 1){	/* it is connecting to C */
			results[64]++;		/* nOH */
			if(aromatic[connAtoms[i][0]] == 1) results[65]++;	/* nOHPh */
			char nai1_H = con_to_H[connAtoms[i][0]];
			short nbi1 = connAtoms[connAtoms[i][0]].size();
			if(nai1_H == 2 && nbi1 == 2) results[66]++;	/* nOHp */
			else if(nai1_H == 1 && nbi1 == 3) results[67]++;	/* nOHs */
			else if(nai1_H == 0 && nbi1 == 4) results[68]++;	/* nOHt */
		}else if(nai == 2){	/* atom-i has two bonds */
			if(isA[connAtoms[i][j]][Molecule::C] == 1){	/* the atom-j is C */
				if(aromatic[connAtoms[i][j]] != 1) co++;	/* C-O (it is not included aromatic) */
				else ph++;
			}
		}
		if(isA[connAtoms[i][j]][Molecule::X] == 1) x++;
		if(isA[connAtoms[i][j]][Molecule::C] == 1) m++;	/* the atom-i has C to connecting */				
	} /* end of FOR-j */

	if(m != 0) results[103]++;	/* nHAcc */
	/* --------------- count components and get results ----------------------- */
	if(co == 2) results[69]++;	/* nROR */
	if(ph == 1 && co == 1) results[70]++;	/* nRORPh */
	if(co >= 1 && x == 1) results[71]++;	/* nROX */
	if(ph >= 1 && x == 1) results[72]++;	/* nROXPh */
}

void Calculator2::S_search(short i, short nai){
	cc = cs = c2s = no = oh = so = s2o = 0;
	for(short j = 0; j < nai; j++){
		if(isA[connAtoms[i][j]][Molecule::C] == 1){ /* atom-j is C */
			if(connBond[i][j] == Bond::SINGLE_BOND)	cs++;	/* C-S */
			if(nai == 1 && connBond[i][0] == Bond::DOUBLE_BOND) c2s++;	/* C=S */
			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
				if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1) cc++;	/* C-C */
			}
		}else if(isA[connAtoms[i][j]][Molecule::O] == 1){ /* atom-j is O */
			if(connBond[i][j] == Bond::SINGLE_BOND){	/* S connecting O with single bond */
				so++;		 /* S-O */
				if(con_to_H[connAtoms[i][j]] == 1) oh++;	 /* OH */						
			}else if(connBond[i][j] == Bond::DOUBLE_BOND) s2o++;		 /* S=O */
		}else if(isA[connAtoms[i][j]][Molecule::N] == 1) no++;/* number of N */

		if(nai == 2){	/* atom-i has two bonds */
			if(isA[connAtoms[i][0]][Molecule::C] == 1){
				if(isA[connAtoms[i][1]][Molecule::C] == 1) results[79]++;	/* nRSR */
				else if(isA[connAtoms[i][1]][Molecule::S] == 1){
					for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
						if(isA[connAtoms[i][j]][Molecule::S] == 1 && isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1) results[80]++;	/* nRSSR */
					}
				}
			}
			if((isA[connAtoms[i][0]][Molecule::C] == 1 && isA[connAtoms[i][1]][Molecule::S] == 1) || (isA[connAtoms[i][1]][Molecule::C] == 1 && isA[connAtoms[i][0]][Molecule::S] == 1)){
				for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
					if(isA[connAtoms[i][j]][Molecule::S] == 1 && isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1) results[80]++;	/* nRSSR */
				}
			}
		}
	} /* end of FOR-j */

	/* --------------- count components and get results ----------------------- */
	if(so == 0 && s2o == 1 && nai == 3) results[74]++;	/* nSO */
	if(so == 0 && s2o == 2 && nai == 3) results[75]++;	/* nSO2 */
	if(so == 1 && s2o == 2 && cs == 1 && nai == 4) results[76]++;	/* nSO3 */
	if(con_to_H[i] == 1 && nai == 1 && connBond[i][0] == Bond::SINGLE_BOND) results[77]++;	/* nSH */
	if(c2s == 1 && cc == 2) results[78]++;	/* nCS */
	if(s2o == 2 && oh == 1 && nai == 4) results[81]++;	/* nSO3H */
	if(s2o == 2 && no == 1 && nai == 4) results[82]++;	/* nSO2N */
}

void Calculator2::P_search(short i, short nai){
	no = n2o = oh = cs = cs1 = cc = co = 0;	
	for(short j = 0; j < nai; j++){
		short aij = connAtoms[i][j];
		char bij = connBond[i][j];
		short naij = connAtoms[aij].size();
		if(isA[connAtoms[i][j]][Molecule::O] == 1){	/* atom-j is O */
			if(connBond[i][j] == Bond::SINGLE_BOND){	/* aotm-i and atom-j connecting with single bond */
				no++;	/* number of -O */
				for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
					if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1) co++;	/* C-O */
				}
			}else if(connBond[i][j] == Bond::DOUBLE_BOND) n2o++;	/* number of =O */
			if(con_to_H[connAtoms[i][j]] == 1) oh++;/* number of OH */
		}else if(isA[connAtoms[i][j]][Molecule::S] == 1){	/* atom-j is S */
			cs++;	/* number of S */
			for(short k = connAtoms[connAtoms[i][j]].size() - 1; k >= 0; k--){
				if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::C] == 1){
					cs1++;
					break;
				}
			}
		}else if(isA[connAtoms[i][j]][Molecule::C] == 1) cc++;	/* number of C */
	}
	/* --------------- count components and get results ----------------------- */
	if(no == 2 && n2o == 1 && cc == 1) results[83]++;	/* nPO3 */
	if(no == 3 && n2o == 1 && co >= 1) results[84]++;	/* nPO4 */
	if(no == 2 && n2o == 1 && cs == 1 && cs1 == 1) results[85]++;	/* nPO3S */
	if(no == 1 && n2o == 1 && cs == 2 && cs1 == 2) results[86]++;	/* nPO2S2 */
	if(no == 1 && n2o == 1 && cs == 1 && cs1 == 1 && cc == 1) results[87]++;	/* nPO2SR */
}

void Calculator2::AliphaticAromatic(){
	try{
		if(mol.isFeaturesReady(::AliphaticAromatic)) return;
		transformCharge();	
		short noHnumAtoms = mol.getAtomsNum(false);

		for(short i = 0; i < noHnumAtoms; i++){	/* compare all atom */
			short nai = connAtoms[i].size(); // nbi == nai

			if(isA[i][Molecule::C] == 1){ /* atomic symbol C */
				oh = c3n = c2n = c2nh = cn = cnn = cn1 = cn2 = cn2h = co = c2o = coc = c2o = cs = 0;
				c2s = cs1 = x = ph = pho = phs = phn2 = chain = chaino = chains = chainn2 = 0;

				/* ---------------------- atom-i has two bonds ----------------------------- */
				if(nai == 2) C_nbi2(i);
				else if(nai == 3) C_nbi3(i);
				cc = c2n = c3n = cn1 = cn2 = x = 0;
				C_ring(i, nai);
			} else if(isA[i][Molecule::N] == 1){ /* atomic symbol N */
				N_search(i, nai);
			} else if(isA[i][Molecule::O] == 1){ /* the atom-i is O */
				O_search(i, nai);
			} else if(isA[i][Molecule::S] == 1){ /* atom-i is S */
				S_search(i, nai);
			} else if(isA[i][Molecule::P] == 1){	/* the atom-i is P */
				P_search(i, nai);
			} /* end of if(phosphorus[i] == 1) */

			if(mol.getAtom(i).getAtomSymbol() == Elements::FLOURINE)	/* atom-i is F */
				results[103]++;	/* nHAcc */
		} /* end of FOR-i*/

		vector<double> fingerprint;
		for(short i = 1; i <= 103; i++)
			fingerprint.push_back(results[i]);
		mol.setFingerprint(fingerprint, ::AliphaticAromatic);
	}catch(exception & ex){
		cout << ex.what() << " thrown in AliphaticAromatic"<< endl;
	}
}
