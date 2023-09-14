#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Calculator3.h"


Calculator3::Calculator3(Molecule &_mol)
:mol(_mol),
con_to_H(mol.getCon_to_H()),
isA(mol.getIsA()),
connBond(mol.getConnectedBondTypes(false)),
connAtoms(mol.getConnectedAtomIDs(false)),
aromatic(mol.getIsAromatic())		
{
	bonds.resize(mol.getAtomsNum(false), 0);
	used.assign(mol.getAtomsNum(false), 0);
	results.resize(111, 0);

	cc = c2 = c3 = co = c2o = c2c = cx = c2x = c3x = car = 0;
	arc = arx = al = ar = al2x = sp3 = oc = oal = oar = oh = 0;
	o1 = o2 = o3 = ro = nx = n2x = sr = s2r = s2o = ssr = x2 = 0;
	x2x = nc2o = ncc = h = m = 0;
	AtomCentredFragments();
}

Calculator3::~Calculator3(void){
}

void Calculator3::AtomCentredFragments(){
	try{
		if(mol.isFeaturesReady(::AtomCentredFragments)) return;

		short noHnumAtoms = mol.getAtomsNum(false);
		for(short i = 0; i < noHnumAtoms; i++){
			bonds[i] = connBond[i].size();
		}
		for(short i = 0; i < noHnumAtoms; i++){	/* compare all atom */
			// reset counts
			cc = c2 = c3 = co = c2o = c2c = cx = c2x = c3x = car = arc = arx = al = ar = al2x = 0;
			sp3 = oc = oal = oar = oh = o1 = o2 = o3 = ro = 0;
			nx = n2x = sr = s2r = s2o = ssr = x2 = x2x = nc2o = ncc = h = m = 0;

			for(short j = connAtoms[i].size() - 1; j >= 0; j--){
				short aij = connAtoms[i][j];
				short naij = connAtoms[aij].size();
				/* ========================== atom-i is a carbon ========================== */
				IisC(i, connBond[i][j], connAtoms[i][j]);
				/* ========================== atom-i is an oxygen ========================== */
				IisO(i, connBond[i][j], connAtoms[i][j]);
				/* ========================== atom-i is a nitrogen ========================== */
				IisN(i, connBond[i][j], connAtoms[i][j]);			
				/* ========================== atom-i is a sulfur ========================== */
				IisS(i, connBond[i][j], connAtoms[i][j]);
			} /* end of FOR-j */

			/* ------------------------- carbon ----------------------- */
			calculate_C_results(i);
			/* ------------------------ Hydrogen ----------------------- */
			calculate_H_results(i);
			/* ------------------------ oxygen ----------------------- */
			calculate_O_results(i);
			/* ------------------------ Nitrogen ----------------------- */
			calculate_N_results(i);		
			/* ------------------------ Sulfur ----------------------- */
			calculate_S_results(i);
		} /* end of FOR-i */

		/* --------------------- printout results(60) to output file ----------------------- */
		vector<double> fingerprint;
		/* output results: C-001 ~ C-026 */
		for(short i = 1; i <= 26; i++)	
			fingerprint.push_back(results[i]);

		/* output results: C-0036 ~ C-041 */
		for(short i = 36; i <= 41; i++)	
			fingerprint.push_back(results[i]);

		/* output results: H-046, H-050 */
		fingerprint.push_back(results[46]);
		fingerprint.push_back(results[50]);

		/* output results: H-052 ~ H-055, O-056 ~ O-060 */
		for(short i = 52; i <= 60; i++)
			fingerprint.push_back(results[i]);

		/* output results: N-066 ~ N-074 */
		for(short i = 66; i <= 74; i++)
			fingerprint.push_back(results[i]);

		/* output results: N-076 ~ N-078 */
		for(short i = 76; i <= 78; i++)	
			fingerprint.push_back(results[i]);

		/* output results: S-106 ~ S-110 */
		for(short i = 106; i <= 110; i++)	
			fingerprint.push_back(results[i]);
		mol.setFingerprint(fingerprint, ::AtomCentredFragments);
	}catch(exception &ex){
		cout << ex.what() << " thrown in AtomCentredFragments" << endl;
	}
}

/* ========================== atom-i is a carbon ========================== */
void Calculator3::IisC(short i, short bij, short aij){	
	if(isA[i][Molecule::C] == 1){
		if(aromatic[i] != 1){	/* the atom-i is a not aromatic C */
			if(bij == Bond::SINGLE_BOND){ /* the bond that between atom-i and aotm-j is singl bond */
				if(isA[aij][Molecule::C] == 1 && aromatic[aij] != 1) cc++; /* C-C */
				else if(isA[aij][Molecule::C] != 1) cx++; /* C-X */
				if(aromatic[aij] == 1) car++; /* C-aromatic */

				/* --------------- Hydrogen: H-046~H-054 ----------------------- */
				if(isA[aij][Molecule::C] == 1){ /* the atom-j is a C */
					if(con_to_H[i] >= 1){ /* has one or more H connected */
						sp3++;	/* number of single bonds; if sp3 == number of bonds, the atom-i is sp3 atomic */
						for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
							if(isA[connAtoms[aij][k]][Molecule::C] != 1) m++;	/* count atoms except C and H */
						}
					}
				}else if(isA[aij][Molecule::O] == 1){ /* there has a single bond between atom-i and atom-j, atom-j is O */
					co++;	/* C-O */
					if(con_to_H[aij] == 1) h++; // C-OH
				}		
			}else if(bij == Bond::DOUBLE_BOND){ /* the bond type is double bond */
				c2++; /* =C */
				if(isA[aij][Molecule::C] == 1) c2c++;
				else if(isA[aij][Molecule::C] != 1) c2x++; // C=X
				else if(isA[aij][Molecule::O] == 1) c2o++;	// C=O
			}else if(bij == Bond::TRIPLE_BOND){ /* the bond type is triple bond */
				c3++; /* #C */
				if(isA[aij][Molecule::C] != 1) c3x++; /* C#X */					
			}
		}else if(bij == Bond::SINGLE_BOND){ /* the bond type is a single bond */
			if(isA[aij][Molecule::C] == 1) arc++;	/* aromatic-C */
			else arx++;	/* aromatic-X */					
		}
	} /* end of if(carbon[i] == 1 && aromatic[i] != 1) */

	if(isA[aij][Molecule::C] == 1 && aromatic[aij] != 1) {/* the atom-j is a C, it is not an aromatic */
		for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
			if(isA[connAtoms[aij][k]][Molecule::C] != 1 && connBond[aij][k] == Bond::DOUBLE_BOND)
				al2x++;	/* Al=X */
		}
	}
}

/* ========================== atom-i is an oxygen ========================== */
void Calculator3::IisO(short i, short bij, short aij){
	if(isA[i][Molecule::O] == 1){	/* the atom-i is an O */
		if(bij == Bond::SINGLE_BOND){
			o1++;	/* -O */
			if(isA[aij][Molecule::C] == 1){ /* the atom-j is C */
				oc++;	/* C-O */
				if(aromatic[aij] != 1){ /* it is not an aromatic C */
					oal++;	/* O-Al */
					short n = 0;	
					for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
						if(connBond[aij][k] == Bond::DOUBLE_BOND ||	connBond[aij][k] == Bond::AROMATIC_BOND) /* double bond, or aromatic bond */
							break;
						else if(isA[connAtoms[aij][k]][Molecule::C] != 1 &&	connAtoms[aij][k] != i)	/* it is alse not the atom-i: is an another heteroatom */
							break;
						else n++;	/* count the number of bonds */
					}
					/* C-OH:	1. the O has only one H and one C to connecting with single bond 
					2. the C without double bond and aromatic bonds 
					3. the C only connected C */
					if(n == connAtoms[aij].size() && bonds[i] == 1 && con_to_H[i] == 1) results[56]++;	/* O-056: alcohol */
					n = 0;
					short p = 0;
					/* if the atom-j has three bonds and atom-i only has one bond */
					if(connAtoms[aij].size() == 3 && bonds[i] == 1){
						for(short k = 2; k >= 0; k--){
							if(connBond[aij][k] == Bond::DOUBLE_BOND){
								if(isA[connAtoms[aij][k]][Molecule::C] == 1) results[57]++;	/* O-057: enol */
								else if(isA[connAtoms[aij][k]][Molecule::O] == 1) n++;
							}else if(connBond[aij][k] == Bond::SINGLE_BOND && isA[connAtoms[aij][k]][Molecule::C] == 1) p++;
						}
					}
					if(n == 1 && p == 1 && bonds[aij] == 3 && bonds[i] == 1 && con_to_H[i] == 1)
						results[57]++;	/* O-057: carboxyl OH */
				} /* end of if(aromatic[aij] != 1) */

				if(aromatic[aij] == 1){ /* the atom-i connected to aromatic */
					oar++;	/* O-aromatic */
					if(con_to_H[i] == 1) results[57]++;	/* O-057: Phenol */							
				}
			} /* end of if(carbon[aij] == 1) */
		}else if(bij == Bond::DOUBLE_BOND) o2++;	/* =O */
		else if(bij == Bond::TRIPLE_BOND) o3++;	/* #O */				
	} /* end of if(oxygen[i] == 1) */

	if(bij == Bond::SINGLE_BOND){
		/* count the next C atom and differentiate it an aromatic or no */
		if(isA[aij][Molecule::C] == 1){
			if(aromatic[aij] != 1)/* the atom-j is a C , it is not an aromatic C */
				al++;	/* the next atom (atom-j) is not an aromatic C */
			else 
				ar++;	/* the next atom (atom-j) is an aromatic C */
			if(connAtoms[aij].size() == 3){
				for(short k = 2; k >= 0; k--){	/* runs over all bonds */
					if(isA[connAtoms[aij][k]][Molecule::O] == 1 && connBond[aij][k] == Bond::DOUBLE_BOND)/* the atom-k is an O, it is a double bond that between atom-j and atom-k */
						nc2o++;	/* -C=O */						
					else if(isA[connAtoms[aij][k]][Molecule::C] == 1 && connBond[aij][k] == Bond::SINGLE_BOND) /* the atom-k is an O, it is a single bond that between atom-j and atom-k */
						ncc++;	/* -C-C */
				}
			}
		}				
	}

	if(isA[aij][Molecule::O] == 1){ /* atom-j is an O */
		if(con_to_H[aij] == 1) oh++;	/* OH */
		if(bij == Bond::SINGLE_BOND){ /* the atom-j connected to atom-i with single bond */
			/* check the third atoms */
			for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
				/* the third atom is a C */
				if(isA[connAtoms[aij][k]][Molecule::C] == 1) ro++;	/* C-O */
			}
		}
	} /* end of if(oxygen[aij] == 1) */
}

/* ========================== atom-i is a nitrogen ========================== */
void Calculator3::IisN(short i, short bij, short aij){	
	if(isA[i][Molecule::N] == 1){
		if(isA[aij][Molecule::C] != 1){ /* the atom-j is not a C */
			if(bij == Bond::SINGLE_BOND){	/* it is connecting to atom-i with single bond */
				nx++;	/* N-X */
				for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
					if(connBond[aij][k] == Bond::DOUBLE_BOND && isA[connAtoms[aij][k]][Molecule::C] != 1) x2x++;	/* -X=X */
				}
			}else if(bij == Bond::DOUBLE_BOND) n2x++;	/* N=X */
		}
	}
}

void Calculator3::IisS(short i, short bij, short aij){
	if(isA[i][Molecule::S] == 1){
		if(isA[aij][Molecule::C] == 1 && bij == Bond::SINGLE_BOND) sr++;	/* S-R */
		if(connBond[i].size() == 2){
			if(isA[aij][Molecule::S] == 1){ /* the atom-j is a S */
				if(used[i] != 1 && used[aij] != 1){ /* the two atoms never used */
					for(short k = connAtoms[aij].size() - 1; k >= 0; k--){
						if(isA[connAtoms[aij][k]][Molecule::C] == 1 && connBond[aij].size() == 1){
								ssr++;	/* S-S-R */
								used[i] = used[aij] = 1;
						}
					}
				}
			}
		} /* end of if(nbi == 2)*/

		if(bij == Bond::DOUBLE_BOND){ /* the atom-i connecting to atom-j with double bond */
			if(isA[aij][Molecule::O] == 1) s2o++;	/* S=O */
			else if(isA[aij][Molecule::C] == 1) s2r++;	/* S=R */					
		}
	} /* end of if(isA[i][Molecule::S] == 1) */
}

/* ------------------------- carbon ----------------------- */
void Calculator3::calculate_C_results(short i){	
	if((cc == 1 && con_to_H[i] == 3) ||	/* AlH3 */
		(car == 1 && con_to_H[i] == 3) ||	/* ArH3 */
		(cc == 0 && con_to_H[i] == 4))		/* CH4 */
		results[1]++;	/* C-001: CH3R/CH4 */

	if((cc == 2 && bonds[i] == 2 && con_to_H[i] == 2) ||	/* Al2CH2 */
		(cc == 1 && car == 1 && bonds[i] == 2 && con_to_H[i] == 2) ||	/* C-CH2-Ar */
		(car == 2 && cc == 0 && bonds[i] == 2 && con_to_H[i] == 2))	/* Ar2CH2 */
		results[2]++;	/* C-002: CH2R2 */

	if(cc + car == 3 && bonds[i] == 3 && con_to_H[i] == 1) results[3]++;	/* C-003: CHR3 */
	if(cc + car == 4 && bonds[i] == 4 && con_to_H[i] == 0) results[4]++;	/* C-004: CR4 */
	if(cx == 1 && bonds[i] == 1 && con_to_H[i] == 3) results[5]++;	/* C-005: CH3X */
	if(cc + car == 1 && cx == 1 && bonds[i] == 2 && con_to_H[i] == 2) results[6]++;	/* C-006: CH2RX */
	if(cx == 2 && con_to_H[i] == 2) results[7]++;	/* C-007: CH2X2 */
	if(cc + car == 2 && cx == 1 && bonds[i] == 3 && con_to_H[i] == 1) results[8]++;	/* C-008: CHR2X */
	if(cc == 1 && cx == 2 && con_to_H[i] == 1) results[9]++;	/* C-009: CHRX2 */
	if(cx == 3 && con_to_H[i] == 1) results[10]++;	/* C-010: CHX3 */
	if(cc + car == 3 && cx == 1 && bonds[i] == 4 && con_to_H[i] == 0) results[11]++;	/* C-011: CR3X */
	if(cc == 2 && cx == 2 && con_to_H[i] == 0) results[12]++;	/* C-012: CR2X2 */
	if(cc == 1 && cx == 3 && con_to_H[i] == 0) results[13]++;	/* C-013: CRX3 */
	if(cx == 4 && con_to_H[i] == 0)	results[14]++;	/* C-014: CX4 */
	if(c2 == 1 && bonds[i] == 1 && con_to_H[i] == 2) results[15]++;	/* C-015: =CH2 */
	if(c2 == 1 && cc + car == 1 && bonds[i] == 2 && con_to_H[i] == 1) results[16]++;	/* C-016: =CHR */
	if(cc == 2 && c2 == 1 && bonds[i] == 3 && con_to_H[i] == 0)	results[17]++;	/* C-017: =CR2 */
	if(cc == 0 && cx == 1 && c2 == 1 && bonds[i] == 2 && con_to_H[i] == 1) results[18]++;	/* C-018: =CHX */
	if(cc + car == 1 && cx == 1 && c2 == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[19]++;	/* C-019: =CRX */
	if(cc == 0 && cx == 2 && c2 == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[20]++;	/* C-020: =CX2 */
	if(c3 == 1 && bonds[i] == 1 && con_to_H[i] == 1) results[21]++;
	if((cc + car == 1 && c3 == 1 && bonds[i] == 2 && con_to_H[i] == 0) ||
		(c2c == 2 && bonds[i] == 2 && con_to_H[i] == 0)) results[22]++;	/* C-022: #CR/R=C=R */
	if(cc == 0 && cx == 1 && c3 == 1 && bonds[i] == 2 && con_to_H[i] == 0) results[23]++;	/* C-023: #CX */
	if(aromatic[i] == 1 && con_to_H[i] == 1) results[24]++;	/* C-024: R--CH--R */
	if(aromatic[i] == 1 && arc == 1) results[25]++;	/* C-025: R--CR--R */
	if(aromatic[i] == 1 && arx == 1) results[26]++;	/* C-026: R--CX--R */
	if(cc == 1 && c2x == 1 && bonds[i] == 2 && con_to_H[i] == 1) results[36]++;	/* C-036: Al-CH=X */
	if(car == 1 && c2x == 1 && bonds[i] == 2 && con_to_H[i] == 1) results[37]++;	/* C-037: Ar-CH=X */
	if(cc == 2 && c2x == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[38]++;	/* C--038: Al-C(=X)-Al */
	if(car == 1 && c2x == 1 && cc == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[39]++;	/* C-039: Ar-C(=X)-R */
	if((cc == 1 && c2x == 1 && cx == 1 && bonds[i] == 3 && con_to_H[i] == 0) ||
		(cc == 1 && c3x == 1 && bonds[i] == 2 && con_to_H[i] == 0) ||
		(c2x == 2 && bonds[i] == 2 && con_to_H[i] == 0))
		results[40]++;	/* C-040: R-C(=X)-X / R-C#X / X=C=X */
	if(cx == 2 && c2x == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[41]++;	/* C-041: X-C(=X)-X */
}

/* ------------------------ Hydrogen ----------------------- */
void Calculator3::calculate_H_results(short i){
	if(isA[i][Molecule::C] != 1) results[50] += con_to_H[i];	/* H-051: H attached to heterotatom */
	if(sp3 == bonds[i]){
		if(m == 0) results[46]++;	/* H-046: H attached to C0(SP3) no X attached to next C */
		else if(m == 1) results[52]++;	/* H-052: H attached to C0(SP3) with 1X attached to next C */
		else if(m == 2) results[53]++;	/* H-053: H attached to C0(SP3) with 2X attached to next C */
		else if(m == 3) results[54]++;	/* H-054: H attached to C0(SP3) with 3X attached to next C */
		else if(m == 4) results[55]++;	/* H-055: H attached to C0(SP3) with 4X attached to next C */
	}
}

/* ------------------------ oxygen ----------------------- */
void Calculator3::calculate_O_results(short i){	
	if(o2 == 1 && bonds[i] == 1 && con_to_H[i] == 0) results[58]++;	/* O-058: O= */
	if(oal == 2 && bonds[i] == 2 && con_to_H[i] == 0) results[59]++;	/* O-059: Al-O-Al */
	if((oal == 1 && oar == 1 && o1 == 2 && bonds[i] == 2 && con_to_H[i] == 0) ||	/* Al-O-Ar */
		(oal == 0 && oar == 2 && o1 == 2 && bonds[i] == 2 && con_to_H[i] == 0) ||	/* Ar-O-Ar */
		(oc == 2 && al2x == 1 && o1 == 2 && bonds[i] == 2 && con_to_H[i] == 0))	/* R-O-C=X */
		results[60]++;	/* O-060: Al-O-Ar / Ar-O-Ar / R-O-C=X */
}

/* ------------------------ Nitrogen ----------------------- */
void Calculator3::calculate_N_results(short i){
	if(isA[i][Molecule::N] == 1){
		if(al == 1 && ar == 0 && bonds[i] == 1 && con_to_H[i] == 2) results[66]++;	/* N-066: Al-NH2 */
		else if(al == 2 && ar == 0 && bonds[i] == 2 && con_to_H[i] == 1) results[67]++;	/* N-067: Al2-NH1 */
		else if(al == 3 && ar == 0 && bonds[i] == 3 && con_to_H[i] == 0) results[68]++;	/* N-068: Al3-N */

		/* the atom-i(N) connected an aromatic with single bond */
		if((al == 0 && ar == 1 && bonds[i] == 1 && con_to_H[i] == 2) ||
			/* the atom-i(N) connected a atom except C and H with single bond */
			(al == 0 && ar == 0 && bonds[i] == 1 && con_to_H[i] == 2))
			results[69]++;	/* N-069: Ar-NH2 / X-NH2 */

		if(al == 1 && ar == 1 && bonds[i] == 2 && con_to_H[i] == 1) results[70]++;	/* N-070: Ar-NH-Al */
		if(al == 2 && ar == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[71]++;	/* N-071: Ar-N-Al2 */
		if((x2x >= 1 && bonds[i] == 3 && con_to_H[i] == 0 && used[i] != 1) ||	/* >N-X=X */
			(nc2o >= 1 && ncc >= 1 && bonds[i] == 3 && con_to_H[i] == 0 && used[i] != 1))	/* RCO-N< */
		{
			results[72]++;	/* N-072: RCO-N< / >N-X=X */
			used[i] = 1;
		}

		if((al == 0 && ar == 2 && bonds[i] == 2 && con_to_H[i] == 1) ||	/* Ar2NH */
			(al == 0 && ar == 3 && bonds[i] == 3 && con_to_H[i] == 0) ||	/* Ar3N */
			(al == 1 && ar == 2 && bonds[i] == 3 && con_to_H[i] == 0))		/* Ar2N-Al */
			results[73]++; /* N=073: Ar2NH / Ar3N / Ar2N-Al */
		short nbi = connBond[i].size();
		if(nbi >= 1){
			if(nbi == 1 && connBond[i][0] == Bond::TRIPLE_BOND && isA[connAtoms[i][0]][Molecule::C] == 1) results[74]++;	/* N-074: R#N */
			/* 1. the atom-i(N) has two bonds: one single bond and one double bond
			2. the double bond connected to a C
			3. the single bond connected an atom that except C and H */
			else if(nbi >= 2){
				if(nbi == 2 && ((connBond[i][0] == Bond::DOUBLE_BOND && connBond[i][1] == Bond::SINGLE_BOND && isA[connAtoms[i][0]][Molecule::C] == 1) ||
					(connBond[i][1] == Bond::DOUBLE_BOND && connBond[i][0] == Bond::SINGLE_BOND && isA[connAtoms[i][0]][Molecule::C] == 1))){
						results[74]++;	/* N-074: R=N- */
				}
			}
		}

		if((ar == 1 && o2 == 2 && bonds[i] == 3 && con_to_H[i] == 0) ||	/* Ar-NO2 */
			(ar == 1 && o2 == 1 && o1 == 1 && bonds[i] == 3 && con_to_H[i] == 0) ||	/* Ar-NO2 */
			(ro >= 1 && o2 == 2 && bonds[i] == 3 && con_to_H[i] == 0) ||	/* RO-NO2 */
			(ro >= 1 && o2 == 1 && o1 == 2 && bonds[i] == 3 && con_to_H[i] == 0))		/* RO-NO2 */
		{
			results[76]++;	/* N-076: Ar-NO2 / RO-NO2 */
		}

		if(al == 1 && o2 == 2 && bonds[i] == 3 && con_to_H[i] == 0) results[77]++;	/* N-077: Al-NO2 */
		if((ar == 1 && n2x == 1 && bonds[i] == 2 && con_to_H[i] == 0) ||	/* Ar-N=X */
			(nx == 1 && n2x == 1 && bonds[i] == 2 && con_to_H[i] == 0))	/* X-N=X */
		{
			results[78]++;	/* N-078: Ar-N=X / X-N=X */
		}
	} /* end of if(isA[i][Molecule::N] = 1) */
}

/* ------------------------ Sulfur ----------------------- */
void Calculator3::calculate_S_results(short i){
	if(isA[i][Molecule::S] == 1){
		if(sr == 1 && bonds[i] == 1 && con_to_H[i] == 1) results[106]++;	/* S-106: R-SH */
		if((sr == 1 && ssr == 1 && bonds[i] == 2 && con_to_H[i] == 0) ||	/* RS-SR */
			(sr == 2 && bonds[i] == 2 && con_to_H[i] == 0))	/* R2S */
		{
			results[107]++;	/* S-107: R2S / RS-SR */
		}

		if(s2r == 1 && bonds[i] == 1 && con_to_H[i] == 0) results[108]++;	/* S-108: R=S */
		if(sr == 2 && s2o == 1 && bonds[i] == 3 && con_to_H[i] == 0) results[109]++;	/* S-109: R-SO-R */
		if(sr == 2 && s2o == 2 && bonds[i] == 4 && con_to_H[i] == 0) results[110]++;	/* S-110: R-SO2-R */
	} /* end of if(isA[i][Molecule::S] == 1) */
}