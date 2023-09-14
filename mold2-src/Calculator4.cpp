#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Calculator4.h"

Calculator4::Calculator4(Molecule & _mol)
:mol(_mol),
connBond(mol.getConnectedBondTypes(false)),
connAtoms(mol.getConnectedAtomIDs(false)),
inRing(mol.getInRing()),
aromatic(mol.getIsAromatic()),
rings(mol.getRings()),
isA(mol.getIsA()),
con_to_H(mol.getCon_to_H()){
	calculateMLogP();
}

Calculator4::~Calculator4(){
}

void Calculator4::calculateMLogP(){
	try{
		if(mol.isFeaturesReady(::MlogP)) return;

		/* -------------------- nuns: count total number of unsaturated bonds ------------------------ */
		short nuns = countOfUnsaturatedBonds() / 2;
		vector<float> atomCount = atomicCount();

		/* ---------------------------------- fcx ------------------------------------- */
		double fcx = atomCount[Molecule::C] + atomCount[Molecule::F] + atomCount[Molecule::Cl] + atomCount[Molecule::Br] + atomCount[Molecule::I];
		double no_nn = atomCount[Molecule::O] + atomCount[Molecule::N];

		double fprx = 0.0;
		double fqn = 0;
		double fncs = 0;
		double famp = 0;
		short nno2 = 0;	

		stepOne(fprx, fqn, fncs, famp, nno2);
		short irng = countIrng();

		short ihb = 0;
		short w = 0;
		short npol = 0;
		short ibetal = 0;
		short ialk = 0;

		/* -------------------------------- ihb ------------------------------ */
		for(int i = rings.size() - 1; i >= 0 ;i--){
			ihbStepOne(i, w, ihb);

			/* -------------------------- npol ----------------------------- */
			calculateNpol(i, w, npol);

			/* -------------------------------- ihb: no-aromatic 5-membered ring or 6-membered ring ------------------------------ */
			checkIhb5_6ring(i, w, ihb);

			/* ---------------------------------- ibetal ------------------------------------- */
			calculateIbetal(i, ibetal);

			/* -------------------------- ialk ----------------------------- */
			calculateIalk(i, ialk);
		}

		/* ------------------------------- calculate the logP ------------------------------------------ */
		double logP = -1.041 + 1.244 * pow(fcx, 0.6) - 1.017 * pow(no_nn, 0.9) + 
			0.406 * fprx - 0.145 * pow(nuns, 0.8) + 0.511 * ihb + 0.268 * npol 
			- 2.215 * famp + 0.912 * ialk - 0.392 * irng - 3.684 * fqn
			+ 0.474 * nno2 + 1.582 * fncs + 0.773 * ibetal;
		mol.setFingerprint(logP, ::MlogP);
	}catch(exception &ex){
		cout << ex.what() << " thrown in Mlogp" << endl;
	}
}

int Calculator4::countOfUnsaturatedBonds(){
	short nuns = 0;
	for(short i = mol.getAtomsNum(false) - 1; i >= 0; i--){
		/* -------------------- nuns: count total number of unsaturated bonds ------------------------ */
		for(short j = connBond[i].size() - 1; j >= 0; j--){
			if(connBond[i][j] == Bond::AROMATIC_BOND){
				nuns++;
				break;
			}else if(connBond[i][j] >= Bond::DOUBLE_BOND){
				nuns++;
			}
		}
	}
	return nuns;
}

vector<float> Calculator4::atomicCount(){
	vector<float> atomCount(9, 0);
	vector<short> & elementCounts = mol.countElements();
	atomCount[Molecule::C] = elementCounts[Elements::CARBON];
	atomCount[Molecule::N] = elementCounts[Elements::NITROGEN];
	atomCount[Molecule::O] = elementCounts[Elements::OXYGEN];
	atomCount[Molecule::P] = elementCounts[Elements::PHOSPHORUS];
	atomCount[Molecule::S] = elementCounts[Elements::SULFUR];
	atomCount[Molecule::F] = elementCounts[Elements::FLOURINE] * 0.5f;
	atomCount[Molecule::Cl] = elementCounts[Elements::CLORINE];
	atomCount[Molecule::Br] = elementCounts[Elements::BROMINE] * 1.5f;
	atomCount[Molecule::I] = elementCounts[Elements::IODINE] * 2.0f;
	return atomCount;
}

void Calculator4::stepOne(double & fprx, double & fqn, double & fncs, double & famp, short & nno2){
	short noHnumAtoms = mol.getAtomsNum(false);
	for(short i = 0; i < noHnumAtoms; i++){
		short nai = connAtoms[i].size();
		if(isA[i][Molecule::N] == 1){ // N
			int x = 0;
			for(short j = nai - 1; j >= 0; j--){
				if(isA[connAtoms[i][j]][Molecule::N] == 1 || isA[connAtoms[i][j]][Molecule::O] == 1){
					if(connBond[i][j] == Bond::SINGLE_BOND) fprx += 1; // NN or NO
				}else if(isA[connAtoms[i][j]][Molecule::O] == 1){ // structure of N-oxide
					fqn += 0.5;	/* it is a N-oxide structure */
					x++;
					break;
				}
			}

			if(x == 0){ // look for quaternary N
				int y = 0;
				for(short j = nai - 1; j >= 0; j--){
					if(connBond[i][j] == Bond::SINGLE_BOND) y++;
				}
				if(y >= 4 && con_to_H[i] == 0) fqn += 1.0; // found the quaternary N
			}

			/* ---------------------------- structure of nitro groups ------------------------ */
			x = 0;	/* count number of -O */
			int y = 0;	/* count number of =O */
			for(short j = nai - 1; j >= 0; j--){
				if(nai == 3){ /* the N has three atoms connection */
					if(isA[connAtoms[i][j]][Molecule::O] == 1){
						if(connBond[i][j] == Bond::SINGLE_BOND) x++;
						else if(connBond[i][j] == Bond::DOUBLE_BOND) y++;
					}
				}
			}
			if((x >= 1 && y == 1) || (y == 2)){ /* -NO2 */
				nno2++;	/* count the number of structure of -NO2 */
			}
		}else if(isA[i][Molecule::O] == 1){ // O
			/* ------------------------------------ fprx ----------------------------------- */
			for(short j = 0; j < nai; j++){
				if(isA[connAtoms[i][j]][Molecule::N] == 1 || isA[connAtoms[i][j]][Molecule::O] == 1){
					if(connBond[i][j] == Bond::SINGLE_BOND) fprx += 1;	/* structure of OO or ON */
				}
			}
		}else if(isA[i][Molecule::C] == 1){  // C
			/* ------------------------------------ fprx ----------------------------------- */
			int h, r, s, w, x, y, z;
			h = r = s = w = x = y = z = 0;
			for(short j = 0; j < nai; j++){
				if(connBond[i][j] == Bond::SINGLE_BOND){ /* found a single bond the belong the C */
					if(isA[connAtoms[i][j]][Molecule::C] == 1) r++;	/* C-C */
					else if(isA[connAtoms[i][j]][Molecule::O] == 1) w++;	/* O-C */
					else if(isA[connAtoms[i][j]][Molecule::N] == 1){ /* the C connection a N with the single bond */
						x++;	/* N-C */
						if(con_to_H[connAtoms[i][j]] == 2) h++;
					}				
				}else if(connBond[i][j] == Bond::DOUBLE_BOND){	/* found a double bond the belong the C */
					if(isA[connAtoms[i][j]][Molecule::S] == 1) s++;	/* S=C */
					else if(isA[connAtoms[i][j]][Molecule::O] == 1) y++;	/* O=C */
					else if(isA[connAtoms[i][j]][Molecule::N] == 1) z++;	/* N=C */
				}
			}
			switch(w + x + y + z){
				case 2:
					fprx += 1;
					break;
				case 3:
					fprx += 3;
					break;
				case 4:
					fprx += 6;
			}

			if(r == 1 && s == 1 && x == 1 && h != 0 && w == 0 && y == 0 && z == 0){
				fprx -= 1;	/* R-C(=S)-N */
			}else if(r == 1 && y == 1 && x == 1 && h != 0 && w == 0 && s == 0 && z == 0){
				fprx -= 1;	/* R-C(=O)-N */
			}

			/* ----------------------------------- fncs: -N=C=S group ---------------------------- */
			if(nai == 2){ /* the atom-i has two bonds */
				x = 0;	/* count the sturcture of -N= */
				y = 0;	/* count the sturcture of =S */				
				for(short j = 0; j < nai; j++){
					if(connBond[i][j] == Bond::DOUBLE_BOND){ /* the atom-i connecting atom-j with double bond */
						if(isA[connAtoms[i][j]][Molecule::N] == 1){
							if((short)connBond[connAtoms[i][j]].size() == 2){
								z = 0;
								for(short k = connBond[connAtoms[i][j]].size() - 1; k >= 0; k--) z += connBond[connAtoms[i][j]][k];
								if(z == 3) x++; /* the sturcture is -N= */
							}
						}else if(isA[connAtoms[i][j]][Molecule::S] == 1) y++;
					}
				}
				if(x == 1 && y == 1) fncs += 1.0;

				/* --------------------------------- fncs: -S-CN ---------------------------------- */
				x = 0;	/* count number of N connecting to C */
				for(short j = 0; j < nai; j++){
					if(connBond[i][j] == Bond::SINGLE_BOND){ /* single bond between atom-i and atom-j */
						if(isA[connAtoms[i][j]][Molecule::S] == 1){
							y = 0;	/* count number of single bonds of the S */
							for(short k = connBond[connAtoms[i][j]].size() - 1; k >= 0; k--){
								if(connBond[connAtoms[i][j]][k] == Bond::SINGLE_BOND) y++; // single bond */
							}
						}
					}
					if(isA[connAtoms[i][j]][Molecule::N] == 1) x++; // found a N connected on C
				}
				if(x == 1 && y == 2) fncs += 0.5; // -S-CN
			} /* end of if(connection_bonds[i][0].num_bonds == 2) */
		}else if(isA[i][Molecule::S] == 1){
			int h, w, x, y, z;
			h = w = x = y = z = 0;
			for(short j = 0; j < nai; j++){
				if(connBond[i][j] == Bond::SINGLE_BOND){	/* found a single bond the belong the S */
					if(isA[connAtoms[i][j]][Molecule::C] == 1) w++; // C-S
					else if(isA[connAtoms[i][j]][Molecule::N] == 1){ /* the S connection a N with the single bond */
						x++;	// N-S
						if(con_to_H[connAtoms[i][j]] == 2) h++;
					}else if(isA[connAtoms[i][j]][Molecule::O] == 1) z++; // O-S
				}else if(connBond[i][j] == Bond::DOUBLE_BOND &&	isA[connAtoms[i][j]][Molecule::S] == 1) y++; // O=S
			}
			if(nai == 4 && w == 1 && x == 1 && y == 2 && h != 0 && z == 0){
				fprx += -1;	/* -S(=O2)-N */
			}else if(x + y + z == 2){
				fprx += 1;
			}else if(x + y + z == 3){
				fprx += 3;
			}else if(x + y + z == 4){
				fprx += 6;
			}
		}else if(isA[i][Molecule::P] == 1){
			int x, y, z;
			x = y = z = 0;
			for(short j = 0; j < nai; j++){
				if(isA[connAtoms[i][j]][Molecule::N] == 1) x++;	/* N-P */
				else if(isA[connAtoms[i][j]][Molecule::O] == 1){ /* the P connected an O */
					if( connBond[i][j] == Bond::SINGLE_BOND) y++;	/* O-P */
				}else if(connBond[i][j] == Bond::DOUBLE_BOND) z++;	/* O=P */
			}
			if(x + y + z == 2){
				fprx += 1;	/* P between N or O */
			}else if(x + y + z == 3){
				fprx += 3;
			}else if(x + y + z == 4){
				fprx += 6;
			}
			if(z != 0 && y != 0){ /* the P has one or more structure as: O-P=O */
				if(y >= z){
					fprx += z;
				}else{
					fprx += y;
				}
			}
		}

		if(isA[i][Molecule::C] == 1 && aromatic[i] != 1){ /* the atom-i is not an aromatic C */
			/* ------------------------------ fprx ------------------------------ */
			if(nai == 3){
				int x, y, z;
				x = y = z = 0;
				for(short j = 0; j < nai; j++){
					if(isA[connAtoms[i][j]][Molecule::O] == 1){
						if(connBond[i][j] == Bond::SINGLE_BOND && inRing[connAtoms[i][j]] != 1) x++; //the atom-j is not a ring atom
						else if(connBond[i][j] == Bond::DOUBLE_BOND) y++; // double bond betwee atom-i and atom-j
					}
				}
				if(x + y == 2) fprx++;
			}

			/* ------------------------------ famp ------------------------------ */
			if((con_to_H[i] >= 1 &&	nai == 2) || (nai == 3)){	/* the C has three bonds */
				int w, x, y, z;	/* count number of N connection on alpha-C */
				w = x = y = z = 0;		/* count number of C connection on alpha-C */
				for(short j = 0; j < nai; j++){
					if(connBond[i][j] == Bond::SINGLE_BOND){	/* atom-i connecting to atom-j with single bond */
						if(isA[connAtoms[i][j]][Molecule::N] == 1 && con_to_H[connAtoms[i][j]] >= 2) w++; /* count number of N connection on alpha-C */
						else if(isA[connAtoms[i][j]][Molecule::C] == 1){	/* the atom-j is a C */
							x++;	/* count number of C connection on alpha-C */
							short naij = connAtoms[connAtoms[i][j]].size();
							if( naij == 3){
								for(short k = 0; k < naij; k++){
									if(isA[connAtoms[connAtoms[i][j]][k]][Molecule::O] == 1){
										if(connBond[connAtoms[i][j]][k] == Bond::DOUBLE_BOND) y++; // O=C
										else if(connBond[connAtoms[i][j]][k] == Bond::SINGLE_BOND && con_to_H[connAtoms[connAtoms[i][j]][k]] == 1) z++;	// count number of OH connection on atom-j
									}
								}
							}
						}
					}
				}
				if(nai == 2){		/* the a-C has two bonds */
					if(w == 1 && x == 1 && y == 1 && z == 1) famp += 1.0;
				}else if(nai == 3){	/* the a-C has three bonds */
					if(w == 1 && x == 2 && y == 1 && z == 1) famp += 1.0;	/* alpha-amino */
				}
			}
		}else if(aromatic[i] == 1){ /* ----------- 2-aminobenzoic acid or 3-aminobenzoic acid or 4-aminobenzoic acid ----------- */
			int x, y, z; //TODO: need furthur confirmation
			z = 0;	/* count number of -NH2 */
			x = 0;	/* count number of =O */ 
			y = 0;	/* count number of -OH */					
			for(short j = 0; j < nai; j++){				
				if(isA[connAtoms[i][j]][Molecule::C] == 1 && aromatic[connAtoms[i][j]] != 1 && connBond[connAtoms[i][j]].size() == 3){
					for(short k = 2; k >= 0; k--){
						if(connBond[connAtoms[i][j]][k] == Bond::DOUBLE_BOND && isA[connAtoms[connAtoms[i][j]][k]][Molecule::O] == 1) x++;
						else if(connBond[connAtoms[i][j]][k] == Bond::SINGLE_BOND && isA[connAtoms[connAtoms[i][j]][k]][Molecule::O] == 1 && con_to_H[connAtoms[connAtoms[i][j]][k]] == 1) y++;
					}
				}else if(isA[connAtoms[i][j]][Molecule::N] == 1 && con_to_H[connAtoms[i][j]] == 2 && 
					connBond[connAtoms[i][j]].size() == 1 && connBond[i][j] == Bond::SINGLE_BOND) z++;
			}
			if(x >= 1 && y >= 1 && z>= 1){ /* found an aminobenzoic acid structure */
				famp += 0.5;	/* aminobenzoic acid */
			}
		} /* end of else if(aromatic[i] == YES) */

		/* ---------------------------------- pyridinecarboxylic acid ------------------------------------- */
		if(isA[i][Molecule::N] == 1 && connBond[i].size() >= 2 && inRing[i] == 1){
			for(short j = rings.size() - 1; j >= 0; j--){  //TODO original code is weird
				int w, x, y, z, r, s;
				w = x = y = z = r = s = 0;
				if(rings[j].size() == 6){	/* the i-ring must be a 6-membered ring */
					for(short k = 5; k >= 0; k--){
						if(rings[j][k] == i){
							for(short m = connAtoms[i].size() - 1; m >= 0; m--){
								if(connBond[i][m] == Bond::DOUBLE_BOND){
									for(short n = 5; n >= 0 ; n--){
										if(connAtoms[i][m] == rings[j][n] && isA[connAtoms[i][m]][Molecule::C] == 1) w++;	// unmber of C=N
									}
								}else if(connBond[i][m] == Bond::SINGLE_BOND){
									for(short n = 5; n >= 0 ; n--){
										if(connAtoms[i][m] == rings[j][n] && isA[connAtoms[i][m]][Molecule::C] == 1) x++;	// number of C-N
									}
								}
							}
						}else if(isA[rings[j][k]][Molecule::C] == 1){
							for(short m = connAtoms[rings[j][k]].size() - 1; m >= 0; m--){
								if(connBond[rings[j][k]][m] == Bond::DOUBLE_BOND){
									for(short n = 5; n >= 0; n--){
										if(connAtoms[rings[j][k]][m] == rings[j][n]) y++;	// number of =C
									}
								}else if(connBond[rings[j][k]][m] == Bond::SINGLE_BOND){
									for(short n = 5; n >= 0 ; n--){
										if(connAtoms[rings[j][k]][m] == rings[j][n]) z++; // number of -C
									}
								}
								if(connAtoms[rings[j][k]].size() == 3){	/* the ring's atom has three connection */
									if(isA[connAtoms[rings[j][k]][m]][Molecule::C] == 1 && inRing[connAtoms[rings[j][k]][m]] != 1){
										for(short n = connAtoms[connAtoms[rings[j][k]][m]].size() - 1; n >= 0; n--){
											if(isA[connAtoms[connAtoms[rings[j][k]][m]][n]][Molecule::O] == 1){
												if( connBond[connAtoms[rings[j][k]][m]][n] == Bond::DOUBLE_BOND) r++; // count number of =O */
												else if(connBond[connAtoms[rings[j][k]][m]][n] == Bond::DOUBLE_BOND && con_to_H[connAtoms[connAtoms[rings[j][k]][m]][n]] == 1) s++;	/* count number of OH */
											}
										}
									}
								}
							}
						}
					}
				}			
				if(w == 1 && x == 1 && y == 3 && z == 3 && r >= 1 && s >= 1){
					famp += 0.5;	/* pyridinecarboxylic acid */
				}
			}
		} /* end of if(nitrogen[i] == YES && ring[i] == YES) */
	} /* end of for(i = 0; i <= connection_table[i][0].num_bonds; i++) */
}

short Calculator4::countIrng(void){ //the book say: it is presence.
	for(int i = rings.size() - 1; i >= 0; i--){
		for(int j = rings[i].size() - 1; j >= 0; j--){
			if(aromatic[rings[i][j]] != 1){/* look for an aromatic ring and count it */
				return 1;	/* count number of no aromatic ring */
			}
		}
	}
	return 0;
}

/* -------------------------- ialk ----------------------------- */
void Calculator4::calculateIalk(short i, short & ialk){	
	int x = 0;	/* count number of C on the ring-i */
	int y = 0;	/* count the double bond of the ring-i */
	short nri = rings[i].size();
	for(short j = 0; j < nri; j++){
		if(isA[rings[i][j]][Molecule::C] == 1){ /* atom-j of ring-i is a C */
			x++;	/* count number of C on the ring-i */
		}
		for(short k = connBond[rings[i][j]].size() - 1; k >= 0; k--){ /* running over all connection atoms of atom-j */
			if(connBond[rings[i][j]][k] == Bond::DOUBLE_BOND){ /* there is a double bond between atom-j and aotm-i */
				for(short m = 0; m < nri; m++){	/* check the atom-k that belong to ring-i or no */
					if(rings[i][m] == connAtoms[rings[i][j]][k]){ /* the atom-k is belong to ring-i */
						y++;	/* count the double bond of the ring-i */
					}
				}
			}
		}
	}
	if(x == nri && y > 0) ialk = 1;
}

/* ---------------------------------- ibetal ------------------------------------- */
void Calculator4::calculateIbetal(short i, short & ibetal){			
	int nri = rings[i].size();
	if(nri == 4){ /* the ring is a 4-membered ring */
		int z = 0;
		for(short j = 3; j >= 0; j--){	/* search all atoms of the i-ring */
			if(isA[rings[i][j]][Molecule::C] == 1){ /* the i-ring's atom-j is a C -- try to find the Beta-C */
				if(connBond[rings[i][j]].size() == 3){ /* the C has three bonds */
					for(short k = 2; k >= 0; k--){ /* search all bonds of the C */
						if(connBond[rings[i][j]][k] == Bond::DOUBLE_BOND){ /* found a double bond from the C */
							if(isA[connAtoms[rings[i][j]][k]][Molecule::O] == 1){ /* found a O=C */
								for(short m = connAtoms[rings[i][j]].size() - 1; m >= 0; m--){
									if(isA[connAtoms[rings[i][j]][m]][Molecule::C] == 1 && connBond[rings[i][j]][m] == Bond::SINGLE_BOND){
										for(short n = 0; n < nri; n++){
											if(rings[i][n] == connAtoms[rings[i][j]][m]){ /* found alpha-C */
												for(short r = connAtoms[rings[i][n]].size() - 1; r >= 0; r--){ /* running over all bonds of alpha-C */
													if(isA[connAtoms[rings[i][n]][r]][Molecule::C] == 1 && connBond[rings[i][n]][r] == Bond::SINGLE_BOND){
														for(short s = 0; s < nri; s++){
															if(rings[i][s] == connAtoms[rings[i][n]][r] && connAtoms[rings[i][n]][r] != rings[i][j]){
																for(short w = connAtoms[rings[i][s]].size() - 1; w >= 0; w--){
																	if(isA[connAtoms[rings[i][s]][w]][Molecule::N] == 1 && connBond[rings[i][s]][w] == Bond::SINGLE_BOND){
																		for(short x = 0; x < nri; x++){
																			if(rings[i][x] == connAtoms[rings[i][s]][w]){
																				for(short y = connAtoms[rings[i][x]].size() - 1; y >= 0; y--){
																					if(connAtoms[rings[i][x]][y] == rings[i][j] && connBond[rings[i][x]][y] == Bond::SINGLE_BOND){
																						ibetal = 1;
																						z = 1;
																						break;
																					}
																				}
																			}
																			if(z == 1) break;
																		}
																	}
																	if(z == 1) break;
																}
															}
															if(z == 1) break;
														}
													}
													if(z == 1) break;
												}
											}
											if(z == 1) break;
										}
									}
									if(z == 1) break;
								}
							}
						}
						if(z == 1) break;
					} /* end of FOR-k */
				}
			}
			if(z == 1) break;
		}
	} /* end of if(my_ring[i][0] == 5) */
}

/* -------------------------------- ihb: no-aromatic 5-membered ring or 6-membered ring ------------------------------ */
void Calculator4::checkIhb5_6ring(short i, short &w, short & ihb){
	short nri = rings[i].size();
	if(w != nri){
		if(nri == 5 || nri == 6){	/* only compare 5-ring and 6-ring */
			w = 0;
			for(short j = nri - 1; j >= 0; j--){
				if(isA[rings[i][j]][Molecule::C] == 1) w++;
			}
			if(w == nri){ /* all atoms of the ring are C */
				for(short j = nri - 1; j >= 0; j--){
					if(connAtoms[rings[i][j]].size() == 3){
						int x = 0;
						/* the intramolecular H-bonds are base on SP2 struecure if there are no on aromatic ring */
						for(short k = connBond[rings[i][j]].size() - 1; k >= 0; k--){
							x += connBond[rings[i][j]][k];
						}
						if(x == 4){	/* sp2: one double bond and two single bonds */
							for(short k = connAtoms[rings[i][j]].size() - 1; k >= 0; k--){
								if(isA[connAtoms[rings[i][j]][k]][Molecule::O] == 1 ||	/* the atom-k is an O */
									isA[connAtoms[rings[i][j]][k]][Molecule::N] == 1){	/* or the atom-k is a N */
										if(con_to_H[connAtoms[rings[i][j]][k]] >= 1){		/* the atom-k has one or more H */
											int y = 0;
											for(short m = 0; m < nri; m++){
												if(rings[i][m] == connAtoms[rings[i][j]][k]){ /* atom-k belong to ring-i */
													y = 1;	/* flag for found the atom is belong to ring-i */
													break;	/* continue search ring-i's atom */
												}
											}
											if(y == 0){	/* the atom-k is not belong to ring-i */
												for(short m = connAtoms[rings[i][j]].size() - 1; m >= 0; m--){ /* search all connecting atoms of my_ring[i][j] */
													if(connAtoms[connAtoms[rings[i][j]][m]].size() == 3 ||	/* found a three connection atom */
														connAtoms[rings[i][j]][m] != connAtoms[rings[i][j]][k]){	/* the atom-m is not the atom-k */
															for(short n = connAtoms[connAtoms[rings[i][j]][m]].size() - 1; n >= 0; n--){ /* running over all connection atoms of atom-m */
																if(isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::O] == 1 ||	/* the atom-m is an O */
																	isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::N] == 1){	/* or the atom-m is a N */
																		if(con_to_H[connAtoms[connAtoms[rings[i][j]][m]][n]] == 0){
																			int z = 0;
																			for(short r = 0; r < nri; r++){
																				if(rings[i][r] == connAtoms[connAtoms[rings[i][j]][m]][n]){ /* the atom-n belong to ring-i */
																					z++;
																					break;
																				}
																			}
																			if(z == 0) /* the atom-m does not belong to ring-i */
																				ihb = 1;
																		}
																}
															}
													}
												}
											}
										} /* end of atom-k has H connection */
										/* the atom-k does not has H connection */
										else if(con_to_H[connAtoms[rings[i][j]][k]] == 0){		/* the atom-k has one or more H */
											short y = 0;
											for(short m = 0; m < nri; m++){
												if(rings[i][m] == connAtoms[rings[i][j]][k]){ /* atom-k belong to ring-i */
													y = 1;	/* flag for found the atom is belong to ring-i */
													break;
												}
											}
											if(y == 0){	/* the atom-k is not belong to ring-i */
												for(short m = connAtoms[rings[i][j]].size() - 1; m >= 0; m--){ /* search all connecting atoms of my_ring[i][j] */
													if(connAtoms[connAtoms[rings[i][j]][m]].size() == 3 ||	/* found a three connection atom */
														connAtoms[rings[i][j]][m] != connAtoms[rings[i][j]][k]){	/* the atom-m is not the atom-k */
															for(short n = connAtoms[connAtoms[rings[i][j]][m]].size() - 1; n >= 0; n--){
																if((isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::O] == 1) ||	/* the atom-m is an O */
																	(isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::N] == 1 &&	/* or the atom-m is a N */
																	con_to_H[connAtoms[connAtoms[rings[i][j]][m]][n]] >= 1)){		/* found the intramolecular H-bond */
																		int z = 0;
																		for(short r = 0; r < nri; r++){
																			if(rings[i][r] == connAtoms[connAtoms[rings[i][j]][m]][n]){ /* the atom-n belong to ring-i */
																				z++;
																				break;
																			}
																		}
																		if(z == 0) /* the atom-m does not belong to ring-i */
																			ihb = 1;
																}
															}
													}
												}
											}
										}
								}
							} /* end of FOR-k */
						}
					} 
				} /* for(j = 1; j < my_ring[i][0]; j++) */
			} //end of IF(all atoms are C of the ring) */
		} /* end of if(5-ring || 6-ring) */
	}
}

/* -------------------------- npol ----------------------------- */
void Calculator4::calculateNpol(short i, short & w, short & npol){
	w = 0;
	short nri = rings[i].size();
	for(short j = 0; j < nri; j++){
		if(aromatic[rings[i][j]] == 1) w++;
	}
	if(w == nri){	/* the ring is an aromatic ring */
		for(short j = 0; j < nri; j++){
			if(connAtoms[rings[i][j]].size() == 3){
				int x, y, z;
				x = y = z = 0;
				for(short k = 2; k >= 0; k--){
					if(isA[connAtoms[rings[i][j]][k]][Molecule::C] != 1){ /* found structure of POL */
						npol++;	/* Ar-X (X: any atom except C nad H) */
						break;	/* jump out from FOR-k -- look for other POL */
					}else{
						if(aromatic[connAtoms[rings[i][j]][k]] != 1){ /* Ar-C */
							for(short m = connAtoms[connAtoms[rings[i][j]][k]].size() - 1; m >= 0; m--){
								if(isA[connAtoms[connAtoms[rings[i][j]][k]][m]][Molecule::C] != 1){ /* found structure of POL */
									if(isA[connAtoms[connAtoms[rings[i][j]][k]][m]][Molecule::O] == 1){
										y++;	/* count number of O */
										if(inRing[connAtoms[connAtoms[rings[i][j]][k]][m]] == 1){ /* the O is member of a ring */
											z++;	/* count number of O that belong a ring */
										}
									}
									x++;	/* flag */
									npol++;	/* Ar-C-X or Ar-C=X or Ar-C#X (: any atom except C or H) */
								}else{ /* the atom-m is a C */
									if(connBond[connAtoms[rings[i][j]][k]][m] == Bond::DOUBLE_BOND){ /* the atom-m connected with double bond */
										if(inRing[connAtoms[connAtoms[rings[i][j]][k]][m]] != 1){  
											x++;
											npol++;	/* Ar-C=C */
										}
									}
								}
							}
						}
					}
					if(x != 0)	/* found the POL structure that connecting to atom-j */
						break;	/* jump out from FOR-k -- for saving the running time and looking for other POL */
				}
				if(y == 2 && z == 0){	/* Ar-COOH */
					npol--;
				}
			} /* end of if(connection_table[my_ring[i][j]][0].num_atoms == 3) */
		} /* end	npol--; of for(j = 1; j <= my_ring[i][0]; j++) */
	} /* end of }all atoms of the ring are aromatic atoms */
}

void Calculator4::ihbStepOne(short i, short & w, short & ihb){
	short nri = rings[i].size();
	for(short j = 0; j < nri; j++){
		if(aromatic[rings[i][j]] == 1){	/* the ring is an aromatic ring */
			if(connAtoms[rings[i][j]].size() == 3){
				for(short k = 2; k >= 0; k--){
					if(isA[connAtoms[rings[i][j]][k]][Molecule::O] == 1){	/* the atom-k is an O */
						if(con_to_H[connAtoms[rings[i][j]][k]] == 1){	/* the atom-k an O and with OH structure */
							for(short m = 2; m >= 0; m--){	/* look for intramolecular H-bonds */
								//TODO: was changed
								if(connAtoms[connAtoms[rings[i][j]][m]].size() == 3 &&	/* found an atom has three connection */
									connAtoms[rings[i][j]][k] != connAtoms[rings[i][j]][m]){		/* the atom-m is a new structure */
										for(short n = connAtoms[connAtoms[rings[i][j]][m]].size() - 1; n >= 0; n--){
											if(isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::N] == 1){	/* the atom-n is a N */
												if(con_to_H[connAtoms[connAtoms[rings[i][j]][m]][n]] == 2){	/* -NH2 */
													ihb = 1;	/* found intramolecular H-bond: OH--NH2 */
												}
											}else if(isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::C] == 1){	/* the atom-n is a C */
												int x, y, z;
												x = y = z = 0;
												for(short r = 0; r < nri; r++){
													if(rings[i][r] == connAtoms[connAtoms[rings[i][j]][m]][n]){ /* the C belong to ring-i */
														x++;
													}
												}
												if(x == 0){
													for(short r = connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]].size() - 1; r >= 0; r--){
														if(isA[connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]][r]][Molecule::O] == 1 && connBond[connAtoms[connAtoms[rings[i][j]][m]][n]][r] == Bond::DOUBLE_BOND){
															y++;
														}else if(isA[connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]][r]][Molecule::C] == 1 && connBond[connAtoms[connAtoms[rings[i][j]][m]][n]][r] == Bond::SINGLE_BOND){
															z++;
														}
													}
													if(y == 1 && z == 2) ihb = 1; /* C(=O)-R2,  found intramolecular H-bond: OH--(-CO-R) */
												}
											}
										}
								}
							}
						}
					} /* end of if(oxygen[connection_table[my_ring[i][j]][k].con_atoms] == YES) */
					else if(isA[connAtoms[rings[i][j]][k]][Molecule::N] == 1){	/* the atom-k is a N */
						if(con_to_H[connAtoms[rings[i][j]][k]] == 2){	/* the atom-k has -NH2 structure */
							for(short m = 2; m >= 0; m--){	/* look for intramolecular H-bonds */
								if(connAtoms[connAtoms[rings[i][j]][m]].size() == 3 &&	/* found an atom has three connection */
									connAtoms[rings[i][j]][k] != connAtoms[rings[i][j]][m]){		/* the atom-m is a new structure */
										for(short n = connAtoms[connAtoms[rings[i][j]][m]].size() - 1; n >= 0; n--){
											if(isA[connAtoms[connAtoms[rings[i][j]][m]][n]][Molecule::C] == 1){
												int x, y, z;
												x = y = z = 0;
												for(short r = 0; r <nri; r++){
													if(rings[i][r] == connAtoms[connAtoms[rings[i][j]][m]][n]) x++;
												}
												if(x == 0){
													for(short r = connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]].size() - 1; r >= 0; r--){
														if(isA[connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]][r]][Molecule::O] == 1){
															if(connBond[connAtoms[connAtoms[rings[i][j]][m]][n]][r] == Bond::SINGLE_BOND)
																if(con_to_H[connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]][r]] == 1) w++;	/* -OH */
															else if(connBond[connAtoms[connAtoms[rings[i][j]][m]][n]][r] == Bond::DOUBLE_BOND) y++;	/* =O */
														}else if(isA[connAtoms[connAtoms[connAtoms[rings[i][j]][m]][n]][r]][Molecule::C] == 1 && connBond[connAtoms[connAtoms[rings[i][j]][m]][n]][r] == Bond::SINGLE_BOND) z++;	/* -C */
													}
													if(w == 1 && y == 1 && z == 1){ /* -C(=O)-R2 */
														ihb = 1;	/* found intramolecular H-bond: NH2--(-COOH) */
													}
												}
											}
										}
								}
							}
						}
					} /* end of else if(nitrogen[connection_table[my_ring[i][j]][k].con_atoms] == YES) */
				} /* end of FOR-k */
			}
		} /* end of if(aromatic[my_ring[i][j]] == YES) */
	}
}
