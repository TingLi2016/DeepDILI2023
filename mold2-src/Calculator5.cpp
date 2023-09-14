#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Elements.h"
#include "Calculator5.h"
#include "Constants.h"

Calculator5::Calculator5(Molecule & _mol)
:mol(_mol){
	numAtoms = mol.getAtomsNum(false);
	doAll();
}

Calculator5::~Calculator5(void){
}

void Calculator5::SumValenceVertex1(){
	try{
		if(mol.isFeaturesReady(::SumValenceVertex1)) return;	
		vector<short> & atomic_num = mol.getAtomicNum(false);
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		double svvc1i = 0.0;
		for(short i = 0; i < numAtoms; i++){	/* running over all atoms */
			short nai = connAtoms[i].size();
			for(short j = nai - 1; j >= 0; j--){	/* running over all first order atoms */
				svvc1i += atomic_num[i] / sqrt(double(nai) * 
					connAtoms[connAtoms[i][j]].size());
			}
		}
		svvc1i /= 2;
		mol.setFingerprint(svvc1i, ::SumValenceVertex1);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SumValenceVertex1" << endl;
	}
}

void Calculator5::ReciprocalDistanceIndex(){
	try{
		if(mol.isFeaturesReady(::ReciprocalDistanceIndex)) return;

		vector<double> rds(numAtoms, 0);
		DistanceMatrix<short> & distanceMatrix = mol.getDistanceMatrix(false);

		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				if(i != j){
					short dis = distanceMatrix.getValueInCell(i, j);
					if(dis > 0)	rds[i] += 1.0 / dis;	/* reciprocal distance sums of row_i */
				}
			}
		}

		double rd1spi = 0.0;
		double srd1spi = 0.0;
		short numBonds = mol.getBondsNum(false);
		for(short i = 0; i < numBonds; i++){
			Bond & bond = mol.getBond(i);
			short a = bond.getFirstAtomNumber();
			short b = bond.getSecondAtomNumber();

			if(sqrt(rds[a] * rds[b]) > 0)
				rd1spi += 1.0 / sqrt(rds[a] * rds[b]);	/* reciprocal distance order-1 sum product index: RD1SPI */

			srd1spi += sqrt(rds[a] * rds[b]);			/* squared reciprocal distance order-1 sum product index: SRD1SPI */
		}
		vector<double> fingerprint(2, 0);
		fingerprint[0] = rd1spi;
		fingerprint[1] = srd1spi;
		mol.setFingerprint(fingerprint, ::ReciprocalDistanceIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in ReciprocalDistanceIndex" << endl;
	}
}

void Calculator5::DistanceInfo(void){ //TODO: this function is weird
	try{
		if(mol.isFeaturesReady(::DistanceInfo)) return;

		short numAtoms = mol.getAtomsNum(false);
		/*------------ calculate sigma of atom-i, slp, atom's connectivity index, and Wiener index ---------- */
		DistanceMatrix<short> & distanceMatrix = mol.getDistanceMatrix(false);
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector<short> &eta = distanceMatrix.getEta();
		int max_e = sum_element(eta);

		int cilp = 0;
		for(short i = 0; i < numAtoms; i++){		
			cilp += eta[i] * connAtoms[i].size();	/* atom's connectivity index in longest path */
		}
		double average_e = double(max_e) / numAtoms;	/* average longest path of the molecule */

		/*------------------------- calculate eccentric ------------------------ */
		double adalp = 0.0;
		for(short i = 0; i < numAtoms; i++){
			adalp += fabs(eta[i] - average_e);
		}
		adalp /= numAtoms;	/* eccentric */

		/*---- calculate average distance degree and average of deviation of distance degree ----- */

		/* calculate average distance degree */
		int wi = mol.getDistanceMatrix(false).sumAllCells();	
		double add = double(wi) / numAtoms;

		/* average of deviation of distance degree */
		double addd = 0.0;
		vector<int> & sigma = mol.getDistanceMatrix(false).getSigma();
		for(short i = 0; i < numAtoms; i++){
			addd += fabs(sigma[i] - add);
		}
		addd /= numAtoms;

		/*---------------------------- shortest path in the molecule ----------------------------- */
		int spm = *min_element(sigma.begin(), sigma.end()); 

		/* --------------------------- shortest path centralization index ------------------------- */
		int spci = wi - numAtoms * spm;

		/*---------------------------- variation ----------------------------- */
		int mv = 1;
		for(short i = 0; i < numAtoms; i++){
			if(abs(sigma[i] - spm) > mv)
				mv = abs(sigma[i] - spm);	/* pick the biggest value of variation up */
		}

		vector<double> fingerprint(8, Constants::ERROR_SIGNAL);
		fingerprint[0] = cilp;
		fingerprint[1] = max_e;
		fingerprint[2] = average_e;
		fingerprint[3] = adalp;
		fingerprint[4] = addd;
		fingerprint[5] = spm;
		fingerprint[6] = spci;
		fingerprint[7] = mv;
		mol.setFingerprint(fingerprint, ::DistanceInfo);
	}catch(exception &ex){
		cout << ex.what() << " thrown in DistanceInfo" << endl;
	}
}

void Calculator5::PetitjeanIndex(void){	
	try{
		if(mol.isFeaturesReady(::PetitjeanIndex)) return;

		vector<short> & eta = mol.getDistanceMatrix(false).getEta();	
		short r = *min_element(eta.begin(), eta.end());
		short d = *max_element(eta.begin(), eta.end());	
		double pi = (d - r);	/* calculate the Petitjean index: PI */
		if(r > 0) pi /= r;
		mol.setFingerprint(pi, ::PetitjeanIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in PetitjeanIndex" << endl;
	}
}

void Calculator5::CentricIndex(void){
	try{
		if(mol.isFeaturesReady(::CentricIndex)) return;

		short numAtoms = mol.getAtomsNum(false);
		vector< set<short> > rt(numAtoms);	
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		for(short i = 0; i < numAtoms; i++){
			rt[i].insert(connAtoms[i].begin(), connAtoms[i].end());
		}

		double slcgi = 0.0;
		double sci = 0;
		short n = numAtoms;

		while(n > 0){	/* running several times */
			short p = 0;
			short m = 0;		
			vector< vector< set<short>::iterator > > mark(numAtoms);
			for(short i = 0; i < numAtoms; i++){ /* running over all atoms */
				int nai = rt[i].size();
				p += nai;	/* count how many vertices in the molecular */
				if(nai == 1){	/* it has a terminal vertex */
					m++;	
					short ai0 = *rt[i].begin();
					set<short>::iterator it = rt[ai0].find(i);
					if(it != rt[ai0].end()){
						mark[ai0].push_back(it);	/* make a mark connection to which atom */
					}
					rt[i].clear();	/* remove the terminal vertex, atom-i has no connection to other atom */
					break;
				} /* end of IF() */
			}

			for(short i = 0; i < numAtoms; i++){
				for(vector< set<short>::iterator >::iterator it = mark[i].begin(); it != mark[i].end(); it++)
					rt[i].erase(*it);/* remove the atom that made mark before */
			}

			sci += m * m;	/* calculate the SCI step by step */

			if(m == 0)	/* has no terminal vertex */
				break;
			else{
				double r = double(m) / numAtoms;
				slcgi += -double(r) * log10(r) / log10(2.0);	/* calculate the structure Lopping centric group index: SLCGI */
			}

			if(p == 1)	/* just one vrtex right now -- centric atom */
				break;
			else 
				if(p == 2){	/* last two atoms -- centric atoms */
					sci--;
					slcgi -= -1.0 / numAtoms * log10(1.0 / numAtoms) / Constants::LOG20;
					break;
				}
				n--;
		}/* end of WHILE() */

		/* calculate the structure Lopping centric group index: SLCGI */
		slcgi += -1.0 / numAtoms * (log10(1.0 / numAtoms)/ Constants::LOG20);

		if(sci != 0)
			sci++;
		else
			slcgi = 0.0;	/* if has no terminal vertex, SLCGI equal to 0 */
		vector<double> fingerprint(2, 0);
		fingerprint[0] = sci;
		fingerprint[1] = slcgi;
		mol.setFingerprint(fingerprint, ::CentricIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in CentricIndex" << endl;
	}
}

void Calculator5::RadialCentricIndex(void){
	try{
		if(mol.isFeaturesReady(::RadialCentricIndex)) return;

		// must get a copy not a reference
		vector<short> eta = mol.getDistanceMatrix(false).getEta();
		/* ---------------------- calculate RCI ---------------------------- */
		double	rci = calculatexDEI(eta);
		mol.setFingerprint(rci, ::RadialCentricIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in RadialCentricIndex" << endl;
	}
}


void Calculator5::DistanceDegree(void){
	try{
		if(mol.isFeaturesReady(::DistanceDegree)) return;

		/* calculate the Rouvray index as 2 times the Wiener W Index: W */
		double w = mol.getDistanceMatrix(false).sumAllCells();

		vector<double> fingerprint(2, 0);
		/* ---------------------- calculate MDDMI ---------------------------- */
		// sigma  will be modified. So get a copy not a reference
		vector<int> sigma = mol.getDistanceMatrix(false).getSigma();
		for(short i = 0; i < numAtoms; i++){
			if(sigma[i] > 0)
				fingerprint[1] += -double(sigma[i]) / w * log10(sigma[i] / w) / Constants::LOG20;
		}

		/* ---------------------- calculate IDDEI ---------------------------- */
		fingerprint[0] = this->calculatexDEI(sigma);
		mol.setFingerprint(fingerprint, ::DistanceDegree);
	}catch(exception &ex){
		cout << ex.what() << " thrown in DistanceDegree" << endl;
	}
}

void Calculator5::KierAtomInfoIndex(void){
	try{
		if(mol.isFeaturesReady(::KierAtomInfoIndex)) return;

		vector<float> valence = mol.getValence();
		/* ------------------------------ A0PII ---------------------------------------- */
		double A0PII = calculatexDEI(valence);

		A0PII *= numAtoms;
		mol.setFingerprint(A0PII, ::KierAtomInfoIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in KierAtomInfoIndex" << endl;
	}
}

template<class T> double Calculator5::calculatexDEI(vector<T> & params){
	/* ------------------- calculate the equivalence calsses: ec[][] ------------ */
	short m = 0;
	vector< vector<short> > ec;
	for(short i = 0; i < numAtoms; i++){
		if(params[i] != 0){	/* does not empty */
			/* working the first element */
			ec.push_back(vector<short>(1, i + 1));

			for(short j = i+1; j < numAtoms; j++){	/* working other elements */
				if(params[i] == params[j]){	/* found same element */
					ec[m].push_back(j + 1);	/* elements of the class */
					params[j] = 0;
				}
			}
			params[i] = 0;
			m++;	/* number of equivalence class */
		}
	}

	/* -------------------------------- xDEI ------------------------------------- */
	double xdei = 0.0;
	short np = ec.size();
	for(short i = 0; i < np; i++){
		double r = double(ec[i].size()) / numAtoms;
		xdei += -r * log10(r) / Constants::LOG20;
	} 
	return xdei;
}

void Calculator5::InfoVertexDegree(void){
	try{
		if(mol.isFeaturesReady(::InfoVertexDegree)) return;

		// bond_type will be modified. So get a copy not a reference
		vector<float> bond_type = mol.getConventionalBondOrder();

		vector<double> fingerprint(2, 0);
		/* -------------------------------- IBI ------------------------------------- */
		double div = numAtoms == 1? 1.0 : 2.0 * mol.getBondsNum(false);
		for(short i = 0; i < numAtoms; i++){
			double r = bond_type[i] / div;
			if(r > 0)
				fingerprint[1] += -r * log10(r) / Constants::LOG20;
		}

		/* -------------------------------- IVDEI ------------------------------------- */
		fingerprint[0] = calculatexDEI(bond_type);

		mol.setFingerprint(fingerprint, ::InfoVertexDegree);
	}catch(exception &ex){
		cout << ex.what() << " thrown in InfoVertexDegree" << endl;
	}
}

void Calculator5::VertexComplexity(void){
	try{
		if(mol.isFeaturesReady(::VertexComplexity)) return;

		/* -------------------- calculate the gfi[i][g] and vertex complexity ---------------------- */
		double vdpci = 0.0;
		double cvdpci = 0.0;
		vector< vector<int> > gfi(numAtoms, vector<int>(numAtoms, 0));
		vector<double> vci(numAtoms, 0);
		vector<double> ui(numAtoms, 0);
		vector<double> vi(numAtoms, 0);
		vector<double> xi(numAtoms, 0);
		vector<double> yi(numAtoms, 0);

		DistanceMatrix<short> & distanceMatrix= mol.getDistanceMatrix(false);	
		/*------------------------- calculate sigma of atom-i ------------------------ */
		vector<int> &sigma = distanceMatrix.getSigma();
		vector<short> &eta = distanceMatrix.getEta();
		double w = distanceMatrix.sumAllCells();

		for(short i = 0; i < numAtoms; i++){
			/* ------- calculate the number of distances from the vertex Vi equal to g: gfi[i][g] ------- */
			for(short j = 0; j < numAtoms; j++){
				gfi[i][distanceMatrix.getValueInCell(i, j)]++;	/* number of distances from the vertex */
			}
			gfi[i][0] = 1;	/* the start value always equal to 1 */


			/* ----------------------------- calculate the vertex distance path count index: VDPCI ---------------------------- */
			for(short j = 0; j <= eta[i]; j++){
				double r = double(gfi[i][j]) / numAtoms;
				if(fabs(r) > Constants::ROUNDOFF)
					vci[i] += -r * log10(r) / Constants::LOG20;	/* vertex complexity */
			}
			vdpci += vci[i];	/* vertex distance path count index: VDPCI */

			/* --------------------------- calculate the complexity vertex distance path count index: CVDPCI ----------------------- */
			for(short j = 0; j < numAtoms; j++){
				double r = double(distanceMatrix.getValueInCell(i, j)) / sigma[i];
				if(fabs(r) > Constants::ROUNDOFF){
					ui[i] += -r * log10(r) / Constants::LOG20;	/* vertex distance complexity */
				}
			}
			cvdpci += ui[i] * sigma[i] / w;	/* calculate the complexity vertex distance path count index: CVDPCI */

			/* ---------------------------- relative vertex distance complexity ------------------------- */
			vi[i] = sigma[i] * log10(double(sigma[i])) / Constants::LOG20 - ui[i];

			/* ----------------- calculate the mean extended local information on distances --------------- */
			for(short j = 1; j <= eta[i]; j++){
				yi[i] += gfi[i][j] * j * log10(double(j)) / Constants::LOG20;	/* the mean extended local information on distances */
			}

			/* ----------------- calculate the extended local information on distances --------------- */
			xi[i] = sigma[i] * log10(double(sigma[i])) / Constants::LOG20 - yi[i];

		} /* end of FOR-i */

		vdpci /= numAtoms;	/* VDPCI */

		if(fabs(cvdpci) > Constants::ROUNDOFF)
			cvdpci = log10(cvdpci) / Constants::LOG20;		/* CVDPCI */
		else
			cvdpci = 0.0;

		vector<double> fingerprint;
		fingerprint.push_back(vdpci);
		fingerprint.push_back(cvdpci);
		short numBonds = mol.getBondsNum(false);	
		double times = numBonds / (numBonds - numAtoms + 2);
		// the vertex distance information index: VDII
		fingerprint.push_back(calculatexVDII(ui, numBonds, times)); 

		/// the relative of vertex distance information index: RVDII
		fingerprint.push_back(calculatexVDII(vi, numBonds, times));

		// the mean of vertex distance information index: MVDII ---------------------------- */
		fingerprint.push_back(calculatexVDII( xi, numBonds, times));

		// the extended of vertex distance information index ---------------------------- */
		fingerprint.push_back(calculatexVDII(yi, numBonds, times));

		mol.setFingerprint(fingerprint, ::VertexComplexity);
	}catch(exception &ex){
		cout << ex.what() << " thrown in VertexComplexity" << endl;
	}
}

double Calculator5::calculatexVDII(vector<double>& x_i, short numBonds, double times){	
	double ret = 0.0;
	for(short i = 0; i < numBonds; i++){
		Bond & bond = mol.getBond(i);
		short k = bond.getFirstAtomNumber();
		short m = bond.getSecondAtomNumber();

		if(fabs(x_i[k] * x_i[m]) > Constants::ROUNDOFF)
			ret += 1.0 / sqrt(x_i[k] * x_i[m]);
	}
	ret *= times;	/* vertex distance information index: VDII */
	return ret;
}


void Calculator5::MaximalEigenvalue(void){	
	try{
		if(mol.isFeaturesReady(::MaximalEigenvalue)) return;

		Array2D<double> a(numAtoms, numAtoms, 0.0);
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		size_t k = 0;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j< numAtoms; j++){
				short dis = distanceMatrix[k++];
				if(dis <= 1)
					a[j][i] = a[i][j] = dis; /* get adjacency matrix and the data start from 1 to N */
			}
		}
		vector<double> ev = getEigenvalue(a); /* function call: calculate the eigenvalues */
		double maxEv = *max_element(ev.begin(), ev.end());
		mol.setFingerprint(maxEv, ::MaximalEigenvalue);
	}catch(exception &ex){
		cout << ex.what() << " thrown in MaximalEigenvalue" << endl;
	}
}

void Calculator5::DistancesOfAtomic(void){
	try{
		if(mol.isFeaturesReady(::DistancesOfAtomic)) return;

		vector<double> ad(36, 0);
		//must get a copy not refer to original array
		vector<short> atomic_num = mol.getAtomicNum(false);

		vector< deque<short> > keyAtoms(8);
		for(short i = 0; i < numAtoms; i++){
			if(atomic_num[i] > 6 && atomic_num[i] < 54){
				switch(atomic_num[i]){
				case 7:
					keyAtoms[0].push_back(i);
					break;
				case 8:
					keyAtoms[1].push_back(i);
					break;
				case 9:
					keyAtoms[4].push_back(i);
					break;
				case 15:
					keyAtoms[3].push_back(i);
					break;
				case 16:
					keyAtoms[2].push_back(i);
					break;
				case 17:
					keyAtoms[5].push_back(i);
					break;
				case 35:
					keyAtoms[6].push_back(i);
					break;
				case 53:
					keyAtoms[7].push_back(i);
					break;
				}
			}
		}

		DistanceMatrix<short> & distanceMatrix = mol.getDistanceMatrix(false);
		short aid[] = {0, 8, 15, 21, 26, 30, 33, 35};
		for( vector<short>::size_type i = 0; i < keyAtoms.size(); i++){ // every key atom groups
			while(keyAtoms[i].size()  > 0){ // for every key atom in group i
				short ai = keyAtoms[i].front();
				for( vector<short>::size_type m = i; m < keyAtoms.size(); m++){
					short id = aid[i] + m - i;
					for( vector<short>::size_type n = 0; n < keyAtoms[m].size(); n++){ 
						short dis = distanceMatrix.getValueInCell(ai, keyAtoms[m][n]);
						if(dis > 2){
							ad[id] += dis;
						}
					}				
				}
				keyAtoms[i].pop_front(); // remove this atom
			}
		}

		mol.setFingerprint(ad, ::DistancesOfAtomic);
	}catch(exception &ex){
		cout << ex.what() << " thrown in DistancesOfAtomic" << endl;
	}
}

void Calculator5::MolecularWalkCount(void){
	try{
		if(mol.isFeaturesReady(::MolecularWalkCount)) return;

		/* ---- calculate molecular walk count of order 1~numAtoms ---- */
		vector< vector<short> > awc(numAtoms, vector<short> (7, 1));
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector<short> one, two;
		for(short i = 0; i < numAtoms; i++){	/* number of atoms */
			awc[i][1] = connAtoms[i].size();	/* walk count order 1 */
			one.assign(connAtoms[i].begin(), connAtoms[i].end());
			for(short j = 2; j <= 6; j++){	/* number of order of walk count order 1 ~ 10 */
				short np = one.size();
				two.clear();
				for(short k = 0; k < np; k++){	/* currect atomic level become from up level */
					/*  number of connection atoms of up level's atoms the currect level's atoms */
					two.insert(two.end(), connAtoms[one[k]].begin(), connAtoms[one[k]].end());
				}
				awc[i][j] = two.size();	/* atomic walk count in order j */
				swap(one, two);
			}

		} /* end of FOR-i */

		/* ---------- claculate the molcular walk count order 0 ~ 10 ---------- */
		short wcoSize = numAtoms > 11 ? numAtoms:11;
		vector<int> wco(wcoSize, 0);
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j <= 6; j++){	/* calculate number of walk count order */
				wco[j] += awc[i][j]; /* molecular walk count as sum of atomic walk walk count */
			}
		}

		/* calculated as half-sum of all atomic walk count (just calculated 10 atomic walk count in the project) */
		double wcm = 0;
		for(short i = 0; i <= 10; i++){
			if(i != 0)
				wco[i] /= 2;	/* molecular walk count is calculated as the half-sum of all atomic walk count */

			wcm += wco[i]; /* total walk count as sum of molecular walk count */
		}


		/* printout results of molecular walk count of order 01~10 to output file */
		vector<double> fingerprint;
		for(short i = 1; i <= 6; i++)
			fingerprint.push_back(wco[i]);
		fingerprint.push_back(wcm);
		mol.setFingerprint(fingerprint, ::MolecularWalkCount);
	}catch(exception &ex){
		cout << ex.what() << " thrown in MolecularWalkCount" << endl;
	}
}

void Calculator5::WalkReturningWalkCount(void){
	try{
		if(mol.isFeaturesReady(::WalkReturningWalkCount)) return;

		vector< vector<int> > srwc(numAtoms, vector<int>(6, 0));
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);

		/* ---- calculate molecular walk-returning walk count of order 1~numAtoms ---- */
		vector<short> one, two;
		for(short i = 0; i < numAtoms; i++){	/* number of atoms */
			srwc[i][0] = 1;	/* walk count order 0 */
			two.clear();
			one.clear();
			for(short j = 0; j < 6; j++){	/* number of walk count order 1 ~ 10 */
				if(j == 0){
					two.assign(connAtoms[i].begin(), connAtoms[i].end());
				}else{
					short na = one.size();
					two.clear();
					for(short k = 0; k < na; k++){	/* currect atomic level become from up level */
						/*  number of connection atoms of up level's atoms the currect level's atoms */
						two.insert(two.end(), connAtoms[one[k]].begin(), connAtoms[one[k]].end());	/* atoms on the level */
					}
				}

				/* ------------------- look for walk-returning atoms ---------------------- */
				short naj = two.size();
				for(short k = 0; k < naj; k++){
					if(two[k] == i)
						srwc[i][j]++;
				}
				swap(one, two);
			}

		} /* end of FOR-i */

		/* ---------- claculate the molcular walk-returning walk count order 0 ~ 10 ---------- */
		vector<double> msrwc(6, 0);
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < 6; j++){	/* calculate number of walk count order */
				msrwc[j] += srwc[i][j]; /* molecular walk-returning walk count as sum of atomic walk walk count */
			}
		}

		mol.setFingerprint(msrwc, ::WalkReturningWalkCount);
	}catch(exception &ex){
		cout << ex.what() << " thrown in WalkReturningWalkCount" << endl;
	}
}

void Calculator5::TopolChargeIndex(void){	
	try{
		if(mol.isFeaturesReady(::TopolChargeIndex)) return;

		//calculate AdjacentMatrix[][] and ReciprocalMatrix[i][j]
		vector< vector<short> > AdjacentMatrix(numAtoms, vector<short>(numAtoms, 0));
		vector< vector<float> > ReciprocalMatrix(numAtoms, vector<float>(numAtoms, 0));
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		size_t k = 0;
		for(short i = 0; i<numAtoms; i++){
			for(short j = i + 1; j<numAtoms; j++){
				short dis = distanceMatrix[k++];
				//copy distance_matrix[i][j].path_count value to AdjacentMatrix[][]
				if(dis == 1)
					AdjacentMatrix[j][i] = AdjacentMatrix[i][j] = 1;
				//copy distance_matrix[i][j].path_count value to ReciprocalMatrix[][]
				if(dis != 0)
					ReciprocalMatrix[j][i] = ReciprocalMatrix[i][j] = 1.0f / (dis * dis);
			}
		}

		//calculate GalvezMatrix[][]
		vector< vector<float> > GalvezMatrix(numAtoms, vector<float>(numAtoms, 0));
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				for(short k = 0; k< numAtoms; k++){
					GalvezMatrix[i][j] += AdjacentMatrix[i][k] * ReciprocalMatrix[j][k];
				}
			}
		}

		//calculate ChargeMatrix[][]
		vector< vector<float> > ChargeMatrix(numAtoms, vector<float>(numAtoms, 0));
		for(short i = 0; i < numAtoms; i++){
			for(short j = 0; j < numAtoms; j++){
				if(i != j)
					ChargeMatrix[i][j] = GalvezMatrix[i][j] - GalvezMatrix[j][i];
				else
					ChargeMatrix[i][j] = GalvezMatrix[i][j];

				if(ChargeMatrix[i][j] < 0)
					ChargeMatrix[i][j] = -ChargeMatrix[i][j];
			}
		}

		//calculate topological charge index of order 1 ~ 10
		vector<double> MTCI(11, 0), MMTCI(11, 0);
		k = 0;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j < numAtoms; j++){
				short dis = distanceMatrix[k++];
				if(dis > 0 && dis < 11)
					MTCI[dis] += ChargeMatrix[i][j]; 
			}
		}

		vector<double> fingerprint;
		for (short i=1; i<11; i++) //print out the MTCI[i] group values 
			fingerprint.push_back(MTCI[i]);

		double SMTMCI = 0.0;
		for (short i=1; i<11; i++){ //print out the results of Mean molecular topological charge index order 1 ~10
			MMTCI[i] = MTCI[i];
			if(numAtoms > 1)
				MMTCI[i] /= (numAtoms - 1);
			SMTMCI += MMTCI[i];
			fingerprint.push_back(MMTCI[i]);
		}
		fingerprint.push_back(SMTMCI);
		mol.setFingerprint(fingerprint, ::TopolChargeIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in TopolChargeIndex" << endl;
	}
}

void Calculator5::SchultzMolecular(void){
	try{
		if(mol.isFeaturesReady(::SchultzMolecular)) return;

		int mti = 0;
		/* printout elements of distance_matrix[][].path_count */
		vector< vector<short> > & connAtoms = mol.getConnectedAtomIDs(false);
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		vector<int> & sigma = mol.getDistanceMatrix(false).getSigma();
		for(short i = 0; i < numAtoms; i++){
			/* calculate the Schultz type Molecular Topological index: SMTI*/
			mti += (connAtoms[i].size() + sigma[i]) * connBond[i].size();
		}
		mol.setFingerprint(mti, ::SchultzMolecular);
	}catch(exception &ex){
		cout << ex.what() << " thrown in SchultzMolecular" << endl;
	}
}

void Calculator5::GutmanMolecular(void){
	try{
		if(mol.isFeaturesReady(::GutmanMolecular)) return;

		int mtdi = 0;
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		size_t k = 0;
		for(short i = 0; i < numAtoms; i++){/* start atoms */
			for(short j = i+1; j < numAtoms; j++)/* endpoint atoms */
				/* calculate the Molecular Topological Distance Index: mtdi */
				mtdi += connBond[i].size() * connBond[j].size() * distanceMatrix[k++];
		}
		mol.setFingerprint(mtdi, ::GutmanMolecular);
	}catch(exception &ex){
		cout << ex.what() << " thrown in GutmanMolecular" << endl;
	}
}

void Calculator5::GutmanMTI_VVD(void){
	try{
		if(mol.isFeaturesReady(::GutmanMTI_VVD)) return;

		// Molecular Topolgical Distance Index of valence vertex degrees: MTDIv
		double mtdiv = 0.0;
		vector<double> mtdiv_valence(numAtoms, 0); //valence electrons that without H
		vector<float> & cbo = mol.getConventionalBondOrder();
		mtdiv_valence.assign(cbo.begin(), cbo.end());
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		for(short i = 0; i < numAtoms; i++){
			/* ---------------------- calculate the conventional bond order ---------------- */
			Atom & atom = mol.getAtom(i);
			char symbol = atom.getAtomSymbol();
			/* ------------------ calculate valence electrons --------------------------- */
			/* valence vertex degree = valence electrons + lone pair electrons */
			if(symbol == Elements::NITROGEN){/* symbol N */
				mtdiv_valence[i] += 2;
			}else if(symbol == Elements::OXYGEN){/* symbol O */
				mtdiv_valence[i] += 4;
			}else if(symbol == Elements::FLOURINE){/* symbol F */
				mtdiv_valence[i] += 6;
			}
			/*valence vertex degree(Si) = valence vertex degree / (its atomic number - valence vertex degree - 1)*/
			else if(symbol == Elements::SILICON){/* symbol Si */
				mtdiv_valence[i] = mtdiv_valence[i] / (13 - mtdiv_valence[i]);
			}
			/* P, S, Cl, Br, I: valence vertex degree = (valence vertex degree + lone pair electrons) /
			(its atomic number - (valence vertex degree + lone pair electrons) - 1)*/
			else if(symbol == Elements::PHOSPHORUS){/* symbol P */
				if(mtdiv_valence[i] > 3)
					mtdiv_valence[i] = mtdiv_valence[i] / (14 - mtdiv_valence[i]);
				else
					mtdiv_valence[i] = (mtdiv_valence[i] + 2) / (12 - mtdiv_valence[i]);
			}else if(symbol == Elements::SULFUR){/* symbol S */
				int n = 0;
				for(short j = 0; j < (short)connBond[i].size(); j++) /* compare bonds' type */
					n += connBond[i][j];

				if(n == 2)
					mtdiv_valence[i] = (mtdiv_valence[i] + 4) / (11 - mtdiv_valence[i]);
				else if(n == 4)
					mtdiv_valence[i] = (mtdiv_valence[i] + 2) / (13 - mtdiv_valence[i]);
				else
					mtdiv_valence[i] = mtdiv_valence[i] / (15 - mtdiv_valence[i]);

			}else if(symbol == Elements::CLORINE){/* symbol Cl */
				mtdiv_valence[i] = (mtdiv_valence[i] + 6) / (10 - mtdiv_valence[i]);
			}else if(symbol == Elements::BROMINE){/* symbol Br */
				mtdiv_valence[i] = (mtdiv_valence[i] + 6) / (28 - mtdiv_valence[i]);
			}else if(symbol == Elements::IODINE){ /* symbol I */
				mtdiv_valence[i] = (mtdiv_valence[i] + 6) / (46 - mtdiv_valence[i]);
			}
		}

		/* ------------------------ calculate the mtdiv ----------------------------------- */
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		size_t k = 0;
		for(short i = 0; i < numAtoms; i++){
			for(short j = i + 1; j < numAtoms; j++){
				mtdiv += mtdiv_valence[i] * mtdiv_valence[j] * distanceMatrix[k++];
			}
		}
		mol.setFingerprint(mtdiv, ::GutmanMTI_VVD);
	}catch(exception &ex){
		cout << ex.what() << " thrown in GutmanMTI_VVD" << endl;
	}
}

void Calculator5::PrincipalQuantum(void){
	try{
		if(mol.isFeaturesReady(::PrincipalQuantum)) return;
		double vepqi = 0; //Valence electrons of principal quantum Index: VEPQI
		for(short i = 0; i < numAtoms; i++){
			vepqi += Constants::val_ele_pri_num[mol.getAtom(i).getAtomSymbol()];
		}
		mol.setFingerprint(vepqi, ::PrincipalQuantum);
	}catch(exception &ex){
		cout << ex.what() << " thrown in PrincipalQuantum" << endl;
	}
}

void Calculator5::molecularTopological(void){
	try{
		if(mol.isFeaturesReady(::MolecularTopological)) return;

		/* calculate the Morlecular size and branching index: MSBI */
		double d1 = 0.0;
		double d2 = 0.0;
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		vector<int> & sigma = mol.getDistanceMatrix(false).getSigma();
		for(short i = 0; i < numAtoms; i++){
			int numBond = connBond[i].size();
			d1 += numBond * (sigma[i] * sigma[i]);
			d2 += numBond * sigma[i];
		}
		double d = d1 / d2;
		double msbi = d > 0 ? sqrt(double(numAtoms)) * log(d) : 0; /* calculate the Morlecular size and branching index: MSBI */

		mol.setFingerprint(msbi, ::MolecularTopological);
	}catch(exception &ex){
		cout << ex.what() << " thrown in MolecularTopological" << endl;
	}
}

void Calculator5::TerminalVertex(void){
	try{
		if(mol.isFeaturesReady(::TerminalVertex)) return;

		/*------------------------ calculate the Index of terminal vertex matrix: ITVM ----------------------*/
		vector<double> product(numAtoms, 1.0);
		int rows = 0;
		int n = 0;
		DistanceMatrix<short> & distanceMatrix = mol.getDistanceMatrix(false);
		vector< vector<char> > & connBond = mol.getConnectedBondTypes(false);
		for(short i = 0; i < numAtoms; i++){
			if(connBond[i].size() == 1){ /* the terminal vertices */
				rows++; /* for check signal terminal vertix */
				for(short j = 0; j< numAtoms; j++){
					short dis = distanceMatrix.getValueInCell(i, j);
					if(dis != 0){ /* using nonzero row elements */
						product[j] *= dis; /* product elements of pendent matrix */
					}
					else
						n = j; /* get the address of 0 element that in the pendent matrix */
				}
			}
		}/* end for_i loop */

		if(rows == 1) /* signal terminal vertix */
			product[n] = 0; /* the element cell change back to original matrix -- from 1 to 0*/

		double s = 0.0;
		double itvm = 0.0;
		for(short i = 0; i < numAtoms; i++)
			s += product[i]; /* collection elements of product[] */

		if(s > 1.0e-9)
			itvm = sqrt(s); /* get Index of terminal vertex matrix: ITVM */

		mol.setFingerprint(itvm, ::TerminalVertex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in TerminalVertex" << endl;
	}
}

// This function calculate WienerIndex and meanWienerIndex
void Calculator5::Wiener_MeanWienerIndex(void){
	try{
		if(mol.isFeaturesReady(::Wiener_MeanWienerIndex)) return;

		/* calculate the Wiener Index: WI */
		int wi = mol.getDistanceMatrix(false).sumAllCells();
		/* calculate Aplwi as mean Wiener index */	 
		double Aplwi = numAtoms == 1 ? wi : double(wi) / (numAtoms * (numAtoms - 1));

		wi /= 2; /* calculate WI as half of total topological distance: Wiener index */
		vector<double> fingerprint(2, 0);
		fingerprint[0] = wi;
		fingerprint[1] = Aplwi;
		mol.setFingerprint(fingerprint, ::Wiener_MeanWienerIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in Wiener_MeanWienerIndex" << endl;
	}	
}

void Calculator5::ReciprocalDistance(void){
	try{
		if(mol.isFeaturesReady(::ReciprocalDistance)) return;

		/* calculate the Reciprocal Distance Wiener-type Index: RIWDM */
		double RIWDM = 0.0;
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		for(int i = distanceMatrix.size() - 1; i >= 0; i--){
			if(distanceMatrix[i] > 0) RIWDM += 1.0 / distanceMatrix[i];
		}
		mol.setFingerprint(RIWDM, ::ReciprocalDistance);
	}catch(exception &ex){
		cout << ex.what() << " thrown in ReciprocalDistance" << endl;
	}
}

void Calculator5::HararyIndex(void){
	try{
		if(mol.isFeaturesReady(::HararyIndex)) return;

		/* calculate the Harary Index */
		double harary = 0.0;
		vector<short> & distanceMatrix = mol.getDistanceMatrix(false).getData();
		for(int i = distanceMatrix.size() - 1; i >= 0; i--){
			if(distanceMatrix[i] > 0)
				harary += 2.0 / (distanceMatrix[i] * (1 + distanceMatrix[i])); /* calculate distance matrix */
		}
		mol.setFingerprint(harary, ::HararyIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in HararyIndex" << endl;
	}
}

void Calculator5::Wiener_ReciprocalWiener_PathIndex(void){
	try{
		if(mol.isFeaturesReady(::Wiener_ReciprocalWiener_PathIndex)) return;

		vector<double> fingerprint(2, 0);
		fingerprint[0] = sum_element(mol.getDistancePath().getData()) / 2.0; /* Wiener-Path index: WPI */

		vector<int> & distancePath = mol.getDistancePath().getData();
		for(int i = distancePath.size() - 1; i >= 0; i--){
			if(distancePath[i] > 0){
				fingerprint[1] += 1.0 / distancePath[i];
			}
		}

		fingerprint[1] /= 2.0; /* reciprocal Wiener-Path index: RWPI */

		mol.setFingerprint(fingerprint, ::Wiener_ReciprocalWiener_PathIndex);
	}catch(exception &ex){
		cout << ex.what() << " thrown in Wiener_ReciprocalWiener_PathIndex" << endl;
	}
}

void Calculator5::doAll(void){
	SumValenceVertex1();
	ReciprocalDistanceIndex();	
	DistanceInfo();
	PetitjeanIndex();
	CentricIndex();
	RadialCentricIndex();
	DistanceDegree();
	KierAtomInfoIndex();
	InfoVertexDegree();
	VertexComplexity();
	MaximalEigenvalue();
	DistancesOfAtomic();
	MolecularWalkCount();
	WalkReturningWalkCount();

	TopolChargeIndex();
	SchultzMolecular();
	GutmanMolecular();
	GutmanMTI_VVD();
	PrincipalQuantum();
	molecularTopological();
	TerminalVertex();
	Wiener_MeanWienerIndex();
	ReciprocalDistance();
	HararyIndex();
	Wiener_ReciprocalWiener_PathIndex();
}