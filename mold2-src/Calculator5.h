#pragma once

class Calculator5
{
	Molecule & mol;
	short numAtoms;
public:
	Calculator5(Molecule & _mol);
	~Calculator5(void);
	void SumValenceVertex1(void);
	void ReciprocalDistanceIndex(void);	
	void DistanceInfo(void);
	void PetitjeanIndex(void);
	void CentricIndex(void);
	void RadialCentricIndex(void);
	void DistanceDegree(void);
	void InfoVertexDegree(void);
	void VertexComplexity(void);
	void MaximalEigenvalue(void);
	void DistancesOfAtomic(void);
	void MolecularWalkCount(void);
	void WalkReturningWalkCount(void);
	void KierAtomInfoIndex(void);
	
	void TopolChargeIndex(void);
	void SchultzMolecular(void);
	void GutmanMolecular(void);
	void GutmanMTI_VVD(void);
	void PrincipalQuantum(void);
	void molecularTopological(void);
	void TerminalVertex(void);
	void Wiener_MeanWienerIndex(void);
	void ReciprocalDistance(void);
	void HararyIndex(void);
	void Wiener_ReciprocalWiener_PathIndex(void);
	void doAll(void);
private:
	double calculatexVDII(vector<double>& x_i, short numBonds, double times);
	template<class T> double calculatexDEI(vector<T> & params);
};
