#pragma once

class Calculator6
{
	Molecule & mol;
	void initializeBurdenMatrix_Lowest(Array2D<double>& burdenmatrix);
	void initializeBurdenMatrix_Highest(Array2D<double>& burdenmatrix);
	short getPolarizabilityValueID(vector< vector<short> > & a2b, char symbol , short i);
	vector<float> calculatePol(void);
public:
	Calculator6(Molecule & _mol);
	~Calculator6(void);

	void polarizability(void);	
	void AtomsAndBonds(void);
	void VanderwaalsVol(void);
	
	void ElectrotopologicalVariation(void);
	void EigenvalueLowest(void);
	void EigenvaluePolLowest(void);
	void EigenvalueHighest(void);
	void EigenvaluePolHighest(void);	
	void Geary_MoranAutoCorrelation(void);
	void SP_number(void);
	void Empirical(void);
	void InformationConent0_5(void);	
	void AutoCorrelation(void);
	void doAll(void);
};
