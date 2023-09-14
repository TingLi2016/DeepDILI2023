#pragma once

const short FeatureLocation [] = {
	0, 18, 121, 123, 135, 142, 150, 152, 156, 157,
	158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
	168, 170, 171, 172, 173, 176, 177, 178, 179, 180,
	181, 184, 185, 186, 193, 199, 230, 231, 233, 234,
	238, 246, 250, 251, 253, 254, 258, 260, 262, 268,
	298, 299, 300, 321, 331, 341, 364, 365, 401, 408,
	414, 446, 478, 510, 531, 555, 563, 587, 595, 610,
	713, 773, 776
};

const short FeatureGroupLength[] = {
	18, 103, 2, 12, 7, 8, 2, 4, 1, 1,
	1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
	2, 2, 1, 1, 3, 1, 1, 7, 1, 1,
	3, 1, 8, 7, 6, 31, 1, 2, 1, 4,
	8, 4, 1, 2,	1, 4, 2, 2, 6, 30,
	1, 1, 21, 10, 10, 23, 1, 36, 7, 6,
	32, 64, 32, 21, 24, 8, 24, 8, 15, 103,
	60, 3, 1
};

const short TOTALFEATURES = 777;
const short TOTALGROUPS = 73;

const short AromaticBondConvert = 0;
const short CountElements = 1;
const short MolecularWeight = 2;
const short AtomsAndBonds = 3;
const short VertexDegree = 4;
const short vanderwaalsVol = 5;
const short polarizability = 6;
const short Zagreb = 7;
const short Zagreb_M2_v = 8;
const short PrincipalQuantum = 9;
const short SchultzMolecular = 10;
const short SchultzMTI_VVD = 11;
const short GutmanMolecular = 12;
const short GutmanMTI_VVD = 13;
const short MolecularTopological = 14;
const short TerminalVertex = 15;
const short WienerIndex = 16;
const short MeanWienerIndex = 17;
const short ReciprocalDistance = 18;
const short HararyIndex = 19;
const short LaplacianMatrix = 20;
const short WienerPathIndex = 21;
const short ReciprocalWienerPathIndex = 22;
const short SecondMoharIndex = 23;
const short MaximumPathIndex = 24;
const short MinMaxPathIndex = 25;
const short AllPathWiener = 26;
const short BondsWeighted = 27;
const short MassWeighted = 28;
const short VanDerWaalsWeighted = 29;
const short ElectronegativityWeighted = 30;
const short PolarizabilityWeighted = 31;
const short BalabanAVDCI = 32;
const short BalabanTypeWeighted = 33;
const short ElectrotopologicalVariation = 34;
const short AllPathConnectivity = 35;
const short SumValenceVertex1 = 36;
const short ReciprocalDistanceIndex = 37;
const short KierAtomInfoIndex = 38;
const short KierPathIndex = 39;
const short DistanceInfo = 40;
const short PathWalk = 41;
const short PetitjeanIndex = 42;
const short CentricIndex = 43;
const short RadialCentricIndex = 44;
const short MeanTotalDistance = 45;
const short DistanceDegree = 46;
const short InfoVertexDegree = 47;
const short VertexComplexity = 48;
const short InformationContent0_5 = 49;
const short MaximalEigenvalue = 50;
const short SpanningTree = 51;
const short SumWeighted = 52;
const short D_D_RingIndex = 53;
const short D_D_RingsAndCircuitsIndex = 54;
const short MolecularPathCount = 55;
const short BalabanShortPathIndex = 56;
const short DistancesOfAtomic = 57;
const short MolecularWalkCount = 58;
const short WalkReturningWalkCount = 59;
const short autocorrelation = 60;
const short GearyAutocorrelation = 61;
const short MoranAutocorrelation = 62;
const short TopolChargeIndex = 63;
const short eigenvalueLowest = 64;
const short eigenvaluePolLowest = 65;
const short eigenvalueHighest = 66;
const short eigenvaluePolHighest = 67;
const short SP_Number = 68;
const short AliphaticAromatic = 69;
const short AtomCentredFragments = 70;
const short Empirical = 71;
const short MlogP = 72;


//combined features
const short Wiener_MeanWienerIndex = 16;
const short Wiener_ReciprocalWiener_PathIndex = 21; // 3
const short BW_MW_VW_EWeighted = 27; // 6
const short BalabanAVDCI_BalabanTypeWeighted = 32; //2
const short Geary_MoranAutoCorrelation = 61; //
