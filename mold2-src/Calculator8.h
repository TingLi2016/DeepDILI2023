#pragma once

template <typename T>
class SymmetricMatrix
{
private:
	vector<T> m_element;
	short m_size;
	size_t k1;
public:
	SymmetricMatrix(short newSize){
		int size = (newSize * (newSize + 1)) >> 1;
		m_element.resize(size);
		for(int i = 0; i < size; i++){
			m_element[i] = 0;
		}
		m_size = newSize;
		k1 = m_size + m_size;
	}

	~SymmetricMatrix(void){
		m_element.clear();
	}

	inline void setValueToCell(short i, short j, T value){
		if(i > j){
			swap(i, j);
		}
		size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
		m_element[index] = value;
	}

	inline T getValueInCell(short i, short j){
		if(i > j){
			swap(i, j);
		}
		size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
		return m_element[index];
	}

	inline void cellValuePlus(short i, short j, T additionValue){
		if(i > j){
			swap(i, j);
		}
		
		size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
		m_element[index] += additionValue;
	}

	inline T cumulativeSum(){
		T sum = 0;		
		for(int i = m_element.size() - 1; i >= 0; i--){			
			sum += m_element[i];
		}
		sum *= 2; // double elements
		size_t x = 0;
		size_t k = m_size;
		size_t maxSize = m_element.size();
		while(x < maxSize){ // minus diagnal elements
			sum -= m_element[x];
			x += k--;			
		}
		return sum;
	}

	void print(ostream & out){
		for(short i = 0; i < m_size; i++){
			for(short j = 0; j < m_size; j++){
				out << this->getValueInCell(i, j) << " ";
			}
			out << endl;
		}
	}
};

class Calculator8 {
	Molecule& mol;
	vector< SymmetricMatrix<float> > dmew;
	vector< SymmetricMatrix<float> > zmvpw;
	vector< vector<short> > pv;
	vector< vector<char> > laplacianMatrix;

	vector< int > vertexDistanceCounts;

	short numAtoms;
	static const char ZW = 0;
	static const char MW = 1;
	static const char VW = 2;	
	static const char PW = 3;

	void calculateWeight(SymmetricMatrix<float> & xw, double &le_ev, double &sum_ev, double &abs_sum_ev);
	void initializePV();
	void initializePW();
	void initializeBW_MW_VW();
	void initializeLaplacianMatrix();
	double valenceElections(char symbol, char charge, short i, double Zagreb_Mx_valence);
	void initializeVertexDistanceCounts();
public:
	Calculator8(Molecule& mol);
	~Calculator8(void);
	void BW_MW_VW_EWeighted();
	void PolarizabilityWeighted();
	void BalabanAVDCI_BalabanTypeWeighted();
	void SumWeighted();
	void BalabanShortPathIndex();
	void SecondMoharIndex();
	void SpanningTree();
	void LaplacianMatrix();
	void Zagreb();
	void MeanTotalDistance();
	void VertexDegree();
	void Zagreb_M2_v();
	void SchultzMTI_VVD();
	void doAll(void );
};
