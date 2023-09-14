#pragma once

template <class T> 
class DistanceMatrix
{
private:
	vector<T> m_element;
	vector<short> eta;
	vector<int> sigma;
	int allSum;
	short m_size;
	size_t k1;
	void calculateSigmaEta();
public:
	DistanceMatrix(void);
	~DistanceMatrix(void);
	bool empty(void);
	void clear(void){
		m_element.clear(); 
		m_size = 0;
		eta.clear();
		sigma.clear();
		allSum = -1;
	}
	void resize(short newSize);
	void setValueToCell(short i, short j, T value);
	T getValueInCell(short i, short j);
	void cellValuePlus(short i, short j, T additionValue);	
	int sumAllCells(void);
	vector<int> & getSigma(void);
	vector<short> & getEta(void);
	int sumCol(short col);
	void print(ostream & out);
	vector<T> & getData(void);
};

template <class T>
vector<T> & DistanceMatrix<T>::getData(void){
	return m_element;
}

template <class T> 
DistanceMatrix<T>::DistanceMatrix(void):m_size(0)
{
	allSum = -1;
}

template <class T> 
DistanceMatrix<T>::~DistanceMatrix(void)
{
}

template <class T> 
int DistanceMatrix<T>::sumAllCells(void){
	if(allSum < 0) calculateSigmaEta();
	return allSum;
}

template <class T> 
void DistanceMatrix<T>::calculateSigmaEta(void){
	sigma.resize(m_size, 0);
	eta.resize(m_size, 0);	

	for(short ix = 0; ix < m_size; ix++){
		size_t k2 = -m_size + ix;
		for(short i = ix - 1; i >= 0; i--){
			int id = (((k1 - i) * (i + 1)) >> 1) + k2;
			sigma[ix] += m_element[id];
			if(m_element[id] > eta[ix]){
				eta[ix] = m_element[id];
			}			
		}
		k2 = (((k1 - ix) * (ix + 1)) >> 1) - m_size + ix + 1;
		for(short i = ix + 1; i < m_size; i++){
			sigma[ix] += m_element[k2];
			if(m_element[k2] > eta[ix]){
				eta[ix] = m_element[k2];
			}
			k2++;
		}
	}
	allSum = sum_element(sigma);
}

template <class T> 
void DistanceMatrix<T>::print(ostream & out){
	for(short i = 0; i < m_size; i++){
		out << endl;
		for(short j = 0; j < m_size; j++){
			out << this->getValueInCell(i, j) << " ";
		}
	}
}

template <class T> 
vector<int> & DistanceMatrix<T>::getSigma(void){
	if(sigma.empty()) calculateSigmaEta();
	return sigma;
}

template <class T> 
vector<short> & DistanceMatrix<T>::getEta(void){
	if(eta.empty()) calculateSigmaEta();
	return eta;
}

template <class T> 
void DistanceMatrix<T>::resize(short newSize)
{
	size_t size = (newSize * (newSize - 1)) >> 1;
	m_element.assign(size, 0);
	m_size = newSize;
	k1 = m_size + m_size - 2;
}

template <class T> 
bool DistanceMatrix<T>::empty(void){
	return m_size == 0;
}

template <class T> 
void DistanceMatrix<T>::setValueToCell(short i, short j, T value)
{
	if(i == j) return; // not save zeros
		
	if(i > j) swap(i, j);
	
	size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
	m_element[index] = value;
}

template <class T> 
T DistanceMatrix<T>::getValueInCell(short i, short j)
{
	if(i == j) return 0; // not save zeros		
	if(i > j) swap(i, j);
	
	size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
	return m_element[index];
}

template <class T> 
void DistanceMatrix<T>::cellValuePlus(short i, short j, T additionValue)
{
	if(i == j) return; // not save zeros		
	if(i > j) swap(i, j);
	size_t index = (((k1 - i) * (i + 1)) >> 1) - (m_size - j);
	m_element[index] += additionValue;
}
