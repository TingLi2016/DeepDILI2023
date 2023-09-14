#include "stdafx.h"

#include "DistancePathMatrix.h"

vector<short> DistancePathMatrix::m_null;

DistancePathMatrix::DistancePathMatrix(void):m_size(0)
{
}

DistancePathMatrix::~DistancePathMatrix(void)
{
}

void DistancePathMatrix::resize(short newSize)
{
	size_t size = (newSize * (newSize - 1)) >> 1;
	m_element.resize(size);
	m_size = newSize;
	k1 = m_size + m_size - 2;
}

void DistancePathMatrix::setPath(short from, short to, vector<short> &path){
	if(from == to) return; // not save zeros
		
	if(from > to) swap(from, to);
	
	size_t index = (((k1 - from) * (from + 1)) >> 1) - (m_size - to);
	m_element[index] = path;
}

vector<short>& DistancePathMatrix::getPath(short from, short to){
	if(from == to) return m_null; // not save zeros
		
	if(from > to) swap(from, to);
	
	size_t index = (((k1 - from) * (from + 1)) >> 1) - (m_size - to);
	return m_element[index];
}

void DistancePathMatrix::pushPathMember(short from, short to, short member){
	if(from == to) return; // not save zeros
		
	if(from > to) swap(from, to);
	
	size_t index = (((k1 - from) * (from + 1)) >> 1) - (m_size - to);
	return m_element[index].push_back(member);
}

short DistancePathMatrix::getPathMember(short from, short to, short ID){
	if(from == to) return -1; // not save zeros
		
	if(from > to) swap(from, to);
	
	size_t index = (((k1 - from) * (from + 1)) >> 1) - (m_size - to);
	return m_element[index][ID];
}

short DistancePathMatrix::getPathSize(short from, short to){
	if(from == to) return 0; // not save zeros
		
	if(from > to) swap(from, to);
	
	size_t index = (((k1 - from) * (from + 1)) >> 1) - (m_size - to);
	return m_element[index].size();
}

bool DistancePathMatrix::empty(void){
	return m_size == 0;
}

void DistancePathMatrix::print(ostream & out){
	for(int i = 0; i < m_size; i++){							
		for(int j = i + 1; j < m_size; j++){
			vector<short> & path = getPath(i, j);
			out << "i=" << i << " j=" << j << endl;
			for( int m = 0; m < (int)path.size(); m++){
				out << path[m] << " ";
			}
			out<<endl;
		}
	}
}


