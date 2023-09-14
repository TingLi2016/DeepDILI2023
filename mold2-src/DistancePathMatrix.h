#pragma once

class DistancePathMatrix
{
	vector< vector<short> > m_element;
	short m_size;
	static vector<short> m_null;
	size_t k1;

public:
	DistancePathMatrix(void);
	~DistancePathMatrix(void);
	void clear(void){m_element.clear(); m_size = 0;}
	bool empty(void);
	void resize(short newSize);
	void setPath(short from, short to, vector<short> &path);
	vector<short>& getPath(short from, short to);
	void pushPathMember(short from, short to, short member);
	short getPathMember(short from, short to, short ID);
	short getPathSize(short from, short to);
	void print(ostream & out);
};
