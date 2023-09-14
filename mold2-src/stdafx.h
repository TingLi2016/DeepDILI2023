// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

//#include "targetver.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <vector>
#include <bitset>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits>
#include <algorithm>

#include "mold2Exception.h"
#include "jama_eig.h"
#include "Features.h"

#define SPACES " \t\r\n"
#define WHITE 0
#define GRAY 1
#define BLACK 2

using namespace std;

void trim(std::string & line);

//template <class T> inline std::string to_string (const T& t);

template <class T> inline std::string to_string (const T & t){
	std::stringstream ss;
	ss << setw(10) << t;
	return ss.str();
}

template <class T> void print(vector<T> & A){
	cout << endl;
	for(int i = 0; i < (int)A.size(); i++)
		cout << A[i] << " ";
}

template <class T> void print(vector< vector<T> > & A){
	for(int i = 0; i < (int)A.size(); i++){
		cout << endl;
		for(int j = 0; j < (int)A[i].size(); j++)
			cout<< A[i][j] << " ";
	}
}


template <class T> inline T sum_element(vector< T > & A){
	T sumV = 0;
	for(int i = A.size() - 1; i >= 0; i--)
		sumV += A[i];
	return sumV;
}

template <class T> inline T sum_element(vector< vector<T> > & A){
	T sumV = 0;
	for(int i = A.size() - 1; i >= 0; i--)
		for(int j = A[0].size() - 1; j >= 0 ; j--)
			sumV += A[i][j];
	return sumV;
}

vector<double> getEigenvalue(Array2D<double> & A);

// TODO: reference additional headers your program requires here
