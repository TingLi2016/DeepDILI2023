// stdafx.cpp : source file that includes just the standard includes
// testSDfile.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"

void trim(std::string & line){
	line.erase(0, line.find_first_not_of(SPACES));
	line.erase(line.find_last_not_of(SPACES) + 1); 
}

vector<double> getEigenvalue(TNT::Array2D<double> & A){
	JAMA::Eigenvalue<double> evEngine(A);
	TNT::Array1D<double> ev;
	evEngine.getRealEigenvalues(ev);
	vector<double> ret(ev.dim());
	for(int i = 0; i < (int)ret.size(); i++)
		ret[i] = ev[i];
	return ret;
}


// TODO: reference any additional headers you need in STDAFX.H
// and not in this file
