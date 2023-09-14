#include "stdafx.h"
#include "Elements.h"

// atomic symbols
const std::string Elements::m_Symbols[] = {
	"H", "He","Li","Be","B", "C", "N", "O", "F", "Ne",
	"Na","Mg","Al","Si","P", "S", "Cl","Ar","K", "Ca",
	"Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn",
	"Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y", "Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
	"Sb","Te","I", "Xe","Cs","Ba","La","Ce","Pr","Nd",
	"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	"Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg",
	"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
	"Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
	"Md","No","Lr"
};	

// const elements
const char Elements::HYDROGEN = 0;
const char Elements::CARBON = 5;
const char Elements::NITROGEN = 6;
const char Elements::OXYGEN = 7;
const char Elements::SILICON = 13;
const char Elements::PHOSPHORUS = 14;
const char Elements::SULFUR = 15;
const char Elements::FLOURINE = 8;
const char Elements::CLORINE = 16;
const char Elements::BROMINE = 34;
const char Elements::IODINE = 52;

//Total elements
const char Elements::TOTALELEMENTS = 103;

//Marker
const char Elements::NOTFOUND = -1;

std::map<std::string, char> Elements::m_SymbolToidMap;

Elements::Elements(void)
{
}

Elements::~Elements(void)
{
}

// Convert atomic order to symbol
string Elements::idToSymbol(char atomicID){
	return m_Symbols[atomicID];
};

// convert symbol to atomic id
char Elements::symbolToId(std::string symbol){
	if(m_SymbolToidMap.empty()){
		for(char i = 0; i < TOTALELEMENTS; i++){			
			m_SymbolToidMap.insert(std::pair<std::string, char>(m_Symbols[i], i));
		}
	}

	std::map<std::string, char>::iterator it = m_SymbolToidMap.find(symbol);
	if(it != m_SymbolToidMap.end()){
		return it->second;
	}else{
		return Elements::NOTFOUND;
	}	
};
