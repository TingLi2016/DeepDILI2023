#pragma once

class Elements
{
public:
	const static char NOTFOUND;
	const static char TOTALELEMENTS;

	const static char HYDROGEN;
	const static char CARBON;
	const static char NITROGEN;	
	const static char OXYGEN;
	const static char SILICON;
	const static char PHOSPHORUS;	
	const static char SULFUR;	
	const static char FLOURINE;	
	const static char CLORINE;	
	const static char BROMINE;	
	const static char IODINE;

	static string idToSymbol(char atomicID);

	static char symbolToId(std::string symbol);

	Elements(void);

	~Elements(void);

private:
	const static std::string m_Symbols[];
	static std::map<std::string, char> m_SymbolToidMap;
};

