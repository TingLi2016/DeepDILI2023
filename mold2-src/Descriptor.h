#pragma once

class Descriptor
{
private:
	unsigned int m_count;
	time_t m_start, m_finish;

	void finishCalculation(ostream &report);
	void printDescriptorHeaderLine(ostream & out);
	void reportNewMolecule(ostream & report);
	void printReportHeader(ostream & report, string inputFile, string& outputFile, string reportFile, string  withHFile, string withoutHFile);
	void reportException(ostream & report, string & msg);
	void printFeatures(ostream & out, string& autoID, Molecule & noHmol);

	void printReleaseInfo(void);
	string generateAutoID(string & preFieldName, string & currentFieldName, int ID);
public:
	Descriptor(void);
	~Descriptor(void);

	void printUsage(ostream & out);
	void DoIt(istream & in, ostream & out);
};
