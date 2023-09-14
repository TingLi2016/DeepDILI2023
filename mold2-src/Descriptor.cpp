#include "stdafx.h"
#include "DistanceMatrix.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "Molecule.h"
#include "Descriptor.h"
#include "SDfile.h"
#include "Features.h"
#include "Calculator5.h"
#include "Calculator6.h"
#include "Calculator3.h"
#include "Calculator4.h"
#include "Calculator2.h"
#include "Calculator7.h"
#include "Calculator8.h"

Descriptor::Descriptor(void)
: m_count(0){}

Descriptor::~Descriptor(void){
}

void Descriptor::printDescriptorHeaderLine(ostream & out) {
  out << "id";
	for(short i = 1; i <= TOTALFEATURES; i++)
		out << "\tD" << setw(3) << setfill('0') << i;
	out << endl;
}


void Descriptor::printFeatures(ostream & out, string& autoID, Molecule & noHmol)
{
  out << noHmol.getHeaderMoleculeId() << "\t";
	noHmol.printFeatures(out);
}

void Descriptor::reportNewMolecule(ostream & report)
{
	report << "Processing molecule " << m_count << "." << endl;
}

void Descriptor::reportException(ostream & report, string & msg)
{
	report << "    " << msg << endl;
}

void Descriptor::printReportHeader(ostream & report, string inputFile, string& outputFile, string reportFile, string  withHFile, string withoutHFile)
{
	report << "The start time is " << ctime(&m_start) << endl;
	report << "Program input file: " << inputFile << endl;
	report << "Output file: " << outputFile << endl;

	report << "Report file: " << reportFile << endl;

	if(withHFile.size() > 0)
		report << "H-added file: " << withHFile << endl;
	if(withoutHFile.size() > 0)
		report << "H-deleted file: " << withoutHFile << endl;
	report << endl << "Processing descriptor names." << endl;
}

void Descriptor::printUsage(ostream & out)
{
	for(short i = 0; i < 70; i++)
		cout << ("-");
	cout << endl << endl;
	cout << "Generate molecular descriptors" << endl;
	cout << "Usage: descriptors -i input.sdf [-o output.txt] [-r report.txt] [-h withH.sdf] [-d withoutH.sdf]" << endl;
	cout << endl;
	cout << "Notes:" << endl;
	cout << "   1. -i file: The input data file must be a SD file: *.sdf." << endl;
	cout << "   2. -o file: Results file [default: dexcriptors.txt]." << endl;
	cout << "   3. -r file: Report file [default: report.txt]."  << endl;
	cout << "   4. -h file: print out hydrogen atoms added SD file [default:no]" << endl;
	cout << "   5. -d file: print out hydrogen atoms deleted SD file [default:no]" << endl;
	cout << "   6. Terms in [] can be omitted." << endl;
	cout << "   7. Input is case senesitive." << endl << endl;
	cout << "National Center for Toxicological Research" << endl;
	cout << "US Food and Drug Administration" << endl;
	cout << endl;
	for(short i = 0; i < 70; i++) cout << "-";
	cout << endl;
	cout << "Press any key ..." << endl;
	cin.ignore();
}

void Descriptor::finishCalculation(ostream & report)
{
	report << "Program finished normally" << endl;
}

void Descriptor::DoIt(istream & in, ostream & out) {
	Molecule mol;
	printDescriptorHeaderLine(out);

	SDfile sdFile;

  m_count = 0;

  if(sdFile.openSDfile(in)) {

    while(sdFile.hasNext()){
      sdFile.getNextMolecule(mol);
      m_count++;

      Calculator2 cal2(mol);
      Calculator3 cal3(mol);
      Calculator4 cal4(mol);
      Calculator5 cal5(mol);
      Calculator6 cal6(mol);
      Calculator7 cal7(mol);
      Calculator8 cal8(mol);

      out << mol.getHeaderMoleculeId();
      mol.printFeatures(out);
    }
  }
  out.flush();
}

void Descriptor::printReleaseInfo(void)
{
	cout << "TODO: add release notice" << endl;
}

string Descriptor::generateAutoID(string & preFieldName, string & currentFieldName, int ID)
{
	string autoID = "";
	if(preFieldName.size() > 0){
		if(preFieldName == currentFieldName){
			autoID = to_string(ID);
		}else{
			autoID = " ";
		}
	}
	return autoID;
}