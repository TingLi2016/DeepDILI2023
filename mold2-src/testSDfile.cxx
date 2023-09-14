// testSDfile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Atom.h"
#include "Bond.h"
#include "DistancePathMatrix.h"
#include "DistanceMatrix.h"
#include "Molecule.h"
#include "Descriptor.h"

int main(int argc, char* argv[])
{
	bool canGo = true;
	Descriptor descriptor;
	if(argc <= 0){
		canGo = false;
	}else{
		map<string, string> params;
		for(int i = 1; i < argc; i += 2){
			params.insert(pair<string, string>(argv[i], argv[i + 1]));
		}
		// params.insert(pair<string, string>("-i", "c:\\Documents and Settings\\zsu.FDA\\My Documents\\Visual Studio 2008\\Projects\\testSDfile\\Release\\LTKB_06022010_test1.sdf"));
		// params.insert(pair<string, string>("-i", "h:\\DescriptorsGenerator\\LTKB_Cmax.sdf"));
		// params.insert(pair<string, string>("-i", "h:\\DescriptorsGenerator\\testSDfile\\LTKB_Cmax.sdf"));
		//params.insert(pair<string, string>("-i", "h:\\DescriptorsGenerator\\testSDfile\\1771.mol"));
		//	params.insert(pair<string, string>("-i", "c:\\Documents and Settings\\zsu.FDA\\My Documents\\Visual Studio 2008\\Projects\\mold2\\Debug\\LTKB_Cmax2.sdf"));
		//	params.insert(pair<string, string>("-i", "e:\\DescriptorsGenerator\\testSDfile\\LTKB_06022010_test1.sdf"));
		// params.insert(pair<string, string>("-i", "h:\\DescriptorsGenerator\\test.mol"));

		map<string, string>::iterator it = params.find("-c");
		if(it != params.end()){
			map<string, string>::iterator it1 = params.find("-os");
			string name;
			if(it1->second == "WIN"){
				name = params.find("-o")->second;
			}else{
				name = params.find("-o")->second.substr(1);
				name.erase(name.size() - 1);
			}
			ofstream out(name.c_str());
			if(out.is_open()){
				descriptor.DoIt(cin, out);
				out.close();
			}
		}else{
			it = params.find("-i");

			if(it == params.end()){
				cout << "No input sdf file" << endl;
				canGo = false;
			}else{
				//check input sdf file
				ifstream in(it->second.c_str());
				if(!in.is_open()){
					cout << "Can not open the input file: " << it->second << endl;
					canGo = false;
				}else{
					in.close();
					string inputFileName = it->second;

					// check output file
					string outputFileName = "output.txt";
					it = params.find("-o");
					if(it != params.end()){
						outputFileName = it->second;
					}
					ofstream out(outputFileName.c_str());
					if(!out.is_open()){
						cout << "Can not open output results file: " << outputFileName << endl;
						canGo = false;
					}else{
						// check report file
						string reportFileName = "report.txt";
						it = params.find("-r");
						if(it != params.end()){
							reportFileName = it->second;
						}
						ofstream rep(reportFileName.c_str());
						if(!rep.is_open()){
							cout << "Can not open report file: " << reportFileName << endl;
							cout << "Report to standard output" << endl;
							reportFileName = "Standard output";
						}
						ostream & report = reportFileName == "Standard output" ? cout:rep;

						// check H-added output file
						string withHFile = "";
						it = params.find("-h");
						if(it != params.end()){
							withHFile = it->second;
							ofstream withH(withHFile.c_str());
							if(!withH.is_open()){
								cout << "Can not open H-added sdf: " << withHFile << endl;
								out.close();
								canGo = false;
							}else{
								withH.close();
							}
						}
						if(canGo){
							// check H-deleted output file
							string withoutHFile = "";
							it = params.find("-d");
							if(it != params.end()){
								withoutHFile = it->second;
								ofstream withoutH(withoutHFile.c_str());
								if(!withoutH.is_open()){
									cout << "Can not open H-deleted sdf: " << withHFile << endl;
									out.close();
									canGo = false;
								}else{
									withoutH.close();
								}
							}
							if(canGo){
								//everything is ok, go!
								descriptor.DoIt(out, report,inputFileName, outputFileName, reportFileName, withHFile, withoutHFile);

								// clean up opened streams
								out.close();
								in.close();
								rep.close();
								cout << "Finished! Press any key ..." << endl;
								cin.ignore();
							}
						}
					}
				}
			}
		}
	}

	if(!canGo){
		// show usage
		descriptor.printUsage(cout);
	}
	return 0;
}

