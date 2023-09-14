// Mold2Debug.cpp : Defines the entry point for the console application.
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
  try
  {
    Descriptor descriptor;
    descriptor.DoIt(cin, cout);
	  return 0;
  }
  catch(mold2Exception& ex)
  {
    cerr << "Error processing molecule \"" << ex.molid() << "\": " << ex.what() << endl;
    return 1;
  }
  catch(std::exception& ex)
  {
    cerr << "An unexpected error occurred during processing." << endl;
    return 2;
  }
}
