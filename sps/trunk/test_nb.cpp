#include "spectrum.h"
#include "batch.h"
#include "filters.h"
#include "projectionutils.h"
#include "PeptideSpectrumMatch.h"
#include "ParameterList.h"

#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;
using namespace specnets;


int main(int argc, char **argv) {
  SpecSet specs;

  specs.loadPklBin(argv[1]);
  cout << "Number of spectra in " << argv[1] << ": " << specs.size() << endl;

  //  specs.SaveSpecSet_mgf(argv[2]);

}
