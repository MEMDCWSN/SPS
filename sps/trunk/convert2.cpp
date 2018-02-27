////////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "PeptideSpectrumMatchSet.h"
#include "Logger.h"

#include "PWizInterface.h"
#include "Specific2.h"
#include "mzxml.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

#define OK    1
#define ERROR 0

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
// Dump specset data (Debug)
////////////////////////////////////////////////////////////////////////////////
void dumpSpecset(SpecSet &specSet)
{
  if (!specSet.size())
  {
    cout << "Specset size is 0" << endl;
    return;
  }

  specnets::Spectrum &spectrum = specSet[0];

  if (!spectrum.size())
  {
    cout << "spectrum size is 0" << endl;
    return;
  }

  cout << "parentMassTol :       " << spectrum.parentMassTol << endl;
  cout << "parentMZ :            " << spectrum.parentMZ << endl;
  cout << "parentMass :          " << spectrum.parentMass << endl;
  cout << "parentCharge :        " << spectrum.parentCharge << endl;
  cout << "scan :                " << spectrum.scan << endl;
  cout << "msLevel :             " << spectrum.msLevel << endl;
  cout << "msFragType :          " << spectrum.msFragType << endl;
  cout << "fileName :            " << spectrum.fileName << endl;
  cout << "resolution :          " << spectrum.resolution << endl;
  cout << "instrument_name :     " << spectrum.instrument_name << endl;
  cout << "ITOL :                " << spectrum.ITOL << endl;
  cout << "ITOLU :               " << spectrum.ITOLU << endl;
  cout << "TOL :                 " << spectrum.TOL << endl;
  cout << "TOLU :                " << spectrum.TOLU << endl;
  cout << "spectrum_quality :    " << spectrum.spectrum_quality << endl;
  cout << "idDist :              " << spectrum.idDist << endl;
  cout << "retention_time :      " << spectrum.retention_time << endl;
  cout << "precursor_intensity : " << spectrum.precursor_intensity << endl;
}
////////////////////////////////////////////////////////////////////////////////
// Process spectra files
////////////////////////////////////////////////////////////////////////////////
int processSpectraFiles(string &outfile, FilenameManager &fm, bool output_scan_index_file)
{
  SpecSet specs;
  string outName;
  vector<TwoValues<unsigned int> > specsInfo; // Scan num (pos.0), msLevel (pos.1)

  // Generate output filename
  if (outfile.length() == 0)
  {
    size_t found = fm.filenameFull.find_last_of(".");
    outName = fm.filenameFull.substr(0, found);
    outName += ".pklbin";
  }
  else
  {
    outName = outfile;
  }

  // load the specset
  int loadOk = 0; //specs.Load(fm.filenameFull.c_str());

  int msl = 2;
  if (fm.extension == "mzxmlall")
  {
    fm.extension = "mzxml";
    msl = 1;
  }
  DEBUG_MSG("Extension: " << fm.extension);

  // load mzxml file
  
  string extension_lower = fm.extension;
  std::transform(extension_lower.begin(),
                 extension_lower.end(),
                 extension_lower.begin(),
                 ::tolower);
  if (extension_lower.compare("mzxml") == 0)
  {
    // auxilizary vector needed
    vector<short> msLevel;
    // hold the return value
    int ret;
    // load method depends on the presence of spectrumscan specification
    try
    {
      DEBUG_MSG("Loading mzXML file with specnets code");
      loadOk = (int)LoadMzxml((char * const )(fm.filenameFull.c_str()),
                              fm.filenameFull,
                              specs,
                              &msLevel,
                              msl);
    }
    catch (...)
    {
      ERROR_MSG("Failed to convert mzXML file using specnets methods");
      loadOk = 0;
    }
  }

  if (loadOk == 0)
  {
    // load other formats
    loadOk = (int)specs.Load(fm.filenameFull.c_str(), fm.extension.c_str());

    // load using pwiz
    if (loadOk == 0)
    {
      DEBUG_MSG("Loading data using pwiz");
      PWizInterface pwiz;
      loadOk = pwiz.loadDataUsingPWiz(fm.filenameFull, specs, msl);

      if (loadOk == -2)
      {
        ERROR_MSG("PWiz libraries are not available.");
        return ERROR;
      }

    }

    if (loadOk == 0)
    {
      stringstream err;
      ERROR_MSG("Error loading file: " << fm.filenameFull);
      return ERROR;
    }
  }

  /*
   cout << "==================================================" << endl;
   for(int i = 0 ; i < specs.size() ; i++) {
   cout << "------------ spectrum " << i << "----------" << endl;
   for(int j = 0 ; j < specs[i].size() ; j++) {
   cout << j << " (" << specs[i][j][0] << " ; " << specs[i][j][1] << ")" << endl;
   }
   }
   */

  if (loadOk)
  {
      if(output_scan_index_file){
          //Saving as output the file name scan and index
          ofstream myfile;
          myfile.open (outName.c_str());
          for(int i = 0; i < specs.size(); i++){
              if(specs[i].scan == 4294967295){
                myfile << fm.filenameFull <<"\t"<<"-1"<<"\t"<<i+1<<endl;
              }
              else{
                myfile << fm.filenameFull <<"\t"<<specs[i].scan<<"\t"<<i+1<<endl;
              }
              
          }
          myfile.close();
      }
      else{
          specs.savePklBin(outName.c_str());
      }
  }
  else
  {
    ERROR_MSG("ERROR parsing file_type.\n");
    ERROR_MSG("Valid types for MS files are: ms2, mgf, mzml, mzxml, mzxmlall, pkl or prms");
    ERROR_MSG("Valid types for identity files are: pepxml, mzidentml");
    return ERROR;
  }

  //dumpSpecset(specs);
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Process peptide files
////////////////////////////////////////////////////////////////////////////////
int processPeptideFiles(string &outfile,
                        FilenameManager &fm,
                        string &modsFile,
                        string &aasFile,
                        string &argv0,
                        bool outputFixed,
                        bool dump)
{
  PeptideSpectrumMatchSet psmSet;
  string outName;

  //cout << "Entering processPeptideFiles()" << endl;

  // Generate output filename
  if (outfile.length() == 0)
  {
    size_t found = fm.filenameFull.find_last_of(".");
    outName = fm.filenameFull.substr(0, found);
    outName += ".psm";
  }
  else
  {
    outName = outfile;
  }

  // load the standard AAs file
  FilenameManager fmStd(argv0.c_str());
  string stdFn = fmStd.path;
  if (stdFn.length())
    stdFn += '/';
  stdFn += "AA_standard.txt";
  AAJumps jumpsStd(1);
  jumpsStd.loadJumps(stdFn.c_str(), false);

  // load the AAs file, if specified
  AAJumps jumps(1);
  AAJumps *AAs = NULL;
  if (aasFile.length())
  {
    jumps.loadJumps(aasFile.c_str(), true);
    AAs = &jumps;
  }

  // define the structure for the search mods
  vector<SearchModsData> searchModsData;

  // load the file
  int loadOk = 0;

  if (loadOk == 0)
  {

    // load using pwiz
    if (loadOk == 0)
    {
      PWizInterface pwiz;
      loadOk = pwiz.openIdent(fm.filenameFull);
      if (loadOk)
      {
        if (dump)
          pwiz.dump();
        pwiz.acquaireMods(searchModsData);
        for (int i = 0; i < searchModsData.size(); i++)
          searchModsData[i].remove(jumpsStd, AAs);
        loadOk = pwiz.loadDataUsingPWiz(psmSet, searchModsData, outputFixed);
      }
    }

    if (loadOk == -2)
    {
      stringstream err;
      cerr << "PWiz libraries are not present." << endl;
      return ERROR;
    }

    if (loadOk == 0)
    {
      stringstream err;
      cerr << "Error loading file: " << fm.filenameFull << endl;
      return ERROR;
    }
  }

  if (modsFile.length() && searchModsData.size())
  {
    ofstream of(modsFile.c_str(), std::ios::out | std::ios::binary);
    if (!of)
    {
      cerr << "Error opening file: " << modsFile << endl;
      return -1;
    }

    // variable mods
    if (searchModsData[0].variable.size())
    {
      stringstream vmm, vmr;
      vmm << "VARIABLE_MODIFICATIONS=";
      vmr << "VARIABLE_MODIFICATIONS_ALLOWED_SITES=";
      for (int i = 0; i < searchModsData[0].variable.size(); i++)
      {
        if (i)
        {
          vmm << ";";
          vmr << ";";
        }
        vmm << searchModsData[0].variable[i].massShift;
        for (int j = 0; j < searchModsData[0].variable[i].residues.size(); j++)
        {
          if (j)
            vmr << ":";
          vmr << searchModsData[0].variable[i].residues[j];
        }
      }
      of << vmm.str() << endl;
      of << vmr.str() << endl;
    }

    // fixed mods
    if (searchModsData[0].fixedMods.size())
    {
      stringstream vmm, vmr;
      vmm << "FIXED_MASS_SHIFTS=";
      vmr << "FIXED_MODIFICATIONS_ALLOWED_SITES=";
      for (int i = 0; i < searchModsData[0].fixedMods.size(); i++)
      {
        if (i)
        {
          vmm << ";";
          vmr << ";";
        }
        vmm << searchModsData[0].fixedMods[i].massShift;
        for (int j = 0; j < searchModsData[0].fixedMods[i].residues.size(); j++)
        {
          if (j)
            vmr << ":";
          vmr << searchModsData[0].fixedMods[i].residues[j];
        }
      }
      of << vmm.str() << endl;
      of << vmr.str() << endl;
    }

    of.close();
  }

  if (loadOk == 1)
  {
    //cout << "Saving..." << endl;
    psmSet.saveToFile(outName.c_str());
  }
  else
  {
    cerr << "ERROR parsing file: " << fm.filenameFull << endl;
    return ERROR;
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Entry point
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  string type, modsFile, aasFile, inFile, outFile, argv0;
  FilenameManager fm;
  bool outputFixed = false;
  bool dump = false;
  
  //Options for Spectrum File Conversion
  bool output_file_scan_index = false;

  argv0 = argv[0];

  Logger::setDefaultLogger(Logger::getLogger(0));

  // add divert for segfault
  addSegFaultDivert2();

  //ERROR_MSG("Testing");

  bool showHelp = false;
  for (size_t i = 0; i < argc; i++)
  {
    string arg(argv[i]);
    if (arg.compare("--help") == 0)
      showHelp = true;
  }

  int i = 1;

  if (argc <= i)
    showHelp = true;

  if (argc > i)
  {
    string arg1(argv[i]);
    if (arg1.compare("--type") == 0)
    {
      i++;
      if (argc <= i)
        showHelp = true;
      else
        type = argv[i++];
    }
  }

  if (argc > i)
  {
    string arg1 = argv[i];
    if (arg1.compare("--mods-file") == 0)
    {
      i++;
      if (argc <= i)
        showHelp = true;
      else
        modsFile = argv[i++];
    }
  }

  if (argc > i)
  {
    string arg1 = argv[i];
    if (arg1.compare("--aas-file") == 0)
    {
      i++;
      if (argc <= i)
        showHelp = true;
      else
        aasFile = argv[i++];
    }
  }

  if (argc > i)
  {
    string arg1 = argv[i];
    if (arg1.compare("--fixed") == 0)
    {
      i++;
      outputFixed = true;
    }
  }

  if (argc > i)
  {
    string arg1 = argv[i];
    if (arg1.compare("--dump") == 0)
    {
      i++;
      dump = true;
    }
  }
  
  //Iterating through all parameters
  for(int k = 0; k < argc; k++)
  {
    string arg1 = argv[k];
    
    //Checking for output_file_scan_index
    if (arg1.compare("--output_file_scan_index") == 0)
    {
      i++;
      output_file_scan_index = true;
    }
    
  }
  

  if (argc <= i)
    showHelp = true;
  else
    inFile = argv[i++];

  if (argc > i)
    outFile = argv[i];

  if (showHelp)
  {
    cerr << "Syntax: " << argv[0]
        << " [--type <input_file_type>] [--mods-file <mods_file_name>] [--aas-file <AAs_file>] [--fixed] [--output_file_scan_index] <input_file_name> [<output_file_name>]"
        << endl;
    cerr
        << "Valid types for MS files are: ms2, mgf, mzml, mzxml, mzxmlall, pkl or prms.\n";
    cerr << "Valid types for identity files are: pepxml, mzidentml\n";
    cerr << "output_file_scan_index -- Outputs scans and indices of a given input spectrum file instead of pklbin\n";
    return -1;
  }

  fm.filenameFull = inFile;
  fm.splitFilename();
  if (type.length())
    fm.extension = type;

  fm.lowerCaseExtension();

  if (fm.extension == "pepxml" || fm.extension == "mzidentml")
    return processPeptideFiles(outFile,
                               fm,
                               modsFile,
                               aasFile,
                               argv0,
                               outputFixed,
                               dump);

  return processSpectraFiles(outFile, fm, output_file_scan_index);
}
////////////////////////////////////////////////////////////////////////////////
