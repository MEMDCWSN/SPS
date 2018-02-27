//#include "batch.h"
#include "spectrum.h"
#include "mzxml.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace specnets;

int main(int argc, char **argv)
{
  SpecSet specs;
  vector<TwoValues<unsigned int> > specsInfo; // Scan num (pos.0), msLevel (pos.1)
  string outName;

  if (argc < 3)
  {
    cerr << "Syntax: " << argv[0]
        << " <file_type> <input_file_name> <output_file_name_prefix>" << endl;
    cerr
        << "        file_type: mgf, ms2, pkl (multiple spectra per file), prms (pepnovo_prm output format)\n";
    return -1;
  }
  for (unsigned int index = 0; argv[1][index] != (char)0; index++)
    argv[1][index] = tolower(argv[1][index]);

  char *outFN;
  unsigned int fnPos;
  if (argc <= 3)
  {
    outFN = (char *)malloc(strlen(argv[2]) + 1);
    strcpy(outFN, argv[2]);
    for (fnPos = strlen(outFN) - 1; fnPos > 0 and outFN[fnPos] != '.'; fnPos--)
      ;
    if (fnPos > 0)
      outFN[fnPos] = 0;
  }
  else
  {
    outFN = (char *)malloc(strlen(argv[3]) + 1);
    strcpy(outFN, argv[3]);
  }

  bool loadOk = false;
  if (strcmp(argv[1], "mgf") == 0)
  {
    if (specs.LoadSpecSet_mgf(argv[2]) <= 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = 2;
    }
    loadOk = true;
  }

  if (strcmp(argv[1], "ms2") == 0)
  {
    if (specs.LoadSpecSet_ms2(argv[2]) <= 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = 2;
    }
    loadOk = true;
  }

  if (strcmp(argv[1], "pkl") == 0)
  {
    if (specs.LoadSpecSet_pkl(argv[2]) <= 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = 2;
    }
    loadOk = true;
  }

  if (strcmp(argv[1], "prms") == 0)
  {
    if (specs.LoadSpecSet_prms(argv[2]) <= 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = 2;
    }
    loadOk = true;
  }

  if (strcmp(argv[1], "prmsv3") == 0)
  {
    if (specs.LoadSpecSet_prmsv3(argv[2]) <= 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = 2;
    }
    loadOk = true;
  }

  vector<short> msLevel;
  if (strcmp(argv[1], "mzxml") == 0)
  {
    if (LoadMzxml(argv[2], specs, &msLevel, 2) == 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = msLevel[i];
    }
    loadOk = true;
    /*  Extra code to extract CID spectra from CID/ETD runs
     SpecSet cid(specs.size());
     unsigned int cidCount=0, specIdx;
     for(specIdx=0; specIdx<specs.size(); specIdx++)
     if(specs[specIdx].msFragType==Spectrum::FragType_CID) cid[cidCount++]=specs[specIdx];
     cid.resize(cidCount);
     outName = string("cid_")+string(outFN)+string(".pklbin");
     cid.SaveSpecSet_pklbin(outName.c_str());
     outName = string("cid_")+string(outFN)+string(".mgf");
     cid.SaveSpecSet_mgf(outName.c_str());
     */
  }

  if (strcmp(argv[1], "mzxmlall") == 0)
  {
    if (LoadMzxml(argv[2], specs, &msLevel, 1) == 0)
    {
      cerr << "ERROR loading " << argv[2] << "!\n";
      return -1;
    }
    specsInfo.resize(specs.size());
    for (unsigned int i = 0; i < specsInfo.size(); i++)
    {
      specsInfo[i][0] = specs[i].scan;
      specsInfo[i][1] = msLevel[i];
    }
    loadOk = true;
  }

  if (loadOk)
  {
    outName = string(outFN) + string(".pklbin");
    specs.savePklBin(outName.c_str());
    //outName = string(outFN) + string(".bin");
    //		if(argc==5) Save_binArray<unsigned int>(argv[4],specsInfo);
    // Save_binArray<unsigned int> (outName.c_str(), specsInfo);
  }
  else
  {
    cerr
        << "ERROR parsing file_type; valid options are ms2, mgf, mzxml, pkl or prms.\n";
    return -1;
  }

  free(outFN);
  return (0);
}
