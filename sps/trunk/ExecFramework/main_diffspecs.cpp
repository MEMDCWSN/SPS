//
//  main_diffspecs - stand alone spectrum file comparator
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  DEBUG_TRACE;
  if (argc != 4 && argc != 5)
  {
    cerr << "Usage: main_diffspecs type specfile1 specfile2 [param_file]" << endl;
    cerr << "       Valid types are: mgf, prms, specset, specpairset, pklbin" << endl;
    return -1;
  }
  
  DEBUG_MSG("Comparing " << argv[2] << " and " << argv[3]);
  
  bool specs_are_different = false;
  
  ParameterList ip;
  DEBUG_TRACE;
  if (argc == 5)
  {
    ip.readFromFile(argv[4]);
  }
  double tolerancePeak = 0.1;
  double toleranceParentMass = 0.1;
  if (ip.exists("TOLERANCE_PEAK")) 
  {
    tolerancePeak = getFloat(ip.getValue("TOLERANCE_PEAK").c_str());
  }

  if (ip.exists("TOLERANCE_PM")) 
  {
    toleranceParentMass = getFloat(ip.getValue("TOLERANCE_PM").c_str());
  }
  DEBUG_VAR(tolerancePeak);
  DEBUG_VAR(toleranceParentMass);

  string type = argv[1];

  if (type == "mgf" || type == "prms" || type == "specset" || type == "pklbin")
  {
    SpecSet spectra1;
    SpecSet spectra2;
    size_t size1 = 0;
    size_t size2 = 0;


    if (type == "mgf")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_mgf(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
      spectra2.LoadSpecSet_mgf(argv[3]);
      size2 = spectra2.size();
      DEBUG_VAR(size2);
    }
    else if (type == "prms")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_prmsv3(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
      spectra2.LoadSpecSet_prmsv3(argv[3]);
      size2 = spectra2.size();
      DEBUG_VAR(size2);
    }
    else if (type == "specset" || type == "pklbin")
    {
      DEBUG_TRACE;
      spectra1.loadPklBin(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
      spectra2.loadPklBin(argv[3]);
      size2 = spectra2.size();
      DEBUG_VAR(size2);
    }

    if (size1 != size2)
    {
      WARN_MSG("Spectrum Sets have different sizes [" << size1 << "] vs [" << size2 << "]");
      specs_are_different = true;
    }
    
    int minSpecSize = size1 < size2 ? size1 : size2;
    for (size_t i = 0; i < minSpecSize; i++)
    {
      //DEBUG_VAR(i);
      size_t peakSize1 = spectra1[i].size();
      //DEBUG_VAR(peakSize1);
      size_t peakSize2 = spectra2[i].size();
      //DEBUG_VAR(peakSize2);

      if (spectra1[i].scan != spectra2[i].scan)
      {
        WARN_MSG("Spectrum [" << i << "] has different scan number [" << 
                 spectra1[i].scan << "] vs [" << spectra2[i].scan << "]");
        specs_are_different = true;
      }
      
      if (spectra1[i].msLevel != spectra2[i].msLevel)
      {
        WARN_MSG("Spectrum [" << i << "] has different MS level [" << 
                 spectra1[i].msLevel << "] vs [" << spectra2[i].msLevel << "]");
        specs_are_different = true;
      }
      
      if (peakSize1 != peakSize2)
      {
        WARN_MSG("Spectrum [" << i << "] has differing number of peaks [" << 
                 peakSize1 << "] vs [" << peakSize2 << "]");
        specs_are_different = true;
      }
      
      // check parent mass tolerance
      double parentMassDiff = fabs(spectra1[i].parentMass - spectra2[i].parentMass);
      if(parentMassDiff > toleranceParentMass) 
      {
        DEBUG_VAR(parentMassDiff);
        DEBUG_VAR(toleranceParentMass);
        WARN_MSG("Spectrum [" << i << "] has differing parent mass [" << setprecision(10) << 
                 spectra1[i].parentMass << "] vs [" << spectra2[i].parentMass << "]");
        specs_are_different = true;
      }
    	
      int minPeakSize = peakSize1 < peakSize2 ? peakSize1 : peakSize2;
      for (size_t j = 0; j < minPeakSize; j++)
      {
        TwoValues<float> peak1 = spectra1[i][j];
        //DEBUG_VAR(peak1.values[0]);
        //DEBUG_VAR(peak1.values[1]);
        TwoValues<float> peak2 = spectra2[i][j];
        //DEBUG_VAR(peak2.values[0]);
        //DEBUG_VAR(peak2.values[1]);

        // Peak value difference
        double minVal1 = fabs(peak1[0] - peak2[0]);
        // Peak intensity difference
        double minVal2 = fabs(peak1[1] - peak2[1]);

        if (minVal1 > tolerancePeak)  
        {
          WARN_MSG("Spectrum [" << setprecision(10) << i << "] Peak [" << j << "] is different [" << 
                   peak1[0] << "," << peak1[1] << "] vs [" <<
                   peak2[0] << "," << peak2[1] << "]" );
          specs_are_different = true;
        }

      }  
      
    }
  }
  else if (type == "specpairset")
  {
    SpectrumPairSet specpairset1;
    SpectrumPairSet specpairset2;
    specpairset1.loadFromBinaryFile(argv[2]);
    specpairset2.loadFromBinaryFile(argv[3]);
    size_t size1 = specpairset1.size();
    DEBUG_VAR(size1);
    size_t size2 = specpairset2.size();
    DEBUG_VAR(size2);
    if (size1 != size2)
    {
      WARN_MSG("Spectrum Pair Sets have different sizes [" << size1 << "] vs [" << size2 << "]");
      specs_are_different = true;
    }

    int minSize = size1 < size2 ? size1 : size2;
    for (size_t i = 0; i < minSize; i++)
    {
      SpectrumPair pair1 = specpairset1[i];
      SpectrumPair pair2 = specpairset2[i];

      if (pair1.spec1 != pair2.spec1)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing first index [" << 
                 pair1.spec1 << "] vs [" << pair2.spec1 << "]");
        specs_are_different = true;
      }
      if (pair1.spec2 != pair2.spec2)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing second index [" << 
                 pair1.spec2 << "] vs [" << pair2.spec2 << "]");
        specs_are_different = true;
      }
      if (pair1.score1 != pair2.score1)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing first score [" << 
                 pair1.score1 << "] vs [" << pair2.score1 << "]");
        specs_are_different = true;
      }
      if (pair1.score2 != pair2.score2)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing second score [" << 
                 pair1.score2 << "] vs [" << pair2.score2 << "]");
        specs_are_different = true;
      }
      if (pair1.shift1 != pair2.shift1)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing first shift [" << 
                 pair1.shift1 << "] vs [" << pair2.shift1 << "]");
        specs_are_different = true;
      }
      if (pair1.shift2 != pair2.shift2)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing second shift [" << 
                 pair1.shift2 << "] vs [" << pair2.shift2 << "]");
        specs_are_different = true;
      }
      if (pair1.specC != pair2.specC)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing specC [" << 
                 pair1.specC << "] vs [" << pair2.specC << "]");
        specs_are_different = true;
      }
      if (pair1.spec2rev != pair2.spec2rev)
      {
        WARN_MSG("Spectrum Pair [" << i << "] has differing reverse flag [" << 
                 pair1.spec2rev << "] vs [" << pair2.spec2rev << "]");
        specs_are_different = true;
      }
    }
    
  }
  else
  {
    cerr << "Usage: main_diffspecs type specfile1 specfile2 [param_file]" << endl;
    cerr << "       Valid types are: mgf, prms, specset, specpairset" << endl;
    return -1;
  }  
  
  if(specs_are_different == true){
      return 1;
  }
  
  return 0;
}

