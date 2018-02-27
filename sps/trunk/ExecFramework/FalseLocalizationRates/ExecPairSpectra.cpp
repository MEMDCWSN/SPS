// Header Include
#include "ExecPairSpectra.h"

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecPairSpectra::ExecPairSpectra(void) :
    ownInput(true), ownOutput(true), m_specSet(0x0), m_peptideResults(0x0)
  {
    m_name = "ExecPairSpectra";
    m_type = "ExecPairSpectra";
  }

  // -------------------------------------------------------------------------
  ExecPairSpectra::ExecPairSpectra(const ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(true), ownOutput(true), m_specSet(0x0),
        m_peptideResults(0x0)
  {
    m_name = "ExecPairSpectra";
    m_type = "ExecPairSpectra";
  }

  // -------------------------------------------------------------------------
  ExecPairSpectra::ExecPairSpectra(const ParameterList & inputParams,
                                   PeptideSpectrumMatchSet * peptideResults,
                                   SpecSet * specSet) :
    ExecBase(inputParams), ownInput(false), ownOutput(false),
        m_specSet(specSet), m_peptideResults(peptideResults)
  {
    m_name = "ExecPairSpectra";
    m_type = "ExecPairSpectra";
  }

  // -------------------------------------------------------------------------
  ExecPairSpectra::~ExecPairSpectra(void)
  {
    if (ownInput)
    {
      if (m_specSet)
      {
        delete m_specSet;
        m_specSet = 0x0;
      }
      if (m_peptideResults)
      {
        delete m_peptideResults;
        m_peptideResults = 0x0;
      }
    }

  }
  // -------------------------------------------------------------------------
  ExecBase * ExecPairSpectra::clone(const ParameterList & inputParams) const
  {
    return new ExecPairSpectra(inputParams);
  }
  // -------------------------------------------------------------------------
  bool filterPsmsByMod(PeptideSpectrumMatchSet &inputSet,
                       PeptideSpectrumMatchSet &outputSet,
                       vector<float> mod,
                       int modCount,
                       int msLevel)
  {
    DEBUG_TRACE;
    vector<psmPtr>::const_iterator it;

    for (it = inputSet.m_psmSet.begin(); it != inputSet.m_psmSet.end(); it++)
    {
      Spectrum * spectrum = (*it)->m_spectrum;

      if (spectrum == NULL)
      {
        WARN_MSG("Unable to find spectrum!");
        continue;
      }

      //make sure ms level matches
      if (msLevel > 0 && spectrum->msLevel != 0 && spectrum->msLevel != msLevel)
      {
        continue;
      }

      vector<float> modifications;
      (*it)->getModifications(modifications);

      if (modifications.size() <= modCount)
      {
        bool modMatch = true;
        for (int i = 0; i < modifications.size(); i++)
        {
          bool add = false;
          for (int j = 0; j < mod.size(); j++)
          {
            if (modifications[i] < mod[j] + .001 && modifications[i] > mod[j]
                - .001)
            {
              add = true;
            }
          }
          modMatch = (modMatch && add);
        }
        if (modMatch)
        {
          psmPtr addPsm(new PeptideSpectrumMatch);
          *addPsm = **it;
          outputSet.m_psmSet.push_back(addPsm);
        }
      }
    }
  }
  // -------------------------------------------------------------------------

  void ExecPairSpectra::separateModUnmod(PeptideSpectrumMatchSet &modifiedSpectra,
                                         SpectrumLibrary &unmodifiedSpectra,
                                         int unmodMsLevel,
                                         bool isUniqueUnmodified)
  {
    vector<psmPtr>::iterator jt;
    for (jt = m_peptideResults->m_psmSet.begin(); jt
        != m_peptideResults->m_psmSet.end(); jt++)
    {
      PeptideSpectrumMatch currResult = **jt;
      if (currResult.m_spectrum == NULL)
      {
        DEBUG_VAR(currResult.m_spectrumFile);
        DEBUG_VAR(currResult.m_scanNum);
        DEBUG_TRACE;
      }
      else
      {
        vector<float> modifications;
        currResult.getModifications(modifications);
        if (modifications.size() == 0)
        {
          if (currResult.m_spectrum->msLevel == unmodMsLevel)
          {
            if (isUniqueUnmodified)
            {
              unmodifiedSpectra.addToLibrary(*(currResult.m_spectrum),
                                             currResult,
                                             currResult.m_score,
                                             false);
            }
            else
            {
              unmodifiedSpectra.addToLibrary(*(currResult.m_spectrum),
                                             currResult);
            }
          }
        }
        else
        {
          bool unmod = true;
          for (int i = 0; i < modifications.size(); i++)
          {
            bool currUnmod = false;

            for (int j = 0; j < m_fixedMods.size(); j++)
            {

              if (modifications[i] == m_fixedMods[j])
              {
                currUnmod = true;
              }
            }
            unmod = unmod && currUnmod;
          }

          if (unmod == true)
          {
            //build unmodified
            stringstream key;
            string unmodified;
            currResult.getUnmodifiedPeptide(unmodified);

            key << unmodified << '_' << currResult.m_charge;

            string newKey = key.str();

            if (isUniqueUnmodified)
            {
              unmodifiedSpectra.addToLibrary(*(currResult.m_spectrum),
                                             currResult,
                                             currResult.m_score,
                                             newKey,
                                             false);
            }
            else
            {
              unmodifiedSpectra.addToLibrary(*(currResult.m_spectrum),
                                             currResult,
                                             newKey);
            }
          }
          else
          {
            psmPtr newUnmod = *jt;
            modifiedSpectra.m_psmSet.push_back(newUnmod);
          }
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  bool compareScore(psmPtr i, psmPtr j)
  {
    return i->m_score > j->m_score;
  }
  // -------------------------------------------------------------------------
  void filterToTopMatch(PeptideSpectrumMatchSet &psms)
  {
    std::tr1::unordered_set < string > peptides;

    sort(psms.m_psmSet.begin(), psms.m_psmSet.end(), compareScore);

    PeptideSpectrumMatchSet tempPsms;

    for (int i = 0; i < psms.size(); i++)
    {
      stringstream key;
      string unmodifiedAnnotation;
      psms.m_psmSet[i]->getUnmodifiedPeptide(unmodifiedAnnotation);

      key << unmodifiedAnnotation << "_" << psms.m_psmSet[i]->m_charge;

      if (peptides.find(key.str()) != peptides.end())
      {
        //ignore
      }
      else
      {
        peptides.insert(key.str());
        tempPsms.m_psmSet.push_back(psms.m_psmSet[i]);
      }
    }
    psms = tempPsms;
  }
  // -------------------------------------------------------------------------
  bool ExecPairSpectra::invoke(void)
  {
    PeptideSpectrumMatchSet modifiedPsms;
    SpectrumLibrary unmodifiedPsms;

    int keepCharge = 0;
    if (m_params.exists("CHARGE"))
    {
      keepCharge = m_params.getValueInt("CHARGE");
    }

    if (m_params.exists("FIXED_POSITION_MODS")) //if we have modifications which we don't allow to float
    {
      string fixedPositionString = m_params.getValue("FIXED_POSITION_MODS");
      vector < string > fixedMods;
      stringSplit(fixedPositionString, fixedMods, ",");
      for (int i = 0; i < fixedMods.size(); i++)
      {
        float fixedMod;
        stringstream ss;
        ss << fixedMods[i];
        ss >> fixedMod;
        m_fixedMods.push_back(fixedMod);
      }
      DEBUG_VAR(m_fixedMods[0]);

    }

    int unmodMsLevel = 2;
    bool isUniqueUnmodified = true; /* by default, we pick the top scoring MS2 spectra for the unmodified hits */

    if (m_params.exists("UNMODIFIED_MS_LEVEL"))
    {
      unmodMsLevel = m_params.getValueInt("UNMODIFIED_MS_LEVEL");
    }

    if (m_params.exists("UNIQUE_UNMODIFIED_PEPTIDES"))
    {
      isUniqueUnmodified = m_params.getValueBool("UNIQUE_UNMODIFIED_PEPTIDES");
    }

    separateModUnmod(modifiedPsms,
                     unmodifiedPsms,
                     unmodMsLevel,
                     isUniqueUnmodified);

    if (m_params.exists("PHOSPHO_MOD") && m_params.getValueBool("PHOSPHO_MOD"))
    {
      if (m_params.exists("MODIFICATIONS") || m_params.exists("MASS_SHIFTS"))
      {
        WARN_MSG("Modifications or mass shifts defined, both will be ignored!");
      }
      m_modifications.push_back(-18.010564686);
      m_modifications.push_back(79.966);

      vector<float> massShifts;
      massShifts.push_back(-18.010564686);
      massShifts.push_back(79.966);
      m_massShifts.push_back(massShifts);
      m_massShifts.push_back(massShifts);
    }
    else
    {
      if (m_params.exists("MODIFICATIONS"))
       {
         string modString = m_params.getValue("MODIFICATIONS");
         vector < string > modifications;
         stringSplit(modString, modifications, ",");
         for (int i = 0; i < modifications.size(); i++)
         {
           float modification;
           stringstream ss;
           ss << modifications[i];
           ss >> modification;
           m_modifications.push_back(modification);
         }
         DEBUG_VAR(m_modifications[0]);
       }

       DEBUG_VAR(m_modifications.size());

       if (m_params.exists("MASS_SHIFTS"))
       {
         string massShiftString = m_params.getValue("MASS_SHIFTS");
         vector < string > massShifts;
         stringSplit(massShiftString, massShifts, ",");
         for (int i = 0; i < massShifts.size(); i++)
         {
           vector<string> currMassShift;
           stringSplit(massShifts[i],currMassShift,":");
           vector<float> currMassShiftFloat;
           for (int j = 0; j < currMassShift.size(); j++)
           {
             float massShift;
             stringstream ss;
             ss << currMassShift[j];
             ss >> massShift;
             currMassShiftFloat.push_back(massShift);
           }
           m_massShifts.push_back(currMassShiftFloat);
         }
         DEBUG_VAR(m_massShifts[0][0]);
       }
       DEBUG_VAR(m_massShifts.size());
    }

    PeptideSpectrumMatchSet filteredPsms;

    int modMsLevel = 2;
    int modCount = 1;

    if (m_params.exists("MODIFIED_MS_LEVEL"))
    {
      modMsLevel = m_params.getValueInt("MODIFIED_MS_LEVEL");
    }

    if (m_params.exists("MAX_MOD_COUNT"))
    {
      modCount = m_params.getValueInt("MAX_MOD_COUNT");
    }

    filterPsmsByMod(modifiedPsms,
                    filteredPsms,
                    m_modifications,
                    modCount,
                    modMsLevel);

    if (m_params.exists("UNIQUE_MODIFIED_PEPTIDES"))
    {
      if (m_params.getValueBool("UNIQUE_MODIFIED_PEPTIDES"))
      {
        filterToTopMatch(filteredPsms);

        DEBUG_VAR(modifiedPsms.size());
      }
    }

    DEBUG_VAR(m_modifications.size());
    DEBUG_VAR(m_massShifts.size());

    DEBUG_VAR(filteredPsms.size());

    PeptideSpectrumMatch psm;

    DEBUG_VAR(unmodifiedPsms.size());

    //delimited file to contain spectrum pair information
    string outputLpInfo = "outputLp.txt";

    string filteredSpectraInfo = "skippedLP.txt";

    if (m_params.exists("OUTPUT_LP_INFORMATION"))
    {
      outputLpInfo = m_params.getValue("OUTPUT_LP_INFORMATION");
    }

    if (m_params.exists("FILTERED_SPECTRA_INFORMATION"))
    {
      filteredSpectraInfo = m_params.getValue("FILTERED_SPECTRA_INFORMATION");
    }

    ofstream outputLpHandle(outputLpInfo.c_str(), ios::binary);

    ofstream filteredSpectraHandle(filteredSpectraInfo.c_str(), ios::binary);

    if (!outputLpHandle.is_open() || !outputLpHandle.good())
    {
      ERROR_MSG("Unable to write to file! " << outputLpInfo);
      return false;
    }

    outputLpHandle
        << "unmodifiedMgf\tmodifiedMgf\tunmodifiedAnnotation\tlpInputFilename\tlpOutputFilename\tallMassMgf\tvariantOutputFile\tfileName\tscanNum"
        << endl;

    filteredSpectraHandle << "#SpectrumFile\tScan#\tAnnotation\tReason" << endl;

    /*get parent ion names because I am lazy. Only used when normalizing the spectrum
     */

    vector < string > parentIonNames;
    m_modelLP.getParentMassIonNames(parentIonNames);

    for (int i = 0; i < filteredPsms.size(); i++)
    {

      if (filteredPsms[i]->m_charge <= 0)
      {
        filteredPsms[i]->setChargeByAnnotation();
      }


      //only include single charge
       if (keepCharge > 0 && filteredPsms[i]->m_charge != keepCharge)
       {
         filteredSpectraHandle << filteredPsms[i]->m_spectrumFile << "\t"
             << filteredPsms[i]->m_scanNum << "\t"
             << filteredPsms[i]->m_annotation << "\t" << "Charge" << endl;
         continue;
       }


      string unmodAnnotation;
      filteredPsms[i]->getUnmodifiedPeptide(unmodAnnotation);
      DEBUG_VAR(unmodAnnotation);
      DEBUG_VAR(filteredPsms[i]->m_annotation);
      vector<float> modifications;
      filteredPsms[i]->getModifications(modifications);
      vector<vector<float> > massShifts(modifications.size());
      for (int modi = 0; modi < m_modifications.size(); modi++)
      {
        for (int modj = 0; modj < modifications.size(); modj++)
        {
          if (m_modifications[modi] < modifications[modj] + .001
              && m_modifications[modi] > modifications[modj] - .001)
          {
            massShifts[modj] = m_massShifts[modi];
          }
        }
      }

      float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

      //output LP
      AAJumps jumps(1);

      if (m_params.exists("FILTER_EXPLAINED_INTENSITY"))
      {
        string ionInclude = "all";

        PeptideSpectrumMatch expPsm = *(filteredPsms[i]);
        DEBUG_VAR(expPsm.m_annotation);
        DEBUG_VAR(expPsm.m_spectrum->size());
        DEBUG_VAR(unmodAnnotation);
        expPsm.m_annotation = unmodAnnotation;
        expPsm.annotate(unmodAnnotation,
                        ionInclude,
                        m_model,
                        0,
                        0,
                        peakTol,
                        jumps);

        DEBUG_TRACE;

        float explainedIntensity =
            SpectrumAnnotStatistics::percentExplainedIntensity(expPsm,
                                                               ionInclude);
        DEBUG_VAR(explainedIntensity);
        if (explainedIntensity
            < m_params.getValueFloat("FILTER_EXPLAINED_INTENSITY"))
        {
          filteredSpectraHandle << filteredPsms[i]->m_spectrumFile << "\t"
              << filteredPsms[i]->m_scanNum << "\t"
              << filteredPsms[i]->m_annotation << "\t" << "explainedIntensity"
              << endl;
          continue;
        }
      }
      vector<psmPtr> unmodifiedPsmMatches;

      if (unmodifiedPsms.getSpectrumLibraryMatches(unmodAnnotation,
                                                   filteredPsms[i]->m_charge,
                                                   unmodifiedPsmMatches))
      {
        DEBUG_MSG("found " << filteredPsms[i]->m_annotation << " "
            << (*unmodifiedPsmMatches[0]).m_annotation << " " << (*unmodifiedPsmMatches[0]).m_charge);

        psmPtr unmodifiedPsm;

        if (!SpectrumLibrary::pickBestHit(filteredPsms[i],
                                          unmodifiedPsmMatches,
                                          m_extractedIons,
                                          m_model,
                                          unmodifiedPsm))
        {
          DEBUG_MSG("Unable to choose best hit!");
          filteredSpectraHandle << filteredPsms[i]->m_spectrumFile << "\t"
              << filteredPsms[i]->m_scanNum << "\t"
              << filteredPsms[i]->m_annotation << "\t" << "noUnmodBestPick"
              << endl;
          continue;
        }

        //build path
        stringstream path;
        if (m_params.exists("OUTPUT_LP_DIRECTORY"))
        {
          stringstream num;
          num << setw(4) << setfill('0');
          num << i;
          path << m_params.getValue("OUTPUT_LP_DIRECTORY") << "/" << num.str()
              << "_" << filteredPsms[i]->m_annotation;
        }
        else
        {
          stringstream num;
          num << setw(4) << setfill('0');
          num << i;
          path << num.str() << "_" << filteredPsms[i]->m_annotation;
        }

        string lpInputFilename = path.str();
        lpInputFilename.append(".lp");
        string lpOutputFilename = path.str();
        lpOutputFilename.append(".txt");
        string mgfOutputFilename = path.str();
        mgfOutputFilename.append(".mgf");
        string lpResultsFilename = path.str();
        lpResultsFilename.append(".txt");
        string variantOutputFile = path.str();
        variantOutputFile.append("_variantOutput.txt");

        FalseLocalizationRates flr(jumps, m_modelLP, peakTol);
        DEBUG_VAR(m_massShifts.size());

        Spectrum * unmodifiedSpectrum = unmodifiedPsm->m_spectrum;
        Spectrum * modifiedSpectrum = filteredPsms[i]->m_spectrum;

        //normalize
        if (m_params.exists("SQRT_NORMALIZE"))
        {
          for (int j = 0; j < unmodifiedSpectrum->size(); j++)
            (*unmodifiedSpectrum)[j][1] = sqrt((*unmodifiedSpectrum)[j][1]);
          for (int j = 0; j < modifiedSpectrum->size(); j++)
            (*modifiedSpectrum)[j][1] = sqrt((*modifiedSpectrum)[j][1]);
        }

        if (m_params.exists("NORMALIZE_SPECTRUM"))
        {
          if (m_params.getValueBool("NORMALIZE_SPECTRUM"))
          {
            PeptideSpectrumMatch unmodifiedTemp;
            PeptideSpectrumMatch modifiedTemp;
            unmodifiedTemp.m_charge = unmodifiedPsm->m_charge;
            modifiedTemp.m_charge = filteredPsms[i]->m_charge;
            unmodifiedTemp.m_spectrum = unmodifiedSpectrum;
            modifiedTemp.m_spectrum = modifiedSpectrum;

            string includeIons;
            stringJoin(includeIons,parentIonNames,",");

            unmodifiedTemp.annotate(unmodifiedPsm->m_annotation,
                                    includeIons,
                                    m_modelLP,
                                    0,
                                    0,
                                    peakTol,
                                    false);

            DEBUG_VAR(filteredPsms[i]->m_annotation);
            modifiedTemp.annotate(filteredPsms[i]->m_annotation,
                                  includeIons,
                                  m_modelLP,
                                  0,
                                  0,
                                  peakTol,
                                  false);

            // remove parent mass ions
            for (int j = unmodifiedSpectrum->size() -1; j >= 0; j--)
            {
              if (unmodifiedTemp.m_peakAnnotations[j].first)
              {
                unmodifiedSpectrum->removePeak(j);
              }
            }

            // remove parent mass ions
            for (int j = modifiedSpectrum->size() -1; j >= 0; j--)
            {
              if (modifiedTemp.m_peakAnnotations[j].first)
              {
                modifiedSpectrum->removePeak(j);
              }
            }

            unmodifiedSpectrum->normalize(1000000);
            modifiedSpectrum->normalize(1000000);
          }
        }

        //only normalize non P peaks
        /*expectedIntensitiesPSM      .annotate(unmodifiedPeptide,
         ,
         m_modelLP,
         0,
         0,
         peakTol,
         false);

         */
        if (m_params.exists("PHOSPHO_MOD") && m_params.getValueBool("PHOSPHO_MOD"))
        {
          flr.generateFlrLPPhospho(*unmodifiedSpectrum,
                                    *modifiedSpectrum,
                                    unmodifiedPsm->m_annotation,
                                    modifications,
                                    massShifts,
                                    lpInputFilename.c_str(),
                                    mgfOutputFilename.c_str(),
                                    variantOutputFile.c_str());
        }
        else
        {
          flr.generateFlrLP(*unmodifiedSpectrum,
                          *modifiedSpectrum,
                          unmodifiedPsm->m_annotation,
                          modifications,
                          massShifts,
                          lpInputFilename.c_str(),
                          mgfOutputFilename.c_str(),
                          variantOutputFile.c_str());
        }

        string unmodifiedMgfFile = path.str();
        unmodifiedMgfFile.append("_unmod.mgf");
        string modMgfFile = path.str();
        modMgfFile.append("_mod.mgf");

        SpecSet specs;
        specs.resize(1);

        specs[0] = *(filteredPsms[i]->m_spectrum);
        DEBUG_VAR(specs[0].size());
        specs.SaveSpecSet_mgf(modMgfFile.c_str());
        specs[0] = *(unmodifiedPsm->m_spectrum);
        DEBUG_VAR(specs[0].size());
        specs.SaveSpecSet_mgf(unmodifiedMgfFile.c_str());

        outputLpHandle << unmodifiedMgfFile << "\t" << modMgfFile << "\t"
            << unmodifiedPsm->m_annotation << "\t" << lpInputFilename << "\t" << lpResultsFilename << "\t"
            << mgfOutputFilename << "\t" << variantOutputFile << "\t"
            << filteredPsms[i]->m_spectrumFile << "\t"
            << filteredPsms[i]->m_scanNum << endl;

      }
      else
      {
        DEBUG_MSG("Not found " << filteredPsms[i]->m_annotation << " " << unmodAnnotation <<
            " " << filteredPsms[i]->m_charge);
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::saveOutputData(void)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::loadInspectSpectrum(const char * specSetListFile,
                                            PeptideSpectrumMatchSet &psms,
                                            SpecSet &outputResults)
  {
    vector < vector<string> > spectraFileList;
    map<string, unsigned int> spectraFileListHeader;
    vector < string > requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(specSetListFile,
                                                "\t",
                                                "#",
                                                spectraFileListHeader,
                                                spectraFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load spectrum files!" << specSetListFile);
      return false;
    }

    vector<int> outputResultIndex(psms.m_psmSet.size(), -1);
    std::tr1::unordered_set < string > passedSpectra;
    outputResults.resize(psms.m_psmSet.size());

    for (int i = 0; i < spectraFileList.size(); i++)
    {
      //get file path
      string path = spectraFileList[i][requiredHeaderIndex[0]];
      DEBUG_VAR(path);

      FilenameManager spectrumFm(path.c_str());

      DEBUG_VAR(spectrumFm.extension);
      DEBUG_VAR(spectrumFm.filename);

      SpecSet currSet(1);

      if (spectrumFm.extension.compare("mzxml") == 0)
      {
        ERROR_MSG("Mzxml format not supported! Please convert to pklbin.");
        return false;
      }

      if (!currSet.Load(path.c_str()))
      {
        ERROR_MSG("Unable to load file! " << path);
        return false;
      }

      currSet.setFilename(spectrumFm.filename);
      DEBUG_VAR(currSet[i].fileName);

      if (spectrumFm.extension.compare("mgf") == 0)
      {
        //ONLY BECAUSE INSPECT DOES NOT DEAL WITH SCAN NUMBERS!
        //REMOVE ME!
        for (int k = 0; k < currSet.size(); k++)
        {
          currSet[k].scan = k+1;
          currSet[k].msLevel = 2;
        }
      }



      currSet.index();

      for (int j = 0; j < psms.m_psmSet.size(); j++)
      {
        psmPtr currPsm = psms.m_psmSet[j];
        string psmFilename;

        FilenameManager psmFm(currPsm->m_spectrumFile.c_str());

        stringstream key;
        key.str("");
        key.clear();
        key << psms.m_psmSet[j]->m_scanNum << "$" << psmFm.filename;

        if (psmFm.filename.compare(spectrumFm.filename) == 0)
        {
          Spectrum * currScan = currSet.getIndex(key.str().c_str());

          if (currScan != NULL)
          {
            outputResults[j] = *currScan;
            outputResultIndex[j] = outputResults.size() - 1;
          }
          else
          {
            DEBUG_VAR(currPsm->m_scanNum);
            DEBUG_TRACE;
          }
        }
      }
    }

    for (int i = 0; i < outputResultIndex.size(); i++)
    {
      if (outputResultIndex[i] > -1)
      {
        psms.m_psmSet[i]->m_spectrum = &(outputResults[i]);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::loadInputData(void)
  {
    if (m_params.exists("OUTPUT_LP_DIRECTORY"))
    {
      if (!mkdirRecurse(m_params.getValue("OUTPUT_LP_DIRECTORY")))
      {
        ERROR_MSG("UNABLE TO GENERATE DIRECTORY");
        return false;
      }
    }

    if (m_peptideResults == NULL)
    {
      DEBUG_TRACE;
      m_peptideResults = new PeptideSpectrumMatchSet; //! Matches each spectrum to its associated AAJump and PSM spectrum.
      m_specSet = new SpecSet; //! Spectra we are looking at
    }

    DEBUG_VAR(m_peptideResults->size());
    DEBUG_VAR(m_specSet->size());

    if (m_peptideResults->size() == 0 && m_specSet->size() == 0)
    {
      if (!loadPsmResults(m_params,
                     *m_peptideResults))
      {
        return false;
      }

      if (m_params.exists("INPUT_SPECTRA_FILE_LIST"))
      {
        if (!loadInspectSpectrum(m_params.getValue("INPUT_SPECTRA_FILE_LIST").c_str(),
                                 *m_peptideResults,
                                 *m_specSet))
        {
          ERROR_MSG("Unable to load spectrum files!" << m_params.getValue("INPUT_SPECTRA_FILE_LIST"));
          return false;
        }
      }
      DEBUG_VAR(m_specSet->size());
    }

    DEBUG_TRACE;
    if (m_params.exists("MS2_SCORING_MODEL_LP"))
    {
      if (!m_modelLP.LoadModel(m_params.getValue("MS2_SCORING_MODEL_LP").c_str()))
      {
        ERROR_MSG("Unable to load model! " << m_params.getValue("MS2_SCORING_MODEL_LP"));
        return false;
      }
    }

    if (m_params.exists("MS2_SCORING_MODEL"))
    {
      if (!m_model.LoadModel(m_params.getValue("MS2_SCORING_MODEL").c_str()))
      {
        ERROR_MSG("Unable to load model! " << m_params.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }

    string inputIons = "b,b-NH3,b-H2O,b2,b2-iso,y,y-NH3,y-H2O,y2,y2-iso,P2-H2O";
    stringSplit(inputIons, m_extractedIons, ",");
    if (m_params.exists("IONS_TO_EXTRACT"))
    {
      inputIons = m_params.getValue("IONS_TO_EXTRACT");
      stringSplit(inputIons, m_extractedIons, ",");
    }

    DEBUG_VAR(m_extractedIons.size());

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::loadOutputData(void)
  {
    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPairSpectra::split(int numSplit)
  {
    //void
  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::merge(void)
  {
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPairSpectra::validateStatisticsParams(std::string & error)
  {
    m_isValid = false;
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    m_isValid = true;
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecPairSpectra::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("MS2_SCORING_MODEL_LP");
    VALIDATE_PARAM_EXIST("MS2_SCORING_MODEL");
    VALIDATE_PARAM_EXIST("OUTPUT_LP_DIRECTORY");
    VALIDATE_PARAM_EXIST("OUTPUT_LP_INFORMATION");

    m_isValid = true;
    return true;
  }

}
