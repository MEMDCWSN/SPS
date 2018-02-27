// Header Include
#include "ExecFlr.h"

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecFlr::ExecFlr(void) :
    ownInput(true), ownOutput(true), m_outputLpLines(0x0),
        m_outputLpHeader(0x0), m_lpRequiredHeader(0x0),
        m_lpRequiredHeaderIndex(0x0), m_model(0x0), m_modifications(0x0),
        m_massShifts(0x0), m_modAllowedSites(0x0), m_variantGroups(0x0),
        m_groupQuantities(0x0), m_indices(0x0), m_theoCosineSimilarity(0x0),
        m_unmodCosineSimilarity(0x0), m_validVariants(0x0)
  {
    m_name = "ExecFlr";
    m_type = "ExecFlr";
  }
  // -------------------------------------------------------------------------
  ExecFlr::ExecFlr(vector<vector<string> > * outputLpLines,
                   map<string, unsigned int> * outputLpHeader,
                   vector<string> * lpRequiredHeader,
                   vector<int> * lpRequiredHeaderIndex,
                   MS2ScoringModel * model,
                   vector<float> * modifications,
                   vector<float> * massShifts,
                   vector<string> * modAllowedSites,
                   ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(false), ownOutput(true),
        m_outputLpLines(outputLpLines), m_outputLpHeader(outputLpHeader),
        m_lpRequiredHeader(lpRequiredHeader),
        m_lpRequiredHeaderIndex(lpRequiredHeaderIndex), m_model(model),
        m_modifications(modifications), m_massShifts(massShifts),
        m_modAllowedSites(modAllowedSites), m_variantGroups(0x0),
        m_groupQuantities(0x0), m_indices(0x0), m_theoCosineSimilarity(0x0),
        m_unmodCosineSimilarity(0x0), m_validVariants(0x0)
  {
    m_name = "ExecFlr";
    m_type = "ExecFlr";
  }
  // -------------------------------------------------------------------------
  ExecFlr::ExecFlr(const ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(true), ownOutput(true),
        m_outputLpLines(0x0), m_outputLpHeader(0x0), m_lpRequiredHeader(0x0),
        m_lpRequiredHeaderIndex(0x0), m_model(0x0), m_modifications(0x0),
        m_massShifts(0x0), m_modAllowedSites(0x0), m_variantGroups(0x0),
        m_groupQuantities(0x0), m_indices(0x0), m_theoCosineSimilarity(0x0),
        m_unmodCosineSimilarity(0x0), m_validVariants(0x0)
  {
    m_name = "ExecFlr";
    m_type = "ExecFlr";
  }

  // -------------------------------------------------------------------------
  ExecFlr::~ExecFlr(void)
  {
    if (ownInput)
    {
      if (m_outputLpLines)
      {
        delete m_outputLpLines;
        m_outputLpLines = 0x0;
      }
      if (m_outputLpHeader)
      {
        delete m_outputLpHeader;
        m_outputLpHeader = 0x0;
      }
      if (m_lpRequiredHeader)
      {
        delete m_lpRequiredHeader;
        m_lpRequiredHeader = 0x0;
      }

      if (m_lpRequiredHeaderIndex)
      {
        delete m_lpRequiredHeaderIndex;
        m_lpRequiredHeaderIndex = 0x0;
      }
      if (m_model)
      {
        delete m_model;
        m_model = 0x0;
      }
      if (m_modifications)
      {
        delete m_modifications;
        m_modifications = 0x0;
      }
      if (m_massShifts)
      {
        delete m_massShifts;
        m_massShifts = 0x0;
      }
      if (m_modAllowedSites)
      {
        delete m_modAllowedSites;
        m_modAllowedSites = 0x0;
      }
    }
    if (ownOutput)
    {
      if (m_indices)
      {
        delete m_indices;
        m_indices = 0x0;
      }
      if (m_groupQuantities)
      {
        delete m_groupQuantities;
        m_groupQuantities = 0x0;
      }
      if (m_variantGroups)
      {
        delete m_variantGroups;
        m_variantGroups = 0x0;
      }
      if (m_theoCosineSimilarity)
      {
        delete m_theoCosineSimilarity;
        m_theoCosineSimilarity = 0x0;
      }
      if (m_unmodCosineSimilarity)
      {
        delete m_unmodCosineSimilarity;
        m_unmodCosineSimilarity = 0x0;
      }
      if (m_validVariants)
      {
        delete m_validVariants;
        m_validVariants = 0x0;
      }
    }
  }
  // -------------------------------------------------------------------------
  ExecBase * ExecFlr::clone(const ParameterList & inputParams) const
  {
    return new ExecFlr(inputParams);
  }
  // -------------------------------------------------------------------------
  bool loadLpSpectra(string &directory,
                     vector<string> &line,
                     vector<int> &headerIndices,
                     Spectrum &unmodSpectrum,
                     Spectrum &modSpectrum,
                     Spectrum &allMassesSpectrum)
  {

    stringstream ss;

    //clear
    ss.str(std::string());
    ss << directory << "/" << line[headerIndices[0]];
    string unmodifiedMgf = ss.str();

    //clear
    ss.str(std::string());
    ss << directory << "/" << line[headerIndices[1]];
    string modifiedMgf = ss.str();

    SpecSet unmodSpecSet;
    SpecSet modSpecSet;

    if (!modSpecSet.LoadSpecSet_mgf(modifiedMgf.c_str()))
    {
      ERROR_MSG("Unable to open " << modifiedMgf);
      return false;
    }

    if (!unmodSpecSet.LoadSpecSet_mgf(unmodifiedMgf.c_str()))
    {
      ERROR_MSG("Unable to open " << unmodifiedMgf);
      return false;
    }

    SpecSet allMassesSpecSet;

    ss.str(std::string());
    ss << directory << "/" << line[headerIndices[4]];
    string allMassFilename = ss.str();

    if (!allMassesSpecSet.LoadSpecSet_mgf(allMassFilename.c_str()))
    {
      WARN_MSG("Unable to read " << allMassFilename);
      return false;
    }
    unmodSpectrum = unmodSpecSet[0];
    modSpectrum = modSpecSet[0];
    allMassesSpectrum = allMassesSpecSet[0];

    DEBUG_VAR(modifiedMgf);
    DEBUG_VAR(unmodifiedMgf);
    DEBUG_VAR(allMassFilename);

    return true;
  }
  // -------------------------------------------------------------------------
  bool loadVariantInformation(string &directory,
                              vector<string> &line,
                              vector<int> &headerIndices,
                              vector<float> &variantQuantities,
                              vector<string> &variants)
  {
    stringstream ss;
    ss.str(std::string());
    ss << directory << "/" << line[headerIndices[3]];
    string lpOutput = ss.str();
    DEBUG_VAR(lpOutput);

    vector<float> peakActivity;
    vector<float> peakError;
    vector<float> tempVariantQuantities;
    vector < string > variantNames;
    double objectiveError;

    if (!LinearEquation::readGlpsolOutput(lpOutput.c_str(),
                                          peakActivity,
                                          peakError,
                                          tempVariantQuantities,
                                          variantNames,
                                          objectiveError))
    {
      WARN_MSG("Unable to read " << lpOutput);
      return false;
    }

    if (tempVariantQuantities.size() == 0)
    {
      WARN_MSG("No variant quantities read!");
      return false;
    }

    ss.str(std::string());
    ss << directory << "/" << line[headerIndices[5]];
    string variantOutput = ss.str();
    DEBUG_VAR(variantOutput);

    vector < vector<string> > variantLines;
    map<string, unsigned int> variantHeader;
    vector < string > variantRequiredHeader;
    variantRequiredHeader.push_back("variantName");
    variantRequiredHeader.push_back("lpVariantName");

    vector<int> variantRequiredHeaderIndex;
    if (!DelimitedTextReader::loadDelimitedFile(variantOutput.c_str(),
                                                "\t",
                                                "#",
                                                variantHeader,
                                                variantLines,
                                                variantRequiredHeader,
                                                variantRequiredHeaderIndex))
    {
      WARN_MSG("Unable to read" << variantOutput);
      return false;
    }
    vector < string > tempVariantNames;
    for (int varIndex = 0; varIndex < variantLines.size(); varIndex++)
    {
      variants.push_back(variantLines[varIndex][variantRequiredHeaderIndex[0]]);
      tempVariantNames.push_back(variantLines[varIndex][variantRequiredHeaderIndex[1]]);
    }

    variantQuantities.resize(tempVariantNames.size());

    DEBUG_VAR(variantNames.size());
    DEBUG_VAR(variantQuantities.size());
    DEBUG_VAR(tempVariantNames.size());
    DEBUG_VAR(tempVariantQuantities.size());

    Timer tm;

    tm.start();
    for (int varIndex = 0; varIndex < tempVariantNames.size(); varIndex++)
    {
      for (int nameIndex = 0; nameIndex < variantNames.size(); nameIndex++)
      {
        if (tempVariantNames[varIndex].compare(variantNames[nameIndex].substr(0,
                                                                              5))
            == 0)
        {
          variantQuantities[varIndex] = tempVariantQuantities[nameIndex];
        }
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecFlr::invoke(void)
  {

    DEBUG_VAR(m_groupQuantities);
    if (m_groupQuantities == 0x0)
    {
      ownOutput = true;
      m_groupQuantities = new vector<vector<float> > ;
      m_variantGroups = new vector<vector<vector<string> > > ;
      m_validVariants = new vector<vector<vector<string> > > ;
      m_indices = new vector<int> ;
      m_theoCosineSimilarity = new vector<float> ;
      m_unmodCosineSimilarity = new vector<float> ;
    }
    DEBUG_VAR(m_groupQuantities);

    if (m_outputLpLines == 0x0 || m_outputLpHeader == 0x0 || m_lpRequiredHeader
        == 0x0 || m_lpRequiredHeaderIndex == 0x0 || m_model == 0x0
        || m_modifications == 0x0 || m_massShifts == 0x0 || m_modAllowedSites
        == 0x0)
    {
      ERROR_MSG("Input pointers not defined!");
      return false;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_TRACE;

    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else
    {
      endBaseIdx = -1;
    }

    DEBUG_VAR(m_params.getValueInt("IDX_START"));
    DEBUG_VAR(m_params.getValueInt("IDX_END"));

    AAJumps jumps(1);

    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

    DEBUG_VAR(peakTol);

    FalseLocalizationRates flr(jumps, *m_model, peakTol);

    bool writeHeader = true;

    string directory = "./";
    if (m_params.exists("OUTPUT_LP_DIRECTORY"))
    {
      directory = m_params.getValue("OUTPUT_LP_DIRECTORY");
    }

    if (endBaseIdx == -1 || endBaseIdx > m_outputLpLines->size() - 1)
    {
      endBaseIdx = m_outputLpLines->size() - 1;
    }

    for (int i = startBaseIdx; i <= endBaseIdx; i++)
    {
      string unmodifiedPeptide =
          (*m_outputLpLines)[i][(*m_lpRequiredHeaderIndex)[2]];
      stringstream ss;
      //clear
      ss.str(std::string());
      ss << directory << "/"
          << (*m_outputLpLines)[i][(*m_lpRequiredHeaderIndex)[0]];
      string unmodifiedMgf = ss.str();

      //load files from LP
      Spectrum unmodSpectrum;
      Spectrum modSpectrum;
      Spectrum allMassesSpectrum;

      if (!loadLpSpectra(directory,
                         (*m_outputLpLines)[i],
                         *m_lpRequiredHeaderIndex,
                         unmodSpectrum,
                         modSpectrum,
                         allMassesSpectrum))
      {
        ERROR_MSG("Unable to load associated LP files!");
        continue;
      }

      vector<float> tempQuantities;
      vector < string > variants;

      if (!loadVariantInformation(directory,
                                  (*m_outputLpLines)[i],
                                  *m_lpRequiredHeaderIndex,
                                  tempQuantities,
                                  variants))
      {
        ERROR_MSG("Unable to load variant information!");
        continue;
      }

      Spectrum outputSpectrum;

      DEBUG_TRACE;

      flr.generateTheoreticalSpectrum(unmodifiedPeptide,
                                      unmodSpectrum,
                                      variants,
                                      tempQuantities,
                                      allMassesSpectrum,
                                      outputSpectrum);

      size_t extFind;
      extFind = unmodifiedMgf.find_last_of(".");
      string filePrefix = unmodifiedMgf.substr(0, extFind);
      filePrefix.append("_theo.mgf");

      SpecSet temp;
      temp.resize(1);
      temp[0] = outputSpectrum;
      if (!temp.SaveSpecSet_mgf(filePrefix.c_str()))
      {
        ERROR_MSG("UNABLE TO OPEN FILE!" << filePrefix);
      }

      DEBUG_TRACE;

      //get mods.
      PeptideSpectrumMatch modified;
      modified.m_annotation = variants[0];

      vector<float> massShifts;
      modified.getModifications(massShifts);
      vector < string > modAllowedSites(massShifts.size());
      vector<float> modifications(massShifts.size());

      DEBUG_VAR(massShifts.size());
      DEBUG_VAR(m_massShifts->size());

      DEBUG_VAR(massShifts[0]);
      DEBUG_VAR((*m_massShifts)[0]);

      for (int modi = 0; modi < m_massShifts->size(); modi++)
      {
        for (int modj = 0; modj < massShifts.size(); modj++)
        {
          if ((*m_massShifts)[modi] < massShifts[modj] + .001
              && (*m_massShifts)[modi] > massShifts[modj] - .001)
          {
            modifications[modj] = (*m_modifications)[modi];
          }
        }
      }

      DEBUG_TRACE;
      DEBUG_VAR(modifications[0]);
      for (int modIdx = 0; modIdx < modifications.size(); modIdx++)
      {
        for (int pModIdx = 0; pModIdx < m_modifications->size(); pModIdx++)
        {
          if (modifications[modIdx] == (*m_modifications)[pModIdx])
          {
            modAllowedSites[modIdx] = (*m_modAllowedSites)[pModIdx];
          }
        }
      }

      vector<float> groupQuantities;
      vector<float> minDistinguishingIntensity;
      vector < vector<string> > variantGroups;

      DEBUG_TRACE;

      flr.groupVariants3(modSpectrum,
                        allMassesSpectrum,
                         variants,
                         tempQuantities,
                         massShifts,
                         variantGroups,
                         groupQuantities,
                         minDistinguishingIntensity,
                         m_params.getValueFloat("MIN_GROUPING_INTENSITY"));

      DEBUG_TRACE;

      vector < vector<string> > validVariants(variantGroups.size());

      for (int varIndex = 0; varIndex < variantGroups.size(); varIndex++)
      {
        vector < string > currVariants;

        flr.groupHasValidVariant(modAllowedSites,
                                 massShifts,
                                 variantGroups[varIndex],
                                 currVariants);
        validVariants[varIndex] = currVariants;
//        DEBUG_VAR(currVariants.size());
      }

      DEBUG_TRACE;

      DEBUG_VAR(groupQuantities.size());
      m_groupQuantities->push_back(groupQuantities);
      m_variantGroups->push_back(variantGroups);
      m_indices->push_back(i);
      m_validVariants->push_back(validVariants);

      //calculate theoretical similarity amd unmod similarity
      float theoSimilarity = 0;
      float unmodModSimilarity = 0;

//      DEBUG_VAR(variantGroups.size());

      for (int varIndex = 0; varIndex < variantGroups.size(); varIndex++)
      {
        if (groupQuantities[varIndex] >= 0)
        {
          string variant = variantGroups[varIndex][0];

          //calculate similarity
          vector < string > ionsToExtract;
          ionsToExtract.resize(11);
          ionsToExtract[0] = "b";
          ionsToExtract[1] = "b-NH3";
          ionsToExtract[2] = "b-H2O";
          ionsToExtract[3] = "b2";
          ionsToExtract[4] = "b2-iso";
          ionsToExtract[5] = "y";
          ionsToExtract[6] = "y-NH3";
          ionsToExtract[7] = "y-H2O";
          ionsToExtract[8] = "y2";
          ionsToExtract[9] = "y2-iso";
          ionsToExtract[10] = "P2-H2O";

          CosineSimilarity c_similarity(ionsToExtract);
          float currTheoSimilarity = c_similarity.similarity(outputSpectrum,
                                                             modSpectrum,
                                                             variant,
                                                             variant,
                                                             *m_model,
                                                             false);

          float currMsSimilarity = c_similarity.similarity(unmodSpectrum,
                                                           modSpectrum,
                                                           unmodifiedPeptide,
                                                           variant,
                                                           *m_model,
                                                           false);
          if (theoSimilarity < currTheoSimilarity)
          {
            theoSimilarity = currTheoSimilarity;
          }

          if (unmodModSimilarity < currMsSimilarity)
          {
            unmodModSimilarity = currMsSimilarity;
          }
        }
      }
      DEBUG_TRACE;
      m_theoCosineSimilarity->push_back(theoSimilarity);
      m_unmodCosineSimilarity->push_back(unmodModSimilarity);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecFlr::saveOutputData(void)
  {
    if (m_groupQuantities == 0x0)
    {
      return false;
    }
    ofstream tempFile;

    if (m_params.exists("OUTPUT_FLR_FILE"))
    {
      float similarityParam = 0;

      if (m_params.exists("MIN_COSINE_SIMILARITY"))
      {
        similarityParam = m_params.getValueFloat("MIN_COSINE_SIMILARITY");
      }

      tempFile.open(m_params.getValue("OUTPUT_FLR_FILE").c_str(), ios::out | ios::trunc | ios::binary);
      if (tempFile.fail())
      {
        ERROR_MSG("Unable to open file! " << m_params.getValue("OUTPUT_FLR_FILE"));
        return false;

      }
      if (tempFile.is_open() && tempFile.good())
      {
        tempFile
            << "groupIndex\tquantity\ttheoSimilarity\tunmodModSimilarity\tvariantGroup\tunmodifiedPeptide\tisValid\tvalidVariants"
            << endl;

        DEBUG_VAR(m_theoCosineSimilarity->size());
        DEBUG_VAR(m_unmodCosineSimilarity->size());
        DEBUG_VAR(m_variantGroups->size());
        DEBUG_VAR(m_groupQuantities->size());
        DEBUG_VAR(m_indices->size());
        DEBUG_VAR(m_validVariants->size());

        for (int i = 0; i < m_indices->size(); i++)
        {
          vector < vector<string> > *variantGroups = &((*m_variantGroups)[i]);
          vector<float> * groupQuantities = &((*m_groupQuantities)[i]);

          DEBUG_VAR((*m_indices)[i]);

          string
              unmodifiedPeptide =
                  (*m_outputLpLines)[(*m_indices)[i]][(*m_lpRequiredHeaderIndex)[2]];

          DEBUG_VAR(variantGroups->size());

          for (int varIndex = 0; varIndex < variantGroups->size(); varIndex++)
          {
            if ((*m_theoCosineSimilarity)[i] >= similarityParam
                && (*groupQuantities)[varIndex] >= 0)
            {
              tempFile << (*m_indices)[i];
              tempFile << "\t" << (*groupQuantities)[varIndex];
              tempFile << "\t" << (*m_theoCosineSimilarity)[i];
              tempFile << "\t" << (*m_unmodCosineSimilarity)[i];
              tempFile << "\t" << (*variantGroups)[varIndex][0];
              for (int group = 1; group < (*variantGroups)[varIndex].size(); group++)
              {
                tempFile << ":" << (*variantGroups)[varIndex][group];
              }
              tempFile << "\t" << unmodifiedPeptide;

              if ((*m_validVariants)[i][varIndex].size() > 0)
              {
                tempFile << "\ttrue\t";
                tempFile << (*m_validVariants)[i][varIndex][0];
                for (int k = 1; k < (*m_validVariants)[i][varIndex].size(); k++)
                {
                  tempFile << ":" << (*m_validVariants)[i][varIndex][k];
                }
              }
              else
              {
                tempFile << "\tfalse\t";
              }
              tempFile << endl;
            }
          }
        }
      }
      else
      {
        ERROR_MSG("Unable to write to file!" << m_params.getValue("OUTPUT_FLR_FILE"));
        return false;
      }
    }
    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFlr::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();
    //SpecSet m_inputSpectra; // the input spectra
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecFlr::loadInputData(void)
  {
    if (m_outputLpLines == NULL)
    {
      m_outputLpLines = new vector<vector<string> > ;
      m_outputLpHeader = new map<string, unsigned int> ;
      m_lpRequiredHeader = new vector<string> ;
      m_lpRequiredHeaderIndex = new vector<int> ;
      m_model = new MS2ScoringModel;
      m_modifications = new vector<float> ;
      m_massShifts = new vector<float> ;
      m_modAllowedSites = new vector<string> ;
      ownInput = true;
    }

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
        m_modifications->push_back(modification);
      }
      DEBUG_VAR((*m_modifications)[0]);
    }

    if (m_params.exists("MASS_SHIFTS"))
    {
      string massShiftString = m_params.getValue("MASS_SHIFTS");
      vector < string > massShifts;
      stringSplit(massShiftString, massShifts, ",");
      for (int i = 0; i < massShifts.size(); i++)
      {
        float massShift;
        stringstream ss;
        ss << massShifts[i];
        ss >> massShift;
        m_massShifts->push_back(massShift);
      }
      DEBUG_VAR((*m_massShifts)[0]);
    }

    if (m_params.exists("MODIFICATIONS_ALLOWED_SITES"))
    {
      string massShiftString = m_params.getValue("MODIFICATIONS_ALLOWED_SITES");
      vector < string > allowedSites;
      stringSplit(massShiftString, allowedSites, ",");
      for (int i = 0; i < allowedSites.size(); i++)
      {
        m_modAllowedSites->push_back(allowedSites[i]);
      }
    }

    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_LP_INFORMATION"))
    {

      m_lpRequiredHeader->push_back("unmodifiedMgf");
      m_lpRequiredHeader->push_back("modifiedMgf");
      m_lpRequiredHeader->push_back("unmodifiedAnnotation");
      m_lpRequiredHeader->push_back("lpOutputFilename");
      m_lpRequiredHeader->push_back("allMassMgf");
      m_lpRequiredHeader->push_back("variantOutputFile");

      if (!DelimitedTextReader::loadDelimitedFile(m_params.getValue("OUTPUT_LP_INFORMATION").c_str(),
                                                  "\t",
                                                  "#",
                                                  *m_outputLpHeader,
                                                  *m_outputLpLines,
                                                  *m_lpRequiredHeader,
                                                  *m_lpRequiredHeaderIndex))
      {
        ERROR_MSG("Unable to load lp information files!" << m_params.getValue("OUTPUT_LP_INFORMATION"));
        return false;
      }
    }

    DEBUG_TRACE;
    if (m_params.exists("MS2_SCORING_MODEL"))
    {
      if (!m_model->LoadModel(m_params.getValue("MS2_SCORING_MODEL").c_str()))
      {
        ERROR_MSG("Unable to load model! " << m_params.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFlr::loadOutputData(void)
  {
    if (m_groupQuantities == 0x0)
    {
      ownOutput = true;
      m_groupQuantities = new vector<vector<float> > ;
      m_variantGroups = new vector<vector<vector<string> > > ;
      m_validVariants = new vector<vector<vector<string> > > ;
      m_indices = new vector<int> ;
      m_theoCosineSimilarity = new vector<float> ;
      m_unmodCosineSimilarity = new vector<float> ;
    }
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;
    vector < string > header;
    vector < vector<string> > lines;

    requiredHeader.push_back("groupIndex"); //0
    requiredHeader.push_back("quantity"); //1
    requiredHeader.push_back("theoSimilarity"); //2
    requiredHeader.push_back("unmodModSimilarity"); //3
    requiredHeader.push_back("variantGroup"); //4
    requiredHeader.push_back("unmodifiedPeptide"); //5
    requiredHeader.push_back("validVariants"); //6

    if (!DelimitedTextReader::loadDelimitedFile(m_params.getValue("OUTPUT_FLR_FILE").c_str(),
                                                "\t",
                                                "#",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load lp information files!" << m_params.getValue("OUTPUT_FLR_FILE"));
      return false;
    }

    int lastGroup = -1;
    int totalGroups = 0;
    for (int i = 0; i < lines.size(); i++)
    {
      int currGroup;
      stringstream ss;
      ss << lines[i][requiredHeaderIndex[0]];
      ss >> currGroup;
      if (currGroup != lastGroup)
      {
        totalGroups++;
      }
    }
    m_variantGroups->resize(totalGroups);
    m_groupQuantities->resize(totalGroups);
    m_validVariants->resize(totalGroups);
    m_indices->resize(totalGroups);
    m_theoCosineSimilarity->resize(totalGroups);
    m_unmodCosineSimilarity->resize(totalGroups);

    lastGroup = -1;
    int groupCount = 0;

    for (int i = 0; i < lines.size(); i++)
    {
      int currGroup;
      stringstream ss;
      ss << lines[i][requiredHeaderIndex[0]];
      ss >> currGroup;
      ss.str("");

      float quantity;
      ss << lines[i][requiredHeaderIndex[1]];
      ss >> quantity;
      ss.str("");

      float theoSimilarity;
      ss << lines[i][requiredHeaderIndex[2]];
      ss >> theoSimilarity;
      ss.str("");

      float unmodModSimilarity;
      ss << lines[i][requiredHeaderIndex[3]];
      ss >> unmodModSimilarity;
      ss.str("");

      (*m_indices)[groupCount] = currGroup;
      (*m_groupQuantities)[groupCount].push_back(quantity);
      (*m_theoCosineSimilarity)[groupCount] = theoSimilarity;
      (*m_unmodCosineSimilarity)[groupCount] = unmodModSimilarity;

      vector<string> variants;
      stringSplit2(lines[i][requiredHeaderIndex[4]],variants,":");
      (*m_variantGroups)[groupCount].push_back(variants);

      vector<string> validVariants;
      stringSplit2(lines[i][requiredHeaderIndex[6]],validVariants,":");
      (*m_validVariants)[groupCount].push_back(validVariants);

      if (currGroup != lastGroup)
      {
        groupCount++;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFlr::split(int numSplit)
  {
    DEBUG_VAR(numSplit);

    if (numSplit < 2)
    {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }
    int numLpLines = m_outputLpLines->size();
    DEBUG_VAR(numLpLines);
    if (numLpLines == 0)
    {
      DEBUG_MSG("Must have at least one spectrum");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);

    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx, (int)m_outputLpLines->size() - 1);
    }
    else
    {
      endBaseIdx = m_outputLpLines->size() - 1;
    }
    DEBUG_VAR(endBaseIdx);

    if (startBaseIdx == endBaseIdx)
    {
      DEBUG_MSG("Must have more than one row");
      return m_subModules;
    }

    DEBUG_TRACE;

    //find remainder;
    int numPerSplit = (int)ceil(((float)endBaseIdx - (float)startBaseIdx + 1)/ (float)numSplit);

    DEBUG_VAR(numPerSplit);

    int startIndex = startBaseIdx;
    int endIndex = startBaseIdx + numPerSplit - 1;

    bool splitAll = false;

    for (int i = 0; i < numSplit; i++)
    {
      if (splitAll)
      {
        continue;
      }
      // We have enough comparisons for the CPU
      // Copy the parameters
      ParameterList childParams(m_params);
      // But set the start and end indices
      char buf[128];
      sprintf(buf, "%d", startIndex);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", endIndex);
      childParams.setValue("IDX_END", buf);

      //DEBUG_MSG("Start [" << startIndex << "] End [" << i << "] Split [" << numSplit << "]");
      ExecBase * theClone = new ExecFlr(m_outputLpLines,
                                        m_outputLpHeader,
                                        m_lpRequiredHeader,
                                        m_lpRequiredHeaderIndex,
                                        m_model,
                                        m_modifications,
                                        m_massShifts,
                                        m_modAllowedSites,
                                        childParams);

      theClone->setName(makeName(m_name, i));

      string directory = "";
      if (m_params.exists("OUTPUT_FLR_DIRECTORY"))
      {
        directory = m_params.getValue("OUTPUT_FLR_DIRECTORY");
      }

      if (directory.length() > 2 && directory[directory.length() - 1] == '/')
      {
        directory = directory.substr(0, directory.length() - 1);
      }

      string baseName = directory + "/" + theClone->getName();
      DEBUG_VAR(baseName);

      theClone->m_params.setValue("OUTPUT_FLR_FILE", baseName + "_flr" + ".txt");

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", i + 1);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      m_subModules.push_back(theClone);
      startIndex = endIndex + 1;
      endIndex = startIndex + numPerSplit - 1;

      if (endIndex > (int)m_outputLpLines->size() - 1)
      {
        endIndex = (int)m_outputLpLines->size() - 1;
      }
      if (endIndex < startIndex)
      {
        splitAll = true;
      }
    }
    DEBUG_VAR(m_subModules.size());
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecFlr::merge(void)
  {
    if (m_groupQuantities == 0x0)
    {
      ownOutput = true;
      m_groupQuantities = new vector<vector<float> > ;
      m_variantGroups = new vector<vector<vector<string> > > ;
      m_validVariants = new vector<vector<vector<string> > > ;
      m_indices = new vector<int> ;
      m_theoCosineSimilarity = new vector<float> ;
      m_unmodCosineSimilarity = new vector<float> ;
    }

    DEBUG_TRACE;

    unsigned int startIdx = (unsigned int)m_params.getValueInt("IDX_START", 0);
    unsigned int endIdx = (unsigned int)m_params.getValueInt("IDX_END", 0);

    // Find out how many total pairs from all the children
    DEBUG_VAR(m_subModules.size());

    // Find out how many total pairs from all the children
    DEBUG_VAR(m_subModules.size());
    int totalResults = 0;

    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecFlr * fpe = (ExecFlr *)m_subModules[i];
      DEBUG_VAR(fpe->m_groupQuantities->size());
      totalResults += fpe->m_groupQuantities->size();
    }

    DEBUG_VAR(totalResults);

    if (totalResults == 0)
    {
      ERROR_MSG("No results! Unable to merge");
      return false;
    }

    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecFlr * fpe = (ExecFlr *)m_subModules[i];
      DEBUG_VAR(fpe);

      for (int resultIdx = 0; resultIdx < (*fpe->m_groupQuantities).size(); resultIdx++)
      {
        vector<float> currQuantities;
        for (int j = 0; j < (*fpe->m_groupQuantities)[resultIdx].size(); j++)
        {
          currQuantities.push_back((*fpe->m_groupQuantities)[resultIdx][j]);
        }

        vector < vector<string> > currVariants;
        for (int j = 0; j < (*fpe->m_variantGroups)[resultIdx].size(); j++)
        {
          vector < string > currVariant;
          for (int varIndex = 0; varIndex
              < (*fpe->m_variantGroups)[resultIdx][j].size(); varIndex++)
          {
            currVariant.push_back((*fpe->m_variantGroups)[resultIdx][j][varIndex]);
          }
          currVariants.push_back(currVariant);
        }

        vector < vector<string> > currValidVariants;
        for (int j = 0; j < (*fpe->m_validVariants)[resultIdx].size(); j++)
        {
          vector < string > currValidVariant;
          for (int varIndex = 0; varIndex
              < (*fpe->m_validVariants)[resultIdx][j].size(); varIndex++)
          {
            currValidVariant.push_back((*fpe->m_validVariants)[resultIdx][j][varIndex]);
          }
          currValidVariants.push_back(currValidVariant);
        }
        //set values
        m_groupQuantities->push_back(currQuantities);
        m_variantGroups->push_back(currVariants);
        m_validVariants->push_back(currValidVariants);
        m_indices->push_back((*fpe->m_indices)[resultIdx]);
        m_theoCosineSimilarity->push_back((*fpe->m_theoCosineSimilarity)[resultIdx]);
        m_unmodCosineSimilarity->push_back((*fpe->m_unmodCosineSimilarity)[resultIdx]);
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecFlr::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("MIN_GROUPING_INTENSITY");
    VALIDATE_PARAM_EXIST("MODIFICATIONS");
    VALIDATE_PARAM_EXIST("MASS_SHIFTS");
    VALIDATE_PARAM_EXIST("MODIFICATIONS_ALLOWED_SITES");

    m_isValid = true;
    return true;
  }

}
