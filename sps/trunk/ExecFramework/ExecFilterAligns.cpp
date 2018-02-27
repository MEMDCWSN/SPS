// Header Include
#include "ExecFilterAligns.h"

// Module Includes
#include "AlignmentUtils.h"
#include "Logger.h"
#include "Filters.h"
#include "FileUtils.h"

// External Includes
#include "utils.h"  // for Save_binArray only
// System Includes
#include "stdlib.h"

using namespace specnets;
using namespace std;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecFilterAligns::ExecFilterAligns(void) :
    m_filteredPairs(0x0), m_spectra(0x0), 
        m_ratios(0x0), m_means(0x0), m_varTerms(0x0), ownInput(true), 
        m_idxKept(0x0), m_pvalues(0x0), ownOutput(true)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::ExecFilterAligns(const ParameterList & inputParams) :
    ExecBase(inputParams), m_filteredPairs(0x0), m_spectra(0x0), 
        m_ratios(0x0), m_means(0x0), m_varTerms(0x0), ownInput(true), 
        m_idxKept(0x0), m_pvalues(0x0), ownOutput(true)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::ExecFilterAligns(const ParameterList & inputParams,
                                     SpectrumPairSet * filteredPairs,
                                     SpecSet * spectra,
                                     PeptideSpectrumMatchSet * filterPsmSet,
                                     vector<TwoValues<float> > * ratios,
                                     vector<TwoValues<float> > * means,
                                     vector<float> * varTerms,
                                     vector<unsigned int> * idxKept,
                                     vector<TwoValues<float> > * pvalues) :
    ExecBase(inputParams), m_filteredPairs(filteredPairs), m_spectra(spectra), 
        m_filterPsmSet(filterPsmSet),
        m_ratios(ratios), m_means(means), m_varTerms(varTerms), ownInput(false), 
        m_idxKept(idxKept), m_pvalues(pvalues), ownOutput(false)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::~ExecFilterAligns(void)
  {
  #if 1
    if (ownInput) {
      if (m_filteredPairs) {
        //delete m_filteredPairs;
      }
      if (m_spectra) {
        //delete m_spectra;
      }
      if (m_ratios) {
        delete m_ratios;
      }
      if (m_means) {
        delete m_means;
      }
      if (m_varTerms) {
        delete m_varTerms;
      }
    } 
    #endif
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFilterAligns::clone(const ParameterList & inputParams) const
  {
    return new ExecFilterAligns(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::invoke(void)
  {
    DEBUG_VAR(m_filteredPairs->size());
    DEBUG_TRACE;
    if (m_params.getValueBool("DEBUG_PARAMS")) {
      stringstream aux;
      m_params.print(aux);
      DEBUG_MSG(aux.str());
    }

    if (m_idxKept == 0x0)
    {
      ownOutput = true;
      m_idxKept = new vector<unsigned int> ();
      m_pvalues = new vector<TwoValues<float> > ();
    }

    float minPValue = m_params.getValueDouble("MAX_PVALUE");
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);
    bool filterTrigs = m_params.getValueBool("FILTER_TRIGS", true);
    bool runAlignGf = m_params.getValueBool("ALIGNGF", false);

    DEBUG_VAR(minPValue);
    DEBUG_VAR(minRatio);
    DEBUG_VAR(pmTol);
    DEBUG_VAR(peakTol);
    DEBUG_VAR(filterTrigs);
    DEBUG_VAR(runAlignGf);

    // Convert varTerms to standard deviations
    for (int pivot = 0; pivot < m_means->size(); pivot++)
    {
      (*m_varTerms)[pivot] = sqrt((*m_varTerms)[pivot] - (*m_means)[pivot][0]
          * (*m_means)[pivot][0]);
    }

    DEBUG_VAR(m_means->size());
    DEBUG_VAR(m_varTerms->size());
    DEBUG_VAR(m_ratios->size());

    if (!runAlignGf) {
      unsigned int numKept = filterAligns(*m_filteredPairs,
                                        *m_idxKept,
                                        *m_pvalues,
                                        *m_means,
                                        *m_varTerms,
                                        *m_ratios,
                                        minPValue,
                                        minRatio,
                                        pmTol,
                                        filterTrigs);

      DEBUG_VAR(numKept);
    } else if (filterTrigs) {
      filterTriangles(*m_filteredPairs, 0, m_filteredPairs->size(), pmTol, *m_idxKept);
      unsigned int idxLast = 0;
      for (unsigned int idxPair = 0; idxPair < m_idxKept->size(); idxPair++) {
        if (idxLast < idxPair) {
          (*m_filteredPairs)[idxLast] = (*m_filteredPairs)[idxPair];
          (*m_idxKept)[idxLast] = (*m_idxKept)[idxPair];
        }
        idxLast++;
      }
      m_filteredPairs->resize(idxLast);
      m_idxKept->resize(idxLast);
    }    
    DEBUG_VAR(m_filteredPairs->size());
    DEBUG_VAR(m_pvalues->size());
    DEBUG_VAR(m_idxKept->size());

    //---------------------------------------------------------------------------
    // FILTER BY #SPEC/CONNECTED UNIQUE PM IF DESIRED
    //---------------------------------------------------------------------------
    int maxVarSpec 			= m_params.getValueInt("MAX_CONNECTED_VAR_SPEC", -1);
    DEBUG_VAR(maxVarSpec);
    if (runAlignGf && maxVarSpec > 0) {
      if (!m_spectra || m_spectra->size() == 0) {
        ERROR_MSG("MAX_CONNECTED_VAR_SPEC specified but no spectra present.");
        return false;
      }
      m_filteredPairs->filter_by_max_spec_per_variant(
              *m_spectra,
              maxVarSpec);
      DEBUG_VAR(m_filteredPairs->size());
    }

    DEBUG_VAR(m_filteredPairs->size());

    
    //---------------------------------------------------------------------------
    // FILTER BY UNIQUE MASS NUMBER IF DESIRED
    //---------------------------------------------------------------------------
    int maxUniqueMassNumber = m_params.getValueInt("MAX_UNIQUE_MASS_NUMBER", -1);
    DEBUG_VAR(maxUniqueMassNumber);
    if (runAlignGf && maxUniqueMassNumber > 0) {
      if (!m_spectra || m_spectra->size() == 0) {
        ERROR_MSG("MAX_UNIQUE_MASS_NUMBER specified but no spectra present.");
        return false;
      }
      m_filteredPairs->filter_by_unique_mass_number_in_component(
              *m_spectra, 
              maxUniqueMassNumber);
      DEBUG_VAR(m_filteredPairs->size());
    }

    DEBUG_VAR(m_filteredPairs->size());
    
    //---------------------------------------------------------------------------
    // FILTER USING PSMs TO SET FDR
    //---------------------------------------------------------------------------
    float edgeFdr = m_params.getValueFloat("FILTERALIGNS_EDGE_FDR");
    DEBUG_VAR(edgeFdr);
    if ( m_params.getValueBool("ALIGNGF", false) && edgeFdr > 0 ) {
      DEBUG_VAR(m_filteredPairs->size());
      if (!m_filteredPairs->filter_by_edge_fdr22(
            *m_spectra, 
            *m_filterPsmSet,
            edgeFdr, 
            pmTol,
            peakTol,
            false)) {
        ERROR_MSG("Unable to filter by edge FDR");
        return false;
      }
      DEBUG_VAR(m_filteredPairs->size());
    }
    
    //---------------------------------------------------------------------------
    // FILTER BY COMPONENT SIZE IF DESIRED
    //---------------------------------------------------------------------------
    int maxComponentSize = m_params.getValueInt("MAX_COMPONENT_SIZE", -1);
    DEBUG_VAR(maxComponentSize);
    if (maxComponentSize > 0) {
      if (!m_filteredPairs->filter_by_component_size(maxComponentSize)) {
        ERROR_MSG("Unable to filter by component size");
        return false;
      }
      DEBUG_VAR(m_filteredPairs->size());
    }


    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::loadInputData(void)
  {
    if (m_filteredPairs == 0) {
      ownInput = true;
      m_filteredPairs = new SpectrumPairSet();
      m_spectra = new SpecSet();
      m_filterPsmSet = new PeptideSpectrumMatchSet();
      m_ratios = new vector<TwoValues<float> > ();
      m_means = new vector<TwoValues<float> > ();
      m_varTerms = new vector<float> ();
    } 
  
    string dataDir;
    if (m_params.exists("INPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("INPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    
    if (m_params.exists("INPUT_ALIGNS")) {
      string fileName = dataDir + m_params.getValue("INPUT_ALIGNS");
      DEBUG_MSG("Loading aligns..." << fileName);
      if (!m_filteredPairs->loadFromBinaryFile(fileName)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
    }

    if (m_params.exists("INPUT_SPECTRA")) {
      string fileName = dataDir + m_params.getValue("INPUT_SPECTRA");
      DEBUG_MSG("Loading: " << fileName);
      if (m_spectra->loadPklBin(fileName.c_str()) < 0) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
    }

    if (m_params.exists("INPUT_RATIOS")) {
      string fileName = dataDir + m_params.getValue("INPUT_RATIOS");
      DEBUG_MSG("Loading ratios..." << fileName);
      if (!Load_binArray(fileName.c_str(), *m_ratios)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
      DEBUG_VAR(m_ratios->size());
    }

    if (m_params.exists("INPUT_MEANS")) {
      string fileName = dataDir + m_params.getValue("INPUT_MEANS");
      DEBUG_MSG("Loading means..." << fileName);
      if (!Load_binArray(fileName.c_str(), *m_means)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
      DEBUG_VAR(m_means->size());
    }

    if (m_params.exists("INPUT_VARIANCE")) {
      string fileName = dataDir + m_params.getValue("INPUT_VARIANCE");
      DEBUG_MSG("Loading variances..." << fileName);
      if (!Load_binArray(fileName.c_str(), *m_varTerms)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
      DEBUG_VAR(m_varTerms->size());
    }

    string psmFilename = m_params.getValue("INPUT_FILTER_PSMS");
    if (!psmFilename.empty()) {
      DEBUG_MSG("Loading PSM file [" << psmFilename << "]...");
      if (!m_filterPsmSet->loadFromFile(m_params.getValue("INPUT_FILTER_PSMS").c_str())) {
        ERROR_MSG("Error reading input filter PSM file.");
        return false;
      }
    } 
    
    psmFilename = m_params.getValue("INPUT_FILTER_MSGFDB_PSMS");
    if (!psmFilename.empty()) {
      DEBUG_MSG("Loading MSGFDB PSM file [" << psmFilename << "]...");
      if (!m_filterPsmSet->loadMSGFDBResultsFile(m_params.getValue("INPUT_FILTER_MSGFDB_PSMS").c_str())) {
        ERROR_MSG("Error reading input filter MSGFDB PSM file.");
        return false;
      }
    } 
    
    psmFilename = m_params.getValue("INPUT_FILTER_MSGFPLUS_PSMS");
    if (!psmFilename.empty()) {
      DEBUG_MSG("Loading MSGF+ PSM file [" << psmFilename << "]...");
      if (!m_filterPsmSet->loadMSGFPlusResultsFile(m_params.getValue("INPUT_FILTER_MSGFPLUS_PSMS").c_str())) {
        ERROR_MSG("Error reading input filter MSGF+ PSM file.");
        return false;
      }
    }

    
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::saveOutputData(void)
  {
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    
    if (m_params.exists("OUTPUT_ALIGNS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_ALIGNS");
      DEBUG_VAR(m_filteredPairs->size());
      DEBUG_MSG("Outputting aligns..." << fileName);
      if (!m_filteredPairs->saveToBinaryFile(fileName)) {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_ALIGNS"));
        return false;
      }
    }
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_PVALUES")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_PVALUES");
      DEBUG_MSG("Outputting pvalues..." << fileName);
      if (!Save_binArray(fileName.c_str(), *m_pvalues)) {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_PVALUES"));
        return false;
      }
    }
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_INDICES")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_INDICES");
      DEBUG_MSG("Outputting indices..." << fileName);
      if (!Save_binArray(fileName.c_str(), *m_idxKept)) {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_INDICES"));
        return false;
      }
    }

    if (m_params.exists("TXT_OUTPUT_ALIGNS")) {
	  string fileName = dataDir + m_params.getValue("TXT_OUTPUT_ALIGNS");
	  DEBUG_VAR(m_filteredPairs->size());
	  DEBUG_MSG("Outputting aligns..." << fileName);
	  if (!m_filteredPairs->saveToTSVFile(fileName)) {
		ERROR_MSG("Could not save: " << m_params.getValue("TXT_OUTPUT_ALIGNS"));
		return false;
	  }
	}

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::saveInputData(std::vector<std::string> & filenames)
  {
    string alignsFilename = m_params.getValue("INPUT_ALIGNS");
    if (!fileExists(alignsFilename))
    {
      m_filteredPairs->saveToBinaryFile(alignsFilename);
      DEBUG_MSG("Saving " << alignsFilename);
    }
    else
    {
      DEBUG_MSG("Not Saving " << alignsFilename << " (already exists)");
    }
    m_params.setValue("INPUT_ALIGNS", alignsFilename);

    string baseFilename = getName();
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(alignsFilename);

    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::loadOutputData(void)
  {
    if (m_idxKept == 0x0) {
      ownOutput = true;
      m_idxKept = new vector<unsigned int> ();
      m_pvalues = new vector<TwoValues<float> > ();
    }
    
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    
    if (m_params.exists("OUTPUT_ALIGNS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_ALIGNS");
      DEBUG_MSG("Loading aligns..." << fileName);
      if (!m_filteredPairs->loadFromBinaryFile(fileName)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PVALUES")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_PVALUES");
      DEBUG_MSG("Loading pvalues..." << fileName);
      if (!Load_binArray(fileName.c_str(), *m_pvalues)) {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_PVALUES"));
        return false;
      }
    }
    if (m_params.exists("OUTPUT_INDICES")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_INDICES");
      DEBUG_MSG("Loading indices..." << fileName);
      if (!Load_binArray(fileName.c_str(), *m_idxKept)) {
        ERROR_MSG("Could not load: " << fileName);
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFilterAligns::split(int numSplit)
  {
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("MAX_PVALUE");
    VALIDATE_PARAM_EXIST("MIN_RATIO");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("FILTER_TRIGS");
    
    m_isValid = true;
    return true;
  }

} // namespace specnets

