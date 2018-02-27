//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecProteosafeToSpecnetsParams.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unistd.h>

using namespace specnets;
using namespace std;

namespace specnets
{

  ExecProteosafeToSpecnetsParams::ExecProteosafeToSpecnetsParams(void)
  {
    AAJumps aajumps(1);
    m_penaltyMatrix = new PenaltyMatrix(aajumps);
    m_name = "ExecProteosafeToSpecnetsParams";
    m_type = "ExecProteosafeToSpecnetsParams";
  }

  ExecProteosafeToSpecnetsParams::ExecProteosafeToSpecnetsParams(const ParameterList & inputParams)
  {
    AAJumps aajumps(1);
    m_penaltyMatrix = new PenaltyMatrix(aajumps);
    m_name = "ExecProteosafeToSpecnetsParams";
    m_type = "ExecProteosafeToSpecnetsParams";
    m_params = inputParams;
  }

  ExecProteosafeToSpecnetsParams::~ExecProteosafeToSpecnetsParams(void)
  {
    delete m_penaltyMatrix;
  }

  ExecBase * ExecProteosafeToSpecnetsParams::clone(const ParameterList & inputParams) const
  {
    return new ExecProteosafeToSpecnetsParams(inputParams);
  }

  bool ExecProteosafeToSpecnetsParams::invoke(void)
  {
    m_params.setValue("RELATIVE_DIR", "1");
    
    if (!m_params.exists("instrument.instrument"))
    {
      ERROR_MSG("Parameter \'instrument.instrument\' not specified!!");
      return false;
    }
    string IT_ID = "ESI-ION-TRAP";
    string FT_ID = "FT-HYBRID";
    string QTOF_ID = "QTOF";
    string Inst_ID = m_params.getValue("instrument.instrument");
    if (Inst_ID == IT_ID) {
      m_params.setValue("DECONV_MS2", "0");
      m_params.setValue("INSTRUMENT_TYPE", "IT");
    } else if (Inst_ID == FT_ID) {
      m_params.setValue("DECONV_MS2", "1");
      m_params.setValue("INSTRUMENT_TYPE", "FT");
    } else if (Inst_ID == QTOF_ID) {
      m_params.setValue("DECONV_MS2", "0");
      m_params.setValue("INSTRUMENT_TYPE", "qtof");
    } else {
      ERROR_MSG("Unrecognized value \'" << Inst_ID << " for parameter \'instrument.instrument\'!!");
      return false;
    }

    if (!m_params.exists("fragmentation.fragmentation")) {
      ERROR_MSG("Parameter \'fragmentation.fragmentation\' not specified!!");
      return false;
    }

    string Mult_ID = "Merge";
    string CID_ID = "CID";
    string HCD_ID = "HCD";
    string ETD_ID = "ETD";
    string Frag_ID = m_params.getValue("fragmentation.fragmentation");
    if (Frag_ID == Mult_ID) {
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("MERGE_SAME_PREC", "1");
      m_params.setValue("PRM_CLUSTERING_NO_OP", "0");
    }
    else if (Frag_ID == CID_ID) {
      m_params.setValue("ACTIVATION", "CID");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("PRM_CLUSTERING_NO_OP", "1");
    }
    else if (Frag_ID == HCD_ID) {
      m_params.setValue("ACTIVATION", "HCD");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("PRM_CLUSTERING_NO_OP", "1");
    }
    else if (Frag_ID == ETD_ID) {
      m_params.setValue("ACTIVATION", "ETD");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("PRM_CLUSTERING_NO_OP", "1");
    } else {
      ERROR_MSG("Unrecognized value \'" << Frag_ID << " for parameter \'fragmentation.fragmentation\'!!");
      return false;
    }

    const string pmTolParam = "tolerance.PM_tolerance";
    float pmTol = m_params.getValueFloat(pmTolParam);
    const string pmTolUnitParam = "tolerance_unit.PM_unit";
    const string ionTolParam = "tolerance.Ion_tolerance";
    const string ionTolUnitParam = "tolerance_unit.Ion_unit";

    if (!m_params.exists(pmTolParam)) {
      ERROR_MSG("Parameter \'" << pmTolParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(pmTolUnitParam)) {
      ERROR_MSG("Parameter \'" << pmTolUnitParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(ionTolParam)) {
      ERROR_MSG("Parameter \'" << ionTolParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(ionTolUnitParam)) {
      ERROR_MSG("Parameter \'" << ionTolUnitParam << "\' not specified!!");
      return false;
    }

    string pmTolUnits = m_params.getValue(pmTolUnitParam);
    const float ppmFac = 1000000.0;

    std::transform(pmTolUnits.begin(),
                   pmTolUnits.end(),
                   pmTolUnits.begin(),
                   ::tolower);

    if (pmTolUnits == "da")  {
      m_params.setValue("TOLERANCE_PM", parseFloat(pmTol, 5));
    }
    else if (pmTolUnits == "ppm") {
      m_params.setValue("TOLERANCE_PM_PPM", parseFloat(pmTol, 5));
      m_params.setValue("TOLERANCE_PM",
                        parseFloat(((pmTol / ppmFac) * 2000.0), 5));
    } else {
      ERROR_MSG("Unknown tolerance type \'" << pmTolUnits << "\' !!");
      return false;
    }

    float ionTol = m_params.getValueFloat(ionTolParam);
    string ionTolUnits = m_params.getValue(ionTolUnitParam);

    std::transform(ionTolUnits.begin(),
                   ionTolUnits.end(),
                   ionTolUnits.begin(),
                   ::tolower);

    if (ionTolUnits == "da") {
      m_params.setValue("TOLERANCE_PEAK", parseFloat(ionTol, 5));
    } else if (ionTolUnits == "ppm") {
      m_params.setValue("TOLERANCE_PEAK_PPM", parseFloat(ionTol, 5));
      m_params.setValue("TOLERANCE_PEAK",
                        parseFloat(((ionTol / ppmFac) * 2000.0), 5));
    } else {
      ERROR_MSG("Unknown tolerance type \'" << ionTolUnits << "\' !!");
      return false;
    }

    if (!m_params.exists("Default.complexity")) {
      ERROR_MSG("Parameter \'Default.complexity\' not specified!!");
      return false;
    }

    const string Simple_ID = "purified simple mixture";
    const string Complex_ID = "complex mixture";
    string complexity_ID = m_params.getValue("Default.complexity");

#if 0
    if (complexity_ID == Simple_ID) {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "1");
      m_params.setValue("PARTIAL_OVERLAPS", "1");
      m_params.setValue("MAX_PVALUE", "0.05");
      m_params.setValue("CLUSTALW_MINSCORE", "500");
    } else if (complexity_ID == Complex_ID) {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "2");
      m_params.setValue("PARTIAL_OVERLAPS", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("CLUSTALW_MINSCORE", "50000");
    }
#endif

    m_params.setValue("SPEC_TYPE_MSMS", "0");

    if (m_params.getValueBool("custom_preprocess.Correct_Charge")) {
      m_params.setValue("GUESS_CHARGE", "yes");
      m_params.setValue("FIX_CHARGE_ZEROS","0");
    } else {
      m_params.setValue("GUESS_CHARGE", "no");
      m_params.setValue("FIX_CHARGE_ZEROS","1");
    }

    // ExecFilterPairs
    m_params.setValue("MIN_SHIFT", "0");
    m_params.setValue("USE_MIN_DIST_57", "1");

    // ExecAssembly
    m_params.setValue("NO_SEQUENCING", "0");
    m_params.setValue("PARALLEL_PATHS", "0");
    
    addParamFromProteoSAFe("DECONV_MS2",
                           "custom_preprocess.Deconvolute_MS",
                           true);
    addParamFromProteoSAFe("MIN_SPECTRUM_QUALITY",
                           "custom_preprocess.Min_Spec_Quality",
                           false);


    bool peformMetaAssembly = m_params.getValueBool("custom_preprocess.Perform_MetaAssembly");
    m_params.setValue("META_ASSEMBLY_NO_OP", peformMetaAssembly ? "0" : "1");

    //Reading File Names
    int mapping_count = 0;
    map<string, string> file_name_mapping;
    std::vector<std::string> search_spectra_names;
    std::vector<std::string> original_names;
    while (1)
    {
      char buf[100];
      sprintf(buf, "upload_file_mapping%i", mapping_count);

      mapping_count++;

      if (!m_params.exists(buf))
        break;

      std::string mapping = m_params.getValue(buf);
      std::string mangled_name = mapping.substr(0, mapping.find("|"));
      std::string original_name = mapping.substr(mapping.find("|") + 1);

      if (mangled_name.find("spec") == string::npos)
        continue;

      file_name_mapping[original_name] = mangled_name;

      std::string path_to_spectra = m_params.getValue("SPECTRA_DIR", "");
      search_spectra_names.push_back(path_to_spectra + "/" + mangled_name);
      original_names.push_back(original_name);
    }

    //Constructing input spectra
    string all_input_spectra_specnetsformat = "";
    string set_filename_input = "";
    for (int i = 0; i < search_spectra_names.size(); i++) {
      all_input_spectra_specnetsformat += search_spectra_names[i];
      all_input_spectra_specnetsformat += ";";
      set_filename_input += original_names[i];
      set_filename_input += ";";
    }

    m_params.setValue("INPUT_SPECS_MS",
                      all_input_spectra_specnetsformat.c_str());

//    m_params.setValue("SET_FILENAMES", set_filename_input.c_str());


    parsePtmSpecifications();
    
    // SET DEFAULT VALUES FOR ALL PARAMETERS THAT HAVE NOT BEEN SET
    addDefaultParameterValues();

    m_params.print(cerr);
    
    // CREATE DUMMY FILES (THIS IS PROBABLY TEMPORARY)
    if (m_params.exists("CREATE_DUMMY_RESULT_SPECNETS")) {
      string filename = m_params.getValue("CREATE_DUMMY_RESULT_SPECNETS");
      ofstream of(filename.c_str(), ios::out);
      if (!of) {
        ERROR_MSG("Unable to create dummy file [" << m_params.getValue("CREATE_DUMMY_RESULT_SPECNETS") << "]");
        return false;
      }
      
      of << "#Scan#	SpectrumFile	Annotation	OrigAnnotation	Protein	dbIndex	numMods	matchOrientation	startMass	Charge	MQScore	p-value	isDecoy	StrictEnvelopeScore	UnstrictEvelopeScore	CompoundName	Organism	FileScanUniqueID	FDR	LibraryName	mzErrorPPM	LibMetaData	Smiles	Inchi	LibSearchSharedPeaks	Abundance	ParentMassDiff	SpecMZ	ExactMass	LibrarySpectrumID" << endl;
      of.close();
    }
    if (m_params.exists("CREATE_DUMMY_RESULT_SPECNETS_DB")) {
      string filename = m_params.getValue("CREATE_DUMMY_RESULT_SPECNETS_DB");
      ofstream of(filename.c_str(), ios::out);
      if (!of) {
        ERROR_MSG("Unable to create dummy file [" << m_params.getValue("CREATE_DUMMY_RESULT_SPECNETS_DB") << "]");
        return false;
      }
      of << "Compound_Name	Ion_Source	Instrument	Compound_Source	PI	Data_Collector	Adduct	Scan	Precursor_MZ	ExactMass	Charge	CAS_Number	Pubmed_ID	Smiles	INCHI	INCHI_AUX	Library_Class	SpectrumID	IonMode	UpdateWorkflowName	LibraryQualityString	TaskID	#Scan#	SpectrumFile	LibraryName	MQScore	Organism	TIC_Query	RT_Query	MZErrorPPM	SharedPeaks	MassDiff	LibMZ	SpecMZ	SpecCharge" << endl;
      of.close();
    }
    if (m_params.exists("CREATE_DUMMY_RESULT_MSPLIT")) {
      string filename = m_params.getValue("CREATE_DUMMY_RESULT_MSPLIT");
      ofstream of(filename.c_str(), ios::out);
      if (!of) {
        ERROR_MSG("Unable to create dummy file [" << m_params.getValue("CREATE_DUMMY_RESULT_MSPLIT") << "]");
        return false;
      }
      of << "#SpectrumFile	Scan#	Annotation	Protein	Charge	cosine(M, A+B)	cosine(M,A)	cosine(A,B)	alpha	res-alpha	#Peak-0.85Intensity	simBias(M,A+B)	simBias(A)	projCos(M,A+B)	projCos(M,A)	meanCos	meanDeltaCos	Precusor(M)	Precursor(A)	spectrumIndex	libIndex1	libIndex2	svm1-score	svm2-score	IsUnique	Proteins" << endl;
      of.close();
    }
    if (m_params.exists("CREATE_DUMMY_RESULT_MSGFDB")) {
      string filename = m_params.getValue("CREATE_DUMMY_RESULT_MSGFDB");
      ofstream of(filename.c_str(), ios::out);
      if (!of) {
        ERROR_MSG("Unable to create dummy file [" << m_params.getValue("CREATE_DUMMY_RESULT_MSGFDB") << "]");
        return false;
      }
      of << "#SpecFile	SpecIndex	Scan#	FragMethod	Precursor	PMError(Da)	Charge	Peptide	Protein	DeNovoScore	MSGFScore	SpecProb	P-value	FDR	PepFDR	IsUnique	Proteins" << endl;
      of.close();
    }
    if (m_params.exists("CREATE_DUMMY_RESULT_MODA")) {
      string filename = m_params.getValue("CREATE_DUMMY_RESULT_MODA");
      ofstream of(filename.c_str(), ios::out);
      if (!of) {
        ERROR_MSG("Unable to create dummy file [" << m_params.getValue("CREATE_DUMMY_RESULT_MODA") << "]");
        return false;
      }
      of << "SpectrumFile	Index	ScanNo	ObservedMW	Charge	CalculatedMW	DeltaMass	Score	Probability	Peptide	Protein	PeptidePosition	ScanOffset	IsUnique	Proteins" << endl;
      of.close();
    }
    
    int numNodes = m_params.getValueInt("NUM_SPLIT_MSGFDB", 1);
    string splitStr;
    for (int i = 0; i < numNodes; i++) {
      splitStr += "spec_";
      char buf[128];
      sprintf(buf, "%d", i);
      splitStr += buf;
      splitStr += ".mgf;";
    }
    m_params.setValue("MERGE_CONVERT_OUTPUT_SPECTRA", splitStr);

    return true;
  }

  bool ExecProteosafeToSpecnetsParams::loadInputData(void)
  {

    return true;
  }

  bool ExecProteosafeToSpecnetsParams::saveOutputData(void)
  {
    string paramOutfile = m_params.getValue("RESULTS_DIR");
    DEBUG_MSG("Saving..." <<  "[" << paramOutfile << "]");
    m_params.writeToFile(paramOutfile);

    string knownModsFile = m_params.getValue("KNOWN_MODS_FILE");
    DEBUG_MSG("Saving..." <<  "[" << knownModsFile << "]");
    m_penaltyMatrix->saveKnownMods(knownModsFile);

    string penaltyFile = m_params.getValue("CLEAVAGE_PENALTY_FILE");
    DEBUG_MSG("Saving..." << "[" << penaltyFile << "]");
    m_penaltyMatrix->saveCleavagePenalties(penaltyFile);

    return true;
  }

  bool ExecProteosafeToSpecnetsParams::saveInputData(std::vector<std::string> & filenames)
  {

    return true;
  }

  bool ExecProteosafeToSpecnetsParams::loadOutputData(void)
  {

    return true;
  }

  std::vector<ExecBase *> const & ExecProteosafeToSpecnetsParams::split(int numSplit)
  {

    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecProteosafeToSpecnetsParams::merge(void)
  {

    return true;
  }

  bool ExecProteosafeToSpecnetsParams::validateParams(std::string & error)
  {
    return true;
  }

  void ExecProteosafeToSpecnetsParams::parsePtmSpecifications(void)
  {
    float pmTol = m_params.getValueFloat("tolerance.PM_tolerance");
    float minAdjustment = m_params.getValueInt("MIN_PARENT_MASS_ADUSTMENT");
    float pmAdjustment = max(pmTol, minAdjustment);
    
    if (m_params.getValueBool("NTERM_PARENT_MASS_ADUSTMENTS")) {
      // These are nterm parent mass adjustment mods for when parent mass is off
      for (float a = 1.0; a <= pmAdjustment; a += 1.0) {
        m_penaltyMatrix->addKnownNtermModification("<", a);
        m_penaltyMatrix->addKnownNtermModification("<", -a);
      }
    }
    //    NOTE: Let's always allow parent mass shifts at cterm of at least 1.0
    for (float a = 1.0; a <= pmAdjustment; a += 1.0) {
      m_penaltyMatrix->addKnownCtermModification(">", a);
      m_penaltyMatrix->addKnownCtermModification(">", -a);
    }
          
    DEBUG_VAR(m_params.getValue("ptm.DEAMIDATION"));
    DEBUG_VAR(m_params.getValueBool("ptm.DEAMIDATION"));
    
    if (m_params.getValueBool("ptm.OXIDATION")) {
      m_penaltyMatrix->addKnownModification("M", 15.994915);
    }
    if (m_params.getValueBool("ptm.LYSINE_METHYLATION")) {
      m_penaltyMatrix->addKnownModification("K", 14.015650);
    }
    if (m_params.getValueBool("ptm.PYROGLUTAMATE_FORMATION")) {
      m_penaltyMatrix->addKnownNtermModification("Q", -17.026549);
    }
    if (m_params.getValueBool("ptm.PHOSPHORYLATION")) {
      m_penaltyMatrix->addKnownModification("STY", 79.966331);
    }
    if (m_params.getValueBool("ptm.NTERM_CARBAMYLATION")) {
      m_penaltyMatrix->addKnownNtermModification("<", 43.005814);
    }
    if (m_params.getValueBool("ptm.NTERM_ACETYLATION")) {
      m_penaltyMatrix->addKnownNtermModification("<", 42.010565);
    }
    if (m_params.getValueBool("ptm.DEAMIDATION")) {
      m_penaltyMatrix->addKnownModification("NQ", 0.984016);
    }

    map<std::string, std::string> groups;
    string group("ptm.custom_PTM");
    m_params.getGroups(groups, group);
    //DEBUG_VAR(groups.size());
    
    string delimiters(",");
    map<std::string, std::string>::iterator itrGroup = groups.begin();
    map<std::string, std::string>::iterator itrGroupEnd = groups.end();
    for ( ; itrGroup != itrGroupEnd; itrGroup++) {
      //DEBUG_VAR(itrGroup->first);
      vector<string> subStrings;
      stringSplit(itrGroup->second, subStrings, delimiters);
      float mass;
      sscanf(subStrings[0].c_str(), "%f", &mass);
      if (subStrings[2] == "opt") {
        //DEBUG_MSG("Mod: " << subStrings[1] << "  " << mass);
        m_penaltyMatrix->addKnownModification(subStrings[1], mass);
      }
      if (subStrings[2] == "opt_nterm") {
        //DEBUG_MSG("Nterm: " << subStrings[1] << "  " << mass);
        m_penaltyMatrix->addKnownNtermModification(subStrings[1], mass);
      }
    }

    if (m_params.getValue("cysteine_protease.protease") == "Trypsin") {
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_START, 'K', -0.5);
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_START, 'R', -0.5);
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_INTERNAL, 'K', -0.5);
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_INTERNAL, 'R', -0.5);
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_END, 'K', -0.5);
      m_penaltyMatrix->addCleavagePenalty(CLEAVAGE_END, 'R', -0.5);
    }

    return;
  }

  void ExecProteosafeToSpecnetsParams::addParamFromProteoSAFe(const string & spsParam,
                                                             const string & proteosafeParam,
                                                             const bool mapYesNoToInt)
  {
    string paramVal = m_params.getValue(proteosafeParam, "default");

    if (paramVal != "default")
    {
      if (mapYesNoToInt && paramVal == "yes")
      {
        m_params.setValue(spsParam, "1");
      }
      else if (mapYesNoToInt && paramVal == "no")
      {
        m_params.setValue(spsParam, "0");
      }
      else if (!mapYesNoToInt)
      {
        m_params.setValue(spsParam, paramVal);
      }
      else
      {
        ERROR_MSG("Invalid parameter combo \'" << proteosafeParam << "=" << paramVal);
      }
    }
  }

  void ExecProteosafeToSpecnetsParams::addDefaultParameterValues(void)
  {
    // Basic parameters
    m_params.addIfDoesntExist("TOLERANCE_PEAK", "0.4");
    m_params.addIfDoesntExist("TOLERANCE_PM", "1.5");
    m_params.addIfDoesntExist("RESOLUTION", "0.1");

    // Preprocessing parameters
    m_params.addIfDoesntExist("CLUSTER_MIN_SIZE", "1");
    m_params.addIfDoesntExist("CLUSTER_MODEL", "LTQ_TRYP");
    m_params.addIfDoesntExist("MSCLUSTER_MIX_PROB ", "0.05");
    m_params.addIfDoesntExist("INSTRUMENT_TYPE", "IT");
    m_params.addIfDoesntExist("MIN_SPECTRUM_QUALITY", "0");
    m_params.addIfDoesntExist("CORRECT_PM", "no");
    m_params.addIfDoesntExist("GUESS_CHARGE", "no");
    m_params.addIfDoesntExist("PEPNOVO_OUTDIR", "spectra");

    // Alignment parameters
    m_params.addIfDoesntExist("AA_DIFF_COUNT", "2");
    m_params.addIfDoesntExist("MIN_SHIFT", "0");
    m_params.addIfDoesntExist("MIN_MOD_MASS", "-100");
    m_params.addIfDoesntExist("MAX_MOD_MASS", "100");
    m_params.addIfDoesntExist("MAX_NUM_MODS", "1");
    m_params.addIfDoesntExist("MIN_RATIO", "0.35");

    m_params.addIfDoesntExist("MAX_PVALUE", "0.045");
    m_params.addIfDoesntExist("MIN_MATCHED_PEAKS", "6"); // Minimum number of matched peaks to consider a spectral alignment
    m_params.addIfDoesntExist("MAX_AA_JUMP", "2");
    m_params.addIfDoesntExist("MIN_OVERLAP_AREA", "0.45");
    m_params.addIfDoesntExist("PENALTY_PTM", "-2000"); // Set to "0" for SpecNets
    m_params.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
    m_params.addIfDoesntExist("FILTER_TRIGS", "no"); // Set to "no" for SpecNets
    m_params.addIfDoesntExist("PARTIAL_OVERLAPS", "1"); // Set to "0" for SpecNets
    m_params.addIfDoesntExist("TAGS_MATCH_FLANK", "1");
    m_params.addIfDoesntExist("TAGS_MATCH_COUNT", "2");

    // Comparative Shotgun Protein Sequencing (CSPS) parameters
    m_params.addIfDoesntExist("CLUSTALW_MINSCORE", "1000000"); // Disabled by default - activate explicitly (i.e., set to ~250) for cSPS projects

    // De novo sequencing parameters
    m_params.addIfDoesntExist("SPSPATH_MIN_NUM_PEAKS", "5");
    m_params.addIfDoesntExist("ADD_ENDPOINTS", "1");
    m_params.addIfDoesntExist("SPSPATH_MIN_NUM_SPECS", "2");
    m_params.addIfDoesntExist("PARALLEL_PATHS", "0");
    m_params.addIfDoesntExist("SPS_MIN_EDGES_TO_COMPONENT", "1");
    m_params.addIfDoesntExist("MIN_METACONTIG_SCORE", "3.0");
    m_params.addIfDoesntExist("MIN_METACONTIG_SIZE", "1");


    // tagsearch/matchma parameters
    m_params.addIfDoesntExist("TAG_LEN", "6");
    m_params.addIfDoesntExist("DOUBLE_AA_JUMPS", "1");
    m_params.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0"); // Set to 2 for SpecNets
    m_params.addIfDoesntExist("MAX_NUM_TAGS", "0");
    m_params.addIfDoesntExist("MAX_NUM_MODS", "2");
    m_params.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7"); // Minimum number of matched peaks between spectrum/database to accept PSM
    m_params.addIfDoesntExist("TAG_MATCH_TOP_SCORING_ONLY", "1");

    // Grid parameters
#if 0    
    m_params.addIfDoesntExist("GRID_TYPE", "sge");
    m_params.addIfDoesntExist("GRID_NUMNODES", "-1");
    m_params.addIfDoesntExist("GRID_NUMCPUS", "1");
    m_params.addIfDoesntExist("GRID_EXE_DIR", "");
    m_params.addIfDoesntExist("GRID_SGE_EXE_DIR", "");
#endif    
  }

}
