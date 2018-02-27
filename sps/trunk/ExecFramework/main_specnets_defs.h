#include <string>

// System Defines
#define TEST_VALID {          \
  DEBUG_VAR(isValid);         \
  if (!isValid) {             \
    ERROR_MSG(errorString);   \
    return false;             \
  }                           \
}

#define TEST_RETURN(__module, __var) {                                                    \
  if(__var.size() == 0) {                                                               \
    ERROR_MSG("invoking " << __module << ": empty return set \"" << #__var << "\".");     \
  /*  return false; */                                                                     \
  }                                                                                   \
}

#define TEST_RETURN_STATUS(__module) {                            \
  DEBUG_VAR( returnStatus);                                     \
  if (!returnStatus) {                                          \
    ERROR_MSG("invoking " << __module << " exited in error.");    \
    return false;                                              \
  }                                                             \
}

#define TEST_SAVE_OUPUT_DATA(__module) {                                      \
  DEBUG_VAR(returnStatus);                                                  \
  if (!returnStatus) {                                                      \
    ERROR_MSG("Saving output data for " << __module << " exited in error.");  \
    return false;                                                           \
  }                                                                         \
}

namespace specnets
{
  using namespace std;
  
  const string DEFAULT_DECONV_MODEL = "model_isoenv.bin";
  const string DEFAULT_DECONV_FILES_LIST = "deconv_files.txt";
  const string DEFAULT_BIN_FILES_FILENAME = "bin_files.txt";
  const string DEFAULT_INPUT_FILES_LIST = "spectra/pklbin_files.txt";
  const string DEFAULT_USER_INPUT_FILES_LIST = "input_index.txt";
  const string DEFAULT_INPUT_FILE_BASE = "specs_ms_";

  const string DEFAULT_INPUT_MAPPING = "input_mapping.bin";

  const string CLUSTER_MSCLUST = "MSCluster";
  const string CLUSTER_PRMS = "PrmClust";

  const string SPECTRA_DIR = "spectra";
  const string SPECTRA_GRID_SCORING_DIR = "spectra/grid_scoring";

  const string SPECS_MS_FILE = "specs_ms.pklbin";
  const string SPECS_MS_MGF_FILE = "specs_ms.mgf";
  const string SPECS_MS_CLUST_FILE = "specs_ms.clust";

  const string SPECS_SCORED_FILE = "specs_scored.pklbin";
  const string SPECS_SCORED_MGF_FILE = "specs_scored.mgf";
  const string SPECS_SCORED_CLUST_FILE = "specs_scored.clust";
  const string SPECS_SCORED_UNCLUST_FILE = "specs_scored_unclustered.pklbin";
  const string SPECS_PEPNOVO_MGF_FILE = "specs_ms_pepnovo.mgf";
  const string SPECS_PEPNOVO_PRMS_FILE = "specs_scored.prms";

  const string ALIGNS_DIR = "aligns";
  const string ALIGNS_SILAC_FILE = "pairs_silac.bin";
  
  enum Stage
  {
    STAGE_BEGIN = 0,
    STAGE_MS2DECONV = 1,
    STAGE_MSCLUSTER = 2,
    STAGE_SCORING = 3,
    STAGE_PRMCLUSTERING = 4,
    STAGE_GENOMS = 5,
    STAGE_FILTERPAIRS = 6,
    STAGE_FILTERALIGNS = 7,
    STAGE_ALIGNMENT = 8,
    STAGE_FILTERSTARPAIRS = 9,
    STAGE_PENALTYGENERATE = 10,
    STAGE_ASSEMBLY = 11,
    STAGE_METAASSEMBLY = 12,
    STAGE_SPECNETS = 13,
    STAGE_CONTIGPROTALIGN = 14,
    STAGE_SPECTAGGENERATE = 15,
    STAGE_SPECPROTALIGN = 16,
    STAGE_PROTPROTALIGN = 17,
    STAGE_HOMOLOGYASSEMBLY = 18,
    STAGE_MERGE = 19,
    STAGE_STATPROTSEQS = 20,
    STAGE_REPORT = 21
  };
  
}


