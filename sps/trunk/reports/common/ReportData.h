////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_DATA_H__
#define __REPORT_DATA_H__
////////////////////////////////////////////////////////////////////////////////
#include "aminoacid.h"
#include "Defines.h"

////////////////////////////////////////////////////////////////////////////////
#define DEFAULT_ROOT_DIRECTORY      "."
#define DEFAULT_TABLES_DIR          "./report"
#define DEFAULT_OUTPUT_DIR          "./report"
#define FILE_CLUSTERMS              "spectra/pklbin_files.txt"
#define FILE_SCAN_FILES             "spectra/bin_files.txt"
#define FILE_CONSENSUS_SPECTRA      "spectra/specs_ms.pklbin"
#define FILE_PROTEINS               "protid.fasta"
#define FILE_CONTIG_INDICES         "spectra/contigs_indices.bin"
#define FILE_HOMGLUE_REF_MIDX       "homology/homglue_ref_midx.pklbin"
#define FILE_HOMGLUE_REF_MP         "homology/homglue_ref_mp.bin"
#define FILE_CONTIGS_MIDX           "homology/contigs_midx.pklbin"
#define FILE_CONTIGS_MP             "homology/contigs_mp.bin"
#define FILE_SPS_SEQS               "assembly/sps_seqs.pklbin"
#define FILE_ABRUIJN                "assembly/component_info.bin"
#define FILE_STARS_SPECTRA          "spectra/stars.pklbin"
#define FILE_CONTIG_SPECTRA         "spectra/contigs.pklbin"
#define FILE_HOMGLUE_MATCHES_MP     "homology/homglue_matches_mp.bin"
#define FILE_HOMGLUE_MATCHES_MIDX   "homology/homglue_matches_midx.pklbin"
#define FILE_HOMGLUE_MATCHES        "homology/homglue_matches.pklbin"
#define FILE_INPUT_SPECTRA_LIST     "spectra/input_index.txt"
#define FILE_CONTIG_NAMES           "homology/ref_sps_names.txt"
#define FILE_INPUT_MAPPING          "spectra/input_mapping.bin"

#define DEFAULT_ANNOTATION_FILE_LOCATION "."

#define DEFAULT_MASS_SHIFT            specnets::AAJumps::massHion
#define DEFAULT_PEAK_MASS_TOLERANCE   0.45
#define DEFAULT_PARENT_MASS_TOLERANCE 1.0
#define DEFAULT_RESOLUTION            0.01
#define DEFAULT_CELLS_PER_LINE        20

#define DEFAULT_TOOL                  1

#define DEFAULT_JOB_NAME          "MyJob"
#define DEFAULT_USER_NAME         "Undefined"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports
{
////////////////////////////////////////////////////////////////////////////////
  /*! \brief ReportData container class

   Helper container class used to aggregate data used by report generater classes.
   Mainly used to reduce the overhead when passing a large number of parameters between
   class methods. Also minimizes changes needed when a new parameter or data member needs
   to be added, changed or removed.

   */
  class ReportData
  {

  public:

    /*! \brief Holder for the executables location.
     */
    string exeDir;

    // Directories and data file names

    /*! \brief Project location.
     */
    string projectDir;
    /*! \brief Directory for outputting the report tables.
     */
    string tablesDir;
    /*! \brief Directory for outputting the report fiels.
     */
    string outDir;
    /*! \brief URL for the server. This string will be embedded in the dynamic report pages.
     */
    string server;

    /*! \brief User password for the reports.
     */
    string pwd;
    string pwd2;

    // directories for report. Used for relocation.
    /*! \brief Project directory in case of relocation. This is used when the projec location directory during
     report generation is different of the project location when the report is under the web server.
     This string should contain the project location (path) under the web server.
     */
    string targetProjectDir;

    /*! \brief Number of cell per line to be used in the tables, in the protein coverage pages.
     */
    int cellsPerLine;

    /*! \brief Data detail dysplay level. Only two levels are currently supported: with and without showing the spectra images.
     */
    int displayLevel;

    /*! \brief Number of CPSs to be used to generate report pages. This is only valid for static reports with multiple pages.
     */
    int cpu;

    /*! \brief Tool used.
     */
    int tool;

    /*! \brief True if the clusters layer does not exists. In this case, the imput spectra maps directly to contigs.
     */
    bool noClusters;

    // allow realign (includes edit boxes) of protein coverage
    /*! \brief True if realign is allowed. This relates to the protein coverage pages. Input boxes for individual
     aminoacids will not be rendered in dynamic reports it this variable is false.
     */
    bool allowRealign;

    /*! \brief allow dynamic reports.
     */
    bool dynamic;

    /*! \brief
     */
    string annotationModelDir;

    /*! \brief
     */
    string annotationModel;

    /*! \brief
     */
    string aasFileDir;

    /*! \brief
     */
    string aasFile;

    /*! \brief
     */
    float massShift;

    /*! \brief
     */
    string annotationModelPrm;

    /*! \brief
     */
    float massShiftPrm;

    /*! \brief
     */
    string annotationModelDirPrm;

    /*! \brief
     */
    float peakMassTol;

    /*! \brief
     */
    float parentMassTol;

    /*! \brief
     */
    float resolution;

    /*! \brief
     */
    int tableID;

    /*! \brief
     */
    int filterField;

    /*! \brief
     */
    string filterData;

    /*! \brief
     */
    int updateField;

    /*! \brief
     */
    string updateData;

    /*! \brief
     */
    int sortColumnIdx;

    /*! \brief
     */
    int sortDirection;

    /*! \brief
     */
    int startRow;

    /*! \brief
     */
    int rows;

    /*! \brief
     */
    ReportData()
    {
      projectDir = DEFAULT_ROOT_DIRECTORY;
      tablesDir = DEFAULT_TABLES_DIR;
      outDir = DEFAULT_OUTPUT_DIR;
      cellsPerLine = DEFAULT_CELLS_PER_LINE;

      annotationModelDir = DEFAULT_ANNOTATION_FILE_LOCATION;
      annotationModel = DEFAULT_ANNOTATION_MODEL;

      massShift = DEFAULT_MASS_SHIFT;
      massShiftPrm = DEFAULT_MASS_SHIFT;

      peakMassTol = DEFAULT_PEAK_MASS_TOLERANCE;
      parentMassTol = DEFAULT_PARENT_MASS_TOLERANCE;
      resolution = DEFAULT_RESOLUTION;

      tableID = -1;
      filterField = -1;
      updateField = -1;
      displayLevel = 100;

      sortColumnIdx = -1;
      sortDirection = 0;

      startRow = 0;
      rows = -1;

      cpu = 1;

      tool = DEFAULT_TOOL;

      noClusters = false;

      allowRealign = true;
      dynamic = true;

    }

  };

////////////////////////////////////////////////////////////////////////////////
  /*! \brief
   */
  class ReportGeneratorData : public ReportData
  {

  public:

    /*! \brief
     */
    string job;

    /*! \brief
     */
    string user;

    /*! \brief
     */
    string filenameClusterMS;

    /*! \brief
     */
    string filenameCluster;

    /*! \brief
     */
    string filenameScanFiles;

    /*! \brief
     */
    string filenameConsensusSpectra;

    /*! \brief
     */
    string filenameProteins;

    /*! \brief
     */
    string filenameContigIndices;

    /*! \brief
     */
    string filenameHomglueRefMidx;

    /*! \brief
     */
    string filenameHomglueRefMp;

    /*! \brief
     */
    string filenameContigsMidx;

    /*! \brief
     */
    string filenameContigsMp;

    /*! \brief
     */
    string filenameSpsSeqs;

    /*! \brief
     */
    string filenameAbruijn;

    /*! \brief
     */
    string filenameStarSpectra;

    /*! \brief
     */
    string filenameContigSpectra;

    /*! \brief
     */
    string filenameHomglueMatchesMp;

    /*! \brief
     */
    string filenameHomglueMatchesMidx;

    /*! \brief
     */
    string filenameHomglueMatches;

    /*! \brief
     */
    string filenameInputSpectraList;

    /*! \brief
     */
    string filenameContigNames;

    /*! \brief
     */
    string filenameInputMapping;

    // constructor
    /*! \brief
     */
    ReportGeneratorData() :
        filenameClusterMS(FILE_CLUSTERMS), filenameScanFiles(FILE_SCAN_FILES), filenameConsensusSpectra(FILE_CONSENSUS_SPECTRA), filenameProteins(FILE_PROTEINS), filenameContigIndices(FILE_CONTIG_INDICES), filenameHomglueRefMidx(FILE_HOMGLUE_REF_MIDX), filenameHomglueRefMp(FILE_HOMGLUE_REF_MP), filenameContigsMidx(FILE_CONTIGS_MIDX), filenameContigsMp(FILE_CONTIGS_MP), filenameSpsSeqs(FILE_SPS_SEQS), filenameAbruijn(FILE_ABRUIJN), filenameStarSpectra(FILE_STARS_SPECTRA), filenameContigSpectra(FILE_CONTIG_SPECTRA), filenameHomglueMatchesMp(FILE_HOMGLUE_MATCHES_MP), filenameHomglueMatchesMidx(FILE_HOMGLUE_MATCHES_MIDX), filenameHomglueMatches(FILE_HOMGLUE_MATCHES), filenameInputSpectraList(FILE_INPUT_SPECTRA_LIST), filenameContigNames(FILE_CONTIG_NAMES), filenameInputMapping(FILE_INPUT_MAPPING),

        job(DEFAULT_JOB_NAME), user(DEFAULT_USER_NAME)

    {
    }
    ;

  };
////////////////////////////////////////////////////////////////////////////////
}
;
// namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
