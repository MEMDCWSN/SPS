////////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_DATA_H__
#define __CLUSTER_DATA_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <list>
#include <map>

#include "SpecSet.h"


////////////////////////////////////////////////////////////////////////////////

#define MS_CLUSTER_MGF_DIR "spectra/out/mgf"
#define MS_CLUSTER_OUT_DIR "spectra/out"
#define MS_CLUSTER_TEMP_DIR "spectra/tmp"
#define MS_CLUSTER_DATA_SUBPATH "spectra/out/clust"
#define MS_CLUSTER_FILE_PREFIX  "clusters_"
#define CLUSTER_DATA_FILENAME   "spectra/clusterData.bin"


using namespace std;
using namespace specnets;

namespace spsReports {

////////////////////////////////////////////////////////////////////////////////
// type ClusterData
//
// Defines the data structure for storing clusters data
//
// Where to load cluster information to
// vector index 1 is cluster #
// Unsigned 1 is spectrum file index
// Unsigned 2 is spectrum index
//
typedef map<unsigned,list<pair<unsigned,unsigned> > > ClusterData_t;

////////////////////////////////////////////////////////////////////////////////
// file clusterData.bin structure
//
// 1 unsinged int -> number of records
// <Records>
//
// Per record:
// 1 unsigned int -> record number
// 1 unsigned int -> number of entries for record
// <Entries>
//
// per entry:
// 1 unsigned int: file index
// 1 unsigned int: spectrum index
//
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ClusterData class

   Defines clusters by assotiating cluster indexes to spectra indexes.

   */
class ClusterData {

  /*! \brief Helper method to sort MsCluster file names
   */
  static bool clusterFilenameSort(const string &a, const string &b);

  /*! \brief Read MsCluster file names, given a directory and a filename prefix
   */
  int  getClusterFileNames(const string &location, const string &prefix, bool sort);

  /*! \brief add a record to summary record list
   */
  void addRecordToSummary(vector<vector<string> > &result, SpecSet &specs_ms, int cluster, int numSpectra, float sumPercursorIntensity, vector<int> &groupCount);

  /*! \brief Check if a given string in filename is equal to any of the substrings in fileGroups
   */
  int checkFileInGroup(string &fileGroup, string &filename);


 public:

  /*! \brief MsCluster clustering data structure
   */
   ClusterData_t data;

  /*! \brief file names of input MsCluster files
   */
   vector<string> fileNames;

  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Default constructor
   */
  ClusterData()  {fileNames.resize(0);};

    //@}

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
   */
  ~ClusterData() {};

  /*! \brief Load data form MsCluster output files.
   */
  int loadMsClusterData(const string &projectDir);

  /*! \brief Load data stored in binary format, but load it as if the data appeared twice
   */
  int loadDataDouble(const string &projectDir);

  /*! \brief Load data stored in binary format
   */
  int loadData(const string &projectDir);

  /*! \brief Output file in binary (.bin) format
   */
  int saveData(const string &projectDir);

  /*! \brief General load method
   */
  int load(const string &projectDir, bool saveBinary = true, bool rebuild = false);

  /*! \brief check if the file is already present. Returns 0 if the file doesn't exist, 1 if it does
   */
  int check(const string &projectDir);

  /*! \brief cluster count method
   */
  int getClusterCount(void)  {return data.size();};

  /*! \brief spectra count methods
   */
  int getSpectraCount(void);

  /*! \brief mapping methods
   */
  int  getConsensusFromInputSpectra(int fileIndex, int inputSpectra);

  /*! \brief Returns the list of pais (fileindex ; cluster index) of spectra that map to a given cluster
   */
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(int consensus);

  /*! \brief output cluster data to cout
   */
  void dump(ostream &sout, bool web);

  /*! \brief write file in csv format
   */
  int writeCsv(vector<string> *inputFiles, vector<vector<vector<int> > > *scanNumbers, string &outFileName, vector<SpecSet> &specs);

  /*! \brief write file in csv format
   */
  int getCsvData(vector<string> *inputFiles, vector<vector<vector<int> > > *scanNumbers, vector<SpecSet> &specs, vector<vector<string> > &data);

  /*! \brief get the summary
   */
  int getSummary(vector<vector<string> > &result, vector<SpecSet> &specs, SpecSet &specs_ms, vector<string> *inputFiles, std::map<string, string> &groups);

};
////////////////////////////////////////////////////////////////////////////////
};  // namespace sps
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
