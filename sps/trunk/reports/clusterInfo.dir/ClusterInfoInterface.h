///////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_INFO_INTERFACE_H__
#define __CLUSTER_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "ClusterData.h"
#include "SpecSet.h"

#include "ClusterSet.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace spsReports;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
  /*! \brief Interface class for Alignplot tool

   Provides command line user interface 

   */
class ClusterInfoInterface {


  // cluster data object
  //ClusterData clusterData;
  ClusterSet clusterData;

  // input files filenames (pklbin)
  vector<string> inputFiles;

  // original files filenames
  vector<string> originalFiles;

  // scan numbers filenames (.bin files)
  vector<string> scanFileNames;

  // the scan numbers
  vector<vector<vector<int> > > scanNumbers;

  //
  string binFiles, projectdir;

  // input filenames of filename containers
  string binFilesNames, inputFilesNames, originalFilesNames;

  // output file and directory
  string outdir, outFileName, outSummaryFilename;

  // Input data (specsets)
  vector<SpecSet> inputData;

  // specs_ms file (pklbin)
  SpecSet specs_ms;



  int  loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles);

    /*! \brief Writes a CSV file.

     Outputs a CSV file, based on an matrix, previously build by the ClusterData object.

     @return 0 if successful; -1 if there was an error
     */
  int writeCsvFile(ParameterList & commandLineParams,
                   vector<string> *inputFiles, 
                   vector<vector<vector<int> > > *scanNumbers, 
                   string &outFileName, 
                   vector<SpecSet> &specs);

    /*! \brief Write the cluster info information file.

     @return 0 if successful; -1 if there was an error
     */
  int WriteClusterSummary(std::map<string, string> &groups, string groupName, vector<string> *inputFiles);

 public:


  // Constructors and destructor
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Default contructor
     */
  ClusterInfoInterface();
    //@}

    //! \name DESTRUCTOR
    //@{
  ~ClusterInfoInterface();
    //@}

  // Option parsing
    /*! \brief Processes the command line params.

     Initializes the default parameter values and parses the command line against the declared parameters.

     @return 0 if successful; -1 if there was an error
     */
  int processOptions(int argc, char **argv);

  // Output help
    /*! \brief Help().

     Outputs program options.

     @return 0 if successful; -1 if there was an error
     */
  int help(ostream &);

  // Output version information
    /*! \brief version().

     Outputs program version information.

     @return 0 if successful; -1 if there was an error
     */
  int version(ostream &);

  // Output error messages
    /*! \brief display error message.

     Outputs an error message and exits program.

     @return 0 if successful; -1 if there was an error
     */
  int error(const string &);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
