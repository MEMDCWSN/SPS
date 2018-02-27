// Header Include
#include "ExecMsCluster.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "utils.h"
#include "ExecMergeConvert.h"

// System Includes
#include <stdio.h>
#include <string.h>

using namespace specnets;
using namespace std;
/*
 const string ExecMsCluster::TMP_DIRECTORY = "spectra/tmp";
 const string ExecMsCluster::OUT_DIRECTORY = "spectra/out";
 const string ExecMsCluster::MGF_DIRECTORY = "spectra/out/mgf";
 const string ExecMsCluster::CLUST_DIRECTORY = "spectra/out/clust";
 */

const string ExecMsCluster::TMP_DIRECTORY = "tmp";
const string ExecMsCluster::OUT_DIRECTORY = "out";
const string ExecMsCluster::MGF_DIRECTORY = "mgf";
const string ExecMsCluster::CLUST_DIRECTORY = "clust";
const string ExecMsCluster::CLUSTERS_FILE_PREFIX = "clusters_";

// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(void) :
    ownInput(true), m_inputSpectraPaths(0x0), ownOutput(true),
        m_outputClusters(0x0), m_outputClusteredSpectra(0x0)
{
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_inputSpectraPaths = new vector<string>;
  m_outputClusters = new ClusterSet;
  m_outputClusteredSpectra = new SpecSet;
}
// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(const ParameterList & params) :
    ExecBase(params), ownInput(true), m_inputSpectraPaths(0x0), ownOutput(true),
        m_outputClusters(0x0), m_outputClusteredSpectra(0x0)
{
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_inputSpectraPaths = new vector<string>;
  m_outputClusters = new ClusterSet;
  m_outputClusteredSpectra = new SpecSet;
}

ExecMsCluster::ExecMsCluster(const ParameterList & inputParams,
                             vector<string> *inputSpectraPaths) :
    ExecBase(inputParams), ownInput(false),
        m_inputSpectraPaths(inputSpectraPaths), ownOutput(true),
        m_outputClusters(0x0), m_outputClusteredSpectra(0x0)
{
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_outputClusters = new ClusterSet;
  m_outputClusteredSpectra = new SpecSet;
}

ExecMsCluster::ExecMsCluster(const ParameterList & inputParams,
                             vector<string> *inputSpectraPaths,
                             ClusterSet *outputClusters,
                             SpecSet *outputClusteredSpectra) :
    ExecBase(inputParams), ownInput(false),
        m_inputSpectraPaths(inputSpectraPaths), ownOutput(false),
        m_outputClusters(outputClusters),
        m_outputClusteredSpectra(outputClusteredSpectra)
{
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
}

// -------------------------------------------------------------------------
ExecMsCluster::~ExecMsCluster(void)
{
  if (ownInput)
  {
    delete m_inputSpectraPaths;
  }
  if (ownOutput)
  {
    delete m_outputClusters;
    delete m_outputClusteredSpectra;
  }
}
// -------------------------------------------------------------------------
ExecBase * ExecMsCluster::clone(const ParameterList & input_params) const
{
  return new ExecMsCluster(input_params);
}
// -------------------------------------------------------------------------
bool ExecMsCluster::invoke(void)
{

  DEBUG_MSG("Entering  ExecMsCluster::invoke()");

  // Get CLUSTER_MIN_SIZE's value
  int clusterMinSize = m_params.getValueInt("CLUSTER_MIN_SIZE");

  // If we are not doing clustering, just merge the input files
  if (clusterMinSize <= 0)
  {
    DEBUG_MSG("CLUSTER_MIN_SIZE is zero, skipping call to MSCluster");
    m_outputClusteredSpectra->resize(0);
    vector<int> locSpecIndicies;
    vector<int> locSpecScans;
    for (int i = 0; i < m_inputSpectraPaths->size(); i++)
    {
      SpecSet locSpecs;
      DEBUG_MSG("Merging spectra from \'" << (*m_inputSpectraPaths)[i] << "\' ...");
      if (!locSpecs.Load((*m_inputSpectraPaths)[i].c_str()))
      {
        ERROR_MSG("Failed to load \'" << (*m_inputSpectraPaths)[i] << "\'");
        return false;
      }
      for (int j = 0; j < locSpecs.size(); j++)
      {
        locSpecs[j].fileIndex = i;
        locSpecIndicies.push_back(j);
        locSpecScans.push_back(locSpecs[j].scan);
      }
      m_outputClusteredSpectra->swapAppendSpecSet(locSpecs, false);
    }
    m_outputClusters->resize(m_outputClusteredSpectra->size());
    for (int i = 0; i < m_outputClusteredSpectra->size(); i++)
    {
      (*m_outputClusters)[i].resize(1);
      (*m_outputClusters)[i].m_index = i;
      (*m_outputClusters)[i][0].m_index = locSpecIndicies[i];
      (*m_outputClusters)[i].m_scan = i + 1;
      (*m_outputClusteredSpectra)[i].scan = i + 1;
      //(*m_outputClusters)[i][0].m_scan = (*m_outputClusteredSpectra)[i].scan;
      (*m_outputClusters)[i][0].m_scan = locSpecScans[i];
      (*m_outputClusters)[i].m_filename =
          (*m_outputClusteredSpectra)[i].fileName;
      (*m_outputClusters)[i][0].m_filename =
          (*m_outputClusteredSpectra)[i].fileName;
      (*m_outputClusters)[i].m_fileIndex =
          (*m_outputClusteredSpectra)[i].fileIndex;
      (*m_outputClusters)[i][0].m_fileIndex =
          (*m_outputClusteredSpectra)[i].fileIndex;
    }
    return true;
  }

  const string baseDir = m_params.getValue("PROJECT_DIR", ".");
  const string tempDir = getPath(baseDir, TMP_DIRECTORY, false);

  if (mkdir_if_not_exist(tempDir.c_str()))
  {
    DEBUG_MSG("Made directory \'" << tempDir << "\'");
  }
  string msClusterInputFileList;
  vector<SpecSet> inputSpectrumInfo;
  formatMsClusterInput(tempDir, inputSpectrumInfo, msClusterInputFileList);

  const string outDir = getPath(baseDir, OUT_DIRECTORY, false);
  if (mkdir_if_not_exist(outDir.c_str()))
  {
    DEBUG_MSG("Made directory \'" << outDir << "\'");
  }

  const string mgfDir = getPath(outDir, MGF_DIRECTORY, false);
  if (mkdir_if_not_exist(mgfDir.c_str()))
  {
    DEBUG_MSG("Made directory \'" << mgfDir << "\'");
  }

  const string clustDir = getPath(outDir, CLUST_DIRECTORY, false);
  if (mkdir_if_not_exist(clustDir.c_str()))
  {
    DEBUG_MSG("Made directory \'" << clustDir << "\'");
  }

  // EXE_DIR
  std::string exeDir = m_params.getValue("EXE_DIR");
  // rtrim.
  rtrim(exeDir);

  // Call MsCluster
  int ret = callMsCluster(exeDir,
                          msClusterInputFileList,
                          tempDir,
                          outDir,
                          baseDir);
  if (ret != 0)
  {
    ERROR_MSG("Error executing MsCluster");
    ERROR_MSG("Check data is centroided and not profile");
    return false;
  }

  vector<string> outMgfDirectoryListing;
  directoryContents(mgfDir, "", ".mgf", outMgfDirectoryListing, true);
  DEBUG_VAR(outMgfDirectoryListing.size());
  m_outputClusteredSpectra->resize(0);
  for (unsigned int i = 0; i < outMgfDirectoryListing.size(); i++)
  {
    //DEBUG_VAR(outMgfDirectoryListing[i]);
    const string specFileOut = getPath(mgfDir,
                                       outMgfDirectoryListing[i],
                                       false);
    SpecSet locSpecsOut;
    DEBUG_VAR(specFileOut);
    if (!locSpecsOut.Load(specFileOut.c_str()))
    {
      ERROR_MSG("Failed to load spectra from \'" << specFileOut << "\'");
      return false;
    }
    m_outputClusteredSpectra->swapAppendSpecSet(locSpecsOut);
  }

  if (!formatMsClusterOutput(clustDir, inputSpectrumInfo))
  {
    return false;
  }

  //Rank filter the clustered spectra if specified
  int rankFiltK = m_params.getValueInt("CLUST_RANK_FILTER", -1);
  if (rankFiltK > 0)
  {
    DEBUG_MSG("Rank-filtering clustered MS/MS spectra with K = " << rankFiltK);
    for (unsigned int i = 0; i < m_outputClusteredSpectra->size(); i++)
    {
      (*m_outputClusteredSpectra)[i].rankFilterPeaks(rankFiltK);
    }
  }

  DEBUG_MSG("Removing folder \'" << clustDir << "\' ...");
  if (removeFolder(clustDir) != 0)
  {
    ERROR_MSG("Failed to remove folder \'" << clustDir << "\' ...");
  }

  DEBUG_MSG("Removing folder \'" << mgfDir << "\' ...");
  if (removeFolder(mgfDir) != 0)
  {
    ERROR_MSG("Failed to remove folder \'" << mgfDir << "\' ...");
  }

  DEBUG_MSG("Removing folder \'" << outDir << "\' ...");
  if (removeFolder(outDir) != 0)
  {
    ERROR_MSG("Failed to remove folder \'" << outDir << "\' ...");
  }

  DEBUG_MSG("Removing folder \'" << tempDir << "\' ...");
  if (removeFolder(tempDir) != 0)
  {
    ERROR_MSG("Failed to remove folder \'" << tempDir << "\'");
  }

  DEBUG_MSG("Exiting  ExecMsCluster::invoke()");

  //m_outputClusters->dump(cerr, false);

  return true;
}

// -------------------------------------------------------------------------
bool ExecMsCluster::loadInputData(void)
{
  if (!m_params.exists("INPUT_SPECS_MS"))
  {
    ERROR_MSG("Missing parameter INPUT_SPECS_MS!!");
    return false;
  }
  splitText(m_params.getValue("INPUT_SPECS_MS").c_str(),
            *m_inputSpectraPaths,
            ";");
  return true;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::saveOutputData(void)
{
  DEBUG_TRACE;

  string outDir = m_params.getValue("PROJECT_DIR");
  
  // Save concatenated SpecSet file
  if (m_params.exists("OUTPUT_SPECTRA")) {
    string fileName = m_params.getValue("OUTPUT_SPECTRA");
    DEBUG_MSG("Saving concatenated SpecSet file to " << fileName );
    if (!ExecMergeConvert::saveSpecsetMultiple(outDir, 
                                               fileName, 
                                               m_outputClusteredSpectra)) {
      ERROR_MSG("Failed to save to [" << outDir << "] [" << fileName << "]");
      return false;
    }
  }

  if (m_params.exists("OUTPUT_CLUSTERS")) {
    string fileName = m_params.getValue("OUTPUT_CLUSTERS");
    fileName = outDir + "/" + fileName;
    DEBUG_MSG("Saving Cluster data to " << fileName );
    if (!m_outputClusters->saveBinaryFile(fileName)) {
      ERROR_MSG("Failed to save to [" << fileName << "]");
      return false;
    }
  }

  return true;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::saveInputData(std::vector<std::string> & filenames)
{
  return false;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::loadOutputData(void)
{
  return false;
}
// -------------------------------------------------------------------------
vector<ExecBase*> const & ExecMsCluster::split(int numSplit)
{
  m_subModules.resize(0);
  return m_subModules;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::merge(void)
{
}
// -------------------------------------------------------------------------
bool ExecMsCluster::validateParams(std::string & error)
{
  m_isValid = false;

  VALIDATE_PARAM_EXIST("EXE_DIR");
  VALIDATE_PARAM_EXIST("PROJECT_DIR");
  VALIDATE_PARAM_EXIST("CLUSTER_MIN_SIZE");
  VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
  VALIDATE_PARAM_EXIST("TOLERANCE_PM");

  m_isValid = true;
  return true;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::formatMsClusterInput(const string &tempDir,
                                         vector<SpecSet> &inputSpectrumInfo,
                                         string &msClusterInputFileList)
{

  if (m_inputSpectraPaths->size() == 0)
  {
    ERROR_MSG("Input cluster file list not specified!!");
    return false;
  }
  string fileLines(""); // List of input files to MSCluster

  // Prepare vector of spectrum headers
  inputSpectrumInfo.resize(m_inputSpectraPaths->size());

  for (unsigned int fIdx = 0; fIdx < m_inputSpectraPaths->size(); fIdx++)
  {
    SpecSet locSpecs;

    // load each input specset one file at a time
    if (!ExecMergeConvert::loadSpecset(m_params.getValue("EXE_DIR"),
                                       (*m_inputSpectraPaths)[fIdx],
                                       &locSpecs))
    {
      return false;
    }

    inputSpectrumInfo[fIdx].resize(locSpecs.size());
    for (unsigned int i = 0; i < locSpecs.size(); i++)
    {
      // Copy header information
      inputSpectrumInfo[fIdx][i].copyNP(locSpecs[i]);
      inputSpectrumInfo[fIdx][i].resize(0);

      // Do some weird stuff so MSCluster does nto crash
      if (locSpecs[i].parentMass > 5000)
      {
        WARN_MSG("MSCluster hack for spectrum " << i << " in " << (*m_inputSpectraPaths)[fIdx]);
        float z = locSpecs[i].parentCharge;
        if (z > 0)
        {
          WARN_MSG("converting to m/z with fake charge 0");
          locSpecs[i].parentMass = (locSpecs[i].parentMass
              + AAJumps::massHion * (z - 1)) / z;
          locSpecs[i].parentCharge = 1;
        }
        else
        {
          WARN_MSG("deleting spectrum");
          locSpecs[i].parentMass = 0;
          locSpecs[i].resize(0);
        }
      }
      else
      {
        inputSpectrumInfo[fIdx][i].parentMass = -1.0;
      }
    }

    // Prepare temporary specset file
    FilenameManager mngr((*m_inputSpectraPaths)[fIdx]);
    string tempSpecFile = getPath(tempDir, mngr.filename, false);
    tempSpecFile += "-format-prec.mgf";

    if (!ExecMergeConvert::saveSpecset(tempSpecFile, &locSpecs, false))
    {
      return false;
    }
    fileLines += tempSpecFile;
    fileLines += "\n";
  }

  msClusterInputFileList = getPath(tempDir, "MSCluster_input_files.txt", false);

  // Write out list of input files to MSCluster
  ofstream myfileOut(msClusterInputFileList.c_str(), ios::binary);
  if (myfileOut.is_open())
  {
    myfileOut << fileLines;
    myfileOut.close();
  }
  else
  {
    ERROR_MSG("Could not write to file " << msClusterInputFileList);
    return false;
  }

  return true;
}

bool ExecMsCluster::formatMsClusterOutput(const string &clustDir,
                                          vector<SpecSet> &inputSpectrumInfo)
{
  vector<string> fileNames;
  // Get msCluster output file names. If error, exit
  if (!directoryContents(clustDir, CLUSTERS_FILE_PREFIX, "", fileNames, false))
  {
    return false;
  }

  DEBUG_TRACE;

  std::sort(fileNames.begin(),
            fileNames.end(),
            ExecMsCluster::clusterFilenameSort);

  m_outputClusters->resize(0);

  unsigned int clusterSpecIdx = 0;
  // Cycle thru the found MsCLuster files
  for (unsigned i = 0; i < fileNames.size(); ++i)
  {
    // File reader object
    ifstream currentFile;
    // compose file name to include project directory
    string aux2 = getPath(clustDir, fileNames[i], false);
    // open the file
    currentFile.open(aux2.c_str(), ifstream::in | ifstream::binary);
    // test for errors
    if (!currentFile)
    {
      ERROR_MSG("Cannot open \'" << aux2 << "\'");
      return false;
    }

    //
    for (unsigned j = 0; currentFile.peek() != EOF; ++j)
    {
      string sline, sname;
      char cplot;
      unsigned nclusteridx, nspectrum3;

      getline(currentFile, sname, '.');
      currentFile >> nclusteridx >> cplot >> nspectrum3;
      getline(currentFile, sname);

      //ClusterData_t::iterator ic = data.insert(data.end(), make_pair(nspectrum, list< pair<unsigned, unsigned> >()));
      Cluster cluster;
      cluster.m_index = clusterSpecIdx;
      //DEBUG_VAR(cluster.m_index);

      while (currentFile && currentFile.peek() != '\n')
      {
        unsigned nclusteridx2, nspectrum2, aux;

        currentFile >> aux >> nclusteridx2 >> nspectrum2;
        getline(currentFile, sline);

        //ic->second.push_back(make_pair(nclusteridx2, nspectrum2));

        //DEBUG_MSG("adding " << nclusteridx2 << "," << nspectrum2);

        // Initialize the cluster with just the file index and local spectrum index
        cluster.add(nspectrum2, -1, nclusteridx2, "");

        //DEBUG_VAR(cluster[cluster.size()-1].m_index);

      }
      //if ( k == 0 ) {
      //  data.erase(ic);
      //} //else

      if (cluster.size() > 0)
      {
        // Add cluster and increment the index counter
        m_outputClusters->push_back(cluster);
        clusterSpecIdx++;
      }
    }
  }

  DEBUG_TRACE;

  if (m_outputClusters->size() != m_outputClusteredSpectra->size())
  {
    ERROR_MSG("Clusters and clustered spectra do not have the same size!!");
    DEBUG_VAR(m_outputClusters->size());
    DEBUG_VAR(m_outputClusteredSpectra->size());
    return false;
  }

  DEBUG_TRACE;

  string outFilename("");
  if (m_params.exists("OUTPUT_SPECTRA"))
  {
    vector<string> files;
    splitText(m_params.getValue("OUTPUT_SPECTRA").c_str(), files, ";");
    FilenameManager mngr(files[0]);
    outFilename = mngr.getFilenameWithExtension();
  }

  int min_clust_size = m_params.getValueInt("CLUSTER_MIN_SIZE", 1);
  //DEBUG_VAR(m_combinedSpecSet->size());
  //DEBUG_VAR(m_inputSpectra.size());
  for (unsigned int i = 0; i < m_outputClusters->size(); i++)
  {
    int oldFileIdx = -1;
    int oldSpecIdx = -1;

    for (unsigned int j = 0; j < (*m_outputClusters)[i].size(); j++)
    {
      // Get the file index
      unsigned int fileIdx = (*m_outputClusters)[i][j].m_fileIndex;
      // Get local spectrum index
      unsigned int specIdx = (*m_outputClusters)[i][j].m_index;

      // Remember first un-clustered spectrum so we can copy its headers to the clustered spectrum
      oldFileIdx = (oldFileIdx < 0) ? fileIdx : oldFileIdx;
      oldSpecIdx = (oldSpecIdx < 0) ? specIdx : oldSpecIdx;

      // Fill in correct scan number from input MS/MS spectra
      (*m_outputClusters)[i][j].m_scan =
          inputSpectrumInfo[fileIdx][specIdx].scan;
      // Fill in correct file name from input MS/MS spectra
      (*m_outputClusters)[i][j].m_filename =
          inputSpectrumInfo[fileIdx][specIdx].fileName;

    }

    //Removing the spectra that have too few spectra in them
    if ((*m_outputClusters)[i].size() < min_clust_size)
    {
      list<int> peakstoremove;
      for (int p = 0; p < (*m_outputClusteredSpectra)[i].size(); p++)
      {
        peakstoremove.push_back(p);
      }

      (*m_outputClusteredSpectra)[i].removePeaks(peakstoremove);
    }

    // Copy headers from original MS/MS b/c they were modified before passing to MSCluster
    double clustParentMZ = (*m_outputClusteredSpectra)[i].parentMZ;

    (*m_outputClusteredSpectra)[i].copyNP(inputSpectrumInfo[oldFileIdx][oldSpecIdx]);
    if ((*m_outputClusteredSpectra)[i].parentMass < 0)
    {
      (*m_outputClusteredSpectra)[i].setParentMZ(clustParentMZ);
    }

    // Set new scan #s to 1-based indices
    (*m_outputClusteredSpectra)[i].scan = i + 1;

    // Fill in file name with output spectra file
    (*m_outputClusteredSpectra)[i].fileName = outFilename;

    (*m_outputClusteredSpectra)[i].fileIndex = 0;

    (*m_outputClusters)[i].m_scan = (*m_outputClusteredSpectra)[i].scan;
    (*m_outputClusters)[i].m_filename = (*m_outputClusteredSpectra)[i].fileName;
    (*m_outputClusters)[i].m_fileIndex =
        (*m_outputClusteredSpectra)[i].fileIndex;
  }

  int removableSize = m_params.getValueInt("REMOVABLE_BAD_CLUSTER_SIZE", 0);
  if( removableSize ){
    DEBUG_MSG("Removing (repeated) bad clusters <=" << removableSize << " in size.");
    removeRedundantClusters(removableSize);
  }

  DEBUG_TRACE;

  return true;
}
////////////////////////////////////////////////////////////////////////////////

bool ExecMsCluster::removeRedundantClusters(int removableSize){

	static double	isolationWindow = 3;//in mz
	static double	miginalWindow   = isolationWindow + 5;//for compensate for MS-Cluster tolerance
	static double	wideFragmenttol = 2;//in Da
	static double 	minCosine 		= 0.5;
	static int	 	minMatchedPeaks = 5;

	SpecSet filteredSpectra(m_outputClusteredSpectra->size());
	for (unsigned int i = 0; i < m_outputClusteredSpectra->size(); i++){
		filteredSpectra[i].copyNP((*m_outputClusteredSpectra)[i]);

		float left = filteredSpectra[i].parentMZ-50, right = filteredSpectra[i].parentMZ+10;
		for(unsigned int k=0; k<(*m_outputClusteredSpectra)[i].size(); k++){
			if( left < (*m_outputClusteredSpectra)[i][k][0] && (*m_outputClusteredSpectra)[i][k][0] < right ) continue;
			filteredSpectra[i].push_back((*m_outputClusteredSpectra)[i][k]);
		}
		filteredSpectra[i].rankFilterPeaks(2);
	}

	vector<set<int> > pairList(filteredSpectra.size());
	for (unsigned int i = 0; i < filteredSpectra.size(); i++){
		for (unsigned int j = i+1; j < filteredSpectra.size(); j++) {

			if( filteredSpectra[j].parentMZ - filteredSpectra[i].parentMZ > miginalWindow ) break;
			if( filteredSpectra[i].parentCharge != filteredSpectra[j].parentCharge ||
					fabs(filteredSpectra[i].parentMZ - filteredSpectra[j].parentMZ) > isolationWindow ) continue;

			unsigned int numMatchedPeaks = 0;
			float score1 = 0, score2 = 0;
			float exactCosine = filteredSpectra[i].scoreMatch( filteredSpectra[j], wideFragmenttol, numMatchedPeaks, score1, score2, true, false, true );

			if( exactCosine < minCosine || numMatchedPeaks < minMatchedPeaks ) continue;
			pairList[i].insert(j);
			pairList[j].insert(i);
		}
	}

	vector<int> toBeRemoved;
	for (unsigned int sgt = 0; sgt < filteredSpectra.size(); sgt++){

		int mySize = (*m_outputClusters)[sgt].size();
		if( mySize > removableSize || pairList[sgt].size() == 0 ) continue; //skip large cluster

		vector<int> pairedClusters;
		for (unsigned int j = 0; j < filteredSpectra.size(); j++) {

			if( filteredSpectra[sgt].parentMZ - filteredSpectra[j].parentMZ > miginalWindow ) continue;
			if( filteredSpectra[j].parentMZ - filteredSpectra[sgt].parentMZ > miginalWindow ) break;

			if( pairList[j].size() == 0 || sgt == j ) continue;

			if( pairList[sgt].find(j) != pairList[sgt].end() ) {
				pairedClusters.push_back(j);
			}
		}

		bool hasTriangle = false;
		int  hasGoodCluster = 0;
		for(int m=0; m<pairedClusters.size(); m++){
			int pair1 = pairedClusters[m];
			for(int n=m+1; n<pairedClusters.size(); n++){
				int pair2 = pairedClusters[n];
				if( pairList[pair1].find(pair2) != pairList[pair1].end() ){
					hasTriangle = true;
					if( (*m_outputClusters)[pair1].size() > mySize || (*m_outputClusters)[pair2].size() > mySize ) hasGoodCluster++;

				//	if( ((*m_outputClusters)[pair1].size() > mySize && (*m_outputClusters)[pair2].size() >= mySize) ||
				//			((*m_outputClusters)[pair1].size() >= mySize && (*m_outputClusters)[pair2].size() > mySize) ) hasGoodCluster++;
				}
			}
		}
		if( hasGoodCluster != 0 && hasTriangle ){
			toBeRemoved.push_back(sgt);
		}
	}

	for(int i=0; i<toBeRemoved.size(); i++){//remove all peaks in the bad cluster
		int rm = toBeRemoved[i];
		list<int> peakstoremove;
		for (int p = 0; p < (*m_outputClusteredSpectra)[rm].size(); p++) {
			peakstoremove.push_back(p);
		}
		(*m_outputClusteredSpectra)[rm].removePeaks(peakstoremove);
	}

	DEBUG_MSG("Of " << filteredSpectra.size() << " clustered spectra, " << toBeRemoved.size() << " were removed." );

	return true;
}


////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
int ExecMsCluster::callMsCluster(const string &exeDir,
                                 const string &msClusterInputFileList,
                                 const string &tempDir,
                                 const string &outDir,
                                 const string &baseDir)
{
  // TOLERANCE_PM
  float pmTol = 0.0;
  if (m_params.exists("TOLERANCE_PM") != 0)
    pmTol = getFloat(m_params.getValue("TOLERANCE_PM").c_str());

  // MsCluster command
  std::string msClusterCommand;

  // Path
  msClusterCommand = exeDir;

  replace(msClusterCommand.begin(), msClusterCommand.end(), '\\', '/');

  // Executable name
  msClusterCommand = getPath(msClusterCommand, "MsCluster_bin", false);

  // --list
  msClusterCommand += " --list ";
  msClusterCommand += msClusterInputFileList;

  // set tmp and out dirs
  msClusterCommand += " --tmp-dir ";  //./spectra/tmp";
  msClusterCommand += tempDir;
  msClusterCommand += " --out-dir ";  //./spectra/out";
  msClusterCommand += outDir;

  // --output-name
  msClusterCommand += " --output-name clusters";
  // --window
  msClusterCommand += " --window ";
  msClusterCommand += parseFloat(pmTol, 7);

  // --fragment-tolerance

  // --model-dir
  msClusterCommand += " --model-dir ";
  msClusterCommand += exeDir;
  msClusterCommand += "/resources/Models_mscluster";

  // --fragment-tolerance
  msClusterCommand += " --fragment-tolerance ";
  msClusterCommand += m_params.getValue("TOLERANCE_PEAK");

  // CLUSTER_MODEL
  if (m_params.exists("CLUSTER_MODEL") != 0)
  {
    msClusterCommand += " --model ";
    msClusterCommand += m_params.getValue("CLUSTER_MODEL");

    if (m_params.getValueBool("GUESS_CHARGE", false))
    {
      msClusterCommand += " --assign-charges ";
    }

    if (m_params.getValueBool("CORRECT_PM", false))
    {
      msClusterCommand += " --correct-pm ";
    }
  }

  // MIN_SPECTRUM_QUALITY
  if (m_params.exists("MIN_SPECTRUM_QUALITY") != 0)
  {
    msClusterCommand += " --sqs ";
    msClusterCommand += m_params.getValue("MIN_SPECTRUM_QUALITY");
  }

  // CLUSTER_PMTOL_PPM
  if (m_params.exists("CLUSTER_PMTOL_PPM") != 0)
  {
    msClusterCommand += " --precursor-ppm ";
    msClusterCommand += m_params.getValue("CLUSTER_PMTOL_PPM");
  }

  if (m_params.exists("MSCLUSTER_MIX_PROB") != 0)
  {
    msClusterCommand += " --mixture-prob ";
    msClusterCommand += m_params.getValue("MSCLUSTER_MIX_PROB");
  }

  msClusterCommand += " --memory-gb 0.35 ";

  // Force MsCluster_bin output into a file
  msClusterCommand += " > ";
  msClusterCommand += getPath(baseDir, "log_mscluster.txt", false);

  //Send to sdtout DEBUG
  DEBUG_MSG("call msCluster: " << msClusterCommand);

  //string libPath = exeDir + "/libs";
  //addEnvironmentVariable(msClusterCommand, "LD_LIBRARY_PATH", libPath);

  return spsSystem(msClusterCommand.c_str());
}

bool ExecMsCluster::clusterFilenameSort(const string &a, const string &b)
{
  if (a.length() == b.length())
    return a.compare(b) < 0;
  return (a.length() < b.length());
}

////////////////////////////////////////////////////////////////////////////////
