////////////////////////////////////////////////////////////////////////////////
#include <dirent.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "ClusterData.h"
#include "ReportDefines.h"
#include "utils.h"

////////////////////////////////////////////////////////////////////////////////
// defines

#define CSV_SEP ';'

////////////////////////////////////////////////////////////////////////////////
// namespaces

using namespace std;
namespace spsReports {

////////////////////////////////////////////////////////////////////////////////
// Sort directory contents
////////////////////////////////////////////////////////////////////////////////
bool ClusterData::clusterFilenameSort(const string &a, const string &b)
{
  if(a.length() == b.length())
    return a.compare(b) < 0;
  return (a.length() < b.length());
}
////////////////////////////////////////////////////////////////////////////////
// Get directory contents
////////////////////////////////////////////////////////////////////////////////
int ClusterData::getClusterFileNames(const string & directory, const string &prefix, bool doSort)
{
  // first, create a pointer to the directory
  DIR *pdir = NULL;
  // Open directory
  pdir = opendir (directory.c_str());

  // Directory listing structure
  struct dirent *pent = NULL;
  // Test for invalid directory
  if (pdir == NULL)  {
    // print an error message and exit the program
    ERROR_MSG("Invalid directory for MsCluster output data: " << directory);
    return -1;
  }

  // Cycle thru all elements in the directory
  while(pent = readdir(pdir)) {
    // Check for errors
    if (pent == NULL) {
      ERROR_MSG("Reading MsCluster output data directory: unexpected error.");
      return -1;
    }

    // if the file matches the prefix, lets added to the clust file list
    string aux = pent->d_name;
    if(aux.compare(0, prefix.length(), prefix) == 0)
      fileNames.push_back(aux);
  }

  // Close the directory
  closedir (pdir);

  // Sort the file list
  if(doSort)
    sort(fileNames.begin(), fileNames.end(), ClusterData::clusterFilenameSort);

  // exit with no error
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Load data form MsCluster output files.
////////////////////////////////////////////////////////////////////////////////
int ClusterData::loadMsClusterData(const string &projectDir)
{
  // Compose MsCluster output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/' && aux[aux.length()-1] != '\\')
    aux += '/';
  aux += MS_CLUSTER_DATA_SUBPATH;
  // Get msCluster output file names. If error, exit
  if(getClusterFileNames(aux, MS_CLUSTER_FILE_PREFIX, true) == -1)
    return -1;

  // Spectrum counter
  unsigned nspectrum = 0;

  // Cycle thru the found MsCLuster files
  for(unsigned i = 0; i < fileNames.size(); ++ i) {
    // File reader object
    ifstream currentFile;
    // compose file name to include project directory
    string aux2 = aux;
    aux2 += '/';
    aux2 += fileNames[i];
    // open the file
    currentFile.open(aux2.c_str(), ifstream::in | ifstream::binary);
    // test for errors
    if(!currentFile) {
      //cerr << "ERROR: cannot open '" << fileNames[i] << endl;
      return -1;
    }

    //
    for (unsigned j = 0; currentFile.peek() != EOF; ++ j)  {
      string sline, sname;
      char cplot;
      unsigned nclusteridx, nspectrum3;

      getline(currentFile, sname, '.');
      currentFile >> nclusteridx >> cplot >> nspectrum3;
      getline(currentFile, sname);

      ClusterData_t::iterator ic = data.insert(data.end(), make_pair(nspectrum, list< pair<unsigned, unsigned> >()));

      unsigned k = 0;

      for (; currentFile && currentFile.peek() != '\n'; ++ k)
      {
        unsigned nclusteridx2, nspectrum2, aux;

        currentFile >> aux >> nclusteridx2 >> nspectrum2;
        getline(currentFile, sline);

        ic->second.push_back(make_pair(nclusteridx2, nspectrum2));
      }
      if ( k == 0 ) {
        data.erase(ic);
      } //else
      nspectrum++;
    }
    nspectrum--;
  }

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Load data stored in binary format
////////////////////////////////////////////////////////////////////////////////
int ClusterData::loadData(const string &projectDir)
{
  // make sure it's empry
  data.clear();

  // compose file name to include project directory
  string aux = projectDir;
  if(aux[aux.length()-1] != '/' && aux[aux.length()-1] != '\\')
    aux += '/';
  aux += CLUSTER_DATA_FILENAME;

  // open file to read from memory
  ifstream file(aux.c_str(), ios::in | ios::binary | ios::ate);
    // if error, say so
  if(!file.is_open()) {
    ERROR_MSG("Cluster data: could not open file " << aux);
    return -1;
  }

  // calculate file size
  file.seekg (0, ios::end);
  unsigned size = file.tellg();
  //cout << "size (load single): " << size << endl;
  // get memory block to load file into
  unsigned *memblock = new unsigned int[size * sizeof(char) / sizeof(unsigned int)];
  // Move to the beggining of file
  file.seekg (0, ios::beg);
  // read it
  file.read((char *)memblock, size);
  // close the file
  file.close();

  if(size < sizeof(unsigned int)) {
    ERROR_MSG("Reading file. Only got " << size << " bytes.");
    return -1;
  }

  // file loading into the structure we desire, from the memory buffer
  int filePointer = 0;
  // get the # of records
  unsigned nRecors = memblock[filePointer++];
  //cout << "nRecors (load single): " << nRecors << endl;
  // parse each record
  for(unsigned i = 0 ; i < nRecors ; i++) {
    // read record number
    unsigned numRecord = memblock[filePointer++];
    // create the record entry
    ClusterData_t::iterator ic = data.insert(data.end(), make_pair(numRecord, list< pair<unsigned, unsigned> >()));
    // read the number of entries for the record
    unsigned numEntries = memblock[filePointer++];
    // pease each entrie
    for(unsigned j = 0 ; j < numEntries ; j++) {
      // read file index
      unsigned fileIdx = memblock[filePointer++];
      // read spectrum index
      unsigned specIdx = memblock[filePointer++];
      // add entry to record
      ic->second.push_back(make_pair(fileIdx, specIdx));
    }
  }

  // delete buffer
  delete[] memblock;

  //return OK
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// Load data stored in binary format, load it to double the number of spectra.
// This is only used by merging routine of CSPS and GenoMS
// Whereas there were N consensus spectra prevously, there will now be 2N spectra.
// the 0th and the Nth spectra are identical, they just have different indices
////////////////////////////////////////////////////////////////////////////////
int ClusterData::loadDataDouble(const string &projectDir)
{
  // make sure it's empry
  data.clear();

  //cout << "cleared data" << endl;
  // compose file name to include project directory
  string aux = projectDir;
  if(aux[aux.length()-1] != '/' && aux[aux.length()-1] != '\\')
    aux += '/';
  aux += CLUSTER_DATA_FILENAME;
  //cout << "cluster file name: " << aux << endl;
  // open file to read from memory
  ifstream file(aux.c_str(), ios::in | ios::binary | ios::ate);
    // if error, say so
  if(!file.is_open()) {
    ERROR_MSG("Cluster data: could not open file " << aux);
    return -1;
  }

  // calculate file size
  file.seekg (0, ios::end);
  unsigned size = file.tellg();
  //cout << "file size: " << size << endl;
  // get memory block to load file into
  unsigned *memblock = new unsigned int[(size) * sizeof(char) / sizeof(unsigned int)];
  //cout << "allcoated memblock" << endl;

  // Move to the beggining of file
  file.seekg (0, ios::beg);
  // read it
  file.read((char *)memblock, size);
  //cout << "read the file 1 time" << endl;
  // close the file
  file.close();

  if(size < sizeof(unsigned int))
    return -1;

  // file loading into the structure we desire, from the memory buffer
  int filePointer = 0;
  // get the # of records
  unsigned nRecors = memblock[filePointer++];
  //cout << "nRecors: " << nRecors << endl;
  // parse each record
  for(unsigned i = 0 ; i < nRecors ; i++) {
    // read record number

    unsigned numRecord = memblock[filePointer++];
    // create the record entry
    ClusterData_t::iterator ic = data.insert(data.end(), make_pair(numRecord, list< pair<unsigned, unsigned> >()));
    // read the number of entries for the record
    unsigned numEntries = memblock[filePointer++];
    // pease each entrie
    for(unsigned j = 0 ; j < numEntries ; j++) {
      // read file index
      unsigned fileIdx = memblock[filePointer++];
      // read spectrum index
      unsigned specIdx = memblock[filePointer++];
      // add entry to record
      ic->second.push_back(make_pair(fileIdx, specIdx));
    }
  }

  //cout << "finished loading first set" << endl;
  // file loading into the structure we desire, from the memory buffer
  filePointer = 0;
  // get the # of records
  nRecors = memblock[filePointer++];

  // parse each record
  for(unsigned i = 0 ; i < nRecors ; i++) {
    // read record number
    unsigned numRecord = memblock[filePointer++] + nRecors;
    // create the record entry
    ClusterData_t::iterator ic = data.insert(data.end(), make_pair(numRecord, list< pair<unsigned, unsigned> >()));
    // read the number of entries for the record
    unsigned numEntries = memblock[filePointer++];
    // pease each entrie
    for(unsigned j = 0 ; j < numEntries ; j++) {
      // read file index
      unsigned fileIdx = memblock[filePointer++];
      // read spectrum index
      unsigned specIdx = memblock[filePointer++];
      // add entry to record
      ic->second.push_back(make_pair(fileIdx, specIdx));
      //cout << "Adding spec " << j << " of " << numEntries << endl;
      //fflush(stdout);
    }
  }

  // delete buffer
  delete[] memblock;

  //return OK
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Output file in binary format
////////////////////////////////////////////////////////////////////////////////
int ClusterData::saveData(const string &projectDir)
{
  // compose file name to include project directory
  string auxStr = projectDir;
  if(auxStr[auxStr.length()-1] != '/' && auxStr[auxStr.length()-1] != '\\')
    auxStr += '/';
  auxStr += CLUSTER_DATA_FILENAME;

  // open file to write to
  ofstream file(auxStr.c_str(), ios::out | ios::binary);
    // if error, say so
  if(!file.is_open()) {
    ERROR_MSG("Cluster data: could not open file for writing: " << auxStr);
    return -1;
  }

  // output file

  // output # of records
  unsigned aux = data.size();
  //cout << "data size when written: " << aux << endl;
  file.write((char*)&aux, sizeof(unsigned));
  // Write records
  for( ClusterData_t::iterator ic = data.begin() ; ic != data.end() ; ic++) {
    // Write record #
    aux = ic->first;
    file.write((char*)&aux, sizeof(unsigned));
    // Write # of entries
    aux = ic->second.size();
    file.write((char*)&aux, sizeof(unsigned));
    // Write entries
    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++) {
      // Write file index
      aux = ic2->first;
      file.write((char*)&aux, sizeof(unsigned));
      // Write spectrum index
      aux = ic2->second;
      file.write((char*)&aux, sizeof(unsigned));
    }
  }
  // close the file
  file.close();

  // exit OK
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// General load method
////////////////////////////////////////////////////////////////////////////////
int ClusterData::load(const string &projectDir, bool saveBinary, bool rebuild)
{
  // first try to load data in binary format
  if(rebuild || loadData(projectDir) == -1) {
    // if cant find it, try to read the MsCluster output files
    if(loadMsClusterData(projectDir) == -1) {
      // if fail, there is an error. Possibly the wrong directory
      return ERROR;
    }
    // If successful, save in binary format
    if(saveBinary)
      if(saveData(projectDir) == -1) {
        // if fail to save, there is an error. Possibly a permissions problem
        return ERROR;
      }
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// check if the file is already present
////////////////////////////////////////////////////////////////////////////////
int ClusterData::check(const string &projectDir)
{
  if(loadData(projectDir) == -1)
    return 0;
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
// Mapping methods
////////////////////////////////////////////////////////////////////////////////
int ClusterData::getConsensusFromInputSpectra(int fileIndex, int inputSpectra)
{
  ClusterData_t::iterator ic;
  for(ic = data.begin() ; ic != data.end() ; ic++) {
    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++)
      if((ic2->first == fileIndex) && (ic2->second == inputSpectra))
        return ic->first;
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
list<pair<unsigned,unsigned> > *ClusterData::getInputSpectraFromConsensus(int consensus)
{
  ClusterData_t::iterator ic;
  ic = data.find(consensus);
  if(ic != data.end())
    return &(ic->second);
  return NULL;
}
////////////////////////////////////////////////////////////////////////////////
int ClusterData::getSpectraCount(void)
{
  int ret = 0;
  ClusterData_t::iterator ic;
  for(ic = data.begin() ; ic != data.end() ; ic++)
    ret += ic->second.size();
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
// Dump contents --> debug
////////////////////////////////////////////////////////////////////////////////
void ClusterData::dump(ostream &sout, bool web)
{
  for(int i = 0 ; i < fileNames.size() ; i++) {
    sout << i << " -- " << fileNames[i];
    if(web)
      sout << "<br>";
    sout << endl;
  }
  ClusterData_t::iterator ic;
  for(ic = data.begin() ; ic != data.end() ; ic++) {
    sout << ic->first << " ---------------- ";
    if(web)
      sout << "<br>";
    sout << endl;
    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++)
      sout << " ; " <<ic2->first << " , " << ic2->second;
      if(web)
        sout << "<br>";
      sout << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ClusterData::addRecordToSummary(vector<vector<string> > &result, SpecSet &specs_ms, int cluster, int numSpectra, float sumPercursorIntensity, vector<int> &groupCount)
{
  float parentMass = 0.0;
  float parentCharge = 0.0;
  float m = 0.0;

  if(specs_ms.size() > cluster) {
    parentMass    = specs_ms[cluster].parentMass;
    parentCharge  = specs_ms[cluster].parentCharge;
    if(parentCharge)
      m           = (parentCharge + parentMass - 1.0) / parentCharge;
    else
      m           = parentMass;
  }

  // store the data
  vector<string> row;
  stringstream aux;

  aux << (cluster + 1);
  row.push_back(aux.str());

  aux.str(std::string());
  aux << numSpectra;
  row.push_back(aux.str());

  aux.str(std::string());
  aux << parentMass;
  row.push_back(aux.str());

  aux.str(std::string());
  aux << parentCharge;
  row.push_back(aux.str());

  aux.str(std::string());
  aux << m;
  row.push_back(aux.str());

  aux.str(std::string());
  aux << sumPercursorIntensity;
  row.push_back(aux.str());

  for(int i = 0 ; i < groupCount.size() ; i++) {
    aux.str(std::string());
    aux << groupCount[i];
    row.push_back(aux.str());
  }


  result.push_back(row);
}
////////////////////////////////////////////////////////////////////////////////
int ClusterData::checkFileInGroup(string &fileGroup, string &filename)
{
  vector<string> res;
  string delim(";");
  stringSplit2(fileGroup, res, delim);

  for(int i = 0 ; i < res.size() ; i++) {
    if(res[i].compare(filename) == 0)
      return i;
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
// 1) cluster index, 1 to numClusters
// 2) number of spectra in cluster
// 3) cluster Spectrum::parentMass
// 4) cluster Spectrum::precursorCharge
// 5) cluster (Spectrum::parentMass + Spectrum::precursorCharge -1) / Spectrum::precursorCharge
// 6) summed Spectrum::precursor_intensity for all spectra in the cluster
int ClusterData::getSummary(vector<vector<string> > &result, vector<SpecSet> &specs, SpecSet &specs_ms, vector<string> *inputFiles, std::map<string, string> &groups)
{
  // Cycle thru records
  for( ClusterData_t::iterator ic = data.begin() ; ic != data.end() ; ic++) {

    // Get record
    int cluster = ic->first;

    int numSpectra                =  0;
    double sumPercursorIntensity  =  0.0;

    // local group holder
    std::vector<int> localGroups;
    localGroups.resize(groups.size());
    for(int i = 0 ; i < groups.size() ; i++)
      localGroups[i] = 0;

    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++) {

      // get the file index
      int fileIndex     = ic2->first;
      int spectrumIndex = ic2->second;

      //get filename
      string filename;
       if(inputFiles)
        if(fileIndex < inputFiles->size())
          filename = (*inputFiles)[fileIndex];

      // get group for this spectrum, and add it
      int idx = 0;
      for(map<std::string, std::string>::iterator it = groups.begin(); it != groups.end(); it++) {
        int ret = checkFileInGroup(it->second, filename);
        if(ret >= 0)
          localGroups[idx]++;
        idx++;
      }

      // add data
      numSpectra++;
      sumPercursorIntensity += specs[fileIndex][spectrumIndex].precursor_intensity;
    }

    addRecordToSummary(result, specs_ms, cluster, numSpectra, sumPercursorIntensity, localGroups);
  }

  // exit
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
int ClusterData::getCsvData(vector<string> *inputFiles, vector<vector<vector<int> > > *scanNumbers, vector<SpecSet> &specs, vector<vector<string> > &ret)
{
  // Cycle thru records
  for( ClusterData_t::iterator ic = data.begin() ; ic != data.end() ; ic++) {
    // Get record
    int cluster = ic->first;
    // Write entries
    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++) {

      // get the filename
      int fileIndex = ic2->first;

      string filename;
      filename = "File index: ";
      filename += parseInt(fileIndex);

      if(inputFiles)
        if(fileIndex < inputFiles->size())
          filename = (*inputFiles)[fileIndex];

      // get spectrum index
      int spectrumIndex = ic2->second;
      string spectrumIndex2 = parseInt(spectrumIndex);

      // Get Scan number
      int scan = -1;
      if(scanNumbers) {
        if(scanNumbers)
          if(fileIndex < scanNumbers->size()) {
            if(spectrumIndex < (*scanNumbers)[fileIndex].size())
              scan = (*scanNumbers)[fileIndex][spectrumIndex][0];
          }
      }

      float parentMass = -1.0;
      int parentCharge = -1;
      float retention_time = 0.0;
      float precursor_intensity;

      if(fileIndex < specs.size())
        if(spectrumIndex < specs[fileIndex].size()) {
          parentMass = specs[fileIndex][spectrumIndex].parentMass;
          parentCharge = specs[fileIndex][spectrumIndex].parentCharge;
          retention_time = specs[fileIndex][spectrumIndex].retention_time;
          precursor_intensity = specs[fileIndex][spectrumIndex].precursor_intensity;
        }

      // store the data
      vector<string> row;
      stringstream aux;

      int cc = cluster + 1;
      aux << cc;
      row.push_back(aux.str());

      row.push_back(filename);

      row.push_back(spectrumIndex2);

      aux.str(std::string());
      aux << scan;
      row.push_back(aux.str());

      aux.str(std::string());
      aux << parentMass;
      row.push_back(aux.str());

      aux.str(std::string());
      aux << parentCharge;
      row.push_back(aux.str());

      aux.str(std::string());
      aux << retention_time;
      row.push_back(aux.str());

      aux.str(std::string());
      aux << precursor_intensity;
      row.push_back(aux.str());

      ret.push_back(row);
    }
  }

  // exit OK
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
int ClusterData::writeCsv(vector<string> *inputFiles, vector<vector<vector<int> > > *scanNumbers, string &outFileName, vector<SpecSet> &specs)
{
  // open file to write to
  ofstream file(outFileName.c_str(), ios::out | ios::binary);
    // if error, say so
  if(!file.is_open()) {
    ERROR_MSG("Cluster data: could not open file for writing: " << outFileName);
    return -1;
  }

  // output file

  // Write records
  for( ClusterData_t::iterator ic = data.begin() ; ic != data.end() ; ic++) {
    // Get record
    int cluster = ic->first;
    // Write entries
    list<pair<unsigned,unsigned> >::iterator ic2;
    for(ic2 = ic->second.begin() ; ic2 != ic->second.end() ; ic2++) {

      // get the filename
      int fileIndex = ic2->first;

      string filename;
      filename = "File index: ";
      filename += parseInt(fileIndex);

      if(inputFiles)
        if(fileIndex < inputFiles->size())
          filename = (*inputFiles)[fileIndex];

      // get spectrum index
      int spectrumIndex = ic2->second;
      string spectrumIndex2 = parseInt(spectrumIndex);

      // Get Scan number
      int scan = -1;
      if(scanNumbers) {
        if(scanNumbers)
          if(fileIndex < scanNumbers->size()) {
            if(spectrumIndex < (*scanNumbers)[fileIndex].size())
              scan = (*scanNumbers)[fileIndex][spectrumIndex][0];
          }
      }

      float parentMass = -1.0;
      int parentCharge = -1;
      float retention_time;
      float precursor_intensity;

      if(fileIndex < specs.size())
        if(spectrumIndex < specs[fileIndex].size()) {
          parentMass = specs[fileIndex][spectrumIndex].parentMass;
          parentCharge = specs[fileIndex][spectrumIndex].parentCharge;
          retention_time = specs[fileIndex][spectrumIndex].retention_time;
          precursor_intensity = specs[fileIndex][spectrumIndex].precursor_intensity;
        }




      // Write to file
      file << cluster << CSV_SEP << filename << CSV_SEP << spectrumIndex2 << CSV_SEP << scan<< CSV_SEP << parentMass << CSV_SEP << parentCharge  << CSV_SEP << retention_time << CSV_SEP << precursor_intensity << endl;
    }
  }
  // close the file
  file.close();

  // exit OK
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
};  //namespace
////////////////////////////////////////////////////////////////////////////////
