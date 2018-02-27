/*
 * ClusterSet.cpp
 *
 *  Created on: Jul 15, 2013
 *      Author: aguthals
 */

#include "ClusterSet.h"
#include "Logger.h"
#include "utils.h"

using namespace std;

namespace specnets
{
  const unsigned short ClusterSet::BIN_VERSION = 1;

  const unsigned short ClusterSet::BIN_SUBVERSION = 1;

  const string ClusterSet::BIN_VERSION_ID = "ClusterSet_binVersion";

  const string ClusterSet::BIN_SUBVERSION_ID = "ClusterSet_binSubVersion";

  bool ClusterSet::saveToBinaryStream(FILE *fp) const
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned int numClusters = size();
    unsigned int count = fwrite(&numClusters, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save # of clusters");
      return false;
    }

    for (unsigned int i = 0; i < numClusters; i++)
    {
      if (!m_clusters[i].saveToBinaryStream(fp))
      {
        return false;
      }
    }

    return true;
  }

  bool ClusterSet::saveBinaryFile(const string& filename) const
  {
    FILE *fp = fopen(filename.c_str(), "wb");
    if (fp == 0)
    {
      ERROR_MSG("Cannot open " << filename);
      return false;
    }

    map<string, unsigned short> versions;
    versions[Cluster::BIN_VERSION_ID] = Cluster::BIN_VERSION;
    versions[Cluster::BIN_SUBVERSION_ID] = Cluster::BIN_SUBVERSION;
    versions[ClusterSet::BIN_VERSION_ID] = ClusterSet::BIN_VERSION;
    versions[ClusterSet::BIN_SUBVERSION_ID] = ClusterSet::BIN_SUBVERSION;

    if (!writeStringMapToBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error saving version info");
      fclose(fp);
      return false;
    }

    if (!saveToBinaryStream(fp))
    {
      ERROR_MSG("Error saving in binary format");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  bool ClusterSet::loadFromBinaryStream(FILE* fp,
                                        map<string, unsigned short>& versions)
  {
    if (fp == 0)
    {
      return false;
    }
    initialize();

    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    unsigned int numClusters;

    unsigned int count = fread(&numClusters, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of clusters");
      return false;
    }

    DEBUG_VAR(numClusters);

    m_clusters.resize(numClusters);

    for (unsigned int i = 0; i < numClusters; i++)
    {
      if (!m_clusters[i].loadFromBinaryStream(fp, versions))
      {
        return false;
      }
    }

    return true;
  }

  bool ClusterSet::loadBinaryFile(const string& filename)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    if (fp == 0)
    {
      ERROR_MSG("Opening " << filename);
      return false;
    }

    map<string, unsigned short> versions;
    if (!readStringMapFromBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error reading version info");
      fclose(fp);
      return false;
    }

    if (!loadFromBinaryStream(fp, versions))
    {
      ERROR_MSG("Error loading in binary format");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Gets # of spectra
  ////////////////////////////////////////////////////////////////////////////////
  int ClusterSet::getSpectraCount(void)
  {
    int ret = 0;
    for (int i = 0; i < m_clusters.size(); i++)
      ret += m_clusters[i].getSpectraCount();
    return ret;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Get consensus spectra index from input files
  ////////////////////////////////////////////////////////////////////////////////
  int ClusterSet::getConsensusFromInputSpectra(int fileIndex, int inputSpectra)
  {
    int ret = -1;
    for (int i = 0; i < m_clusters.size(); i++)
      if ((ret = m_clusters[i].getConsensusFromInputSpectra(fileIndex,
                                                            inputSpectra))
          != -1)
        return ret;
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Gets input spectra list given a cluster index
  ////////////////////////////////////////////////////////////////////////////////
  list<pair<unsigned, unsigned> > *ClusterSet::getInputSpectraFromConsensus(int consensus)
  {
    //DEBUG_VAR(m_clusters.size());

    list<pair<unsigned, unsigned> > *ret = NULL;
    for (int i = 0; i < m_clusters.size(); i++)
      if (ret = m_clusters[i].getInputSpectraFromConsensus(consensus))
        return ret;
    return NULL;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Dump cluster information
  ////////////////////////////////////////////////////////////////////////////////
  void ClusterSet::dump(ostream &sout, bool web)
  {
    for (int i = 0; i < m_clusters.size(); i++)
      m_clusters[i].dump(sout, web);
    sout << "--------------------------------";
  }

  /**
   * Returns the cluster information in a table format
   * @param inputFiles The input spectra file names
   * @param scanNumbers A vector of vectors containing the scan numbers.
   * @param specs The actual input spectra
   * @param data the 2D matrix containing the data
   * @return 1 if there was no error.
   */
  ////////////////////////////////////////////////////////////////////////////////
  // Get cluster information in a tabular format (to be used by clusterInfo)
  ////////////////////////////////////////////////////////////////////////////////
  int ClusterSet::getCsvData(vector<string> *inputFiles,
                             vector<vector<vector<int> > > *scanNumbers,
                             vector<SpecSet> &specs,
                             OutputTable &data,
                             const bool writeHeader)
  {
    if (writeHeader && m_clusters.size() > 0)
    {
      m_clusters[0].getCsvDataHeader(data);
    }
    for (int i = 0; i < m_clusters.size(); i++){
      if(i%1000 == 0)
        cout<<"Getting Data for cluster " <<i<<" of "<<m_clusters.size()<<endl;
      m_clusters[i].getCsvData(data, specs, inputFiles, scanNumbers);
    }

    return 1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Helper method that adds a record to the summary list (to be used by clusterInfo)
  ////////////////////////////////////////////////////////////////////////////////
  void ClusterSet::addRecordToSummary(OutputTable &result,
                                      SpecSet &specs_ms,
                                      int cluster,
                                      int numSpectra,
                                      float sumPercursorIntensity,
                                      vector<int> &groupCount)
  {
    float parentMass = 0.0;
    float parentCharge = 0.0;
    float m = 0.0;

    if (specs_ms.size() > cluster)
    {
      parentMass = specs_ms[cluster].parentMass;
      parentCharge = specs_ms[cluster].parentCharge;
      if (parentCharge)
        m = (parentCharge + parentMass - 1.0) / parentCharge;
      else
        m = parentMass;
    }

    // store the data
    vector<pair<string, bool> > row;
    stringstream aux;

    aux << (cluster + 1);
    row.push_back(pair<string, bool>(aux.str(), false));

    aux.str(std::string());
    aux << numSpectra;
    row.push_back(pair<string, bool>(aux.str(), false));

    aux.str(std::string());
    aux << parentMass;
    row.push_back(pair<string, bool>(aux.str(), false));

    aux.str(std::string());
    aux << parentCharge;
    row.push_back(pair<string, bool>(aux.str(), false));

    aux.str(std::string());
    aux << m;
    row.push_back(pair<string, bool>(aux.str(), false));

    aux.str(std::string());
    aux << sumPercursorIntensity;
    row.push_back(pair<string, bool>(aux.str(), false));

    for (int i = 0; i < groupCount.size(); i++)
    {
      aux.str(std::string());
      aux << groupCount[i];
      row.push_back(pair<string, bool>(aux.str(), false));
    }

    result.values.push_back(row);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Get summary data in a 2D array of strings (to be used with ClusterInfo)
  ////////////////////////////////////////////////////////////////////////////////

  int ClusterSet::getSummary(OutputTable &result,
                             vector<SpecSet> &specs,
                             SpecSet &specs_ms,
                             vector<string> *inputFiles,
                             std::map<string, string> &groups)
  {
    //Creating groups to files mapping
    std::map<string, std::set<string> > group_to_files_mapping;
    for (map<std::string, std::string>::iterator it = groups.begin();
          it != groups.end(); it++)
    {
        vector<string> res;
        string delim(";");
        stringSplit2(it->second, res, delim);
        
        std::set<string> file_list_set;
        
        for(int i = 0; i < res.size(); i++){
            file_list_set.insert(res[i]);
        }
        
        group_to_files_mapping[it->first] = file_list_set;
    }
    
    // Cycle thru records
    for (int i = 0; i < m_clusters.size(); i++)
    {

      int numSpectra = 0;
      double sumPercursorIntensity = 0.0;

      int cluster = m_clusters[i].getCluster();

      // local group holder
      std::vector<int> localGroups;
      localGroups.resize(groups.size());
      for (int j = 0; j < groups.size(); j++)
        localGroups[j] = 0;

      m_clusters[i].getSummary(localGroups,
                               specs,
                               specs_ms,
                               group_to_files_mapping,
                               numSpectra,
                               sumPercursorIntensity,
                               inputFiles);

      addRecordToSummary(result,
                         specs_ms,
                         cluster,
                         numSpectra,
                         sumPercursorIntensity,
                         localGroups);
    }

    // exit
    return 1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Find a cluster object, given it's index
  ////////////////////////////////////////////////////////////////////////////////

  Cluster *ClusterSet::findCluster(const int clusterIdx)
  {
    for (int i = 0; i < m_clusters.size(); i++)
      if (clusterIdx == m_clusters[i].getCluster())
        return &m_clusters[i];
    return NULL;
  }

  void ClusterSet::getScanToFileMapping(map<int, list<pair<int, string> > >& clusterInfo)
  {
    clusterInfo.clear();
    for (unsigned int i = 0; i < size(); i++)
    {
      list<pair<int, string> > clust;

      for (unsigned int j = 0; j < m_clusters[i].size(); j++)
      {
        clust.push_back(pair<int, string>(m_clusters[i][j].m_scan,
                                          m_clusters[i][j].m_filename));
      }
      clusterInfo[m_clusters[i].m_scan] = clust;
    }
  }

////////////////////////////////////////////////////////////////////////////////

}
