/*
 * ClusterSet.h
 *
 *  Created on: Jul 15, 2013
 *      Author: aguthals
 */

#ifndef CLUSTERSET_H_
#define CLUSTERSET_H_

#include <vector>
#include <map>

#include "Cluster.h"
#include "SpecSet.h"
#include "utils.h"

#define __USE_CLUSTER_SET__

using namespace std;

namespace specnets
{
  class ClusterSet
  {
  public:

    // binary format version
    static const unsigned short BIN_VERSION;

    // binary format sub-version
    static const unsigned short BIN_SUBVERSION;

    // string identifier for version
    static const string BIN_VERSION_ID;

    // string identifier for sub-version
    static const string BIN_SUBVERSION_ID;

    /**
     * Default constructor, constructs an empty set of clusters
     */
    ClusterSet() :
        m_clusters(0)
    {
    }

    // file names of input MsCluster files
    vector<string> fileNames;

    /**
     * Allocates room for a number of clusters
     * @param numClusters number of clusters to make room for
     */
    ClusterSet(unsigned int numClusters) :
        m_clusters(numClusters)
    {
    }

    /**
     * Clears out the set of clusters and makes room for a new set
     * @param numClusters number of clusters to make room for
     * @return
     */
    inline void initialize(unsigned int numClusters = 0)
    {
      m_clusters.resize(0);
      if (numClusters > 0)
      {
        m_clusters.resize(numClusters);
      }
    }

    /**
     * Returns the number of clusters in this set
     * @return number of clusters
     */
    inline unsigned int size() const
    {
      return m_clusters.size();
    }

    /**
     * Resizes the set of clusters
     */
    inline void resize(unsigned int numSpecs)
    {
      m_clusters.resize(numSpecs);
    }

    /**
     * Adds a single cluster to the end of the set
     * @param newCluster
     * @return new size of the set
     */
    inline unsigned int push_back(const Cluster& newCluster)
    {
      m_clusters.push_back(newCluster);
      return m_clusters.size();
    }

    /**
     * Index operator for accessing clusters
     */
    inline Cluster &operator[](unsigned int i)
    {
      return m_clusters[i];
    }

    /**
     * Const index operator for accessing clusters
     */
    inline const Cluster &operator[](unsigned int i) const
    {
      return m_clusters[i];
    }

    inline ClusterSet& operator=(const ClusterSet& other)
    {
      if (this == &other)
        return *this;

      m_clusters = other.m_clusters;
      return (*this);
    }

    /**
     * Saves all class data fields to an open file in the latest binary format
     * @param fp
     * @return true if all fields were successfully saved, false if not
     */
    virtual bool saveToBinaryStream(FILE* fp) const;

    /**
     * Saves the class structure and all associated data in binary format
     * @param filename
     * @return true if file was saved successfully, false if not.
     */
    virtual bool saveBinaryFile(const string& filename) const;

    /**
     * Loads all class data fields from an open file in the specified binary format
     * @param fp
     * @param versions encodes BIN_VERSION and BIN_SUBVERSION, which should be referenced
     *                 by BIN_VERSION_ID and BIN_SUBVERSION_ID, respectively.
     * @return true if all fields were successfully loaded, false if not
     */
    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    /**
     * Loads the class structure and all associated data in binary format
     * @param filename
     * @return true if file was loaded successfully, false if not.
     */
    virtual bool loadBinaryFile(const string& filename);

    /**
     * Returns the total number of clusters
     * @return The number of clusters. -1 if there was an error.
     */
    virtual int getClusterCount(void)
    {
      return m_clusters.size();
    }
    ;

    /**
     * Returns the total number of spectra
     * @return The number of spectra. -1 if there was an error.
     */
    virtual int getSpectraCount(void);

    /**
     * Returns the cluster index from the input spectra index plus file index
     * @param fileIndex The index of the file containing the spectra
     * @param inputSpectra The spectra index
     * @return The cluster index. -1 if it was not found.
     */
    virtual int getConsensusFromInputSpectra(int fileIndex, int inputSpectra);

    /**
     * Returns the spectra assotiated with the cluster.
     * @param The cluster index
     * @return A list containing <file index ; spectrum index> pairs, which identify a spectra
     */
    virtual list<pair<unsigned, unsigned> >
    *getInputSpectraFromConsensus(int consensus);

    /**
     * Dumps the clusters info
     * @param sout The output stream
     * @param web If true, adds HTML tags
     */
    virtual void dump(ostream &sout, bool web);

    /**
     * Returns the cluster information in a table format
     * @param inputFiles The input spectra file names
     * @param scanNumbers A vector of vectors containing the scan numbers.
     * @param specs The actual input spectra
     * @param data the 2D matrix containing the data
     * @return 1 if there was no error.
     */
    virtual int getCsvData(vector<string> *inputFiles,
                           vector<vector<vector<int> > > *scanNumbers,
                           vector<SpecSet> &specs,
                           OutputTable &data,
                           const bool writeHeader = true);

    /**
     * Returns the clusters summary in a 2D string matrix for ClusterInfo
     * @param result the 2D matrix containing the summary data
     * @param specs The input spectra, grouped by specset
     * @param specs_ms The clusters spectra
     * @param inputFiles The input spectra file names
     * @param groups User defined group information.
     * @return 1 if there was no error.
     */
    virtual int getSummary(OutputTable &result,
                           vector<SpecSet> &specs,
                           SpecSet &specs_ms,
                           vector<string> *inputFiles,
                           std::map<string, string> &groups);

    void getScanToFileMapping(map<int, list<pair<int, string> > >& clusterInfo);

  protected:

    /**
     * Returns Adds a record item to the summary structure (clusterInfo)
     * @param result the 2D matrix containing the summary data
     * @param specs_ms The clustered spectra
     * @param cluster The clusters index
     * @param numSpectra The spectra count
     * @param sumPercursorIntensity Sum of percursor intensity.
     * @param groupCount The number of groups
     */
    void addRecordToSummary(OutputTable &result,
                            SpecSet &specs_ms,
                            int cluster,
                            int numSpectra,
                            float sumPercursorIntensity,
                            vector<int> &groupCount);

    Cluster *findCluster(const int clusterIdx);

    // Set of clustered spectra
    vector<Cluster> m_clusters;
  };
}

#endif /* CLUSTERSET_H_ */
