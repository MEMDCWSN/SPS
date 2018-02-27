/*
 * Cluster.h
 *
 *  Created on: Jul 15, 2013
 *      Author: aguthals
 */
///////////////////////////////////////////////////////////////////////////////
#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <stdio.h>

#include "SpecSet.h"
#include "OutputTable.h"

///////////////////////////////////////////////////////////////////////////////
using namespace std;

///////////////////////////////////////////////////////////////////////////////
namespace specnets
{

  ///////////////////////////////////////////////////////////////////////////////
  //
  ///////////////////////////////////////////////////////////////////////////////
  class SpectrumItem
  {
  public:

    int m_scan; // Scan number of spectrum
    int m_index; // Index number of spectrum
    int m_fileIndex; // Index of the source file in the input spectra file list
    string m_filename; // Source filename of spectrum

    /**
     * Default constructor, sets invalid scan, index, and filename
     */
    SpectrumItem() :
        m_scan(-1), m_index(-1), m_filename(""), m_fileIndex(-1)
    {
    }
    ;

    /**
     * Allocates with an input scan, index, and filename
     * @param scan
     * @param index
     * @param fn
     */
    SpectrumItem(int index, int scan, const int fi, const string &fn) :
        m_scan(scan), m_index(index), m_fileIndex(fi), m_filename(fn)
    {
    }

    /**
     * Re-initializes by setting the scan, indexm and filename
     * @param clusterIdx
     * @param scan
     * @param filename
     */
    inline void initialize(int clusterIdx = -1, int scan = -1, int fileIndex =
                               -1,
                           string filename = "")
    {
      m_scan = scan;
      m_index = clusterIdx;
      m_filename = filename;
      m_fileIndex = fileIndex;
    }
  };

  ///////////////////////////////////////////////////////////////////////////////
  // SpectrumItem data fields of Cluster reference the clustered spectrum
  ///////////////////////////////////////////////////////////////////////////////
  class Cluster : public SpectrumItem
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

    // NO NEED FOR "m_clusterIdx", THIS IS ALREADY STORED IN SpectrumItem::m_index
    // cluster index
    // int m_clusterIdx;

    Cluster() :
        m_clusteredSpectra(0), SpectrumItem()
    {
    }

    /**
     * Destructor
     */
    virtual ~Cluster()
    {
    }

    /**
     * Resets all indices to -1 and clears the set of spectra
     * @param clusterIdx new index of this cluster
     * @param clusterScan new scan of this cluster
     * @param clusterFilename filename of this cluster
     * @param numSpecs number of spectra to make room for
     * @return
     */
    inline void initialize(int clusterIdx = -1,
                           int clusterScan = -1,
                           int fileIndex = -1,
                           string clusterFilename = "",
                           unsigned int numSpecs = 0)
    {
      SpectrumItem::initialize(clusterIdx,
                               clusterScan,
                               fileIndex,
                               clusterFilename);
      m_clusteredSpectra.resize(0);
      if (numSpecs > 0)
      {
        m_clusteredSpectra.resize(numSpecs);
      }
    }

    inline Cluster& operator=(const Cluster& other)
    {
      if (this == &other)
        return *this;

      SpectrumItem::operator =((const SpectrumItem&)other);
      m_clusteredSpectra = other.m_clusteredSpectra;
      return *this;
    }

    /**
     * Adds a spectrum to this cluster
     * @param other
     * @return
     */
    inline virtual void push_back(const SpectrumItem & other)
    {
      m_clusteredSpectra.push_back(other);
    }

    /**
     * Adds a spectrum to this cluster
     * @param scan
     * @param index
     * @param fn
     * @return
     */
    inline virtual void add(int specIndex,
                            int scan,
                            const int fileIndex,
                            const string &fn = "")
    {
      SpectrumItem si(specIndex, scan, fileIndex, fn);
      m_clusteredSpectra.push_back(si);
    }

    /**
     * Adds a spectrum to this cluster while leaving its scan #s and filenames un-set
     * @param index
     * @return
     */
    virtual void add(int specIndex);

    /**
     * Returns the number of spectra in this cluster
     */
    inline unsigned int size() const
    {
      return m_clusteredSpectra.size();
    }

    /**
     * Resizes the set of clustered spectra
     */
    inline void resize(unsigned int numSpecs)
    {
      m_clusteredSpectra.resize(numSpecs);
    }

    /**
     * Index operator for accessing clustered spectra
     */
    inline SpectrumItem &operator[](unsigned int i)
    {
      return m_clusteredSpectra[i];
    }

    /**
     * Const index operator for accessing clustered spectra
     */
    inline const SpectrumItem &operator[](unsigned int i) const
    {
      return m_clusteredSpectra[i];
    }

    /**
     * Saves all class data fields to an open file in the latest binary format
     * @param fp
     * @return true if all fields were successfully saved, false if not
     */
    virtual bool saveToBinaryStream(FILE* fp) const;

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
     * Returns the total number of spectra
     * @return The number of spectra. -1 if there was an error.
     */
    virtual int getSpectraCount(void)
    {
      return m_clusteredSpectra.size();
    }
    ;

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
     * Dumps the cluster info
     * @param sout The output stream
     * @param web If true, adds HTML tags
     */
    virtual void dump(ostream &sout, bool web);

    /**
     * Gets a summary information item
     * @param localGroups Group information
     * @param specs specsets of input spectra
     * @param specs_ms specset of clustered spectra
     * @param groups Groups names
     */
    virtual int getSummary(std::vector<int> &localGroups,
                           vector<SpecSet> &specs,
                           SpecSet &specs_ms,
                           std::map<string, std::set<string> > &groups,
                           int &numSpectra,
                           double &sumPercursorIntensity,
                           vector<string> *inputFiles);

    /**
     * Returns the cluster information in a table format
     * @param data the 2D matrix containing the data
     * @param specs The actual input spectra
     * @return 1 if there was no error.
     */
    virtual int getCsvData(OutputTable &data,
                           vector<SpecSet> &specs,
                           vector<string> *inputFiles,
                           vector<vector<vector<int> > > *scanNumbers);

    /**
     * Writes header for cluster information table
     * @return 1 if there was no error.
     */
    virtual int getCsvDataHeader(OutputTable &data);

    inline int getCluster(void)
    {
      return m_index;
    }
    ;

  protected:

    /**
     * Returns Checks if a filename belongs to a group
     * @param fileGroup The filegroups
     * @param filename The filename
     * @return The goup index, -1 if not found
     */
    int checkFileInGroup(std::set<string> &fileGroup, string &filename);

    // All spectra that are in this cluster
    vector<SpectrumItem> m_clusteredSpectra;

  };
///////////////////////////////////////////////////////////////////////////////
}
///////////////////////////////////////////////////////////////////////////////
#endif /* CLUSTER_H_ */
///////////////////////////////////////////////////////////////////////////////
