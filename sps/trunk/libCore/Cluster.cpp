/*
 * Cluster.cpp
 *
 *  Created on: Jul 15, 2013
 *      Author: aguthals
 */

#include "Cluster.h"
#include "Logger.h"
#include "utils.h"

namespace specnets
{

  const unsigned short Cluster::BIN_VERSION = 1;

  const unsigned short Cluster::BIN_SUBVERSION = 1;

  const string Cluster::BIN_VERSION_ID = "Cluster_binVersion";

  const string Cluster::BIN_SUBVERSION_ID = "Cluster_binSubVersion";

  ///////////////////////////////////////////////////////////////////////////////
  void Cluster::add(int specIndex)
  {
    int scan = -1;
    string aux;
    int fileIndex = -1;
    add(specIndex, -1, -1, aux);
  }

  ///////////////////////////////////////////////////////////////////////////////
  int Cluster::getConsensusFromInputSpectra(int fileIndex, int inputSpectra)
  {
    for (int i = 0; i < m_clusteredSpectra.size(); i++)
      if ((m_clusteredSpectra[i].m_fileIndex == fileIndex)
          && (m_clusteredSpectra[i].m_index == inputSpectra))
        return m_index;
    return -1;
  }
  ////////////////////////////////////////////////////////////////////////////////
  list<pair<unsigned, unsigned> > *Cluster::getInputSpectraFromConsensus(int consensus)
  {
    //DEBUG_VAR(m_index);
    //DEBUG_VAR(consensus);

    if (consensus != m_index)
      return NULL;

    list<pair<unsigned, unsigned> > *ret = new list<pair<unsigned, unsigned> >;

    for (int i = 0; i < m_clusteredSpectra.size(); i++)
      ret->push_back(make_pair<unsigned, unsigned>(m_clusteredSpectra[i].m_fileIndex,
                                                   m_clusteredSpectra[i].m_index));
    return ret;
  }
  ////////////////////////////////////////////////////////////////////////////////
  void Cluster::dump(ostream &sout, bool web)
  {
    sout << "--------------------------------" << endl;
    sout << "** Cluster " << m_index << endl;
    sout << "--------------------------------" << endl;
    for (int i = 0; i < m_clusteredSpectra.size(); i++)
    {
      sout << "(" << m_clusteredSpectra[i].m_fileIndex << ";"
          << m_clusteredSpectra[i].m_filename << m_clusteredSpectra[i].m_index
          << ";" << m_clusteredSpectra[i].m_scan << ")" << endl;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  int Cluster::checkFileInGroup(std::set<string> &fileGroup, string &filename)
  {
    if(fileGroup.find(filename) != fileGroup.end()){
        return 1;
    }
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////////

  int Cluster::getSummary(std::vector<int> &localGroups,
                          vector<SpecSet> &specs,
                          SpecSet &specs_ms,
                          std::map<string, std::set<string> > &groups,
                          int &numSpectra,
                          double &sumPercursorIntensity,
                          vector<string> *inputFiles)
  {
    for (int i = 0; i < m_clusteredSpectra.size(); i++)
    {
      // get the file index
      int fileIndex = m_clusteredSpectra[i].m_fileIndex;
      int spectrumIndex = m_clusteredSpectra[i].m_index;

      //get filename
      string filename = m_clusteredSpectra[i].m_filename;
      ;

      if (inputFiles)
        if (fileIndex < inputFiles->size())
          filename = (*inputFiles)[fileIndex];

      // get group for this spectrum, and add it
      int idx = 0;
      for (map<std::string, std::set<string> >::iterator it = groups.begin();
          it != groups.end(); it++)
      {
        int ret = checkFileInGroup(it->second, filename);
        if (ret >= 0){
          localGroups[idx]++;
        }
        idx++;
      }

      // add data
      numSpectra++;
      if (fileIndex >= 0 && fileIndex < specs.size() &&
          spectrumIndex < specs[fileIndex].size())
        sumPercursorIntensity +=
            specs[fileIndex][spectrumIndex].precursor_intensity;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////

  int Cluster::getCsvData(OutputTable &data,
                          vector<SpecSet> &specs,
                          vector<string> *inputFiles,
                          vector<vector<vector<int> > > *scanNumbers)
  {
    // Get record
    int cluster = m_index;

    for (int i = 0; i < m_clusteredSpectra.size(); i++)
    {

      // get spectrum index
      int spectrumIndex = m_clusteredSpectra[i].m_index;
      string spectrumIndex2 = parseInt(spectrumIndex);

      // get file index
      int fileIndex = m_clusteredSpectra[i].m_fileIndex;

      // Get Scan number
      int scan = m_clusteredSpectra[i].m_scan;

      // Get Scan number
      if (scanNumbers && scan == -1)
      {
        if (scanNumbers)
          if (fileIndex < scanNumbers->size())
          {
            if (spectrumIndex < (*scanNumbers)[fileIndex].size())
              scan = (*scanNumbers)[fileIndex][spectrumIndex][0];
          }
      }

      // get the filename
      string filename;
      if (m_clusteredSpectra[i].m_filename.length())
        filename = m_clusteredSpectra[i].m_filename;
      else
      {
        filename = "File index: ";
        filename += parseInt(m_clusteredSpectra[i].m_fileIndex);
      }

      if (inputFiles)
        if (fileIndex >= 0 && fileIndex < inputFiles->size())
          filename = (*inputFiles)[fileIndex];

      float parentMass = -1.0;
      int parentCharge = -1;
      float retention_time = -1.0;
      float precursor_intensity = -1.0;

      if (fileIndex >= 0 && fileIndex < specs.size())
        if (spectrumIndex < specs[fileIndex].size())
        {
          parentMass = specs[fileIndex][spectrumIndex].parentMass;
          parentCharge = specs[fileIndex][spectrumIndex].parentCharge;
          retention_time = specs[fileIndex][spectrumIndex].retention_time;
          precursor_intensity =
              specs[fileIndex][spectrumIndex].precursor_intensity;
        }

      // store the data
      vector<pair<string, bool> > row;
      stringstream aux;

      int cc = cluster + 1;
      aux << cc;
      row.push_back(pair<string, bool>(aux.str(), false));

      row.push_back(pair<string, bool>(filename, false));

      row.push_back(pair<string, bool>(spectrumIndex2, false));

      aux.str(std::string());
      aux << scan;
      row.push_back(pair<string, bool>(aux.str(), false));

      aux.str(std::string());
      aux << parentMass;
      row.push_back(pair<string, bool>(aux.str(), false));

      aux.str(std::string());
      aux << parentCharge;
      row.push_back(pair<string, bool>(aux.str(), false));

      aux.str(std::string());
      aux << retention_time;
      row.push_back(pair<string, bool>(aux.str(), false));

      aux.str(std::string());
      aux << precursor_intensity;
      row.push_back(pair<string, bool>(aux.str(), false));

      data.values.push_back(row);
    }

    return 1;
  }

  int Cluster::getCsvDataHeader(OutputTable &data)
  {
    vector<pair<string, bool> > row;
    row.push_back(pair<string, bool>("#ClusterIdx", false));
    row.push_back(pair<string, bool>("#Filename", false));
    row.push_back(pair<string, bool>("#SpecIdx", false));
    row.push_back(pair<string, bool>("#Scan", false));
    row.push_back(pair<string, bool>("#ParentMass", false));
    row.push_back(pair<string, bool>("#Charge", false));
    row.push_back(pair<string, bool>("#RetTime", false));
    row.push_back(pair<string, bool>("#PrecIntensity", false));
    data.values.push_back(row);
  }
  ////////////////////////////////////////////////////////////////////////////////

  bool Cluster::saveToBinaryStream(FILE *fp) const
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned int count;
    count = fwrite(&m_scan, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save scan #");
      return false;
    }
    count = fwrite(&m_index, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save index #");
      return false;
    }
    count = fwrite(&m_fileIndex, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save file index #");
      return false;
    }
    vector<string> fName(1);
    fName[0] = m_filename;
    if (!writeStringsToBinaryStream(fp, fName))
    {
      ERROR_MSG("Failed to save filename");
      return false;
    }

    const unsigned int numItems = size();

    count = fwrite(&numItems, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save # of spectra");
      return false;
    }

    if (numItems == 0)
    {
      return true;
    }

    int* indexBuf = (int*)malloc(sizeof(int) * numItems);
    for (unsigned int i = 0; i < numItems; i++)
    {
      indexBuf[i] = m_clusteredSpectra[i].m_scan;
    }
    count = fwrite(indexBuf, sizeof(int), numItems, fp);

    if (count == 0)
    {
      ERROR_MSG("Failed to save clustered scan #s");
      free(indexBuf);
      return false;
    }

    for (unsigned int i = 0; i < numItems; i++)
    {
      indexBuf[i] = m_clusteredSpectra[i].m_index;
    }
    count = fwrite(indexBuf, sizeof(int), numItems, fp);

    if (count == 0)
    {
      ERROR_MSG("Failed to save clustered indices");
      free(indexBuf);
      return false;
    }

    for (unsigned int i = 0; i < numItems; i++)
    {
      indexBuf[i] = m_clusteredSpectra[i].m_fileIndex;
    }
    count = fwrite(indexBuf, sizeof(int), numItems, fp);
    free(indexBuf);
    if (count == 0)
    {
      ERROR_MSG("Failed to save clustered filename indices");
      return false;
    }

    vector<string> fNames(numItems);
    for (unsigned int i = 0; i < numItems; i++)
    {
      fNames[i] = m_clusteredSpectra[i].m_filename;
    }
    if (!writeStringsToBinaryStream(fp, fNames))
    {
      ERROR_MSG("Failed to save clustered filenames");
      return false;
    }

    return true;
  }

  bool Cluster::loadFromBinaryStream(FILE* fp,
                                     map<string, unsigned short>& versions)
  {
    if (fp == 0)
    {
      return false;
    }
    initialize(-1);

    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    unsigned int count;

    count = fread(&m_scan, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read scan #");
      return false;
    }

    count = fread(&m_index, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read index #");
      return false;
    }

    count = fread(&m_fileIndex, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read file index #");
      return false;
    }

    vector<string> fName;
    if ((!readStringsFromBinaryStream(fp, fName)) || fName.size() != 1)
    {
      ERROR_MSG("Failed to read filename");
      return false;
    }
    m_filename = fName[0];

    unsigned int numItems = 0;
    count = fread(&numItems, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of clustered spectra");
      return false;
    }
    m_clusteredSpectra.resize(numItems);

    if (numItems == 0)
    {
      return true;
    }

    int* indexBuf = (int*)malloc(sizeof(int) * numItems);
    count = fread(indexBuf, sizeof(int), numItems, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read scan #s");
      free(indexBuf);
      return false;
    }
    for (unsigned int i = 0; i < numItems; i++)
    {
      m_clusteredSpectra[i].m_scan = indexBuf[i];
    }

    count = fread(indexBuf, sizeof(int), numItems, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read indices");
      free(indexBuf);
      return false;
    }
    for (unsigned int i = 0; i < numItems; i++)
    {
      m_clusteredSpectra[i].m_index = indexBuf[i];
    }

    count = fread(indexBuf, sizeof(int), numItems, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read filename indices");
      free(indexBuf);
      return false;
    }
    for (unsigned int i = 0; i < numItems; i++)
    {
      m_clusteredSpectra[i].m_fileIndex = indexBuf[i];
    }
    free(indexBuf);

    vector<string> fNames;
    if ((!readStringsFromBinaryStream(fp, fNames)) || fNames.size() != numItems)
    {
      ERROR_MSG("Failed to read filenames");
      return false;
    }
    for (unsigned int i = 0; i < numItems; i++)
    {
      m_clusteredSpectra[i].m_filename = fNames[i];
    }

    return true;
  }
}
