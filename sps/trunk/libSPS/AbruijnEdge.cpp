/*
 * AbruijnEdge.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#include "AbruijnEdge.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short AbruijnEdge::BIN_VERSION = 1;
  const unsigned short AbruijnEdge::BIN_SUBVERSION = 1;

  const string AbruijnEdge::BIN_VERSION_ID = "AbruijnEdge_binVersion";
  const string AbruijnEdge::BIN_SUBVERSION_ID = "AbruijnEdge_binSubVersion";

  AbruijnEdge* AbruijnEdge::castEdgePtr(Edge* edgePtr)
  {
    AbruijnEdge* temp = dynamic_cast<AbruijnEdge*> (edgePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Edge \'" << edgePtr->toString() << "\' to an AbruijnEdge");

      abort();
    }
    return temp;
  }

  AbruijnEdge::AbruijnEdge(void) :
    Edge(), m_mass(0), m_rWeight(0), m_fLabel(""), m_jumpLabel(""),
        m_rLabel(""), m_isCorrect(false), m_isAnnotated(false),
        m_modified(false), m_specIDs(), m_expandOnly(false)
  {
  }

  AbruijnEdge::AbruijnEdge(const AbruijnEdge& other) :
    Edge(), m_mass(0), m_rWeight(0), m_fLabel(""), m_jumpLabel(""),
        m_rLabel(""), m_isCorrect(false), m_isAnnotated(false),
        m_modified(false), m_specIDs(), m_expandOnly(false)
  {
    this->operator =(other);
  }

  AbruijnEdge::AbruijnEdge(const list<AbruijnEdge*>& edgeList,
                           const bool& lookReverse) :
    Edge(), m_mass(0), m_rWeight(0), m_fLabel(""), m_jumpLabel(""),
        m_rLabel(""), m_isCorrect(false), m_isAnnotated(false),
        m_modified(false), m_specIDs(), m_expandOnly(false)
  {
    unsigned int numMerged = mergeEdges(edgeList, lookReverse);
  }

  string AbruijnEdge::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << m_fLabel << "," << m_rLabel << " (" << weight << "/" << m_rWeight
        << ")";
    return out.str();
    /*
     out << weight << " \'" << m_fLabel << "\'";
     if (m_jumpLabel != m_fLabel)
     {
     out << "\\\'" << m_jumpLabel << "\'";
     }
     out << " | \'" << m_rLabel << "\'";
     return out.str();*/
  }

  string AbruijnEdge::toString(void) const
  {
    ostringstream out;
    out << Edge::toString() << " : " << getGraphvizLabel();
    return out.str();
  }

  AbruijnEdge & AbruijnEdge::operator=(const AbruijnEdge &other)
  {
    if (this == &other)
    {
      return *this;
    }

    Edge::operator=((const Edge&)other);
    m_mass = other.getMass();
    m_fLabel = other.getFLabel();
    m_jumpLabel = other.getJumpLabel();
    m_rLabel = other.getRLabel();
    m_isCorrect = other.m_isCorrect;
    m_isAnnotated = other.m_isAnnotated;
    m_specIDs = other.m_specIDs;
    m_rWeight = other.m_rWeight;

    return *this;
  }

  void AbruijnEdge::copy(Edge& otherEdge)
  {
    this->operator =(*castEdgePtr(&otherEdge));
  }

  unsigned int AbruijnEdge::mergeEdges(const list<AbruijnEdge*>& edgeList,
                                       const bool& lookReverse,
                                       list<AbruijnEdge*>* usedEdges)
  {
    if (usedEdges != 0)
    {
      usedEdges->clear();
    }

    if (edgeList.size() == 0)
    {
      return 0;
    }
    std::set<string> mergedSpecIDs;
    list<AbruijnEdge*> sortedEdgeList(edgeList);
    sortedEdgeList.sort(AbruijnEdge::MergeSizeCompare(lookReverse));
    double totMass = 0;
    unsigned int numUniqueEdges = 1;

    list<AbruijnEdge*>::const_iterator eIt = sortedEdgeList.begin();
    this->copy(**eIt);
    if (usedEdges != 0)
    {
      usedEdges->push_back(*eIt);
    }
    mergedSpecIDs.insert(m_specIDs.begin(), m_specIDs.end());
    totMass += m_mass;
    eIt++;

    for (; eIt != sortedEdgeList.end(); eIt++)
    {
      bool foundSameSpec = false;
      for (list<string>::const_iterator sIt = (*eIt)->m_specIDs.begin(); sIt
          != (*eIt)->m_specIDs.end(); sIt++)
      {
        if (mergedSpecIDs.count(*sIt) > 0)
        {
          foundSameSpec = true;
          break;
        }
      }
      if (foundSameSpec)
      {
        continue;
      }
      if (usedEdges != 0)
      {
        usedEdges->push_back(*eIt);
      }

      mergedSpecIDs.insert((*eIt)->m_specIDs.begin(), (*eIt)->m_specIDs.end());

      if ((*eIt)->m_fLabel.length() > m_fLabel.length())
      {
        m_fLabel = (*eIt)->m_fLabel;
      }
      if ((*eIt)->m_rLabel.length() > m_rLabel.length())
      {
        m_rLabel = (*eIt)->m_rLabel;
      }

      m_jumpLabel = (*eIt)->getJumpLabel();
      totMass += (*eIt)->getMass();
      weight += (*eIt)->getWeight();
      m_rWeight += (*eIt)->getRWeight();
      numUniqueEdges++;
    }
    if (numUniqueEdges > 0)
    {
      setMass(totMass / ((double)numUniqueEdges));
      loadIDs(mergedSpecIDs.begin(), mergedSpecIDs.end());
    }
    return numUniqueEdges;
  }

  bool AbruijnEdge::saveToBinaryStream(FILE* fp) const
  {
    if (!Edge::saveToBinaryStream(fp))
    {
      return false;
    }
    unsigned int count;

    count = fwrite(&m_mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the AbruijnEdge mass");
      return false;
    }

    count = fwrite(&m_isCorrect, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the AbruijnEdge correctness");
      return false;
    }

    count = fwrite(&m_isAnnotated, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the AbruijnEdge annotated-ness");
      return false;
    }

    count = fwrite(&m_modified, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the AbruijnEdge modified flag");
      return false;
    }

    vector<string> labels(3);
    labels[0] = m_fLabel;
    labels[1] = m_jumpLabel;
    labels[2] = m_rLabel;

    if (!writeStringsToBinaryStream(fp, labels))
    {
      ERROR_MSG("Could not write the AbruijnEdge labels");
      return false;
    }

    return true;
  }

  bool AbruijnEdge::loadFromBinaryStream(FILE* fp,
                                         map<string, unsigned short>& versions)
  {
    if (!Edge::loadFromBinaryStream(fp, versions))
    {
      return false;
    }
    unsigned int count;

    unsigned short edgeVersion = versions[BIN_VERSION_ID];
    unsigned short edgeSubVersion = versions[BIN_SUBVERSION_ID];

    count = fread(&m_mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the AbruijnEdge mass");
      return false;
    }

    count = fread(&m_isCorrect, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the AbruijnEdge correctness");
      return false;
    }

    count = fread(&m_isAnnotated, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the AbruijnEdge annotated-ness");
      return false;
    }

    count = fread(&m_modified, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the AbruijnEdge modified flag");
      return false;
    }

    vector<string> labels;
    if ((!readStringsFromBinaryStream(fp, labels)) || labels.size() != 3)
    {
      ERROR_MSG("Could not read the AbruijnEdge labels");
      DEBUG_VAR(labels.size());
      return false;
    }
    m_fLabel = labels[0];
    m_jumpLabel = labels[1];
    m_rLabel = labels[2];

    return true;
  }
}

