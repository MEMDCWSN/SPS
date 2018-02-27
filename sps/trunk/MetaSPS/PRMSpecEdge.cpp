/*
 * PRMSpecEdge.cpp
 *
 *  Created on: Jun 23, 2011
 *      Author: aguthals
 */

#include "PRMSpecEdge.h"

namespace specnets
{

  void PRMSpecEdge::appendEdges(PRMSpecEdge* edge1,
                                PRMSpecEdge* edge2,
                                PRMSpecEdge* outEdge1Plus2)
  {

    if (edge1->spec2rev)
    {
      ERROR_MSG("Cannot append an edge to a reversed edge "<<
          "(can be done but not yet supported). edge1 cannot "<<
          "be reversed but edge2 can");
      abort();
    }

    int node1, node2, node3;

    if (edge1->spec1 == edge2->spec1)
    {
      node1 = edge1->spec2;
      node2 = edge1->spec1;
      node3 = edge2->spec2;
    }
    else if (edge1->spec1 == edge2->spec2)
    {
      node1 = edge1->spec2;
      node2 = edge1->spec1;
      node3 = edge2->spec1;
    }
    else if (edge1->spec2 == edge2->spec1)
    {
      node1 = edge1->spec1;
      node2 = edge1->spec2;
      node3 = edge2->spec2;
    }
    else if (edge1->spec2 == edge2->spec2)
    {
      node1 = edge1->spec1;
      node2 = edge1->spec2;
      node3 = edge2->spec1;
    }
    else
    {
      ERROR_MSG("Cannot append edges if indices do not match");
      DEBUG_VAR(edge1->spec1);
      DEBUG_VAR(edge1->spec2);
      DEBUG_VAR(edge2->spec1);
      DEBUG_VAR(edge2->spec2);
      abort();
    }

    float forwardShift1 = edge1->getShift(node1, node2);
    float reversedShift1 = edge1->getReversedShift(node1, node2);

    float forwardShift2;
    float reversedShift2;

    if (edge2->spec2rev && edge2->spec2 == node2)
    {
      forwardShift2 = edge2->getReversedShift(node2, node3);
      reversedShift2 = edge2->getShift(node2, node3);
    }
    else
    {
      forwardShift2 = edge2->getShift(node2, node3);
      reversedShift2 = edge2->getReversedShift(node2, node3);
    }

    outEdge1Plus2->spec1 = node1;
    outEdge1Plus2->spec2 = node3;
    outEdge1Plus2->shift1 = forwardShift1 + forwardShift2;
    outEdge1Plus2->shift2 = reversedShift1 + reversedShift2;
    outEdge1Plus2->spec2rev = edge2->spec2rev;
  }

  PRMSpecEdge& PRMSpecEdge::operator=(const PRMSpecEdge &other)
  {
    index = other.index;
    SpectrumPair::operator =((SpectrumPair&)other);
  }

  PRMSpecEdge& PRMSpecEdge::operator=(const SpectrumPair &other)
  {
    SpectrumPair::operator =((SpectrumPair&)other);
  }

  void PRMSpecEdge::initializeEdge(const int idx,
                                   const SpectrumPair& contigPair,
                                   const SpecSet& contigs)
  {
    (*this) = contigPair;
    index = idx;
    shift2 = contigs[spec1].parentMass - (shift1 + contigs[spec2].parentMass);
  }

  float PRMSpecEdge::getShift(int idx1, int idx2) const
  {
    if (idx1 == spec1 && idx2 == spec2)
    {
      return shift1;
    }
    else if (idx2 == spec1 && idx1 == spec2)
    {
      return 0.0 - shift1;
    }
    else
    {
      ERROR_MSG("Invalid indices " << idx1 << " and " << idx2
          << ", looking for " << spec1 << " and " << spec2);
      abort();
      return 0;
    }
  }

  float PRMSpecEdge::getReversedShift(int idx1, int idx2) const
  {
    if (idx1 == spec1 && idx2 == spec2)
    {
      return shift2;
    }
    else if (idx2 == spec1 && idx1 == spec2)
    {
      return 0.0 - shift2;
    }
    else
    {
      ERROR_MSG("Invalid indices " << idx1 << " and " << idx2
          << ", looking for " << spec1 << " and " << spec2);
      abort();
      return 0;
    }
  }

  void PRMSpecEdge::reverse(bool reverse1, bool reverse2)
  {
    if (reverse1 && reverse2)
    {
      float temp = shift1;
      shift1 = shift2;
      shift2 = temp;
    }
    else if (reverse1 && (!reverse2))
    {
      spec2rev = !spec2rev;
      float temp = shift1;
      shift1 = shift2;
      shift2 = temp;
    }
    else if ((!reverse1) && reverse2)
    {
      spec2rev = !spec2rev;
    }

  }

  void PRMSpecEdge::reverse(int idx)
  {
    if (idx != spec1 && idx != spec2)
    {
      ERROR_MSG("Invalid index " << idx << ", looking for " << spec1 << " or " << spec2);
      abort();
      return;
    }

    if (idx == spec1)
    {
      reverse(true, false);
    }
    else
    {
      reverse(false, true);
    }
  }
}
