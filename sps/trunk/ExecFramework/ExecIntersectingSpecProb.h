/*
 * ExecIntersectingSpecProb.h
 *
 *  Created on: Nov 11, 2013
 *      Author: aguthals
 */

#ifndef EXECINTERSECTINGSPECPROB_H_
#define EXECINTERSECTINGSPECPROB_H_

#include "ExecBase.h"
#include "SpecSet.h"
#include "PeptideSpectrumMatchSet.h"
#include "ClusterSet.h"
#include "SpectrumPairSet.h"
#include "MassErrorModel.h"

using namespace std;

namespace specnets
{
  class ExecIntersectingSpecProb : public ExecBase
  {
  public:

    static void generatePSMAligns(const ParameterList &inputParams,
                                  const AAJumps &inputJumps,
                                  const SpecSet &inputSpectra,
                                  const ClusterSet &inputClusters,
                                  SpectrumPairSet &outputPairs);

    //static const int DEBUG_IDX1;
    //static const int DEBUG_IDX2;

    ExecIntersectingSpecProb(void);

    ExecIntersectingSpecProb(const ParameterList & inputParams);

    ExecIntersectingSpecProb(const ParameterList &inputParams,
                              SpecSet *inputSpectra,
                              PeptideSpectrumMatchSet *inputPSMs,
                              SpectrumPairSet *inputPairs,
                              ClusterSet* inputClusters,
                              MassErrorModel *inputErrorModelTarget,
                              MassErrorModel *inputErrorModelDecoy,
                              SpecSet *inputMsSpectra,
                              AAJumps *jumps);

    ExecIntersectingSpecProb(const ParameterList &inputParams,
                              SpecSet *inputSpectra,
                              PeptideSpectrumMatchSet *inputPSMs,
                              SpectrumPairSet *inputPairs,
                              ClusterSet* inputClusters,
                              MassErrorModel *inputErrorModelTarget,
                              MassErrorModel *inputErrorModelDecoy,
                              SpecSet *inputMsSpectra,
                              AAJumps *jumps,
                              SpectrumPairSet *outputPairs,
                              vector<vector<double> > *outputProbs);

    virtual ~ExecIntersectingSpecProb(void);

    virtual ExecBase * clone(const ParameterList & inputParams) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

  protected:

    bool ownInput;
    bool ownOutput;

    SpecSet *m_inputSpectra;
    PeptideSpectrumMatchSet *m_inputPSMs;
    SpectrumPairSet *m_inputPairs;
    ClusterSet *m_inputClusters;
    AAJumps *m_jumps;
    MassErrorModel *m_inputErrorModelTarget;
    MassErrorModel *m_inputErrorModelDecoy;
    SpecSet *m_inputMsSpectra;

    SpectrumPairSet *m_outputPairs;
    vector<vector<double> > *m_outputProbs;
  };
}

#endif /* EXECINTERSECTINGSPECPROB_H_ */
