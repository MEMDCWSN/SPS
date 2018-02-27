/*
 * ExecSpectralProbability.h
 *
 *  Created on: Nov 5, 2013
 *      Author: aguthals
 */

#ifndef EXECSPECTRALPROBABILITY_H_
#define EXECSPECTRALPROBABILITY_H_

#include "ExecBase.h"
#include "SpecSet.h"
#include "PeptideSpectrumMatchSet.h"
#include "ClusterSet.h"
#include "MassErrorModel.h"

namespace specnets
{
  class ExecSpectralProbability : public ExecBase
  {
  public:

    static const int DEBUG_SCAN;

    ExecSpectralProbability(void);

    ExecSpectralProbability(const ParameterList & inputParams);

    ExecSpectralProbability(const ParameterList &inputParams,
                            vector<SpecSet*> *inputSpectra,
                            PeptideSpectrumMatchSet *inputPSMs,
                            MassErrorModel *inputErrorModelTarget,
                            MassErrorModel *inputErrorModelDecoy,
                            SpecSet *inputMsSpectra,
                            AAJumps *jumps);

    ExecSpectralProbability(const ParameterList &inputParams,
                            vector<SpecSet*> *inputSpectra,
                            PeptideSpectrumMatchSet *inputPSMs,
                            MassErrorModel *inputErrorModelTarget,
                            MassErrorModel *inputErrorModelDecoy,
                            SpecSet *inputMsSpectra,
                            AAJumps *jumps,
                            SpecSet *outputSpectra,
                            ClusterSet *outputClusters,
                            PeptideSpectrumMatchSet *outputPSMs);

    virtual ~ExecSpectralProbability(void);

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

    void removeDecoyPSMsMatchingTarget();

    bool ownInput;
    bool ownOutput;

    bool decoysRemoved;

    vector<SpecSet*> *m_inputSpectra;
    PeptideSpectrumMatchSet *m_inputPSMs;
    MassErrorModel *m_inputErrorModelTarget;
    MassErrorModel *m_inputErrorModelDecoy;
    SpecSet *m_inputMsSpectra;

    SpecSet *m_outputSpectra;
    PeptideSpectrumMatchSet *m_outputPSMs;
    ClusterSet *m_outputClusters;

    AAJumps *m_jumps;
  };
}

#endif /* EXECSPECTRALPROBABILITY_H_ */
