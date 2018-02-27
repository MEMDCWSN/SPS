/*
 * ExecTrainMassError.h
 *
 *  Created on: Mar 6, 2014
 *      Author: aguthals
 */

#ifndef EXECTRAINMASSERROR_H_
#define EXECTRAINMASSERROR_H_

#include "ExecBase.h"
#include "SpecSet.h"
#include "PeptideSpectrumMatchSet.h"
#include "MassErrorModel.h"

using namespace std;

namespace specnets
{
  class ExecTrainMassError : public ExecBase
  {
  public:

    static const int MIN_ALLOWED_PSMS;

    ExecTrainMassError(void);

    ExecTrainMassError(const ParameterList & inputParams);

    ExecTrainMassError(const ParameterList & inputParams,
                       SpecSet * inputSpectra,
                       SpecSet * inputSpectraPRM,
                       PeptideSpectrumMatchSet * inputPSMs);

    ExecTrainMassError(const ParameterList & inputParams,
                       SpecSet * inputSpectra,
                       SpecSet * inputSpectraPRM,
                       PeptideSpectrumMatchSet * inputPSMs,
                       MassErrorModel * outputModelTarget,
                       MassErrorModel * outputModelDecoy);

    virtual ~ExecTrainMassError(void);

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
    SpecSet *m_inputSpectraPRM;
    PeptideSpectrumMatchSet *m_inputPSMs;

    MassErrorModel *m_outputModelTarget;
    MassErrorModel *m_outputModelDecoy;
  };
}

#endif /* EXECTRAINMASSERROR_H_ */
