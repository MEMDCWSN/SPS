/*
 * ExecFilterContigPairs.h
 *
 *  Created on: Dec 13, 2010
 *      Author: aguthals
 */

#ifndef EXECFILTERCONTIGPAIRS_H_
#define EXECFILTERCONTIGPAIRS_H_

// Module Includes
#include "ExecBase.h"
#include "ExecMergeConvert.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes
#include "spectrum.h"
#include "SpecSet.h"
#include "MetaSPS/prm_alignment.h"
#include "SpectrumPairSet.h"
#include "utils.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecFilterContigPairs : public ExecBase
  {
  public:
    ExecFilterContigPairs(void);

    ExecFilterContigPairs(const ParameterList & inputParams);

    ExecFilterContigPairs(const ParameterList & inputParams,
                          SpecSet * inputContigs,
                          SpectrumPairSet * outputAlignments);

    ExecFilterContigPairs(const ParameterList & inputParams,
                          SpecSet * inputContigs);

    virtual ~ExecFilterContigPairs(void);

    virtual ExecBase * clone(const ParameterList & input_params) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

  private:
    SpecSet * m_inputContigs; //! The set of input contigs to be aligned
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

    SpectrumPairSet * m_outputAlignments; //! Output contig/contig alignments
    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
  };
}

#endif /* EXECFILTERCONTIGPAIRS_H_ */
