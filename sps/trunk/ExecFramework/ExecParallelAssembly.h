/*
 * ExecParallelAssembly.h
 *
 *  Created on: Mar 26, 2012
 *      Author: aguthals
 */

#ifndef EXECPARALLELASSEMBLY_H_
#define EXECPARALLELASSEMBLY_H_

// Module Includes
#include "ExecBase.h"
#include "ExecMergeConvert.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes
#include "spectrum.h"
#include "SpecSet.h"
#include "AbruijnGraph2.h"
#include "abruijn.h"
#include "clusters.h"
#include "utils.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecParallelAssembly : public ExecBase
  {
  public:
    ExecParallelAssembly(void);

    ExecParallelAssembly(const ParameterList & inputParams);

    ExecParallelAssembly(const ParameterList & inputParams,
                         abinfo_t * inputAbruijn,
                         SpecSet * inputStars,
                         Clusters * outputContigShifts,
                         abinfo_t * outputAbruijn);

    ExecParallelAssembly(const ParameterList & inputParams,
                         abinfo_t * inputAbruijn,
                         SpecSet * inputStars);

    virtual ~ExecParallelAssembly(void);

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

    abinfo_t * m_inputAbruijn; //! Input abinfo
    SpecSet * m_stars; //! Input stars assembled by input contigs
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

    Clusters * m_contigShifts; //! Sets of shifts/endpoints for spectra in assembled contigs
    abinfo_t * m_outputAbruijn; //! Specifies which spectrum peaks were assembled into which abruijn vertices (see abruijn.h)
    bool ownOutput; //! Does this object "own" the output data pointers (and hence has to free them)
  };
}

#endif /* EXECPARALLELASSEMBLY_H_ */
