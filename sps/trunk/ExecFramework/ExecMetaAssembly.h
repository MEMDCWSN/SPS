/*
 * ExecMetaAssembly.h
 *
 *  Created on: Mar 26, 2012
 *      Author: aguthals
 */

#ifndef EXECMETAASSEMBLY_H_
#define EXECMETAASSEMBLY_H_

// Module Includes
#include "ExecBase.h"
#include "ExecMergeConvert.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "ExecParallelAssembly.h"

// External Includes
#include "spectrum.h"
#include "SpecSet.h"
#include "MetaSPS/ContigNetwork.h"
#include "MetaSPS/CombineContigs.h"
#include "SpectrumPairSet.h"
#include "utils.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecMetaAssembly : public ExecBase
  {
  public:
    ExecMetaAssembly(void);

    ExecMetaAssembly(const ParameterList & inputParams);

    ExecMetaAssembly(const ParameterList & inputParams,
                     SpecSet * inputContigs,
                     SpectrumPairSet * inputPairs,
                     abinfo_t * inputAbruijn,
                     SpecSet * inputStars,
                     Clusters * outputContigShifts,
                     abinfo_t * outputAbruijn);

    ExecMetaAssembly(const ParameterList & inputParams,
                     SpecSet * inputContigs,
                     SpectrumPairSet * inputPairs,
                     abinfo_t * inputAbruijn,
                     SpecSet * inputStars);

    virtual ~ExecMetaAssembly(void);

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
    SpecSet * m_spectra; //! Input contig spectra
    SpectrumPairSet * m_spectrumPairs; //! Input set of spectral pairs
    abinfo_t * m_inputAbruijn; //! Input abinfo
    SpecSet * m_stars; //! Input stars assembled by input contigs
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

    SpecSet* m_overlaps;
    vector<vector<int> > * m_prot_match;
    DB_fasta* m_fasta;

    Clusters * m_contigShifts; //! Sets of shifts/endpoints for spectra in assembled contigs
    abinfo_t * m_outputAbruijn; //! Specifies which spectrum peaks were assembled into which abruijn vertices (see abruijn.h)
    bool ownOutput; //! Does this object "own" the output data pointers (and hence has to free them)
  };
}

#endif /* EXECMETAASSEMBLY_H_ */
