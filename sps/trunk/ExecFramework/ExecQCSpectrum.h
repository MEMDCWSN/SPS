/*
 * ExecQCSpectrum.h
 *
 *  Created on: Aug 29, 2013
 *      Author: aguthals
 */

#ifndef EXECQCSPECTRUM_H_
#define EXECQCSPECTRUM_H_

// Module Includes
#include "ExecBase.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "SpectrumAnnotStatistics.h"
#include "SpectrumAnnotParameterList.h"
#include "spectrum_scoring.h"
#include "ClusterSet.h"

// External Includes
#include "utils.h"
#include "ExecMergeConvert.h"

namespace specnets
{
  class ExecQCSpectrum : public ExecBase
  {
  public:

    static void loadMS2ScoringModels(const string &modelDir,
                                     MS2ScoringModelSet &modelSet);

    ExecQCSpectrum();

    ExecQCSpectrum(const ParameterList & inputParams);

    ExecQCSpectrum(const ParameterList & inputParams,
                   MS2ScoringModelSet *ms2ScoringModels,
                   SpecSet *inputSpectra,
                   ClusterSet *inputClusters,
                   PeptideSpectrumMatchSet *inputPSMs,
                   SpectrumAnnotParameterList * statsParams,
                   PeptideSpectrumMatchSet *outputPSMs,
                   std::vector<vector<float> > * spectraStats,
                   std::vector<string> * spectraStatsHeader);

    virtual ~ExecQCSpectrum(void);

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

    MS2ScoringModelSet *m_ms2ScoringModels; // All scoring models
    SpecSet *m_inputSpectra;
    ClusterSet *m_inputClusters;
    PeptideSpectrumMatchSet *m_inputPSMs;

    SpectrumAnnotParameterList * m_statsParams; //! The statistics that we are considering

    PeptideSpectrumMatchSet *m_outputPSMs;

    std::vector<vector<float> > * m_spectraStats; //! Vector of vectors containing peptide statistics
    std::vector<string> * m_spectraStatsHeader; //!Vector containing m_spectraStats header information


    bool ownInput;
    bool ownOutput;

    /*
     void prepareStatsSpectra(const PeptideSpectrumMatchSet &psms, const map<
     int, list<pair<int, string> > > &scanInfo, OutputTable &statsTable);

     void prepareStatsSummary(const PeptideSpectrumMatchSet &psms, const map<
     int, list<pair<int, string> > > &scanInfo, OutputTable &statsTable);
     */

  };
}

#endif /* EXECQCSPECTRUM_H_ */
