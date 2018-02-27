/*
 * ExecQualityControl.h
 *
 *  Created on: Jul 8, 2013
 *      Author: aguthals
 */

#ifndef EXECQUALITYCONTROL_H_
#define EXECQUALITYCONTROL_H_

// Module Includes
#include "ExecBase.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "ExecReportSPSStats.h"
#include "ExecStatistics.h"
#include "spectrum_scoring.h"

// External Includes
#include "utils.h"
#include "ExecMergeConvert.h"

namespace specnets
{
  class ExecQualityControl : public ExecReportSPSStats
  {
  public:

    static void loadMS2ScoringModels(const string &modelDir,
                                     MS2ScoringModelSet &modelSet);

    ExecQualityControl();

    ExecQualityControl(const ParameterList & inputParams);

    virtual ~ExecQualityControl(void);

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

    // All files loaded from SPS project path
    ParameterList *m_paramsSPS; // Project parameters (sps.params)
    SpecSet *m_spectraRawMS2; // Un-clustered input MS/MS spectra
    SpecSet *m_spectraClustMS2; // Clustered MS/MS spectra (specs_ms.pklbin)
    SpecSet *m_spectraScored; // Scored PRM spectra (specs_scored.pklbin)
    MS2ScoringModelSet *m_ms2ScoringModels; // All scoring models

    /*
     * Remaining SPS files are declared in ExecReportSPSStats.h
     *
     * SpecSet* m_spectraStars; // star spectra assembled into contigs (stars.pklbin)
     * SpecSet* m_contigs; // input contigs (sps_seqs.pklbin)
     * abinfo_t* m_abinfo; // abinfo detailing what was assembled in each contig (component_info.bin)
     * SpecSet* m_contigOverlaps; // matchma results overlapping contigs with target proteins
     * vector<vector<int> >* m_contigProtMatch; // matchma results assigning contigs to target proteins
     * DB_fasta* m_fasta; // target proteins
     * PeptideSpectrumMatchSet* m_spectraClustIDs; // database search results for clustered spectra
     * PeptideSpectrumMatchSet *m_spectraRawIDs; // database search results for un-clustered spectra
     * vector<string>* m_targetProts; // subset of targeted FASTA database proteins
     * SpecSet* m_proteinSpectra; // Spectral representation of FASTA database
     * vector<string>* m_peptides; // Un-modified peptide ID strings for matching to protein
     *
     * // Maps 1-based spectrum index to scan number and filename of raw MS/MS spectrum
     * map<int, list<pair<int, string> > > m_rawScanInfo;
     * // Maps 1-based spectrum index of clustered spectra to scan number and filename of raw MS/MS spectrum
     * map<int, list<pair<int, string> > > m_clustScanInfo;
     *
     * MappedSpecnets* m_mappedProj; // SPS project mapped to target proteins
     */

    // DECLARE ADDITIONAL STATISTICS MODULES HERE
    OutputTable *m_statsSpectraRawMS2;
    OutputTable *m_statsSpectraClustMS2;
    OutputTable *m_statsSpectraScored;
    OutputTable *m_statsSpectraStar;
    MappedContigSetTable *m_statsContigs;

    OutputTable *m_statsSummaryRawMS2;
    OutputTable *m_statsSummaryClustMS2;
    OutputTable *m_statsSummaryScored;
    OutputTable *m_statsSummaryStar;
    MappedSPSStatTable *m_statsSummaryContigs;

    void prepareStatsSpectra(const PeptideSpectrumMatchSet &psms, const map<
        int, list<pair<int, string> > > &scanInfo, OutputTable &statsTable);

    void prepareStatsSummary(const PeptideSpectrumMatchSet &psms, const map<
        int, list<pair<int, string> > > &scanInfo, OutputTable &statsTable);

  };
}

#endif /* EXECQUALITYCONTROL_H_ */
