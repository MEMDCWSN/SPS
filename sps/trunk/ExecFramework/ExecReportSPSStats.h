/*
 * ExecReportSPSStats.h
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#ifndef EXECREPORTSPSSTATS_H_
#define EXECREPORTSPSSTATS_H_

// Module Includes
#include "ExecBase.h"
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes
#include "db_fasta.h"
#include "PeptideSpectrumMatchSet.h"
#include "MappedSPSStatTable.h"
#include "MappedContigStatTable.h"
#include "MappedContigSetTable.h"
#include "utils.h"
#include "ExecMergeConvert.h"
#include "ClusterData.h"

// System Includes
#include <string>
#include <vector>

namespace specnets
{
  class ExecReportSPSStats: public ExecBase
  {
  public:
    ExecReportSPSStats(void);

    ExecReportSPSStats(const ParameterList & inputParams);

    virtual ~ExecReportSPSStats(void);

    virtual ExecBase * clone(const ParameterList & input_params) const;

    virtual bool invoke(void);

    virtual bool loadInputData(void);

    virtual bool saveOutputData(void);

    virtual bool saveInputData(std::vector<std::string> & filenames);

    virtual bool loadOutputData(void);

    virtual std::vector<ExecBase *> const & split(int numSplit);

    virtual bool merge(void);

    virtual bool validateParams(std::string & error);

    void FilterFastaProteins(set<int>& target_proteins,
                             vector<string>& put_proteins);

    void FilterSpecIds(vector<string>& put_peptides);

    void ChopContigEnds(int endsChop);

  protected:
    SpecSet* m_contigs; // input contigs
    SpecSet* m_spectraStars; // star spectra assembled into contigs
    abinfo_t* m_abinfo; // abinfo detailing what was assembled in each contig
    SpecSet* m_contigOverlaps; // matchma results overlapping contigs with target proteins
    vector<vector<int> >* m_contigProtMatch; // matchma results assigning contigs to target proteins
    DB_fasta* m_fasta; // target proteins
    MS2ScoringModel* m_model; // scoring model for annotating spectra with b/y ions
    PeptideSpectrumMatchSet* m_spectraClustIDs; // database search results for clustered spectra
    PeptideSpectrumMatchSet *m_spectraRawIDs; // database search results for un-clustered spectra
    vector<string>* m_targetProts; // subset of targeted FASTA database proteins
    SpecSet* m_proteinSpectra; // Spectral representation of FASTA database
    vector<string>* m_peptides; // Un-modified peptide ID strings for matching to protein

    // Maps 1-based spectrum index to scan number and filename of raw MS/MS spectrum
    map<int, list<pair<int, string> > > m_rawScanInfo;
    // Maps 1-based spectrum index of clustered spectra to scan number and filename of raw MS/MS spectrum
    map<int, list<pair<int, string> > > m_clustScanInfo;


    MappedSpecnets* m_mappedProj; // SPS project mapped to target proteins
    bool ownOutput;
    bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

  };

}

#endif /* EXECREPORTSPSSTATS_H_ */
