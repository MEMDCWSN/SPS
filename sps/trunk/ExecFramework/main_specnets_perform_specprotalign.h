#include "abruijn.h"
#include "SpecSet.h"

#include <list>
#include <string>
#include <vector>

namespace specnets
{
  class ParameterList;
  class SpectrumPairSet;
  class SpecSet;
  class PenaltyMatrix;
  class Clusters;
  class PeptideSpectrumMatchSet;
  class PenaltyMatrix;
  class DB_fasta;

  bool performPenaltyGen(ParameterList & ip,
                         SpectrumPairSet & filteredPairs,
                         SpecSet & prmSpectra,
                         PenaltyMatrix & penaltyMatrixMods,
                         map<int, map<int, float> > & scanSpecificPenalties);
                         
  bool performContigProtAlign(ParameterList & ip,
                              SpecSet & contigSpectra,
                              DB_fasta & db,
                              DB_fasta & dbDecoy,
                              PenaltyMatrix & penaltyMatrixBlosum,
                              PenaltyMatrix & penaltyMatrixMods,
                              map<int, map<int, float> > & scanSpecificPenalties,
                              abinfo_t & contigAbinfo,
                              PeptideSpectrumMatchSet & filterPsmSet,
                              SpecSet & matchedSpectraAll,
                              SpecSet & matchedSpectra,
                              PeptideSpectrumMatchSet & psmSet,
                              PeptideSpectrumMatchSet & psmSetDecoy,
                              PeptideSpectrumMatchSet & psmSetFdr,
                              bool gridExecutionFlag,
                              bool resume);

  bool performSpecTagGen(ParameterList & ip,
                         abinfo_t & contigAbinfo,
                         SpecSet & contigSpectra,
                         SpecSet & starSpectra,
                         SpecSet & matchedContigs,
                         PeptideSpectrumMatchSet & psmContigTags,
                         PeptideSpectrumMatchSet & psmContigTagsDecoy);

  bool performSpecProtAlign(ParameterList & ip,
                            SpecSet & contigSpectra,
                            SpecSet & prmSpectra,
                            DB_fasta & db,
                            DB_fasta & dbDecoy,
                            PeptideSpectrumMatchSet & psmTag,
                            PeptideSpectrumMatchSet & psmTagDecoy,
                            PenaltyMatrix & penaltyMatrixBlosum,
                            PenaltyMatrix & penaltyMatrixMods,
                            map<int, map<int, float> > & scanSpecificPenalties,
                            abinfo_t & contigAbinfo,
                            PeptideSpectrumMatchSet & filterPsmSet,
                            SpecSet & matchedSpectraAll,
                            SpecSet & matchedSpectra,
                            PeptideSpectrumMatchSet & psmSet,
                            PeptideSpectrumMatchSet & psmSetDecoy,
                            PeptideSpectrumMatchSet & psmSetFdr,
                            bool gridExecutionFlag,
                            bool resume);

  void makeSpectrumTagsFromContig(ParameterList & ip,
                               vector<
                                  vector<sps::tuple<unsigned int, float, bool> > > & 
                                      outputAssembledShifts,
                               ofstream & ofs);
}

