#include "ExecSpecTagGen.h"

#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"

#define DEBUG_TAG_GEN 0


using namespace std;
using namespace specnets;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecSpecTagGen::ExecSpecTagGen(void) :
    m_contigAbinfo(0x0), m_contigSpectra(0x0), m_starSpectra(0x0), 
    m_matchedContigs(0x0), m_psmContigTags(0x0), m_psmContigTagsDecoy(0x0), 
    ownInput(true), ownOutput(true)
  {
    DEBUG_TRACE
    m_name = "ExecSpecTagGen";
    m_type = "ExecSpecTagGen";
  }

  // -------------------------------------------------------------------------
  ExecSpecTagGen::ExecSpecTagGen(const ParameterList & inputParams) :
    ExecBase(inputParams), 
    m_contigAbinfo(0x0), m_contigSpectra(0x0), m_starSpectra(0x0), 
    m_matchedContigs(0x0), m_psmContigTags(0x0), m_psmContigTagsDecoy(0x0), 
    ownInput(true), ownOutput(true)
  {
    DEBUG_TRACE
    m_name = "ExecSpecTagGen";
    m_type = "ExecSpecTagGen";
  }

  // -------------------------------------------------------------------------
  ExecSpecTagGen::ExecSpecTagGen(const ParameterList & inputParams,
                                 abinfo_t * contigAbinfo,
                                 SpecSet * contigSpectra,
                                 SpecSet * starSpectra,
                                 SpecSet * matchedContigs,
                                 PeptideSpectrumMatchSet * psmContigTags,
                                 PeptideSpectrumMatchSet * psmContigTagsDecoy) :
    ExecBase(inputParams), 
        m_contigAbinfo(contigAbinfo), 
        m_contigSpectra(contigSpectra), 
        m_matchedContigs(matchedContigs),
        m_starSpectra(starSpectra), 
        m_psmContigTags(psmContigTags), 
        m_psmContigTagsDecoy(psmContigTagsDecoy), 
        ownInput(false), ownOutput(false)
  {
    DEBUG_TRACE
    m_name = "ExecSpecTagGen";
    m_type = "ExecSpecTagGen";
  }

  // -------------------------------------------------------------------------
  ExecSpecTagGen::~ExecSpecTagGen(void)
  {
    if (ownInput)
    {
      if (m_contigAbinfo) {
        delete m_contigAbinfo;
      }
      if (m_contigSpectra) {
        delete m_contigSpectra;
      }
      if (m_starSpectra) {
        delete m_starSpectra;
      }
      if (m_matchedContigs) {
        delete m_matchedContigs;
      }
      if (m_psmContigTags) {
        delete m_psmContigTags;
      }
      if (m_psmContigTagsDecoy) {
        delete m_psmContigTagsDecoy;
      }
    }
    if (ownOutput)
    {
      // No output
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecSpecTagGen::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecTagGen(inputParams);
  }

  //---------------------------------------------------------------------------
  void ExecSpecTagGen::makeSpectrumTagsFromContig(
        vector<
            vector<sps::tuple<unsigned int, float, bool> > > & outputAssembledShifts,
        ofstream & ofs)
  {
    SpecSet & matchedContigs = *m_matchedContigs;  // for ease of use
    
    int debugIndex = -1;

    if (m_params.exists("DEBUG_SINGLECONTIG")) {
      debugIndex = m_params.getValueInt("DEBUG_SINGLECONTIG") - 1;
      DEBUG_VAR(debugIndex);
    }

    if (DEBUG_TAG_GEN) {
      DEBUG_VAR(outputAssembledShifts.size());
      for (int i = 0; i < outputAssembledShifts.size(); i++)  {
        DEBUG_MSG(i << "   " << outputAssembledShifts[i].size());
        for (int j = 0; j < outputAssembledShifts[i].size(); j++)  {
          int spectrumIndex = outputAssembledShifts[i][j].m0;
          DEBUG_MSG("        " << spectrumIndex);
        }
      }
    }

    for (int i = 0; i < matchedContigs.size(); i++)
    {
      if (DEBUG_TAG_GEN) DEBUG_VAR(i);
      if (DEBUG_TAG_GEN) DEBUG_VAR(matchedContigs[i].psmList.size());
      if (matchedContigs[i].psmList.size() == 0)
      {
        continue;
      }

      list<psmPtr>::iterator itr = matchedContigs[i].psmList.begin();
      list<psmPtr>::iterator itr_end = matchedContigs[i].psmList.end();
      for (; itr != itr_end; itr++)
      {
        psmPtr origPsm = *itr;
        int contigIndex = origPsm->m_scanNum - 1;
        if (DEBUG_TAG_GEN) DEBUG_MSG("   " << contigIndex);

        if (debugIndex != -1 &&  contigIndex != debugIndex) {
          continue;
        }

        if (DEBUG_TAG_GEN && contigIndex == debugIndex) {
          DEBUG_VAR(origPsm->m_scanNum);
          DEBUG_VAR(origPsm->m_matchOrientation);
          DEBUG_VAR(contigIndex);
          DEBUG_VAR(outputAssembledShifts[contigIndex].size());
        }

        for (int j = 0; j < outputAssembledShifts[contigIndex].size(); j++)
        {
          int spectrumIndex = outputAssembledShifts[contigIndex][j].m0;

          // LARS: This is a hack because there appears to be a problem with spectrum 0
          if (spectrumIndex == 0) {
            continue;
          }
          int spectrumMassShift = (int)outputAssembledShifts[contigIndex][j].m1;

          if (DEBUG_TAG_GEN && contigIndex == debugIndex)
            DEBUG_MSG("  " << contigIndex << "  " << spectrumIndex << "  " << spectrumMassShift);

          psmPtr ptrSpectrumPsm(new PeptideSpectrumMatch);
          *ptrSpectrumPsm = *origPsm;

          if (DEBUG_TAG_GEN && contigIndex == debugIndex)
          {
            DEBUG_VAR(ptrSpectrumPsm->m_startMass);
            DEBUG_VAR(matchedContigs[i].parentMass);
            DEBUG_VAR((*m_starSpectra)[spectrumIndex].parentMass);
          }

          if (origPsm->m_matchOrientation == 1)
          {
            // If the contig was reversed then we need to "reverse" the starting location for the spectrum within the contig
            ptrSpectrumPsm->m_startMass += matchedContigs[i].parentMass
                - spectrumMassShift - (*m_starSpectra)[spectrumIndex].parentMass;
          }
          else
          {
            ptrSpectrumPsm->m_startMass += spectrumMassShift;
          }
          ptrSpectrumPsm->m_scanNum = (*m_starSpectra)[spectrumIndex].scan;

          if (DEBUG_TAG_GEN && contigIndex == debugIndex)
          {
            DEBUG_VAR(ptrSpectrumPsm->m_startMass);
            DEBUG_VAR(ptrSpectrumPsm->m_scanNum);
          }

          // Set the spectrum pointer to this spectrum
          ptrSpectrumPsm->m_spectrum = &((*m_starSpectra)[spectrumIndex]);

          // Is the spectrum reversed wrt the contig?
          //   If so we have to reverse it!
          if (outputAssembledShifts[contigIndex][j].m2 == 1)
          {
            if (DEBUG_TAG_GEN && contigIndex == debugIndex)
              DEBUG_MSG("REVERSING THE REVERSAL FLAG");
            ptrSpectrumPsm->m_matchOrientation = !origPsm->m_matchOrientation;
          }
          ptrSpectrumPsm->saveToFile(ofs);

          if (DEBUG_TAG_GEN && contigIndex == debugIndex)
            DEBUG_TRACE;

        } // for (int j = 0; j < outputAssembledShifts[scanIndex].size(); j++) {
      } // for ( ; itr != itr_end; itr++) {
    } // for (int i = 0; i < matchedContigs.size(); i++) {

    return;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::invoke(void)
  {
    DEBUG_TRACE;
    // Get the starting masses of the each spectrum in each contig
    vector<vector<sps::tuple<unsigned int, float, bool> > > outputAssembledShifts;
    DEBUG_VAR(m_contigAbinfo->size());

    getAssembledShifts(*m_contigSpectra,
                       *m_contigAbinfo,
                       outputAssembledShifts);
    DEBUG_VAR(outputAssembledShifts.size());
    DEBUG_VAR(m_starSpectra->size());

    ofstream ofsTarget(m_params.getValue("OUTPUT_CONTIG_STAR_TAGS_TGT").c_str(), ios::out);
    if (!ofsTarget)
    {
      ERROR_MSG("Error opening [" << 
                   m_params.getValue("OUTPUT_CONTIG_STAR_TAGS_TGT").c_str() << "]");
      exit(-101);
    }
    PeptideSpectrumMatch::saveHeaderToFile(ofsTarget);

    m_matchedContigs->clearPsms();  // Clear any PSMs
    m_psmContigTags->addSpectra(m_matchedContigs); // Associate target PSMs
    DEBUG_VAR(m_psmContigTags->size());

    makeSpectrumTagsFromContig(outputAssembledShifts,
                               ofsTarget);
    ofsTarget.close();

    ofstream ofsDecoy(m_params.getValue("OUTPUT_CONTIG_STAR_TAGS_DEC").c_str(), ios::out);
    if (!ofsDecoy)
    {
      ERROR_MSG("Error opening [" << 
                   m_params.getValue("OUTPUT_CONTIG_STAR_TAGS_DEC").c_str() << "]");
      exit(-101);
    }
    PeptideSpectrumMatch::saveHeaderToFile(ofsDecoy);

    m_matchedContigs->clearPsms();  // Clear any PSMs
    m_psmContigTagsDecoy->addSpectra(m_matchedContigs); // Associate decoy PSMs
    DEBUG_VAR(m_psmContigTagsDecoy->size());

    makeSpectrumTagsFromContig(outputAssembledShifts,
                               ofsDecoy);
    ofsDecoy.close();
    
    m_matchedContigs->clearPsms();  // Leave it clear (safest)

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::loadInputData(void)
  {
    if (ownInput)
    {
      if (!m_contigAbinfo) {
        m_contigAbinfo = new abinfo_t;
      }
      if (!m_contigSpectra) {
        m_contigSpectra = new SpecSet;
      }
      if (!m_starSpectra) {
        m_starSpectra = new SpecSet;
      }
      if (!m_matchedContigs) {
        m_matchedContigs = new SpecSet;
      }
      if (!m_psmContigTags) {
        m_psmContigTags = new PeptideSpectrumMatchSet;
      }
      if (!m_psmContigTagsDecoy) {
        m_psmContigTagsDecoy = new PeptideSpectrumMatchSet;
      }
    }
    if (ownOutput)
    {
      // No output
    }

    if (m_params.exists("INPUT_ABINFO")) {
      DEBUG_MSG("Loading Abinfo.. [" << m_params.exists("INPUT_ABINFO") << "]");
      if (!Load_abinfo(m_params.getValue("INPUT_ABINFO").c_str(), 
                       *m_contigAbinfo)) {
        ERROR_MSG("Error reading contig abinfo file [" << 
                   m_params.getValue("INPUT_ABINFO").c_str() << "]");
        return false;
      }
    }

    if (m_params.exists("INPUT_CONTIGS")) {
      DEBUG_MSG("Loading Contigs.. [" << m_params.exists("INPUT_CONTIGS") << "]");
      if (m_contigSpectra->loadPklBin(m_params.getValue("INPUT_CONTIGS").c_str()) <= 0) {
        ERROR_MSG("Error reading contig shifts file [" << 
                   m_params.getValue("INPUT_CONTIGS").c_str() << "]");
        return false;
      }
    }
    
    if (m_params.exists("INPUT_STARS_PKLBIN")) {
      DEBUG_MSG("Loading Stars Pklbin.. [" << m_params.exists("INPUT_STARS_PKLBIN") << "]");
      if (m_starSpectra->loadPklBin(
              m_params.getValue("INPUT_STARS_PKLBIN").c_str()) < 0) {
        ERROR_MSG("Error reading input star spectra file [" << 
                   m_params.getValue("INPUT_STARS_PKLBIN").c_str() << "]");
        return false;
      }
    }    

    if (m_params.exists("INPUT_MATCHED_CONTIGS_PKLBIN")) {
      DEBUG_MSG("Loading Matched Contigs Pklbin.. [" << m_params.exists("INPUT_MATCHED_CONTIGS_PKLBIN") << "]");
      if (m_matchedContigs->loadPklBin(
              m_params.getValue("INPUT_MATCHED_CONTIGS_PKLBIN").c_str()) < 0) {
        ERROR_MSG("Error reading input matched contigs file [" << 
                   m_params.getValue("INPUT_MATCHED_CONTIGS_PKLBIN").c_str() << "]");
        return false;
      }
    }    

    if (m_params.exists("INPUT_CONTIG_TAGS_TGT")) {
      DEBUG_MSG("Loading Target Contag Tags.. [" << m_params.exists("INPUT_CONTIG_TAGS_TGT") << "]");
      if (!m_psmContigTags->loadFromFile(
               m_params.getValue("INPUT_CONTIG_TAGS_TGT").c_str())) {
        ERROR_MSG("Error reading input contig target tags [" << 
                   m_params.getValue("INPUT_CONTIG_TAGS_TGT").c_str() << "]");
        return false;
      }
      DEBUG_VAR(m_psmContigTags->size());
    }
        
    if (m_params.exists("INPUT_CONTIG_TAGS_DEC")) {
      DEBUG_MSG("Loading Decoy Contag Tags.. [" << m_params.exists("INPUT_CONTIG_TAGS_DEC") << "]");
      if (!m_psmContigTagsDecoy->loadFromFile(
              m_params.getValue("INPUT_CONTIG_TAGS_DEC").c_str())) {
        ERROR_MSG("Error reading input contig decoy tags [" << 
                   m_params.getValue("INPUT_CONTIG_TAGS_DEC").c_str() << "]");
        return false;
      }
      DEBUG_VAR(m_psmContigTagsDecoy->size());
    }
        
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::saveOutputData(void)
  {
    // Output data is written inside invoke
    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::saveInputData(std::vector<std::string> & filenames)
  {
    return true;
  }
  

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::loadOutputData(void)
  {
    return true;
  }
 
  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecSpecTagGen::split(int numSplit)
  {
    //EMPTY
  }

  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::merge(void)
  {
    //EMPTY
  }
 
  // -------------------------------------------------------------------------
  bool ExecSpecTagGen::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("INPUT_ABINFO");
    VALIDATE_PARAM_EXIST("INPUT_CONTIGS");
    VALIDATE_PARAM_EXIST("INPUT_STARS_PKLBIN");
    VALIDATE_PARAM_EXIST("INPUT_MATCHED_CONTIGS_PKLBIN");
    VALIDATE_PARAM_EXIST("OUTPUT_CONTIG_STAR_TAGS_TGT");
    VALIDATE_PARAM_EXIST("OUTPUT_CONTIG_STAR_TAGS_DEC");

    m_isValid = true;
    return true;
  }

}

