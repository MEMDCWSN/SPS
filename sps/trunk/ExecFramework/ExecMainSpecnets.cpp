// Header Includes
#include "ExecMainSpecnets.h"

// Module Includes
#include "ExecTagSearch.h"
#include "ExecSpecNetsPropagation.h"
#include "ExecStatistics.h"
#include "ExecSvmStatistics.h"
#include "ExecFdrPeptide.h"
#include "ExecHomologyAssembly.h"
#include "Logger.h"
#include "ClusterData.h"    // Needed when using Adrian's clustering
using namespace spsReports; // Needed when using Adrian's clustering

// SpecNets Includes
#include "tags.h"
#include "tuple.h"
#include "utils.h"

using namespace std;
using namespace specnets;

namespace specnets
{

  ExecMainSpecnets::ExecMainSpecnets(void) :
    m_specsMsSpectra(0x0), m_specsScoredSpectra(0x0), m_starSpectra(0x0),
        m_pairs(0x0), m_model(0x0), m_db(0x0), m_penaltyMatrixBlosum(0x0),
        m_penaltyMatrixMods(0x0), ownInput(true), m_psms(0x0), m_origPsms(0x0),
        m_psmSpectra(0x0), m_psms_midx(0x0), m_psms_mp(0x0),
        m_snets_contigs(0x0), m_snets_midx(0x0), m_snets_mp(0x0),
        ownOutput(true)
  {
    m_name = "ExecMainSpecnets";
    m_type = "ExecMainSpecnets";
  }

  // -------------------------------------------------------------------------

  ExecMainSpecnets::ExecMainSpecnets(const ParameterList & inputParams) :
    ExecBase(inputParams), m_specsMsSpectra(0x0), m_specsScoredSpectra(0x0),
        m_starSpectra(0x0), m_pairs(0x0), m_model(0x0), m_db(0x0), m_penaltyMatrixBlosum(0x0),
        m_penaltyMatrixMods(0x0), ownInput(true), m_psms(0x0), m_psmSpectra(0x0), m_psms_midx(0x0),
        m_psms_mp(0x0), m_snets_contigs(0x0), m_snets_midx(0x0),
        m_snets_mp(0x0), ownOutput(true)
  {
    m_name = "ExecMainSpecnets";
    m_type = "ExecMainSpecnets";
  }

  // -------------------------------------------------------------------------

  ExecMainSpecnets::ExecMainSpecnets(const ParameterList & inputParams,
                                     SpecSet * msSpectra,
                                     SpecSet * scoredSpectra,
                                     SpecSet * starSpectra,
                                     SpectrumPairSet * pairs,
                                     MS2ScoringModel * model,
                                     DB_fasta * db,
                                     PenaltyMatrix * penaltyMatrixBlosum,
                                     PenaltyMatrix * penaltyMatrixMods,
                                     PeptideSpectrumMatchSet * psms,
                                     PeptideSpectrumMatchSet * origPsms,
                                     SpecSet * psms_spectra,
                                     SpecSet * psms_midx,
                                     vector<vector<int> > * psms_mp,
                                     SpecSet * snets_contigs,
                                     SpecSet * snets_midx,
                                     vector<vector<int> > * snets_mp) :
    ExecBase(inputParams),

    m_specsMsSpectra(msSpectra), m_specsScoredSpectra(scoredSpectra),
        m_starSpectra(starSpectra), m_pairs(pairs), m_model(model), m_db(db), m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods),
        ownInput(false), m_psms(psms), m_origPsms(origPsms), m_psmSpectra(psms_spectra),
        m_psms_midx(psms_midx), m_psms_mp(psms_mp),
        m_snets_contigs(snets_contigs), m_snets_midx(snets_midx),
        m_snets_mp(snets_mp), ownOutput(false)
  {
    m_name = "ExecMainSpecnets";
    m_type = "ExecMainSpecnets";
  }

  // -------------------------------------------------------------------------

  ExecMainSpecnets::~ExecMainSpecnets(void)
  {
    if (ownInput)
    {
      if (m_penaltyMatrixMods)
      {
        delete m_penaltyMatrixMods;
        m_penaltyMatrixMods = 0x0;
      }
      if (m_penaltyMatrixBlosum)
      {
        delete m_penaltyMatrixBlosum;
        m_penaltyMatrixBlosum = 0x0;
      }
      if (m_specsMsSpectra)
      {
        delete m_specsMsSpectra;
        m_specsMsSpectra = 0x0;
      }
      if (m_specsScoredSpectra)
      {
        delete m_specsScoredSpectra;
        m_specsScoredSpectra = 0x0;
      }
      if (m_starSpectra)
      {
        delete m_starSpectra;
        m_starSpectra = 0x0;
      }
      if (m_origPsms)
      {
        delete m_origPsms;
        m_origPsms = 0x0;
      }
      if (m_pairs)
      {
        delete m_pairs;
        m_pairs = 0x0;
      }
      if (m_model)
      {
        delete m_model;
        m_model = 0x0;
      }
      if (m_db)
      {
        delete m_db;
        m_db = 0x0;
      }
    }
    if (ownOutput)
    {
      if (m_psms)
      {
        delete m_psms;
        m_psms = 0x0;
      }
      if (m_psmSpectra)
      {
        delete m_psmSpectra;
        m_psmSpectra = 0x0;
      }
      if (m_snets_contigs)
      {
        delete m_snets_contigs;
        m_snets_contigs = 0x0;
      }
      if (m_snets_midx)
      {
        delete m_snets_midx;
        m_snets_midx = 0x0;
      }
      if (m_snets_mp)
      {
        delete m_snets_mp;
        m_snets_mp = 0x0;
      }
    }

  }

  // -------------------------------------------------------------------------

  ExecBase * ExecMainSpecnets::clone(const ParameterList & inputParams) const
  {
    return new ExecMainSpecnets(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::annotateNodes(DB_fasta *db,
                                       PeptideSpectrumMatchSet * psms,
                                       SpecSet * psmSpectra,
                                       ExecSvmStatistics * svm,
                                       bool storeMatchInfo,
                                       AAJumps &aminoacids)
  {

    SpecSet * psms_midx; //! Per-PSM indices of matched spectrum/protein masses
    vector < vector<int> > *psms_mp; //! Per-PSM: best-matched protein (-1 if none), # mods, prefix/suffix match (cols 0-2, respectively)
    SpectrumPairSet * propagation_pairs; //! Set of pairs used for propagations
    SpecSet * propagation_midx; //! Per-propagation indices of matched spectrum/spectrum masses (in m_pepsAsSpecs)
    vector < vector<int> > *propagation_info; //! Per-spectrum propagation/PSM info (see ExecSpecNetsPropagation->m_projInfo)

    if (storeMatchInfo)
    {
      psmSpectra = m_psmSpectra;
      propagation_pairs = &m_propagation_pairs;
      propagation_info = &m_propagation_info;
    }
    else
    {
      psmSpectra = new SpecSet;
      propagation_pairs = new SpectrumPairSet;
      propagation_info = new vector<vector<int> > ;
    }

    // Initialize SVM used in ExecSpecNetsPropagation
    PeptideSpectrumMatchSet psmsSVM;
    DEBUG_TRACE;

    // ---------------------------------------------
    //   ExecSpecNetsPropagation
    // ---------------------------------------------
    ParameterList paramsPropagation(m_params);
    paramsPropagation.setValue("OUTPUT_SPECS_MATCHIDX",
                               "specnets/snets_midx.pklbin");
    paramsPropagation.setValue("OUTPUT_ALIGNS", "specnets/snets_pairs.bin");
    paramsPropagation.setValue("OUTPUT_SPECS_PROJ",
                               "specnets/snets_specs.pklbin");
    paramsPropagation.setValue("OUTPUT_ANNOTINFO", "specnets/snets_annots.bin");
    DEBUG_TRACE;
    ExecSpecNetsPropagation snetsPropagation(paramsPropagation,
                                             m_specsMsSpectra,
                                             m_specsScoredSpectra,
                                             m_starSpectra,
                                             & aminoacids,
                                             svm,
                                             m_model,
                                             m_pairs,
                                             db,
                                             m_penaltyMatrixMods,
                                             m_penaltyMatrixBlosum,
                                             psms,
                                             propagation_pairs,
                                             propagation_info,
                                             psmSpectra);
    DEBUG_TRACE;
    if (not snetsPropagation.invoke())
      return false;
    DEBUG_TRACE;
    if (not snetsPropagation.saveOutputData())
      return false;
    DEBUG_TRACE;
/*
cerr << "ExecMainSpecnets::annotateNodes() -- >> Outputting psmSpectra mp/midx info:\n";
for (unsigned int i = 0; i < psmSpectra->size(); i++)
{
	cerr << " -- -- spectrum " << i << ", scan = " << (*psmSpectra)[i].scan
			<< ", seq = ";
	if (!(*psmSpectra)[i].psmList.empty())
	{
		psmPtr psm = (*psmSpectra)[i].psmList.front();
		cerr << psm->m_annotation << ", mp = " << psm->m_dbIndex << ", midx:\n";
		for (unsigned int j = 0; j < psm->m_matchedPeaks.size(); j++)
			cerr << psm->m_matchedPeaks[j][0] << " " << psm->m_matchedPeaks[j][1]
			                                                                   << "\n";
	}
	else
		cerr << "\n";
}
cerr.flush();
*/
    psms->getPSMSet(psmSpectra);
    psms->saveToFile("specnets/snets_psms_pre-specprob.txt", true);
    DEBUG_TRACE;

    m_psms_midx->savePklBin("specnets/snets_midx_preFDR.pklbin");
    Save_binArray("specnets/snets_mp_preFDR.bin", *m_psms_mp);
/*
cerr<<" -- >> SVM scores:\n";
for(unsigned int i=0; i<m_psms->size(); i++) {
	cerr<<" -- -- psm "<<i<<", spectrum parent mass = "<<(*m_psms)[i]->m_spectrum->parentMass
	<<", spectrum scan = "<<(*m_psms)[i]->m_spectrum->scan << " parent mz " << (*m_psms)[i]->m_spectrum->parentMZ
	<<", seq = "<<(*m_psms)[i]->m_annotation<<", score = "<<(*m_psms)[i]->m_score
	<<" (decoy = "<<(*m_psms)[i]->m_isDecoy<<")\n";
}
cerr.flush();
*/
    if (not storeMatchInfo)
    {
      delete psmSpectra;
      delete psms_midx;
      delete psms_mp;
      delete propagation_pairs;
      delete propagation_midx;
      delete propagation_info;
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::invoke(void)
  {
    unsigned int idxSpec;
    float peakTol = m_params.getValueDouble("TOLERANCE_PEAK", 0.45);
    bool useSVM = m_params.getValueInt("SPECNETS_USE_SVM", 1);

    if (!m_starSpectra)
    {
      return false;
    }

    if (ownOutput)
    {
      m_psms = new PeptideSpectrumMatchSet;
      m_psms_midx = new SpecSet;
      m_psms_mp = new vector<vector<int> > ;

      m_snets_contigs = new SpecSet;
      m_snets_midx = new SpecSet;
      m_snets_mp = new vector<vector<int> > ;

      m_psmSpectra = new SpecSet;
    }

    //reindex spectra
    for (unsigned int i = 0; i < m_starSpectra->size(); i++)
    {
     (*m_starSpectra)[i].scan = i+1;
     (*m_specsMsSpectra)[i].scan = i+1;
     (*m_specsScoredSpectra)[i].scan = i+1;
    }

    // Temporary workaround for messed up scan numbers in ExecPrmScoring
/*
    for (unsigned int i = 0; i < m_starSpectra->size(); i++)
    {
      cerr << "Spectrum " << i << ", scan = " << (*m_starSpectra)[i].scan << ":\n";
      (*m_starSpectra)[i].output(cerr);
    }
*/
    // Propagate target annotations
    DEBUG_TRACE;
//    m_psms->loadInspectResultsFile(m_params.getValue("INSPECT_PSMS").c_str(), true);
/*
cerr<<"Prior INSPECT_PSMS:\n";
for(unsigned int i=0; i<m_origPsms->size(); i++) {
	psmPtr dbgPsm = (*m_origPsms)[i];
	cerr<<"  i = "<<i<<", m_scanNum = "<<dbgPsm->m_scanNum
			<<", m_annotation = "<<dbgPsm->m_annotation
			<<", m_dbIndex = "<<dbgPsm->m_dbIndex
			<<", m_spectrum = "<<dbgPsm->m_spectrum
			<<"\n";
}
*/
    m_origPsms->resize(0);
    if(m_params.exists("INSPECT_PSMS")) {
        m_origPsms->loadInspectResultsFile(m_params.getValue("INSPECT_PSMS").c_str(), m_params.getValueBool("SCAN_ZERO_INDEX",1));
    }
	DEBUG_VAR(m_params.getValue("INPUT_SPEC_IDS"));
    if(m_params.exists("INPUT_SPEC_IDS")) {
    	// Load MSGFDB results from svn/sps/qc/test_data_silac/MSGFDB_res_01FDR.tsv
    	DEBUG_MSG("Loading MS-GFDB IDs from " << m_params.getValue("INPUT_SPEC_IDS"));
    	if (!m_origPsms->loadMSGFDBResultsFile(m_params.getValue("INPUT_SPEC_IDS").c_str(), false))
    	{
    		return false;
    	}
    	DEBUG_VAR(m_origPsms->size());

    	// Stuff below is needed only for projects using Adrian's clustering (e.g., triples/SILAC)

    	// Specify path to project directory as INPUT_CLUSTERS_DIR (ie "svn/sps/qc/test_sps_silac")
    	ClusterData clusterInfo;
    	DEBUG_MSG("Loading clusters from " << m_params.getValue("INPUT_CLUSTERS_DIR"));
    	if (!clusterInfo.loadData(m_params.getValue("INPUT_CLUSTERS_DIR")))
    	{
    		ERROR_MSG("Failed to load clusters from " << m_params.getValue("INPUT_CLUSTERS_DIR"));
    		return false;
    	}
    	map<int, list<int> > clusterInfoMap;
    	for (map<unsigned, list<pair<unsigned, unsigned> > >::iterator clustIt =
    			clusterInfo.data.begin(); clustIt != clusterInfo.data.end(); clustIt++)
    	{
    		list<int> scans;

    		// scans are always 1-based, these indices are 0-based
    		for (list<pair<unsigned, unsigned> >::iterator scanIt =
    				clustIt->second.begin(); scanIt != clustIt->second.end(); scanIt++)
    		{
    			scans.push_back(((int)scanIt->second) + 1);
    		}
    		clusterInfoMap[((int)clustIt->first) + 1] = scans;
    	}
    	// this method will assign each cluster the "consensus" PSM if one of its un-clustered spectra was identified
    	m_origPsms->cluster(clusterInfoMap);

    	stringstream str;
    	str<<m_params.getValue("INPUT_SPEC_IDS")<<"2";
    	m_origPsms->saveToFile(str.str().c_str(), true);
//    	exit(0);

    	DEBUG_VAR(m_origPsms->size());
    }
/*
cerr<<"(re-)Loaded INSPECT_PSMS:\n";
for(unsigned int i=0; i<m_origPsms->size(); i++) {
	psmPtr dbgPsm = (*m_origPsms)[i];
	cerr<<"  i = "<<i<<", m_scanNum = "<<dbgPsm->m_scanNum
			<<", m_annotation = "<<dbgPsm->m_annotation
			<<", m_dbIndex = "<<dbgPsm->m_dbIndex
			<<", m_spectrum = "<<dbgPsm->m_spectrum
			<<"\n";
}
*/
//	DEBUG_TRACE;
//	m_origPsms->addSpectra(m_starSpectra,true);
    DEBUG_TRACE;

    // ---------------------------------------------
    //   SVM calculation
    // ---------------------------------------------
    PeptideSpectrumMatchSet psmsSVM;
    DEBUG_TRACE;
    ExecSvmStatistics svmModule(m_params, m_origPsms, m_psms);
    DEBUG_TRACE;
    string error;
    if (!svmModule.validateStatisticsParams(error))
    {
      DEBUG_MSG(error);
      return false;
    }

    DEBUG_TRACE;
/*  Code below is crashing with decoys for test_sps_silac_nofilt
     if (!svmModule.invokeWithStatistics(m_specsMsSpectra,
                                        m_specsScoredSpectra,
                                        m_starSpectra,
                                        m_model))
    {
      DEBUG_MSG("Unable to generate svm scores!");
      return false;
    }
*/
    // Just copy input PSMs directly to PSMs used for propagation
    (*m_psms) = (*m_origPsms);

    DEBUG_TRACE;
/*
cerr << " -- >> SVM scores 1:\n";
for (unsigned int i = 0; i < m_psms->size(); i++)
{
	if ((*m_psms)[i]->m_spectrum)
	{
		cerr << " -- -- psm " << i << ", spectrum parent mass = "
				<< (*m_psms)[i]->m_spectrum->parentMass << ", spectrum scan = "
				<< (*m_psms)[i]->m_spectrum->scan << ", seq = "
				<< (*m_psms)[i]->m_annotation << ", score = "
				<< (*m_psms)[i]->m_score << " (decoy = " << (*m_psms)[i]->m_isDecoy <<")";
		if (not (*m_psms)[i]->m_spectrum->psmList.empty()) {
			cerr <<", m_spectrum score = "<<(*m_psms)[i]->m_spectrum->psmList.front()->m_score;
		}
		cerr << ")\n";
	}
}
cerr.flush();
*/
    cerr << "FIX NEEDED - set matched protein indices when loading INSPECT_PSMS into PeptideSpectrumMatchSet\n";
    ERROR_MSG("FIX NEEDED - set matched protein indices when loading INSPECT_PSMS into PeptideSpectrumMatchSet\n");
    for (unsigned int i = 0; i < m_psms->size(); i++)
      m_psms->m_psmSet[i]->m_dbIndex = 0;

    DEBUG_TRACE;

    //Generate output Spectrum object.
    m_psmSpectra->resize(m_starSpectra->size());
    Spectrum prmSpectrum;
    float prmMixRatio = 0.5;
    for (unsigned int idxSpec = 0; idxSpec < m_starSpectra->size(); idxSpec++)
    {
    	(*m_psmSpectra)[idxSpec] = (*m_starSpectra)[idxSpec];
    	(*m_psmSpectra)[idxSpec].psmList.clear();

    	// Add attenuated PRM scores to assist with site localization during propagation
    	prmSpectrum = (*m_specsScoredSpectra)[idxSpec];
    	for(unsigned int idxPeak=0; idxPeak<prmSpectrum.size(); idxPeak++)
    		prmSpectrum[idxPeak][1] *= prmMixRatio;
    	(*m_psmSpectra)[idxSpec].mergePeakList(prmSpectrum, (Spectrum *) 0);

    	// GF-specific unit-rounding of peak scores and masses
    	for(unsigned int i=0; i<(*m_psmSpectra)[idxSpec].size(); i++) {  // GF uses rounded masses/scores
    		(*m_psmSpectra)[idxSpec][i][0] = round(0.9995*(*m_psmSpectra)[idxSpec][i][0]);
    		(*m_psmSpectra)[idxSpec][i][1] = round((*m_psmSpectra)[idxSpec][i][1]);
    	}
    	(*m_psmSpectra)[idxSpec].parentMass = round(0.9995*(*m_psmSpectra)[idxSpec].parentMass);
    	(*m_psmSpectra)[idxSpec].setResolution(1,true);  // Merge peaks with the same rounded mass
    }

//    m_psms->addSpectra(m_psmSpectra, true);

    DEBUG_TRACE;


    // Temporary hard-coded settings for Spectrum::getGFProbability() while
    //   waiting for code to get these from inputParams
    AAJumps aminoacids(1); // Used to calculate spectral probabilities
    for(unsigned int i=0; i<aminoacids.size(); i++) {
    	aminoacids.masses[i] = round(aminoacids.masses[i]);
    }
    aminoacids.forceTolerance(0,1,true); // Enforce 1 Da resolution

    DEBUG_TRACE;
    if(useSVM) {
        if (not annotateNodes(m_db, m_psms, m_psmSpectra, &svmModule, true, aminoacids))
        {
          DEBUG_MSG("Error executing ExecMainSpecnets::annotateNodes()");
          return false;
        }
    } else {
        if (not annotateNodes(m_db, m_psms, m_psmSpectra, (ExecSvmStatistics *) 0, true, aminoacids))
        {
          DEBUG_MSG("Error executing ExecMainSpecnets::annotateNodes()");
          return false;
        }

        // Replace match scores with spectral probabilities

        // Temporary hard-coded settings for Spectrum::getGFProbability() while
        //   waiting for code to get these from inputParams
/*        AAJumps aminoacids(1); // Used to calculate spectral probabilities
        for(unsigned int i=0; i<aminoacids.size(); i++) {
        	aminoacids.masses[i] = round(aminoacids.masses[i]);
        }
        aminoacids.forceTolerance(0,1,true); // Enforce 1 Da resolution
*/
        vector<pair<unsigned int,bool> > ntermMods;
        ntermMods.resize(2);
        ntermMods[0].first = 42;   ntermMods[0].second = false;  // N-term acetylation
        ntermMods[1].first = 111;  ntermMods[1].second = true;   // N-term PyroGlu
        vector<unsigned int> mods;
        mods.resize(2);
        mods[0] = 136;  // Lys+8
        mods[1] = 166;  // R+10
        // mods[0] = 147;  // Oxidized Met
        // mods[1] = 115;  // Deamidation on N
        // mods[2] = 129;  // Deamidation on Q
        for (unsigned int idxSpec = 0; idxSpec < m_psmSpectra->size(); idxSpec++)
        {
        	if(not (*m_psmSpectra)[idxSpec].psmList.empty()) {
            	float matchScore;
            	bool isReversed;
            	Spectrum tmpSpec = (*m_starSpectra)[idxSpec];

            	// Add attenuated PRM scores to assist with site localization during propagation
            	prmSpectrum = (*m_specsScoredSpectra)[idxSpec];
            	for(unsigned int idxPeak=0; idxPeak<prmSpectrum.size(); idxPeak++)
            		prmSpectrum[idxPeak][1] *= prmMixRatio;
            	tmpSpec.mergePeakList(prmSpectrum, (Spectrum *) 0);

            	for(unsigned int i=0; i<tmpSpec.size(); i++) {  // GF uses rounded masses/scores
            		tmpSpec[i][0] = round(0.9995*tmpSpec[i][0]);
            		tmpSpec[i][1] = round(tmpSpec[i][1]);
            	}
            	tmpSpec.filterLowMassPeaks(50);
            	tmpSpec.filterHighMassPeaks(tmpSpec.parentMass-AAJumps::massMH-50);
            	tmpSpec.parentMass = round(0.9995*tmpSpec.parentMass);
            	tmpSpec.setResolution(1,true);  // Merge peaks with the same rounded mass

            	psmPtr psm = (*m_psmSpectra)[idxSpec].psmList.front();
            	matchScore = MatchSpecToPeptide(tmpSpec,
            			                        psm->m_annotation.c_str(),
        			                            peakTol, 0, false, &isReversed, &aminoacids);
                if(isReversed) {
                	tmpSpec.reverse(0.0);
                	tmpSpec.filterLowMassPeaks(50);
                	tmpSpec.filterHighMassPeaks(tmpSpec.parentMass-AAJumps::massMH-50);
                	tmpSpec.setResolution(1,true);  // Merge peaks with the same rounded mass

                	// matchScore needs to be recalculated because some PRMs may fall
                	//  out of range after reversing the spectrum (e.g., if at mass 57 before reversing)
                	matchScore = MatchSpecToPeptide(tmpSpec,
                			                        psm->m_annotation.c_str(),
            			                            peakTol, 0, false, &isReversed, &aminoacids);
                }

            	vector<float> pepMasses;
            	aminoacids.getPRMMasses(psm->m_annotation.c_str(), pepMasses);

cerr<<"(*m_psmSpectra)["<<idxSpec<<"].scan = "<<(*m_psmSpectra)[idxSpec].scan<<", annotation = "<<psm->m_annotation
  	<<", score = "<<psm->m_score<<", matchScore = "<<matchScore;
    			psm->m_score = tmpSpec.getGFProbability(matchScore,
    					                                (unsigned int) round(pepMasses[pepMasses.size()-1]),
    					                                false,
    					                                &aminoacids,
    					                                2,
    					                                &mods,
    					                                &ntermMods);
cerr<<", specprob = "<<psm->m_score;
    			psm->m_score = -log10( psm->m_score );
cerr<<", -log10(specprob) = "<<psm->m_score<<"\n";

				DEBUG_VAR(idxSpec);
				DEBUG_VAR(psm->m_score);

/*
       	(*m_starSpectra)[idxSpec].getGFProbability(0,&aminoacids,0,&mods,&ntermMods);

       	for(list<psmPtr>::iterator iter = (*m_psmSpectra)[idxSpec].psmList.begin(); iter != (*m_psmSpectra)[idxSpec].psmList.end() ; iter++) {
       		matchScore = MatchSpecToPeptide((*m_starSpectra)[idxSpec],
       			                            (*iter)->m_annotation.c_str(),
       			                            peakTol, 0, false);
cerr<<"(*m_psmSpectra)["<<idxSpec<<"].scan = "<<(*m_psmSpectra)[idxSpec].scan<<", annotation = "<<(*iter)->m_annotation
  	<<", score = "<<(*iter)->m_score<<", matchScore = "<<matchScore;
//				(*iter)->m_score = -log10( (*m_psmSpectra)[idxSpec].getGFProbability((*iter)->m_score) );
		(*iter)->m_score = -log10( (*m_psmSpectra)[idxSpec].getGFProbability(matchScore) );
cerr<<", -log10(specprob) = "<<(*iter)->m_score<<"\n";
       	}
*/
        	}
        }

        m_psms->getPSMSet(m_psmSpectra);

        DEBUG_TRACE;
        list<psmPtr>::iterator psmIter;
        for(unsigned int specIdx=0; specIdx<m_psmSpectra->size(); specIdx++) {
        	psmIter = (*m_psmSpectra)[specIdx].psmList.begin();
        	while(psmIter!=(*m_psmSpectra)[specIdx].psmList.end()) {
        		if((*psmIter)->m_score<8)
        			psmIter = (*m_psmSpectra)[specIdx].psmList.erase(psmIter);
        		else psmIter++;
        	}
        }
        PeptideSpectrumMatchSet tmpPSMs;
        DEBUG_TRACE;
        tmpPSMs.getPSMSet(m_psmSpectra);
        DEBUG_MSG("Opening output PTM frequency matrix file: " << m_params.getValue("OUTPUT_SNETS_PTM_MATRIX"));
        string fname("min_score_8_");
        fname += m_params.getValue("OUTPUT_SNETS_PTM_MATRIX");
        m_psms->saveModMatrix(fname.c_str());

        DEBUG_TRACE;

    }
    DEBUG_TRACE;


    /*
     // Load + propagate decoy annotations
     PeptideSpectrumMatchSet psmsDecoy;
     DEBUG_TRACE;
     m_psms->loadInspectResultsFile(m_params.getValue("INSPECT_PSMS_DECOY").c_str(), true);
     DEBUG_TRACE;
     if(not annotateNodes(m_db, &psmsDecoy, true)) {
     DEBUG_MSG("Error executing ExecMainSpecnets::annotateNodes()");
     return false;
     }
     DEBUG_TRACE;
     */

    return true;

#ifdef JUNE_BUG_FIX
    // ---------------------------------------------
    //   SVM calculation
    // ---------------------------------------------
    PeptideSpectrumMatchSet psmsSVM;
    DEBUG_TRACE;
    ExecSvmStatistics svmModule(m_params,
        m_psms,
        & psmsSVM);
    DEBUG_TRACE;
    string error;
    if (!svmModule.validateStatisticsParams(error))
    {
      DEBUG_MSG(error);
      return false;
    }
    DEBUG_TRACE;

    cerr<<" -- >> SVM scores 1:\n";
    for(unsigned int i=0; i<m_psms->size(); i++)
    {
      cerr<<" -- -- psm "<<i<<", spectrum parent mass = "<<(*m_psms)[i]->m_spectrum->parentMass
      <<", spectrum scan = "<<(*m_psms)[i]->m_spectrum->scan
      <<", seq = "<<(*m_psms)[i]->m_annotation<<", score = "<<(*m_psms)[i]->m_score
      <<" (decoy = "<<(*m_psms)[i]->m_isDecoy<<")\n";
    }
    cerr.flush();

    DEBUG_TRACE;
    if (!svmModule.invokeWithStatistics(m_specsMsSpectra,
            m_specsScoredSpectra,
            m_starSpectra,
            m_model))
    {
      DEBUG_MSG("Unable to generate svm scores!");
      return false;
    }
    DEBUG_TRACE;

    cerr<<" -- >> SVM scores 2:\n";
    for(unsigned int i=0; i<m_psms->size(); i++)
    {
      cerr<<" -- -- psm "<<i<<", spectrum parent mass = "<<(*m_psms)[i]->m_spectrum->parentMass
      <<", spectrum scan = "<<(*m_psms)[i]->m_spectrum->scan
      <<", seq = "<<(*m_psms)[i]->m_annotation<<", score = "<<(*m_psms)[i]->m_score<<", score = "<<psmsSVM[i]->m_score
      <<" (decoy = "<<(*m_psms)[i]->m_isDecoy<<")\n";
    }
    cerr.flush();

    DEBUG_VAR(m_specsMsSpectra->size());
    DEBUG_VAR(m_specsScoredSpectra->size());
    DEBUG_VAR(m_starSpectra->size());
    DEBUG_VAR(m_psms->size());

    Spectrum ms = (*m_specsMsSpectra)[0];
    Spectrum scored = (*m_specsScoredSpectra)[0];
    Spectrum star = (*m_starSpectra)[0];

    PeptideSpectrumMatch currPsm = *((*m_psms)[0]);

    //NOTE This is an example to show how getSvmScore works. This isn't necessary in ExecMainSpecnets
    //DELETE ME
    if (!svmModule.getSvmScore(&ms,
            &scored,
            &star,
            currPsm,
            m_model))
    {
      DEBUG_MSG("Unable to generate svm scores!");
      return false;
    }
    DEBUG_VAR(currPsm.m_score);
    //END DELETE ME

    cerr<<" -- >> SVM scores 1:\n";
    for(unsigned int i=0; i<m_psms->size(); i++)
    {
      cerr<<" -- -- psm "<<i<<", spectrum parent mass = "<<(*m_psms)[i]->m_spectrum->parentMass
      <<", spectrum scan = "<<(*m_psms)[i]->m_spectrum->scan
      <<", seq = "<<(*m_psms)[i]->m_annotation<<", score = "<<(*m_psms)[i]->m_score
      <<" (decoy = "<<(*m_psms)[i]->m_isDecoy<<")\n";
    }
    cerr.flush();

    DEBUG_TRACE;

    ParameterList fdrParams;
    fdrParams.addIfExists(m_params, "TOLERANCE_PEAK");
    fdrParams.addIfExists(m_params, "INPUT_RESULTS_TYPE");
    //    m_psms->clear();  // Reuse as container for FDR-filtered PSMs
    m_psms->m_psmSet.resize(0);
    DEBUG_TRACE;
    ExecFdrPeptide fdrModule(fdrParams, &psmsSVM, m_psms);
    DEBUG_TRACE;
    if (!fdrModule.invoke())
    {
      DEBUG_MSG("Unable to generate fdr results!");
      return false;
    }
    DEBUG_TRACE;

    // ---------------------------------------------
    //   Reset midx/mp for annotations rejected by FDR
    // ---------------------------------------------
    for(idxSpec=0; idxSpec<m_psmSpectra->size(); idxSpec++)
    (*m_psmSpectra)[idxSpec].psmList.clear();
    DEBUG_TRACE;
    m_psms->addSpectra(m_psmSpectra);
    DEBUG_TRACE;

    for(idxSpec=0; idxSpec<m_psmSpectra->size(); idxSpec++)
    {
      if((*m_psmSpectra)[idxSpec].psmList.empty())
      {
        cerr<<"(after FDR) no annotation for spectrum "<<idxSpec<<"\n";
        (*m_psms_midx)[idxSpec].resize(0);
        (*m_psms_mp)[idxSpec][0]=-1;
        (*m_psms_mp)[idxSpec][1]=0;
        (*m_psms_mp)[idxSpec][2]=0;
      }
    }
    DEBUG_TRACE;

    m_psms_midx->SaveSpecSet_pklbin("specnets/snets_midx.pklbin");
    Save_binArray("specnets/snets_mp.bin", *m_psms_mp);
#endif

    /*
    // ---------------------------------------------
    //   ExecHomologyAssembly
    // ---------------------------------------------
    DEBUG_TRACE;
    ParameterList homologyAssemblyParams(m_params);
    homologyAssemblyParams.setValue("INPUT_SPECNETS_MATCHED","specnets/snets_midx.pklbin");
    homologyAssemblyParams.setValue("INPUT_SPECNETS_ALIGNS","specnets/snets_pairs.bin");
    homologyAssemblyParams.setValue("INPUT_SPECNETS_ALIGNS_ALL","aligns/pairs_stars.bin");

    homologyAssemblyParams.setValue("SPEC_TYPE_MSMS", "0");
    homologyAssemblyParams.setValue("GRAPH_TYPE", "2");
    homologyAssemblyParams.setValue("MIN_CONTIG_SET", "1");
    homologyAssemblyParams.setValue("EDGE_SCORE_TYPE", "1");
    homologyAssemblyParams.setValue("OUTPUT_SPECS",
                                    "homology/homglue_matches.pklbin");
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MIDX",
                                    "homology/homglue_ref_midx.pklbin");
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MP",
                                    "homology/homglue_ref_mp.bin");
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MIDX",
                                    "homology/homglue_matches_midx.pklbin");
    homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MP",
                                    "homology/homglue_matches_mp.bin");
    homologyAssemblyParams.setValue("OUTPUT_CSV", "homology/specnets.csv");
    homologyAssemblyParams.setValue("OUTPUT_CSPS_ABRUIJN",
                                    "homology/specnets_abruijn.bin");
    homologyAssemblyParams.addIfExists(m_params, "RESOLUTION");
    homologyAssemblyParams.addIfExists(m_params, "MAX_MOD_MASS");
    homologyAssemblyParams.addIfExists(m_params, "TOLERANCE_PEAK");
    homologyAssemblyParams.addIfExists(m_params, "TOLERANCE_PM");
    homologyAssemblyParams.addIfExists(m_params, "MAX_AA_JUMP");
    homologyAssemblyParams.addIfExists(m_params, "PENALTY_PTM");
    homologyAssemblyParams.addIfExists(m_params, "PENALTY_SAME_VERTEX");
    homologyAssemblyParams.writeToFile("debug_homology.params");

    ExecHomologyAssembly moduleMapSpecnets(homologyAssemblyParams,
                                           m_psmSpectra,
                                           m_db,
                                           //    		m_psms_midx,
                                           //    		m_psms_mp,
                                           m_snets_contigs,
                                           m_snets_midx,
                                           m_snets_mp);

    DEBUG_TRACE;
    if (not moduleMapSpecnets.invoke())
      return false;
    DEBUG_TRACE;

    if (not moduleMapSpecnets.saveOutputData())
      return false;
    DEBUG_TRACE;

    // Output indices of spectral networks matched to the database (all at this stage)
    vector<unsigned int> matchedSpecNetsIndices;
    matchedSpecNetsIndices.resize(m_snets_contigs->size());
    for(unsigned int i=0; i<matchedSpecNetsIndices.size(); i++)
    	matchedSpecNetsIndices[i] = i;
    Save_binArray("spectra/contigs_indices.bin", matchedSpecNetsIndices);

    return true;
    */
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::loadInputData(void)
  {
    return false;

/*
	if (ownInput)
    {
      m_starSpectra = new SpecSet;
      m_specsMsSpectra = new SpecSet;
      m_specsScoredSpectra = new SpecSet;
      m_pairs = new SpectrumPairSet;
      m_inspectResults = new PeptideSpectrumMatchSet;
      m_specnetsResults = new PeptideSpectrumMatchSet;
      m_model = new MS2ScoringModel;
      m_networkPairs = new SpectrumPairSet;
    }

    if (m_params.exists("STARS_INPUT_PKLBIN"))
    {
      if (m_starSpectra->LoadSpecSet_pklbin(m_params.getValue("STARS_INPUT_PKLBIN").c_str())
          <= 0)
      {
        DEBUG_MSG("Could not load " << m_params.getValue("STARS_INPUT_PKLBIN"));
        return false;
      }

     }

     DEBUG_MSG("STARS LOADED");

    if (m_params.exists("SPECS_SCORED_INPUT_PKLBIN"))
    {
      if (m_specsScoredSpectra->LoadSpecSet_pklbin(m_params.getValue("SPECS_SCORED_INPUT_PKLBIN").c_str())
          <= 0)
      {
        DEBUG_MSG("Could not load " << m_params.getValue("SPECS_SCORED_INPUT_PKLBIN"));
        return false;
      }
    }

     DEBUG_MSG("SPECS SCORED LOADED");

    if (m_params.exists("SPECS_MS_INPUT_PKLBIN"))
    {
      if (m_specsMsSpectra->LoadSpecSet_pklbin(m_params.getValue("SPECS_MS_INPUT_PKLBIN").c_str())
          <= 0)
      {
        DEBUG_MSG("Could not load " << m_params.getValue("SPECS_MS_INPUT_PKLBIN"));
        return false;
      }
    }

     DEBUG_MSG("MS LOADED");

    if (m_params.exists("MS2_SCORING_MODEL"))
    {
      if (!m_model->LoadModel(m_params.getValue("MS2_SCORING_MODEL").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }
    DEBUG_MSG("SCORING MODEL LOADED");

     if (m_params.exists("INPUT_SPECNETS_RESULTS"))
     {
     if (m_params.getValue("INPUT_RESULTS_TYPE").compare("inspect") == 0)
     {
     if (!m_specnetsResults->loadInspectResultsFile(m_params.getValue("INPUT_SPECNETS_RESULTS").c_str(),
     m_params.getValueBool("SCAN_ZERO_INDEX")))
     {
     DEBUG_MSG("Could not load" << m_params.getValue("INPUT_SPECNETS_RESULTS"));
     return false;
     }
     }
     else
     {
     if (!m_specnetsResults->loadSpecnetsResultsFile(m_params.getValue("INPUT_SPECNETS_RESULTS").c_str(),
     m_params.getValueBool("SCAN_ZERO_INDEX")))
     {
     DEBUG_MSG("Could not load" << m_params.getValue("INPUT_SPECNETS_RESULTS"));
     return false;
     }
     }
     }

     DEBUG_MSG("SPECNETS RESULTS LOADED");

     DEBUG_VAR(m_specnetsResults->size());
     */
    /*
     if (ownInput) {
     if (!m_spectra)
     m_spectra = new SpecSet;
     if (!m_pairs)
     m_pairs = new SpectrumPairSet;
     }
     m_spectra->resize(0);
     m_pairs->resize(0);

     if (ownOutput) {
     if (!m_projSpectra)
     m_projSpectra = new SpecSet;
     if (!m_projPeaks)
     m_projPeaks = new SpecSet;
     if (!m_projPairs)
     m_projPairs = new SpectrumPairSet;
     if (!m_projInfo)
     m_projInfo = new vector<vector<int> > ;
     }
     m_projSpectra->resize(0);
     m_projPeaks->resize(0);
     m_projPairs->resize(0);
     m_projInfo->resize(0);

     if (m_params.exists("AMINO_ACID_MASSES")) {
     AAJumps tmpJumps(-1);
     tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
     }

     if (!m_params.exists("INPUT_SPECS_PKLBIN")) {
     ERROR_MSG("Parameters are incomplete. INPUT_SPECS_PKLBIN is missing.");
     return false;
     }
     else if (m_spectra->LoadSpecSet_pklbin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
     <= 0 or m_spectra->size() == 0) {
     ERROR_MSG("Error reading input spectra from " <<  m_params.getValue("INPUT_SPECS_PKLBIN"));
     return false;
     }

     if (!m_params.exists("INPUT_ALIGNS")) {
     ERROR_MSG("Parameters are incomplete. INPUT_ALIGNS is missing.");
     return false;
     }
     else if ( !m_pairs->loadFromBinaryFile(m_params.getValue("INPUT_ALIGNS").c_str()) ) {
     ERROR_MSG("Error reading set of pairs from " << m_params.getValue("INPUT_ALIGNS"));
     return false;
     }
     */
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::saveOutputData(void)
  {
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_RAW_RESULTS"))
    {
      DEBUG_MSG("Opening output raw scores file...");
      m_psms->saveToFile(m_params.getValue("OUTPUT_RAW_RESULTS").c_str());
    }

    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_SNETS_PTM_MATRIX"))
    {
      DEBUG_MSG("Opening output PTM frequency matrix file: " << m_params.getValue("OUTPUT_SNETS_PTM_MATRIX"));
      m_psms->saveModMatrix(m_params.getValue("OUTPUT_SNETS_PTM_MATRIX").c_str());
    }
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecMainSpecnets::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecMainSpecnets::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("SPECS_MS_STATISTICS_CONFIG");
    VALIDATE_PARAM_EXIST("SPECS_SCORED_STATISTICS_CONFIG");
    VALIDATE_PARAM_EXIST("STARS_STATISTICS_CONFIG");

    m_isValid = true;
    return true;
  }

} // namespace specnets

