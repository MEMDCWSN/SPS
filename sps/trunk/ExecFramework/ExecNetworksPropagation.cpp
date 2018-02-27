/*
 * ExecNetworksPropagation.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: isna
 */

#include "ExecNetworksPropagation.h"
#include "ClusterSet.h"

namespace specnets {

	static int max_step_propagation = 30;
	// -------------------------------------------------------------------------
	ExecNetworksPropagation::ExecNetworksPropagation(void) :
    	      ExecBase(), m_prmSpectra(0x0), m_spectralPairs(0x0), m_psmResults(0x0),
    	      	  m_seedFDR(0), m_edgeFDR(0), m_propagationResults(0x0), ownInput(true), ownOutput(true)
	{
		m_name = "ExecNetworksPropagation";
		m_type = "ExecNetworksPropagation";
	}

	// -------------------------------------------------------------------------
	ExecNetworksPropagation::ExecNetworksPropagation(const ParameterList & inputParams) :
	    ExecBase(inputParams), m_prmSpectra(0x0), m_spectralPairs(0x0), m_psmResults(0x0),
	    	m_seedFDR(0), m_edgeFDR(0), m_propagationResults(0x0), ownInput(true), ownOutput(true)
	{
		m_name = "ExecNetworksPropagation";
		m_type = "ExecNetworksPropagation";
	}

	// -------------------------------------------------------------------------
	ExecNetworksPropagation::~ExecNetworksPropagation(void) {
		if( ownInput ){
			delete m_prmSpectra;
			delete m_spectralPairs;
			delete m_psmResults;
		}

		if( ownOutput ){
			if( m_propagationResults != 0x0 ){
				for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
					if( (*vi) != 0 ) delete (*vi);
				}
				delete m_propagationResults;
			}
		}

	}

	// -------------------------------------------------------------------------
	ExecBase * ExecNetworksPropagation::clone(const ParameterList & inputParams) const {
		return new ExecNetworksPropagation(inputParams);
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::invoke(void) {

		float minSeedScore 		= m_params.getValueDouble("MIN_SEED_SCORE", 0);
		bool  removeMultiplexed	= m_params.getValueBool("RM_MULTIPLEX_PSM", true);

		float overallFDR 		= m_params.getValueDouble("PROPAGATION_FDR", 0.01);
		float edgeFDR_Damping 	= m_params.getValueDouble("EDGE_FDR_DAMPING", 0.8);
		float seedFDR_Damping 	= m_params.getValueDouble("SEED_FDR_DAMPING", 0.6);

		if( m_prmSpectra->size()*2 < m_spectralPairs->size() ) {
			seedFDR_Damping = 0.3;
		}//*/

		float edgeFDR		= overallFDR*edgeFDR_Damping;
		float seedFDR 		= overallFDR*seedFDR_Damping;

		max_step_propagation = m_params.getValueDouble("PROPAGATION_MAX_STEP", 5);

		float pmTol 		= m_params.getValueDouble("TOLERANCE_PM", 1);
		float fragTol 		= m_params.getValueDouble("TOLERANCE_PEAK", 0.5);

		DEBUG_VAR(pmTol);
		DEBUG_VAR(fragTol);
		DEBUG_VAR(seedFDR);
		DEBUG_VAR(edgeFDR);
		DEBUG_VAR(minSeedScore);

		//filter seeds
		PeptideSpectrumMatchSet selectedSeeds;
		m_seedFDR = selectSeeds(selectedSeeds, seedFDR, minSeedScore, removeMultiplexed);

		//filter edges
		if( m_params.exists("ANNOTATED_PAIRS") )
			m_spectralPairs->filter_by_edge_fdr33(*m_prmSpectra, selectedSeeds, edgeFDR, pmTol, fragTol, m_params.getValue("ANNOTATED_PAIRS"), true);
		else m_spectralPairs->filter_by_edge_fdr22(*m_prmSpectra, selectedSeeds, edgeFDR, pmTol, fragTol, true);

		m_edgeFDR = m_spectralPairs->getFDR();
		DEBUG_MSG( "#Selected EDGES: " << m_spectralPairs->size() << ", EDGE FDR: " << m_edgeFDR );

		//propagate seeds
		m_propagationResults = new vector<NetworkedPSM*>(m_prmSpectra->size());
		propagateAnnotations(selectedSeeds, fragTol);

		//control FDR
		selectIdentificationsByFDR(overallFDR);

		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::loadInputData(void) {

		if( m_prmSpectra == 0x0 ) {
		  ownInput = true;
		  m_prmSpectra 		= new SpecSet;
		  m_spectralPairs 	= new SpectrumPairSet;
		  m_psmResults 		= new PeptideSpectrumMatchSet;
		}

		if (m_params.exists("INPUT_PRM_SPECS")){
			m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_SPECS").c_str());
			DEBUG_MSG( "#PRM Spectra: " << m_prmSpectra->size() );
		}

		if (m_params.exists("INPUT_EDGES")){
			m_spectralPairs->loadFromBinaryFile(m_params.getValue("INPUT_EDGES").c_str());
			DEBUG_MSG( "#INPUT EDGES: " << m_spectralPairs->size() );
		}

		if (m_params.exists("INPUT_SEEDS")){
			m_psmResults->loadMSGFDBResultsFile(m_params.getValue("INPUT_SEEDS").c_str());
			DEBUG_MSG( "#INPUT PSMS: " << m_psmResults->size() );
		}
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::saveOutputData(void) {

		map<int, string> proteins;
		for (vector<psmPtr>::iterator psmIterator = m_psmResults->m_psmSet.begin();
						psmIterator != m_psmResults->m_psmSet.end(); psmIterator++) {
			proteins[(*psmIterator)->m_scanNum-1]= (*psmIterator)->m_protein;
		}

		if( m_params.exists("SPEC_ID_OUTPUT") ){//raw spectrum level identification

			ClusterSet clust_info;
			clust_info.loadBinaryFile(m_params.getValue("SPECS_CLUST", ""));
			string specdir = m_params.getValue("SPECS_DIR", "");
			string pklinfo = m_params.getValue("SPECS_PKLBIN", "");

			DEBUG_VAR(clust_info.fileNames.size());
			DEBUG_VAR(pklinfo);
			DEBUG_VAR(specdir);

			vector<string> pklbinList;
			char buf[512];
			FILE *in = fopen(pklinfo.c_str(), "r");
			while( fgets(buf, 512, in) != NULL ) {
				for(int i=strlen(buf)-1; i>0; i--){
					if( buf[i]=='\r' || buf[i]=='\n' ) buf[i]=0;
					if( isalpha(buf[i]) ) break;
				}
				char pklbin[512];
				sprintf(pklbin, "%s%s",  specdir.c_str(), strrchr(buf, '/') );
				pklbinList.push_back( string(pklbin) );
			}
			fclose(in);

			vector<vector<int> >   specCS(pklbinList.size());
			vector<vector<float> > specPM(pklbinList.size());

			for(int pk=0; pk<pklbinList.size(); pk++){
				SpecSet specs;
				specs.loadPklBin(pklbinList[pk].c_str());
				if( specs.size() != 0 ){
					for(int i=0; i<specs.size(); i++){
						specCS[pk].push_back(specs[i].parentCharge);
						specPM[pk].push_back(specs[i].parentMass);
					}
				}
			}

			string outputpath = m_params.getValue("SPEC_ID_OUTPUT", "");
			DEBUG_MSG("The spectrum identification results were written to " + outputpath);

			FILE *txtfp = fopen(outputpath.c_str(), "w");
			fprintf(txtfp, "Cluster\tcCharge\tcMW\t");
			fprintf(txtfp, "Filename\tSpecIndex\tScan\tCharge\tMW\tAnnotation\tIDLayer\tProtein\n");

			int index = 0;
			for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {

				if( (*vi) != 0 ) {//has ID

					Cluster cluster = clust_info[index];
					if( index != cluster.m_index ) {
						DEBUG_MSG("Cluster index ERROR!");
					}

					string prot = proteins[(*vi)->getSeedIndex()];
					for(int i=0; i<cluster.size(); i++){

						fprintf(txtfp, "%d\t%d\t%.4f\t",
								(*m_prmSpectra)[index].scan, (*m_prmSpectra)[index].parentCharge, (*m_prmSpectra)[index].parentMass );

						fprintf(txtfp, "%s\t%d\t%d\t%d\t%.4f\t", cluster[i].m_filename.c_str(), cluster[i].m_index+1, cluster[i].m_scan,
								specCS[cluster[i].m_fileIndex][cluster[i].m_index],  specPM[cluster[i].m_fileIndex][cluster[i].m_index]);

						fprintf(txtfp, "%s\t%d\t%s\t", (*vi)->getAnnotation().c_str(), (*vi)->getLayer(), prot.c_str());
						fprintf(txtfp, "\n");
					}
				}
				index++;
			}
			fclose(txtfp);
		}

		if( m_params.exists("PROPAGATION_OUTPUT") ){//cluster level identification


			vector<int> specCompID;
			vector<int> compSize;
			m_spectralPairs->getSpecComponentID(*m_prmSpectra, specCompID, compSize);

			vector<int> idedInNet(compSize.size());
			vector<int> propaGroupSize(m_propagationResults->size());

			int index = 0;
			for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
				if( (*vi) != 0 ) {
					propaGroupSize[(*vi)->getSeedIndex()]++;
					if( specCompID[index] != -1 ){
						idedInNet[specCompID[index]]++;
					}
				}
				index++;
			}

			string outputpath = m_params.getValue("PROPAGATION_OUTPUT", "");
			DEBUG_MSG("The propagation results were written to " + outputpath);

			FILE *txtfp = fopen(outputpath.c_str(), "w");
			fprintf(txtfp, "SpecIndex\tScanNo\tCharge\tMW\t");
			fprintf(txtfp, "Layer\tMatchScore\tAnnotation\tNumMods\tSeedSpecIndex\tParentSpecIndex\tMassDiffFromParent\tPropagatedAlignScore\tProtein\t");
			fprintf(txtfp, "GroupSize\tNetworkID\tNetworkSize\t#NetID\n");

			index = 0;
			int sinNet = compSize.size()+1;
			for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
				if( (*vi) != 0 ) {
					fprintf(txtfp, "%d\t", index+1);
					fprintf(txtfp, "%d\t%d\t%.4f\t", (*m_prmSpectra)[index].scan, (*m_prmSpectra)[index].parentCharge, (*m_prmSpectra)[index].parentMass );
					fprintf(txtfp, "%s\t%s\t%d\t", (*vi)->toString().c_str(), proteins[(*vi)->getSeedIndex()].c_str(), propaGroupSize[index]-1 );

					int compID= specCompID[index];
					if( compID == - 1 ) fprintf(txtfp, "%d\t%d\t%d", sinNet++, 1, 1);
					else fprintf(txtfp, "%d\t%d\t%d", compID+1, compSize[compID], idedInNet[compID] );
					fprintf(txtfp, "\n");
				}
				index++;
			}
			fclose(txtfp);
		}

		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::saveInputData(std::vector<std::string> & filenames) {
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::loadOutputData(void) {
		return true;
	}

	// -------------------------------------------------------------------------
	vector<ExecBase*> const & ExecNetworksPropagation::split(int numSplit) {
		return m_subModules;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::merge(void) {
		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecNetworksPropagation::validateParams(std::string & error) {

	  m_isValid = false;
	  VALIDATE_PARAM_EXIST("INPUT_PRM_SPECS");
	  VALIDATE_PARAM_EXIST("INPUT_EDGES");
	  VALIDATE_PARAM_EXIST("INPUT_SEEDS");

	  m_isValid = true;
	  return true;
	}

	// -------------------------------------------------------------------------
	float ExecNetworksPropagation::selectSeeds(PeptideSpectrumMatchSet & selectedSeeds,
			                                   float seedFDR, float minScoreforSeeds, bool removeMultiplexed) {

		set<int> presentScans, redundantScans;
		if( removeMultiplexed ){
			for (vector<psmPtr>::iterator psmIterator = m_psmResults->m_psmSet.begin();
					psmIterator != m_psmResults->m_psmSet.end(); psmIterator++) {

				if( presentScans.find((*psmIterator)->m_scanNum) != presentScans.end() ) {
					redundantScans.insert((*psmIterator)->m_scanNum);
				}
				else presentScans.insert((*psmIterator)->m_scanNum);
			}
			presentScans.clear();
		}
		DEBUG_VAR(redundantScans.size());

		selectedSeeds.resize(0);
		set<string> uniquePept;
		float actualFDR = 0;
		for (vector<psmPtr>::iterator psmIterator = m_psmResults->m_psmSet.begin();
				psmIterator != m_psmResults->m_psmSet.end(); psmIterator++) {

			//removing possibly co-eluted spectra
		//	if( redundantScans.find((*psmIterator)->m_scanNum) != redundantScans.end() ) continue;

			//removing different CS spectra
			if( (*psmIterator)->m_charge != (*m_prmSpectra)[(*psmIterator)->m_scanNum-1].parentCharge ) continue;

			//removing badly matching spectra
			if( seedFDR < (*psmIterator)->m_fdr ||
					(*psmIterator)->m_compScore < minScoreforSeeds || -log10((*psmIterator)->m_score) < minScoreforSeeds ) continue;

			if( presentScans.find((*psmIterator)->m_scanNum) == presentScans.end() )
				presentScans.insert((*psmIterator)->m_scanNum);
			else continue;//one PSM per scan

			if( actualFDR < (*psmIterator)->m_fdr ) actualFDR = (*psmIterator)->m_fdr;
			selectedSeeds.m_psmSet.push_back(*psmIterator);
			uniquePept.insert((*psmIterator)->m_origAnnotation);
		}
		DEBUG_MSG( "#Selected SEEDS: " << selectedSeeds.size() << ", SEED FDR: " << actualFDR << ", #Unique Pepts: " << uniquePept.size() );

		return actualFDR;
	}

	void ExecNetworksPropagation::propagateAnnotations(PeptideSpectrumMatchSet & seeds, float fragTol, float minEdgeScore)
	{

		int dataSize = m_prmSpectra->size();
		vector<int> 			identified(dataSize); //id level by propagation, seed level=1
		vector<float> 			openScore(dataSize);

		for (vector<psmPtr>::iterator psmIterator = seeds.m_psmSet.begin(); psmIterator != seeds.m_psmSet.end(); psmIterator++) {
			identified[(*psmIterator)->m_scanNum-1]  = 1; //zero-based
			openScore[(*psmIterator)->m_scanNum-1]   = 1000;
			(*m_propagationResults)[(*psmIterator)->m_scanNum-1] = new NetworkedPSM((*psmIterator)->m_origAnnotation, (*psmIterator)->m_scanNum-1);
		}

		vector<vector<pair<int, float> > > * neighbors= new vector<vector<pair<int, float> > >(dataSize);
		m_spectralPairs->sort_descending_by_score();
		for (int i = 0; i < m_spectralPairs->size(); i++) {
		  int spec1 = (*m_spectralPairs)[i].spec1, spec2 = (*m_spectralPairs)[i].spec2;
		  float score = min((*m_spectralPairs)[i].score1, (*m_spectralPairs)[i].score2);
		  if( score < minEdgeScore ) continue;

		  if( identified[spec1] == 0 ) {
			  (*neighbors)[spec1].push_back( pair<int, float>(spec2, score) );
			  if( identified[spec2] == 1 && openScore[spec1] < score ) openScore[spec1] = score;
		  }
		  if( identified[spec2] == 0 ) {
			  (*neighbors)[spec2].push_back( pair<int, float>(spec1, score) );
			  if( identified[spec1] == 1 && openScore[spec2] < score ) openScore[spec2] = score;
		  }
		}
		m_spectralPairs->sort_pairs_by_index();

		for(int i=0; i<dataSize; i++){
		  if( identified[i] != 0 || (*neighbors)[i].size() == 0 ) continue;
		  if( identified[(*neighbors)[i][0].first] != 0 )
			  depthFirstPropagation( *neighbors, (*neighbors)[i][0].first, i, identified, openScore, fragTol );
		}

		set<int> readyToOpen;
		for(int i=0; i<dataSize; i++){
		  if( identified[i] != 0 || (*neighbors)[i].size() == 0 ) continue;
		  if( openScore[i] != 0 ) readyToOpen.insert(i);
		}

		while( 0 < readyToOpen.size() ){

		  float maxScore    = 0;
		  int 	specToOpen  = -1;

		  for( set<int>::iterator it = readyToOpen.begin(); it != readyToOpen.end(); ++it ){
			  if( maxScore < openScore[*it] ){
				  maxScore = openScore[*it];
				  specToOpen = *it;
			  }
			  else if( maxScore == openScore[*it] && *it < specToOpen ){
				  specToOpen = *it;
			  }
		  }
		  if( specToOpen == -1 ) break;

		  int from= -1, to= specToOpen;
		  vector<pair<int, float> > nbors = (*neighbors)[to];
		  for(int i=0; i<nbors.size(); i++){
			  if( identified[nbors[i].first] != 0 ) {
				  from = nbors[i].first;
				  break;
			  }
		  }

		  depthFirstPropagation( *neighbors, from, to, identified, openScore, readyToOpen, fragTol );
		}

		delete neighbors;
	}

	void ExecNetworksPropagation::depthFirstPropagation(vector<vector<pair<int, float> > > & neighbors,
													    int from,
													    int to,
													    vector<int> & identified,
													    vector<float> & openScore,
													    float fragTol)
	{
		identified[to] = identified[from] + 1;
	//	openScore[to] = (openScore[from] < openScore[to])? openScore[from] : openScore[to];
		(*m_propagationResults)[to] = (*m_propagationResults)[from]->propagateToNeighborSpectrum(identified[from],
							from, (*m_prmSpectra)[to].parentMass-(*m_prmSpectra)[from].parentMass, openScore[to], (*m_prmSpectra)[to], fragTol);

		if( identified[to] > max_step_propagation ) {
			for(int i=0; i<neighbors[to].size(); i++){
				pair<int, float> neigh = neighbors[to][i];
				for( vector<pair<int, float> >::iterator it=neighbors[neigh.first].begin();
						it!=neighbors[neigh.first].end(); it++ ){
					if( (*it).first == to ) {
						neighbors[neigh.first].erase(it);
						break;
					}
				}
			}
			return;
		}

		for(int i=0; i<neighbors[to].size(); i++){
			pair<int, float> neigh = neighbors[to][i];
			if( identified[neigh.first] != 0 ) continue;
			if( neighbors[neigh.first][0].first == to ) {
				openScore[neigh.first] = neigh.second;
				depthFirstPropagation( neighbors, to, neigh.first, identified, openScore, fragTol );
			}
			else {
				if( openScore[neigh.first] < neigh.second ) openScore[neigh.first] = neigh.second;
			}
		}
	}

	void ExecNetworksPropagation::depthFirstPropagation(vector<vector<pair<int, float> > > & neighbors,
											  		    int from,
													    int to,
													    vector<int> & identified,
													    vector<float> & openScore,
													    set<int> & readyToOpen,
													    float fragTol)
	{
		identified[to] = identified[from] + 1;
	//	openScore[to] = (openScore[from] < openScore[to])? openScore[from] : openScore[to];
		(*m_propagationResults)[to] = (*m_propagationResults)[from]->propagateToNeighborSpectrum(identified[from],
							from, (*m_prmSpectra)[to].parentMass-(*m_prmSpectra)[from].parentMass, openScore[to], (*m_prmSpectra)[to], fragTol);

		readyToOpen.erase(to);
		if( identified[to] > max_step_propagation ) {
			for(int i=0; i<neighbors[to].size(); i++){
				pair<int, float> neigh = neighbors[to][i];
				for( vector<pair<int, float> >::iterator it=neighbors[neigh.first].begin();
						it!=neighbors[neigh.first].end(); it++ ){
					if( (*it).first == to ) {
						neighbors[neigh.first].erase(it);
						break;
					}
				}
			}
			return;
		}

		for(int i=0; i<neighbors[to].size(); i++){
			pair<int, float> neigh = neighbors[to][i];
			if( identified[neigh.first] != 0 ) continue;
			if( neighbors[neigh.first][0].first == to ) {
				openScore[neigh.first] = neigh.second;
				depthFirstPropagation( neighbors, to, neigh.first, identified, openScore, readyToOpen, fragTol );
			}
			else {
				if( openScore[neigh.first] < neigh.second ) {
					openScore[neigh.first] = neigh.second;
					readyToOpen.insert(neigh.first);
				}
			}
		}
	}

	int ExecNetworksPropagation::selectIdentificationsByFDR(float specifiedFDR)
	{
		vector<int> numOfID(max_step_propagation+1);

		for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
			if( (*vi) != 0 && (*vi)->getLayer() <= max_step_propagation ) {
				numOfID[(*vi)->getLayer()]++;
			}
		}

		int    totalID = numOfID[0]; //#seedes
		double stepFDR = m_seedFDR;
		double putativeFalse = totalID * stepFDR;

		int maxPropaStep = 0, putty = 0;

		for(int i=1; i<numOfID.size(); i++){

			if( numOfID[i] == 0 ) break;

			stepFDR = stepFDR + m_edgeFDR - stepFDR*m_edgeFDR;
			double addedFalse = numOfID[i]*stepFDR;

			if( specifiedFDR < ( (putativeFalse+addedFalse)/(totalID+numOfID[i]) ) ){
				maxPropaStep = i-1;
				putty = (int)((totalID*specifiedFDR-putativeFalse)/(stepFDR-specifiedFDR));
				putativeFalse += putty*stepFDR;
				totalID += putty;
				break;
			}
			putativeFalse += addedFalse;
			totalID += numOfID[i];
			maxPropaStep = i;
		}

		float aggfer = putativeFalse/totalID;
		DEBUG_MSG("Propagation results: " << totalID << " identifications at aggregate FDR " << aggfer << ", max propagation step: " << maxPropaStep);

		float putty_thr = 1000000;
		if( putty != 0 ) {
			vector<float> scores;
			for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
				if( (*vi) != 0 && (*vi)->getLayer() == maxPropaStep+1 ) {
					scores.push_back((*vi)->getAlignmentScore());
				}
			}
			sort(scores.rbegin(), scores.rend());
			putty_thr = scores[putty];
		}

		for(vector<NetworkedPSM*>::iterator vi=m_propagationResults->begin(); vi!=m_propagationResults->end(); ++vi) {
			if( (*vi) != 0 && maxPropaStep < (*vi)->getLayer() ) {
				if( (*vi)->getLayer() == maxPropaStep+1 && (*vi)->getAlignmentScore() > putty_thr ) continue;

				delete (*vi);
				(*vi) = 0;
			}
		}
		return maxPropaStep;
	}

} /* namespace specnets */





