
#include "Logger.h"
#include "ExecQCSpectralPairs.h"

#include "PeptideSpectrumMatchSet.h"
#include "PeptideSpectrumMatch.h"
#include "SpectrumPairSet.h"
#include "projectionutils.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <string.h>
#include <math.h>

namespace specnets
{
	
	bool PSM::compareBySpectrumIncex(const PSM& m1, const PSM& m2) { return m1.index < m2.index; }

	static double proton = 1.007276;
	static vector<PSM> getPSMListfromMSGFDB(const string & fname, double msgfScoreCut, double fdrCut);
	static set<string> getSPSPair(const string & filename, std::map<std::string, SpectrumPair> & pair_data_map);
	static int round(double a) {
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}

	// -------------------------------------------------------------------------
	ExecQCSpectralPairs::ExecQCSpectralPairs(void) {
		m_name = "ExecQCSpectralPairs";
		m_type = "ExecQCSpectralPairs";
	}

	// -------------------------------------------------------------------------
	ExecQCSpectralPairs::ExecQCSpectralPairs(const ParameterList & inputParams) :
	    ExecBase(inputParams)
	{
		m_name = "ExecQCSpectralPairs";
		m_type = "ExecQCSpectralPairs";
	}


	// -------------------------------------------------------------------------
	ExecQCSpectralPairs::~ExecQCSpectralPairs(void) {
	}

	// -------------------------------------------------------------------------
	ExecBase * ExecQCSpectralPairs::clone(const ParameterList & inputParams) const {
		return new ExecQCSpectralPairs(inputParams);
	}
	
	void ExecQCSpectralPairs::generate_precision_recall_curve(  set<string> spsPair, 
                                                                    vector<PSM> pList,
                                                                    bool        c_chargeSeperate, 
                                                                    double      c_minOverlap,
                                                                    double  c_pmTolerance,
                                                                    double  c_fdrCut,
                                                                    double  c_msgfScoreCut,
                                                                    std::map<std::string, SpectrumPair>  pair_data_map){
            std::cout<<"Cosine\tPrecision\tRecall"<<std::endl;
            for(float i = 0.f; i < 1.f; i+= 0.10){
                set<string> filtered_spsPair;
                
                for (std::set<string>::iterator it=spsPair.begin(); it!=spsPair.end(); ++it){
                    //std::cout << ' ' << *it <<std::endl;
                    float score1 = pair_data_map[*it].score1;
                    if(score1 > i){
                        filtered_spsPair.insert(*it);
                    }
                    //std::cout<<score1<<std::endl;
                }
                std::cout<<i<<"\t";
                generate_precision_recall(filtered_spsPair, pList, c_chargeSeperate, c_minOverlap, c_pmTolerance, c_fdrCut, c_msgfScoreCut);
            }
            
            
            
            
            
        }
	
	void ExecQCSpectralPairs::generate_precision_recall(set<string> spsPair, 
                                                            vector<PSM> pList,
                                                            bool        c_chargeSeperate, 
                                                            double      c_minOverlap,
                                                            double  c_pmTolerance,
                                                            double  c_fdrCut,
                                                            double  c_msgfScoreCut
                                                            ){

                 
            
            int gold_dist0[2] = {0, 0}; //1: correct, 0: incorrect
            int gold_dist1[2] = {0, 0};
            int gold_part[2][11];

            int sps_dist0[2] = {0, 0}; //1: correct, 0: incorrect
            int sps_dist1[2] = {0, 0};
            int sps_part[2][11];
            for(int i=0; i<11; i++){
                    gold_part[0][i] = 0;
                    gold_part[1][i] = 0;
                    sps_part[0][i] = 0;
                    sps_part[1][i] = 0;
            }
            
            int msgfpair = 0, spspair = 0;
            int isize = pList.size();
            
            
            for(int i=0; i<isize; i++){
                    PSM one = pList[i];
                    for(int k=i+1; k<isize; k++){
                            PSM two = pList[k];

                            bool diffSeq = one.isDiffSeq(two);
                            if( c_chargeSeperate && (one.getChargeState() != two.getChargeState()) ) continue;

                            double align[] = {0, 0};
                            one.getPeptOverlap(two, align);

                            int distance = (int)align[0];
                            double overlap = align[1];

                            if( overlap < c_minOverlap ) continue;

                            char tppair[500];
                            sprintf(tppair, "%d_%d", one.getSpectrumIndex(), two.getSpectrumIndex());
                            int spsdected = ( spsPair.find(string(tppair)) == spsPair.end() )? 0 : 1;

                            bool considered = false, correctly = false, identified = false;
                            if(false){
                            //if( fabs(one.getObservedMW()-two.getObservedMW()) < c_pmTolerance ){
                                    if( overlap == 1 ){
                                            if( distance == 0 ) {
												gold_dist0[1]++;
												correctly = true;
												considered = true;
                                            }
                                            else if( distance == 1 ) {
												gold_dist1[1]++;
												correctly = true;
												considered = true;
											}
                                            else if( diffSeq ) {
												gold_dist0[0]++;
												considered = true;
                                            }

                                            if( spsdected == 1 ) {
												identified = true;
												if( distance == 0 ) sps_dist0[1]++;
												else if( distance == 1 ) sps_dist1[1]++;
												else if( diffSeq ) sps_dist0[0]++;
                                            }
                                    }
                            }
                            else{
                                    if( overlap == 1 ) {
                                            if( distance == 1 ) {
                                                    gold_dist1[1]++;
                                                    correctly = true;
                                                    considered = true;
                                            }
                                            else if( diffSeq ) {
                                                    gold_dist1[0]++;
                                                    considered = true;
                                            }
                                            else if( distance == 0 ) {
												gold_dist0[1]++;
												correctly = true;
												considered = true;
											}

                                            if( spsdected == 1 ) {
                                                    identified = true;
                                                    if( distance == 1 ) sps_dist1[1]++;
                                                    else if( diffSeq ) sps_dist1[0]++;
                                                    else if( distance == 0 ) sps_dist0[1]++;
                                            }
                                    }
                                    else {
                                            int bin = (int)(overlap*10);

                                            if( distance == 0 ) {
                                                    gold_part[1][bin]++;
                                                    correctly = true;
                                                    considered = true;
                                            }
                                            else if( distance == -1 ) {
                                                    gold_part[0][bin]++;
                                                    considered = true;
                                            }

                                            if( spsdected == 1 ){
                                                    identified = true;
                                                    if( distance== 0 ) sps_part[1][bin]++;
                                                    else if( distance == -1 ) sps_part[0][bin]++;
                                            }
                                    }
                            }

                            if( !considered ) continue;

                            if( correctly || identified ) {
                            }

                    }
            }
            
            /*printf( "\nStatistics for SPS pairs --------------------------------------\n" );
            printf( " for dist. 0, sensitivity %d / %d = %.2f\n", sps_dist0[1], gold_dist0[1], ((float)sps_dist0[1]/gold_dist0[1]*100) );
            printf( "                 accuracy %d / %d = %.2f | Incorrect: %d\n", sps_dist0[1], (sps_dist0[0]+sps_dist0[1]),
                            ((float)sps_dist0[1]/(sps_dist0[0]+sps_dist0[1])*100), sps_dist0[0]);

            printf( " for dist. 1, sensitivity %d / %d = %.2f\n", sps_dist1[1], gold_dist1[1], ((float)sps_dist1[1]/gold_dist1[1]*100) );
            printf( "                 accuracy %d / %d = %.2f | Incorrect: %d\n", sps_dist1[1], (sps_dist1[0]+sps_dist1[1]),
                            ((float)sps_dist1[1]/(sps_dist1[0]+sps_dist1[1])*100), sps_dist1[0] );*/
                            
            std::cout<<((float)sps_dist1[1]/(sps_dist1[0]+sps_dist1[1])*100)<<"\t"<<((float)sps_dist1[1]/gold_dist1[1]*100)<<std::endl;
            
        }

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::invoke(void) {

		printf( "ExecQCSpectralPairs -------------------------------------------\n\n" );

		string msgf_file = 	m_params.getValue("MSGFDB_PSM_FILE", "");
		string sps_pair = 	m_params.getValue("SPS_PAIR_FILE", "");
		string out_file = 	m_params.getValue("OUTPUT_FILE", "");
//		string fn_file = 	m_params.getValue("OUTPUT_FILE_FN", "");

		bool 	c_chargeSeperate 	= (m_params.getValueInt("CS_DISTINCT", 1))? true : false;
		double 	c_minOverlap 		= m_params.getValueDouble("MIN_OVERLAP", 0.5);
		double 	c_pmTolerance 		= m_params.getValueDouble("TOLERANCE_PM", 2.5);
		double 	c_fdrCut 			= m_params.getValueDouble("MIN_FDR", 0.01);
		double 	c_msgfScoreCut 		= m_params.getValueDouble("MIN_MSGFSCORE", 10);
                
        bool    c_generatecurve        = (m_params.getValueInt("GENERATECURVE", 1))? true : false;

        PeptideSpectrumMatchSet m_psmResults;
        m_psmResults.loadMSGFDBResultsFile(msgf_file.c_str());
        if( m_psmResults.size() == 0 ) { return false; }
        cout << "Input: " << m_psmResults.size() << " PSMs, ";

        SpectrumPairSet m_spectralPairs;
        m_spectralPairs.loadFromBinaryFile(sps_pair.c_str());
        cout << m_spectralPairs.size() << " SPS Pairs"<< endl << endl;
        if( m_spectralPairs.size() == 0 ) { return false; }

		DEBUG_VAR(c_chargeSeperate);
		DEBUG_VAR(c_minOverlap);

		if(c_generatecurve){
			m_spectralPairs.cosine_precision_recall_curve(m_psmResults,
			        									  c_minOverlap,
			        									  c_chargeSeperate,
			        									  true);
			return true;
		}

        m_spectralPairs.cosine_precision_recall(m_psmResults,
        										out_file,
        										c_minOverlap,
        										c_chargeSeperate,
        										true);

/*
		vector<PSM> pList = getPSMListfromMSGFDB(msgf_file, c_msgfScoreCut, c_fdrCut);
		if( pList.size() == 0 ) { return false; }
		printf("Input: %d PSMs, ", pList.size());
                
        std::map<std::string, SpectrumPair>  pair_data_map;

		set<string> spsPair = getSPSPair(sps_pair, pair_data_map);
		if( spsPair.size() == 0 ) { return false; }

		printf("%d SPS Pairs\n\n", spsPair.size());
                
                if(c_generatecurve){
                    generate_precision_recall_curve(spsPair, 
                                                    pList,
                                                    c_chargeSeperate, 
                                                    c_minOverlap,
                                                    c_pmTolerance,
                                                    c_fdrCut,
                                                    c_msgfScoreCut,
                                                    pair_data_map);
                    return true;
                }

		int gold_dist0[2] = {0, 0}; //1: correct, 0: incorrect
		int gold_dist1[2] = {0, 0};
		int gold_part[2][11];

		int sps_dist0[2] = {0, 0}; //1: correct, 0: incorrect
		int sps_dist1[2] = {0, 0};
		int sps_part[2][11];
		for(int i=0; i<11; i++){
			gold_part[0][i] = 0;
			gold_part[1][i] = 0;
			sps_part[0][i] = 0;
			sps_part[1][i] = 0;
		}

		FILE* out = fopen(out_file.c_str(), "w");
		fprintf(out, "overlap\tdistance\tspec1\tpept1\tcs1\tmsgfscore1\t-log(specprob1)\tprotein1\tspec2\tpept2\tcs2\tmsgfscore2\t-log(specprob2)\tprotein2\tcorrect\tspspair\tcosine\tmzdiff\n" );

		bool fnprint = false;
		FILE* fn_out = fopen(fn_file.c_str(), "w");
		if( fn_out != 0 ){
			fnprint = true;
			printf("generated false-negative pairs\n");
			fprintf(fn_out, "overlap\tspec1\tpept1\tcs1\tmsgfscore1\t-log(specprob1)\tspec2\tpept2\tcs2\tmsgfscore2\t-log(specprob2)\n" );
		}

		int msgfpair = 0, spspair = 0;
		int isize = pList.size();
		for(int i=0; i<isize; i++){
			PSM one = pList[i];
			for(int k=i+1; k<isize; k++){
				PSM two = pList[k];

				bool diffSeq = one.isDiffSeq(two);
				if( c_chargeSeperate && (one.getChargeState() != two.getChargeState()) ) continue;

				double align[] = {0, 0};
				one.getPeptOverlap(two, align);

				int distance = (int)align[0];
				double overlap = align[1];

				if( overlap < c_minOverlap ) continue;

				char tppair[500];
				sprintf(tppair, "%d_%d", one.getSpectrumIndex(), two.getSpectrumIndex());
				int spsdected = ( spsPair.find(string(tppair)) == spsPair.end() )? 0 : 1;

				bool considered = false, correctly = false, identified = false;
				if( fabs(one.getObservedMW()-two.getObservedMW()) < c_pmTolerance ){
					if( overlap == 1 ){
						if( distance == 0 ) {
							gold_dist0[1]++;
							correctly = true;
							considered = true;
						}
						else if( distance == 1 ) {
							gold_dist1[1]++;
							correctly = true;
							considered = true;
						}
						else if( diffSeq ) {
							gold_dist0[0]++;
							considered = true;
						}

						if( spsdected == 1 ) {
							identified = true;
							if( distance == 0 ) sps_dist0[1]++;
							else if( distance == 1 ) sps_dist1[1]++;
							else if( diffSeq ) 	sps_dist0[0]++;
						}
					}
				}
				else{
					if( overlap == 1 ) {
						if( distance == 1 ) {
							gold_dist1[1]++;
							correctly = true;
							considered = true;
						}
						else if( diffSeq ) {
							gold_dist1[0]++;
							considered = true;
						}
						else if( distance == 0 ) {
							gold_dist0[1]++;
							correctly = true;
							considered = true;
						}

						if( spsdected == 1 ) {
							identified = true;
							if( distance == 1 ) sps_dist1[1]++;
							else if( diffSeq ) 	sps_dist1[0]++;
							else if( distance == 0 ) sps_dist0[1]++;
						}
					}
					else {
						int bin = (int)(overlap*10);

						if( distance == 0 ) {
							gold_part[1][bin]++;
							correctly = true;
							considered = true;
						}
						else if( distance == -1 ) {
							gold_part[0][bin]++;
							considered = true;
						}

						if( spsdected == 1 ){
							identified = true;
							if( distance== 0 ) sps_part[1][bin]++;
							else if( distance == -1 ) sps_part[0][bin]++;
						}
					}
				}

				if( !considered ) continue;

				if( correctly || identified ) {
					fprintf(out, "%.4f\t%d\t", overlap, distance);
                                        
                                        string peptide1 = one.toString();
                                        string peptide2 = two.toString();
                                        
					fprintf(out, "%s\t%s\t", peptide1.c_str(), peptide2.c_str());
					if( correctly ) fprintf(out, "1\t");
					else fprintf(out, "0\t");
					fprintf(out, "%d\t%f\t%f\n", spsdected, pair_data_map[tppair].score1, pair_data_map[tppair].shift1);
				}
				else if( fnprint ){
					fprintf(fn_out, "%.4f\t", overlap);
					fprintf(fn_out, "%s\t%s\n", one.toString().c_str(), two.toString().c_str());
				}
			}
		}

		printf( "MSGF pairs ----------------------------------------------------\n" );
		printf( " for distance 0, correct: %d / incorrect: %d\n", gold_dist0[1], gold_dist0[0] );
		printf( " for distance 1, correct: %d / incorrect: %d\n", gold_dist1[1], gold_dist1[0] );
		printf( " for prefix/suffix,\n" );
		for(int i=9; i>=c_minOverlap*10; i--){
			printf(	" - overlap %.1f,  correct: %d / incorrect: %d\n", (float)i/10, gold_part[1][i], gold_part[0][i] );
		}

		printf( "\nStatistics for SPS pairs --------------------------------------\n" );
		printf( " for dist. 0, sensitivity %d / %d = %.2f\n", sps_dist0[1], gold_dist0[1], ((float)sps_dist0[1]/gold_dist0[1]*100) );
		printf( "                precision %d / %d = %.2f | Incorrect: %d\n", sps_dist0[1], (sps_dist0[0]+sps_dist0[1]),
				((float)sps_dist0[1]/(sps_dist0[0]+sps_dist0[1])*100), sps_dist0[0]);

		printf( " for dist. 1, sensitivity %d / %d = %.2f\n", sps_dist1[1], gold_dist1[1], ((float)sps_dist1[1]/gold_dist1[1]*100) );
		printf( "                precision %d / %d = %.2f | Incorrect: %d\n", sps_dist1[1], (sps_dist1[0]+sps_dist1[1]),
				((float)sps_dist1[1]/(sps_dist1[0]+sps_dist1[1])*100), sps_dist1[0] );

		printf( "\n for prefix/suffix,\n" );
		for(int i=9; i>=c_minOverlap*10; i--){
			printf( " - overlap %.1f, sensitivity %d / %d = %.2f\n", (float)i/10, sps_part[1][i], gold_part[1][i], ((float)sps_part[1][i]/gold_part[1][i]*100) );
			printf( "                  precision %d / %d = %.2f | Incorrect: %d\n", sps_part[1][i], (sps_part[0][i]+sps_part[1][i]),
					((float)sps_part[1][i]/(sps_part[0][i]+sps_part[1][i])*100), sps_part[0][i] );
		}

		printf( "---------------------------------------------------------------\n" );

		fclose(out);
		if( fnprint ) fclose(fn_out);//*/

		return true;
	}

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::loadInputData(void) {
	}
	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::saveOutputData(void) {
	}

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::saveInputData(std::vector<std::string> & filenames) {
	}

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::loadOutputData(void) {
	}

	// -------------------------------------------------------------------------
	vector<ExecBase*> const & ExecQCSpectralPairs::split(int numSplit) {
	}

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::merge(void) {
	}

	// -------------------------------------------------------------------------
	bool ExecQCSpectralPairs::validateParams(std::string & error) {
		m_isValid = false;

	  VALIDATE_PARAM_EXIST("MSGFDB_PSM_FILE");
	  VALIDATE_PARAM_EXIST("SPS_PAIR_FILE");
	  VALIDATE_PARAM_EXIST("OUTPUT_FILE");

	  m_isValid = true;
	  return true;
	}

	vector<PSM> getPSMListfromMSGFDB(const string & fname, double msgfScoreCut, double fdrCut) { // msgfdb results file path
		vector<PSM> pList;
                
                //Loading specset using real code and not the shit below
                PeptideSpectrumMatchSet psms;
                psms.loadMSGFDBResultsFile(fname.c_str());
                
                DEBUG_MSG("Loaded " << psms.size()<<" PSMs");
                
                for(int i = 0; i < psms.size(); i++){
                    if(psms[i]->m_annotation.find("!") != string::npos){
                        //DEBUG_MSG("MIXTURE SKIP");
                        continue;
                    }
                    
                    string massaged_annotation = specnets_annotation_to_msgf(msp_annotation_to_specnets(psms[i]->m_annotation));
                    //cout<<psms[i]->m_annotation<<"\t"<<massaged_annotation<<endl;
                    
                    if( psms[i]->m_score < msgfScoreCut ) continue;
                    if( psms[i]->m_fdr > fdrCut ) continue;
                    
                    PSM cur(psms[i]->m_scanNum, massaged_annotation, psms[i]->m_charge, psms[i]->m_mz, psms[i]->m_protein);
                    cur.setScore(psms[i]->m_score, psms[i]->m_pValue, psms[i]->m_fdr);
                    pList.push_back( cur );
                    
                    //cout<<massaged_annotation<<"\t"<<psms[i]->m_charge<<"\t"<<psms[i]->m_score<<"\t"<<psms[i]->m_pValue<<"\t"<<psms[i]->m_fdr<<endl;
                    //cout<<psms[i]->m_scanNum<<"\t"<<massaged_annotation<<"\t"<<psms[i]->m_charge<<"\t"<<psms[i]->m_mz<<endl;
                }
                
                sort(pList.begin(), pList.end(), PSM::compareBySpectrumIncex);
                
                return pList;
                
/*
		char buf[4096];
		FILE* in = fopen(fname.c_str(), "r");
		if (in == 0) {
			ERROR_MSG("Can not open: " << fname);
			return pList;
		}

		fgets(buf, 4096, in);//MS-GFDB header
		while( fgets(buf, 4096, in) != NULL ){
			char tp_pept[1024]="";
			int tp_specIndex=0, tp_cs=0, tp_msgfscr=0;
			float tp_pmz=0, tp_sprob=0, tp_fdr=0;

			char *tok = NULL;
			if( (tok = strtok(buf, "\t")) == NULL ) continue; //#SpecFile
			int column = 1;
			while( tok = strtok(NULL, "\t") ){
				switch(column){
				case 1: tp_specIndex= atoi(tok); 			break;
			//	case 2: Scan#
			//	case 3: FragMethod
				case 4: tp_pmz 		= atof(tok); 			break;
			//	case 5: PMError
				case 6: tp_cs 		= atoi(tok); 			break;
				case 7: strcpy(tp_pept, tok); 				break;
			//	case 8: Protein
			//	case 9: DeNovoscore
				case 10: tp_msgfscr = atoi(tok); 			break;
				case 11: tp_sprob	= -log10(atof(tok));	break;
			//	case 12: P-value
				case 13: tp_fdr 	= atof(tok); 			break;

				default: break;
				}
				if( column++ == 13 ) break;
			}
			if( column != 14 ) continue;
			if( tp_msgfscr < msgfScoreCut ) continue;
			if( tp_fdr > fdrCut ) continue;

			PSM cur(tp_specIndex, string(tp_pept), tp_cs, tp_pmz);
			cur.setScore(tp_msgfscr, tp_sprob, tp_fdr);
                        
                        cout<<tp_pept<<"\t"<<tp_cs<<"\t"<<tp_msgfscr<<"\t"<<tp_fdr<<endl;

			pList.push_back( cur );
		}
		fclose(in);
		sort(pList.begin(), pList.end(), PSM::compareBySpectrumIncex);

		return pList;*/
	}

	set<string> getSPSPair(const string & filename, std::map<std::string, SpectrumPair> & pair_data_map) { // ./aligns/ pair*.bin file path
		set<string> pair_set;
		SpectrumPairSet pairs;
		pairs.loadFromBinaryFile(filename.c_str());
		unsigned int numEntries = pairs.size();
		for (unsigned int i = 0; i < numEntries; i++){
			char temp[50];
			sprintf(temp, "%d_%d", pairs[i].spec1+1, pairs[i].spec2+1);
                        pair_data_map[string(temp)] = pairs[i];
			pair_set.insert( string(temp) );
		}
		return pair_set;
	}

	PSM::~PSM() {
		// TODO Auto-generated destructor stub
	}

	PSM::PSM(int specindex, string annopept, int cs, double obpmz, string input_protein){
		index = specindex;
		rawAnno = annopept;
		charge = cs;
		obsvMW = (obpmz-proton) * cs;
                protein = input_protein;

		string cvt = annopept;
		for(int i=0; i<cvt.length(); i++){
			if( cvt[i] == 'I' ) cvt[i] = 'L';
			if( cvt[i] == 'Q' ) cvt[i] = 'K';
		}

		stripSeq = "";
		string 			delta;
		vector<int> 	msite;
		vector<string> 	msdelta;

		bool prevAA = true;
		for(int i=2; i<cvt.length()-2; i++){
			if( isalpha(cvt[i]) ){
				stripSeq += cvt[i];
				if(!prevAA){
					msdelta.push_back(delta);
					delta.clear();
				}
				prevAA = true;
			}
			else if(prevAA) {
				modCount++;
				if( stripSeq.length() == 0 ) msite.push_back(0);
				else msite.push_back(stripSeq.length()-1);
				prevAA = false;
				delta += cvt[i];
			}
			else if(!prevAA) {
				delta += cvt[i];
			}
		}
		if(!prevAA){
			msdelta.push_back(delta);
		}

		mods = (int *)new int[stripSeq.length()];
		for(int i=0; i<stripSeq.length(); i++) mods[i]=0;
		if( modCount > 0 ){
			for(int k=0; k<msite.size(); k++){
				mods[msite[k]] += round(atof(msdelta[k].c_str())) ;
			}
		}
	}

	void PSM::getPeptOverlap(PSM& x, double align[2]){

		align[0] = 0; // distance
		if( rawAnno.compare(x.rawAnno) == 0 ) {
			align[1] = 1; // %overlap
			return;
		}

		PSM smallpsm = *this, largepsm = x;
		if( stripSeq.length() > x.stripSeq.length() ){
			smallpsm = x;
			largepsm = *this;
		}

		int shift = 0, smlen=smallpsm.stripSeq.length(), lalen=largepsm.stripSeq.length();
		align[1] = (double)smlen/lalen;

		if( smlen == lalen ){
			for(int i=0; i<smlen; i++){
				if( smallpsm.stripSeq[i] != largepsm.stripSeq[i] || smallpsm.mods[i] != largepsm.mods[i] ) shift++;
			}
		//	if( shift < 2 )
			{
				align[0] = shift;
				return;
			}
		}
		else {
			if( largepsm.stripSeq.find(smallpsm.stripSeq) == 0 ){
				for(int i=0; i<smlen; i++){
					if( smallpsm.mods[i] != largepsm.mods[i] ) shift++;
				}
			//	if( shift == 0 )
				{
					align[0] = shift;
					return;
				}
			}

			shift = 0;
			int pos = lalen - smlen;
			if( largepsm.stripSeq.find(smallpsm.stripSeq, pos) == pos ){
				for(int i=0; i<smlen; i++){
					if( smallpsm.mods[i] != largepsm.mods[i+pos] ) shift++;
				}
			//	if( shift == 0 )
				{
					align[0] = shift;
					return;
				}
			}
		}

		align[0] = -1; // different peptides of different lengths
	}

	void PSM::setScore(float scr, float sprob, float fdr){
		score = scr;
		sFDR = fdr;
		specProb = sprob;
	}

	bool PSM::isDiffSeq(PSM& x){
		return stripSeq.compare(x.stripSeq) != 0;
	}

	string PSM::toString(){
		char temp[500];
                string specnetted_annot;
                specnets::PeptideSpectrumMatch::inspectToSpecNets(rawAnno, specnetted_annot);
		sprintf(temp, "%d\t%s\t%d\t%d\t%f\t%s", index, specnetted_annot.c_str(), charge, (int)score, specProb, protein.c_str());
		return string(temp);
	}

}
