#include "PeptideSpectrumMatchSet.h"
#include "SpectralLibrary.h"
#include "projectionutils.h"
#include "spectrum_window_filter.h"
#include "spectrum.h"
#include "time.h"

//System Includes
#include <map>
#include <vector>
#include <algorithm>

namespace specnets
{

    int SpectralLibrary::createlibrary( float envelope_score_filter,
                                        float pvalue_filter,
                                        MS2ScoringModel &model,
                                        vector<string> &ionsToExtract,
                                        string allIons,
                                        string aminoacidexclusions,
                                        vector<int> charge_filter,
                                        bool filter_best_rep){
        return createlibrary(envelope_score_filter,
                             pvalue_filter,
                             -1,
                             model,
                             ionsToExtract,
                             allIons,
                             aminoacidexclusions,
                             charge_filter,
                             filter_best_rep);

    }

    int SpectralLibrary::createlibrary( float envelope_score_filter,
                                        float pvalue_filter,
                                        float fdr_filter,
                                        MS2ScoringModel &model,
                                        vector<string> &ionsToExtract,
                                        string allIons,
                                        string aminoacidexclusions,
                                        vector<int> charge_filter,
                                        bool filter_best_rep){

        DEBUG_MSG("ORIGINAL LIB SIZE: "<< specs.size());

        //Creating annotation ends
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                string annot = (*it)->m_annotation;
                string new_annot = create_annotation_ends(annot);       //D add ".*" and "*."
                (*it)->m_annotation = new_annot;
            }
        }

        cout<<"After Creating Annotation Ends: "<<specs.size()<<endl;       //D specs = a set of spectra

        filter_no_psm();
        cout<<"Filtering No PSM: "<<specs.size()<<endl;

        if(aminoacidexclusions.length() > 0){
            filter_aminoacids(aminoacidexclusions);
            filter_no_psm();
            cout<<"Filter Amino Acids: "<<specs.size()<<endl;
        }


        if(envelope_score_filter > -1){
            filter_envelope_filter(envelope_score_filter);
            filter_no_psm();
            cout<<"Envelope Filter: "<<specs.size()<<endl;
        }


        if(pvalue_filter > -1){
            filter_pvalue(pvalue_filter);
            filter_no_psm();
            cout<<"P Value Filter: "<<specs.size()<<endl;
        }



        if(fdr_filter > -1){
            filter_fdr(fdr_filter);
            filter_no_psm();
            cout<<"FDR Filter of "<<fdr_filter<<": "<<specs.size()<<endl;
        }


        if(false){
            filter_multiple_interpretations();
            filter_no_psm();
            cout<<"Filtering Multiple Interpretations: "<<specs.size()<<endl;
        }

        if(charge_filter.size() != 0){
            filter_forcharge(charge_filter);
            cout << endl << "charge_filter.size() = " << charge_filter.size();
            cout << "charge_filter[0] = " << charge_filter[0] << endl;
            cout<<"Charge Filter: "<<specs.size()<<endl;
        }

        if(false){      //D  never filter by this criterion
            vector<int> accepted_fragmentation;
            accepted_fragmentation.push_back(Spectrum::FragType_CID);
            filter_forfragmentation(accepted_fragmentation);
            cout<<"Fragmentation Filter: "<<specs.size()<<endl;
        }

        if(filter_best_rep){
            filter_best_representative(model, ionsToExtract, allIons);      //D get a consensus spectrum for each group of spectra of the same annotation and charge
            cout<<"After Best Rep: "<<specs.size()<<endl;
        }

        bool filter_min_peaks = false;
        if(filter_min_peaks){       //D  never filter by this criterion
            //filter_simple_spectra(model, ionsToExtract, allIons);
            //filter_no_psm();
            cout<<"After Spectra Simplicity: "<<specs.size()<<endl;
        }


        return 0;
    }

    int SpectralLibrary::projection(string target_annotation,
                                    MS2ScoringModel model,
                                    vector<string> ions_to_extract,
                                    string allIons,
                                    Spectrum & outputspectrum){

        //Making sure the annotation looks nice
        target_annotation = create_annotation_ends(target_annotation);

        SpecSet projection_spectra;
        vector<float> projections_set_cosine_depression;
        get_possible_projections(target_annotation,
                                 model,
                                 ions_to_extract,
                                 allIons,
                                 projection_spectra,
                                 projections_set_cosine_depression);

        if(projection_spectra.size() == 0) return -1;

        vector<Spectrum *> spectrum_group;
        for(int spectrum_idx = 0; spectrum_idx < projection_spectra.size(); spectrum_idx++)
            spectrum_group.push_back(&projection_spectra[spectrum_idx]);

        get_consensus(spectrum_group, model, ions_to_extract, allIons, outputspectrum);

        outputspectrum.psmList.clear();
        psmPtr psm(new PeptideSpectrumMatch());
        outputspectrum.psmList.push_back(psm);
        psm->m_annotation = target_annotation;
        psm->m_spectrum = &outputspectrum;


        return 0;
    }


    int SpectralLibrary::search_target_decoy(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             int scoring_method){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;


        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;      //D greater than the mass tolerance
            if(query_spec.parentCharge != 0 && charge != 0 && query_spec.parentCharge != charge) continue;      //D what is query spec?

            float sim = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);

            float dot_bias = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]);
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);

            if(tr1::get<7>(scores_tuple_decoy[i])->m_organism.size() > 0){
                psm->m_organism.push_back(tr1::get<7>(scores_tuple_decoy[i])->m_organism[tr1::get<7>(scores_tuple_decoy[i])->m_organism.size()-1]);
            }

            if(tr1::get<7>(scores_tuple_decoy[i])->m_compound_name.size() > 0){
                psm->m_compound_name.push_back(tr1::get<7>(scores_tuple_decoy[i])->m_compound_name[tr1::get<7>(scores_tuple_decoy[i])->m_compound_name.size()-1]);
            }

            search_results_decoy.push_back(psm);
        }


        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);

        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;

            //DEBUG_MSG(query_spec.parentMZ<<"\t"<<library_mass<<"\t"<<charge<<"\t"<<query_spec.parentCharge<<"\t"<<target_library_ptr[library_idx]->psmList.front()->m_compound_name[0]);

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && charge != 0 && query_spec.parentCharge != charge) continue;
            float sim = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<7>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front();

            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]);
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);

            if(tr1::get<7>(scores_tuple[i])->m_organism.size() > 0){
                psm->m_organism.push_back(tr1::get<7>(scores_tuple[i])->m_organism[tr1::get<7>(scores_tuple[i])->m_organism.size()-1]);
            }

            if(tr1::get<7>(scores_tuple[i])->m_compound_name.size() > 0){
                psm->m_compound_name.push_back(tr1::get<7>(scores_tuple[i])->m_compound_name[tr1::get<7>(scores_tuple[i])->m_compound_name.size()-1]);
            }

            search_results.push_back(psm);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }


        if(target_top_scoring > decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<endl;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_spectrumFile  = search_results[0]->m_spectrum->fileName;
            output_psm->m_organism      = search_results[0]->m_organism;
            output_psm->m_compound_name = search_results[0]->m_compound_name;
            //output_psm->m_organism.insert(output_psm->m_organism.end(), search_results[0]->m_organism.begin(), search_results[0]->m_organism.end());
            //output_psm->m_compound_name.insert(output_psm->m_compound_name.end(), search_results[0]->m_compound_name.begin(), search_results[0]->m_compound_name.end());

            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = false;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_spectrumFile  = search_results_decoy[0]->m_spectrum->fileName;
            output_psm->m_organism      = search_results_decoy[0]->m_organism;
            output_psm->m_compound_name = search_results_decoy[0]->m_compound_name;
            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = true;

            return 0;
        }

        return -1;
    }


    /*! \brief Default Spectral Library Search
        Simply returns top_hits using full spectrum cosine scoring
     */

    int SpectralLibrary::search(Spectrum &query_spec, PeptideSpectrumMatchSet &output_psms, float parentmass_tolerance, float ion_tolerance, int top_hits, int analog_search, float score_threshold, int library_search_quality, int minimum_shared_peaks, bool forceExact){
        //AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        //psmPtr psm(new PeptideSpectrumMatch);

        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            //masses.clear();
            //aajumps.getPRMMasses(specs[library_idx].psmList.front()->m_annotation.c_str(), masses);


            charge = specs[library_idx].parentCharge;
            //float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
            float library_mass = specs[library_idx].parentMZ;
            float query_mass = query_spec.parentMZ;

            float mz_error = fabs(query_mass - library_mass);
            float mass_error = fabs(specs[library_idx].parentMass - query_spec.parentMass);


            if(mz_error > parentmass_tolerance and analog_search == 0) continue;

            //if(mass_error > analog_search and analog_search > 0) continue;

            if(specs[library_idx].psmList.size() > 0){
                if(specs[library_idx].psmList.front()->m_library_quality > library_search_quality){
                    continue;
                }
            }


            float ms_error_ppm = mz_error/library_mass * 1000000;

            unsigned int shared_peaks = 0;


            //float sim = full_spectrum_similarity(specs[library_idx], query_spec, shared_peaks);
            float score1, score2;
            float sim = specs[library_idx].scoreMatch(query_spec, ion_tolerance, shared_peaks, score1, score2, forceExact);

            if(sim != sim) sim = 0.f;

            //cout<<"DEBUG_BACKGROUND:\t"<<query_mass<<"\t"<<library_mass<<"\t"<<mz_error<<"\t"<<ms_error_ppm<<"\t"<<query_spec.scan<<"\t"<<specs[library_idx].scan<<"\t"<<abs((int)((int)specs[library_idx].scan - (int)query_spec.scan))<<"\t"<<sim<<"\t"<<shared_peaks<<endl;


            if(sim < score_threshold) continue;

            if(shared_peaks < minimum_shared_peaks) continue;

            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = &query_spec;
            psm->m_score = sim;
            psm->m_dbIndex = specs[library_idx].psmList.front()->m_dbIndex;
            psm->m_annotation = (string)specs[library_idx].psmList.front()->m_annotation;
            psm->m_spectrumFile = query_spec.fileName;
            psm->m_organism = specs[library_idx].psmList.front()->m_organism;
            psm->m_compound_name = specs[library_idx].psmList.front()->m_compound_name;
            psm->m_scanNum = query_spec.scan;
            psm->m_library_name = get_only_filename(specs[library_idx].psmList.front()->m_spectrumFile);
            psm->m_spectrumID = specs[library_idx].psmList.front()->m_spectrumID;
            psm->m_mz_error_ppm = ms_error_ppm;
            psm->m_shared_peaks = shared_peaks;
            psm->m_mz = query_spec.parentMZ;


            //Metadata temp for Metabolomic searches
            psm->m_fdr = specs[library_idx].ITOL;                                       //RT of the library
            psm->m_pValue = query_spec.retention_time;                                  //RT of the query
            psm->m_strict_envelope_score = specs[library_idx].getTotalIonCurrent();     //TIC of library
            psm->m_unstrict_envelope_score = query_spec.getTotalIonCurrent();           //TIC of query
            psm->m_startMass = specs[library_idx].parentMZ;                             //MZ of match
            psm->m_parentmass_difference = mass_error;                                  //Mass Difference
            stringstream ss;
            ss << specs[library_idx].scan;
            psm->m_protein = specs[library_idx].psmList.front()->m_notes + ":" + ss.str();               //filename



            if(specs[library_idx].psmList.front()->m_charge == 0){
                psm->m_charge = query_spec.parentCharge;
            }
            else{
                psm->m_charge = specs[library_idx].psmList.front()->m_charge;
            }

            search_results.push_back(psm);


            /*
            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = &specs[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = library_idx;
            tr1::get<3>(similarity_tuple) = (string)specs[library_idx].psmList.front()->m_annotation;
            scores_tuple.push_back(similarity_tuple);
            */
        }

        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(search_results.m_psmSet.begin(), search_results.m_psmSet.end(), search_results_comparator_psmPtr);



        for(int i = 0; i < min((int)search_results.size(), top_hits); i++){
            output_psms.push_back(search_results[i]);
        }

        return 0;

        /*sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);

        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_scanNum = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            search_results.push_back(psm);
        }*/

        if(search_results.size() > 0){
            //float match_score = 0.6*search_results[0]->m_score;
            //if(search_results.size() == 1) match_score += 0.4 * search_results[0]->m_score;
            //else match_score += 0.4 * search_results[1]->m_score;

            //output_psm->m_spectrum      = search_results[0]->m_spectrum;
            //output_psm->m_score         = search_results[0]->m_score;
            //output_psm->m_scanNum       = search_results[0]->m_scanNum;
            //output_psm->m_annotation    = search_results[0]->m_annotation;
            //output_psm->m_score         = match_score;

            return 0;
        }

        /*for(int score_idx = 0; score_idx < min((int)scores_tuple.size(), 1); score_idx++){
            cout<<tr1::get<0>(scores_tuple[score_idx])->psmList.front()->m_annotation<<"\t"<<tr1::get<1>(scores_tuple[score_idx])<<"\t0BasedIdx\t"<<tr1::get<2>(scores_tuple[score_idx])<<endl;
        }*/

        return -1;
    }

    /*! \brief Loads an mgf file and adds it to current specset

     */
    unsigned int SpectralLibrary::LoadSpecSet_additionalmgf(const char * filename){
        SpecSet tempspecs;
        tempspecs.LoadSpecSet_mgf(filename);

        DEBUG_MSG("Loaded Spec MGF");

        for(int specidx = 0; specidx < tempspecs.size(); specidx++){
            specs.push_back(tempspecs[specidx]);        //D add one
            list<psmPtr>::iterator it;
            for ( it=specs[specs.size()-1].psmList.begin() ; it != specs[specs.size()-1].psmList.end(); it++ ){
                (*it)->m_spectrum = &(specs[specs.size()-1]);
            }
        }

        return 0;
    }

    void SpectralLibrary::filter_no_psm(){
        SpectralLibrary temp_lib;       //D temporary library
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            if(specs[spectrum_idx].psmList.size() >= 1){        //D psmList a list of sequences matched to the spectrum?
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());      //D copy the half-interval [ temp_lib.specs.begin(), temp_lib.specs.begin() )
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    /*! \brief Filters out PSMs that are not the best scoring PValue

     */
    void SpectralLibrary::filter_multiple_interpretations(){
        SpectralLibrary temp_lib;
        //Filtering out spectra with multiple interpretations, or 0 interpretations
        //D 0 interpretations???
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            if(specs[spectrum_idx].psmList.size() > 1){
                specs[spectrum_idx].psmList.clear();
            }
        }
    }

    /*! \brief Filters out spectra that exceed the fdr threshold

     */
    void SpectralLibrary::filter_fdr(float fdr_filter){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            //DEBUG_MSG("MYFDR\t"<<specs[spectrum_idx].psmList.front()->m_fdr);
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_fdr > fdr_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
        }
    }

    /*! \brief Filters out spectra that exceed the envelope score       //D what is the score?

     */
    void SpectralLibrary::filter_envelope_filter(float envelope_score_filter){
        SpectralLibrary temp_lib;
        //Filtering based on envelope scores
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //cout<<"Envelope: "<<specs[spectrum_idx].psmList.front()->m_strict_envelope_score<<endl;

            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_strict_envelope_score > envelope_score_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }

            /*

            if(specs[spectrum_idx].psmList.front()->m_strict_envelope_score < envelope_score_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }

    /*! \brief Filters out spectra that exceed the pvalue threshold

     */
    void SpectralLibrary::filter_pvalue(float pvalue_filter){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues




        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_pValue > pvalue_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            /*if(specs[spectrum_idx].psmList.front()->m_pValue < pvalue_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }


    void SpectralLibrary::filter_forcharge(vector<int> charge){     //D only keep spectra whose charge is in vector<int> charge
        DEBUG_MSG("void SpectralLibrary::filter_forcharge(vector<int> charge){" << endl);
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){

            if(specs[spectrum_idx].parentCharge == 0){
                specs[spectrum_idx].parentCharge = specs[spectrum_idx].psmList.front()->m_charge;
            }

            int spec_charge = specs[spectrum_idx].psmList.front()->m_charge;
            if(spec_charge <= 0) spec_charge = specs[spectrum_idx].parentCharge;
            //D original code: if(spec_charge == -1) spec_charge = specs[spectrum_idx].parentCharge;

            DEBUG_MSG("spec_charge = " << spec_charge);
            bool valid_charge = false;
            for(int charge_idx = 0; charge_idx < charge.size(); charge_idx++){
                if(spec_charge == charge[charge_idx]){
                    valid_charge = true;
                    break;
                }
            }
            if(valid_charge){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    void SpectralLibrary::filter_forfragmentation(vector<int> accepted_fragmentation){
        SpectralLibrary temp_lib;
        //Filtering based on fragmentation
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            bool valid_fragmentation = false;
            int spec_fragmentation = specs[spectrum_idx].msFragType;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(spec_fragmentation == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }
            if(valid_fragmentation){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }


    void SpectralLibrary::filter_aminoacids(string aminoacidexclusions){        //D a string of exclusion aminoacids
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //string annotation = specs[spectrum_idx].psmList.front()->m_annotation;

            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_annotation.find_first_of(aminoacidexclusions) != string::npos){
                        specs[spectrum_idx].psmList.erase(it);      //D erase every annotations of the psmList which contain any exclusion aminoacid
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            //if(annotation.find_first_of(aminoacidexclusions) == string::npos){
            //    temp_lib.specs.push_back(specs[spectrum_idx]);
            //}
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }


    /*! \brief Finds spectrum that best represents each peptide

     */
    void SpectralLibrary::filter_best_representative(   MS2ScoringModel &model,
                                                        vector<string> &ionsToExtract,
                                                        string allIons){
        //Grouping up peptides with the same annotation
        map<string, vector<Spectrum *> > same_annotation_clusters;
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string annotation = specs[spectrum_idx].psmList.front()->m_annotation;      //D only get the front annotation
            annotation += ('0' + specs[spectrum_idx].parentCharge);       //D consider both annotation and charge
            same_annotation_clusters[annotation].push_back(&specs[spectrum_idx]);
        }

        cout<<"Unique Annotations:"<<same_annotation_clusters.size()<<endl;

        //Choosing representative peptide
        SpectralLibrary temp_lib;
        map<string, vector<Spectrum *> >::iterator it;
        for ( it=same_annotation_clusters.begin() ; it != same_annotation_clusters.end(); it++ ){
            Spectrum consensus_spectrum;
            get_consensus((*it).second, model, ionsToExtract, allIons, consensus_spectrum);     //D calculate a consensus spectrum for this group of spectra?
            temp_lib.specs.push_back(consensus_spectrum);
        }

        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();

        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }

    void SpectralLibrary::get_consensus( vector<Spectrum *> spectrum_group,
                                MS2ScoringModel model,
                                vector<string> ionsToExtract,
                                string allIons,
                                Spectrum &consensus_spectrum){
        //Setting psms to be correct
        for(int spec_idx = 0; spec_idx < spectrum_group.size(); spec_idx++){
            spectrum_group[spec_idx]->psmList.front()->m_spectrum = spectrum_group[spec_idx];
        }


        float best_score = 0.f;
        int best_idx = 0;
        for(int cluster_idx1 = 0; cluster_idx1 < spectrum_group.size(); cluster_idx1++){
            float spec_score = 0.f;
            for(int cluster_idx2 = 0; cluster_idx2 < spectrum_group.size(); cluster_idx2++){
                if(cluster_idx1 == cluster_idx2) continue;
                spec_score += spectrum_similarity(  spectrum_group[cluster_idx1]->psmList.front(),
                                                    spectrum_group[cluster_idx2]->psmList.front(),
                                                    spectrum_group[cluster_idx2]->psmList.front()->m_annotation.length() - 4,       //D -4 because of "*." and ".*"???
                                                    model,
                                                    ionsToExtract,
                                                    allIons);
            }
            if (spec_score > best_score){
                best_idx = cluster_idx1;
                best_score = spec_score;
            }
        }
        consensus_spectrum = *(spectrum_group[best_idx]);

    }

    int SpectralLibrary::get_possible_projections(string target_annotation, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpecSet &projection_specs, vector<float> &projections_set_cosine_depression){
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string library_annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            if(getStringDifference(target_annotation, library_annotation) == 1){
                Spectrum new_synth_spec;

                psmPtr psm = specs[spectrum_idx].psmList.front();
                string source_annotation = psm->m_annotation;

                vector<pair<float, float> > ion_mass_intensity_pair_vector;

                string annotation = cleanAnnotationEnds(source_annotation);
                psm->annotate(annotation,allIons,model,0,0,.45);

                extractIons(psm,annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);

                int diff_location = getDifferenceIndex(annotation, target_annotation) - 2;

                //cout<<"Dif Location: "<<diff_location<<endl;
                //cout<<"Source Annotation: "<<annotation<<endl;

                //Calculating the expected cosine depression
                map<string, float> aa_substitution_map = getAminoAcidSubstitutionLookup();
                //cout<<annotation[diff_location+2]<<"\t"<<target_annotation[diff_location+2]<<endl;
                float cosine_depression = getSubstitutionCosineDepression(annotation[diff_location+2], target_annotation[diff_location+2], aa_substitution_map);
                //cout<<"Cosine Depression: "<<cosine_depression<<endl;

                //If we have a cosine depression of less than a value, we skip this projection
                if(cosine_depression < 0.5){
                    //cout<<"Transformation Not Supported"<<endl;
                    //continue;
                }

                AAJumps aajumps(1);
                vector<float> masses;
                //cout<<"AAJUMPS annot: "<<annotation<<"\t"<<target_annotation<<endl;
                aajumps.getPRMMasses(annotation.c_str(), masses);

                int charge = 2;
                map<char, float> aamasses = getAminoAcidLookup();
                float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float library_mass = (getMass(annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                aajumps.getPRMMasses(target_annotation.c_str(), masses);
                float target_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float target_mass = (getMass(target_annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;

                float mass_difference = library_mass - target_mass;

                //cout<<"Library Mass: "<<library_mass<<"\t"<<"Target mass: "<<target_mass<<"\t"<<mass_difference<<endl;
                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                    //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                }


                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    if(ions_to_extract[i/(annotation.length()-4)].find('b') != -1){
                        if( (i%(annotation.length()-4) + 1) < diff_location + 1) //For B ions don't need to edit
                            continue;

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('y') != -1){
                        if( (i%(annotation.length()-4) + 1) < ((annotation.length()-4) - diff_location) ){ //For Y ions don't need to edit
                            //cout<<"Contineu: "<<(i%(annotation.length()-4) + 1)<<endl;
                            continue;
                        }

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('P') != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                }

                new_synth_spec.resize(0);
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
                int new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].first > 0.f){
                        new_peaklist_size++;
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }

                new_synth_spec.resize(new_peaklist_size);
                int new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_synth_spec[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        new_synth_spec[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }

                //adding projection to this
                projections_set_cosine_depression.push_back(cosine_depression);
                psmPtr target_psm(new PeptideSpectrumMatch);
                target_psm->m_annotation = target_annotation;
                new_synth_spec.psmList.push_back(target_psm);
                new_synth_spec.parentMZ = target_mass;
                new_synth_spec.parentCharge = specs[spectrum_idx].parentCharge;
                projection_specs.push_back(new_synth_spec);
                target_psm->m_spectrum = &projection_specs[projection_specs.size()-1];
            }
        }
        return 0;
    }


    vector<vector<float> > SpectralLibrary::find_global_average_spectrum(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        psmPtr psm(new PeptideSpectrumMatch);
        vector<vector<float> > averages;

        for(int sizes = 0; sizes < 100; sizes++){
            vector<float> average;
            int psm_count = 0;
            for(int i = 0; i < specs.size(); i++){
                if(specs[i].psmList.front()->m_annotation.length()-4 == sizes){
                    psm_count++;
                    //cout<<"I: "<<i<<endl;
                    //cout<<m_psmSet[i]->m_annotation<<endl;
                    Spectrum tempspec = specs[i];
                    //Making Max to 1000
                    preprocess_spectrum_intensities_max_intensity(&tempspec, 1000.f);
                    //Applying SQRT
                    preprocess_spectrum_intensities(&tempspec, 0, 1);

                    psm->m_annotation = specs[i].psmList.front()->m_annotation;
                    psm->m_spectrum = &tempspec;

                    psm->annotate(psm->m_annotation,allIons,model,0,0,.45);

                    vector<float> ion_mass_intensity_pair_vector;

                    extractIons(psm,psm->m_annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);
                    norm_vector(ion_mass_intensity_pair_vector);

                    if(average.size() == 0){
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average.push_back(ion_mass_intensity_pair_vector[j]);
                        }
                        continue;
                    }
                    else{
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average[j] += (ion_mass_intensity_pair_vector[j]);
                        }
                    }
                }
            }
            cout<<sizes<<"\t";
            for(int j = 0; j < average.size(); j++){
                average[j] = average[j]/psm_count;
                cout<<average[j]<<"\t";
            }
            cout<<endl;
            averages.push_back(average);
        }

        /*
        average_spectrum->peakList.clear();
        sort( average.begin(), average.end(), mass_intensity_pair_mass_comp);
        int new_peaklist_size = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].first > 0.f){
                new_peaklist_size++;
                //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
            }
        }

        average_spectrum->peakList.resize(new_peaklist_size);
        int new_peaklist_idx = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].second > 0.f){
                average_spectrum->peakList[new_peaklist_idx][0] = average[p].first;
                average_spectrum->peakList[new_peaklist_idx][1] = average[p].second;
                new_peaklist_idx++;
            }
        }*/

        return averages;
    }


    int longest_common_substring(string A, string B)
    {
        int lA = A.length();
        int lB = B.length();
        if ((lA == 0) or (lB == 0))
            return 0;

        int F[lA+1][lB+1];
        int i, j;

        for (i = 0; i <= lA; i++)
            F[i][0] = 0;
        for (j = 0; j <= lB; j++)
            F[0][j] = 0;

        for (i = 1; i <= lA; i++)
            for (j = 1; j <= lB; j++)
                if (A[i-1] == B[j-1])
                    F[i][j] = F[i-1][j-1] + 1;
                else
                    F[i][j] = max(F[i-1][j], F[i][j-1]);

        return F[lA][lB];
    }


    /*D
    //D decoy generating based on longest common subsequence
    SpectralLibrary SpectralLibrary::create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        unsigned int seed = 0;
        srand(seed);

        vector<string> library_peptides;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            annotation = create_annotation_ends(annotation);        //D add .* and *.
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            stripped_annotation += ('0' + specs[library_idx].parentCharge);     //D sequence + parentCharge
            library_peptides.push_back(stripped_annotation);        //D push into library_peptides not the original specs set
        }

        sort(library_peptides.begin(), library_peptides.end());

        SpectralLibrary decoy_library;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;

            //Randomize the annotation
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            string annotation_orig_stripped = stripped_annotation;

            string random_annotation;

            bool valid_decoy_found = false;
            int best_l = 1000000000;        //D for longest common substrings
            string best_random_annotation;
            for(int random_retries = 0 ; random_retries < 100; random_retries++){
                random_annotation = create_decoy_peptide(annotation, specs[library_idx].parentCharge);
                //D cout << random_annotation << endl;

                string search_random_annotation = random_annotation;
                search_random_annotation += ('0'  + specs[library_idx].parentCharge);

                if(binary_search(library_peptides.begin(), library_peptides.end(), search_random_annotation)){
                    continue;   //try again, already in library
                }
                else{
                    valid_decoy_found = true;
                    int l = longest_common_substring(annotation_orig_stripped, random_annotation);
                    //D cout << "l = " << l << endl;
                    if (l < random_annotation.length()/2)       //D accept immediately if these two annotation "share" less than 50%
                    {
                        best_l = l;
                        best_random_annotation = random_annotation;
                        break;
                    }
                    else
                    {
                        if (l < best_l)
                        {
                            best_l = l;
                            best_random_annotation = random_annotation;
                        }
                    }
                }
            }       //D end for retries

            if (!valid_decoy_found){
                cout<<"No Valid Decoy Found"<<endl;
                continue;
            }

            random_annotation = best_random_annotation;

            cout<<decoy_library.size()<<"\t"<<annotation_orig_stripped<<"\t"<<random_annotation << "\t" << best_l << endl;      //D index << original annotation << decoy annotation

            int peptide_length = getpeptideLength(annotation_orig_stripped);

            //Extracting the ions
            vector<pair <float, float> > ion_mass_intensity_pair_vector;
            specs[library_idx].psmList.front()->annotate(annotation,allIons,model,0,0,0.45);
            extractIons(specs[library_idx].psmList.front(), peptide_length, model, ions_to_extract, ion_mass_intensity_pair_vector, 0, 0);


            //D Collecting the peaks that are not annotated in specs[library_idx]
            vector<pair <float, float> > unannotated_peaks;
            for(int peak_idx = 0; peak_idx < specs[library_idx].size(); peak_idx++){
                bool annotated = false;
                for(int annotated_idx = 0; annotated_idx < ion_mass_intensity_pair_vector.size(); annotated_idx++){
                    if(specs[library_idx][peak_idx][0] == ion_mass_intensity_pair_vector[annotated_idx].first &&
                        specs[library_idx][peak_idx][1] == ion_mass_intensity_pair_vector[annotated_idx].second){
                        annotated = true;
                        break;
                    }       //D the peak_idx is annotated
                }

                if(!annotated){
                    pair<float, float> unannotated_peak;
                    unannotated_peak.first = specs[library_idx][peak_idx][0];
                    unannotated_peak.second = specs[library_idx][peak_idx][1];
                    unannotated_peaks.push_back(unannotated_peak);
                }
            }


            vector<string> original_prefix_array;
            vector<string> original_suffix_array;
            vector<string> random_prefix_array;
            vector<string> random_suffix_array;
            generate_prefix_suffix_peptide(annotation_orig_stripped, original_prefix_array, original_suffix_array);     //D original vs random annotation
            generate_prefix_suffix_peptide(random_annotation, random_prefix_array, random_suffix_array);


            AAJumps aajumps(1);

            for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                if( ions_to_extract[i/(peptide_length)].find('b') != -1 ||
                    ions_to_extract[i/(peptide_length)].find('a') != -1){       //D for prefix

                    int charge = 2;
                    vector<float> masses;


                    string orig_prefix = original_prefix_array[i % (peptide_length)];
                    string random_prefix = random_prefix_array[i % (peptide_length)];

                    aajumps.getPRMMasses(create_annotation_ends(orig_prefix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_prefix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;


                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }


                //D for suffix
                if(ions_to_extract[i/(peptide_length)].find('y') != -1){        //D a, b and y?
                    int charge = 2;
                    vector<float> masses;


                    string orig_suffix = original_suffix_array[i % (peptide_length)];
                    string random_suffix = random_suffix_array[i % (peptide_length)];

                    aajumps.getPRMMasses(create_annotation_ends(orig_suffix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_suffix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;


                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){       //D charge 2?

                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }

                if(ions_to_extract[i/(peptide_length)].find('P') != -1){
                    ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first;      //D what???
                }
            }



            Spectrum* new_synth_spec = new Spectrum();      //D create a new decoy spectrum

            //Adding back in the unannotated peaks
            //ion_mass_intensity_pair_vector.clear();
            ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());

            //Adding these peaks to a spectrum
            new_synth_spec->resize(0);
            sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
            int new_peaklist_size = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){
                    new_peaklist_size++;        //D count the number of created peaks
                    //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                }
            }


            new_synth_spec->resize(new_peaklist_size);
            int new_peaklist_idx = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){     //D only consider peaks with intensity value greater than 0
                    (*new_synth_spec)[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                    (*new_synth_spec)[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                    new_peaklist_idx++;
                }
            }


            //adding projection to this
            psmPtr decoy_psm(new PeptideSpectrumMatch);
            decoy_psm->m_annotation = create_annotation_ends(random_annotation);
            decoy_psm->m_charge = specs[library_idx].psmList.front()->m_charge;
            //D careful to consider. The following statement was commented
            //decoy_psm->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];

            new_synth_spec->psmList.push_back(decoy_psm);
            new_synth_spec->parentCharge = specs[library_idx].parentCharge;
            new_synth_spec->parentMass = specs[library_idx].parentMass;
            new_synth_spec->parentMZ = specs[library_idx].parentMZ;
            decoy_library.specs.push_back(*new_synth_spec);
            //D decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];
            decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = new_synth_spec;     //D m_spectrum -> a spectrum pointer
            //delete new_synth_spec;
            DEBUG_MSG("The decoy spectrum's size = " << decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum->size() << endl);

        }

        for(int i = 0; i < decoy_library.specs.size(); i++){
            DEBUG_MSG("hot check! " << decoy_library.specs[i].psmList.front()->m_spectrum->size());
        }
        return decoy_library;
    }
    D*/


    //D decoy generating based on cosine distance
    /*D
    SpectralLibrary SpectralLibrary::create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        unsigned int seed = 0;
        srand(seed);

        vector<string> library_peptides;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            annotation = create_annotation_ends(annotation);        //D add .* and *.
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            stripped_annotation += ('0' + specs[library_idx].parentCharge);     //D sequence + parentCharge
            library_peptides.push_back(stripped_annotation);        //D push into library_peptides not the original specs set
        }

        sort(library_peptides.begin(), library_peptides.end());

        SpectralLibrary decoy_library;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;

            //Randomize the annotation
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            string annotation_orig_stripped = stripped_annotation;

            string random_annotation;

            bool valid_decoy_found = false;
            for(int random_retries = 0 ; random_retries < 100; random_retries++){
                random_annotation = create_decoy_peptide(annotation, specs[library_idx].parentCharge);
                //D cout << random_annotation << endl;

                string search_random_annotation = random_annotation;
                search_random_annotation += ('0'  + specs[library_idx].parentCharge);

                if(binary_search(library_peptides.begin(), library_peptides.end(), search_random_annotation)){
                    continue;   //try again, already in library
                }

                cout<<library_idx<<"\t"<<decoy_library.size()<<"\t"<<annotation_orig_stripped<<"\t"<<random_annotation << endl;      //D index << original annotation << decoy annotation

                int peptide_length = getpeptideLength(annotation_orig_stripped);

                //Extracting the ions
                vector<pair <float, float> > ion_mass_intensity_pair_vector;
                specs[library_idx].psmList.front()->annotate(annotation,allIons,model,0,0,0.45);
                extractIons(specs[library_idx].psmList.front(), peptide_length, model, ions_to_extract, ion_mass_intensity_pair_vector, 0, 0);

                //D Collecting the peaks that are not annotated in specs[library_idx]
                vector<pair <float, float> > unannotated_peaks;
                for(int peak_idx = 0; peak_idx < specs[library_idx].size(); peak_idx++){
                    bool annotated = false;
                    for(int annotated_idx = 0; annotated_idx < ion_mass_intensity_pair_vector.size(); annotated_idx++){
                        if(specs[library_idx][peak_idx][0] == ion_mass_intensity_pair_vector[annotated_idx].first &&
                            specs[library_idx][peak_idx][1] == ion_mass_intensity_pair_vector[annotated_idx].second){
                            annotated = true;
                            break;
                        }       //D the peak_idx is annotated
                    }

                    if(!annotated){
                        pair<float, float> unannotated_peak;
                        unannotated_peak.first = specs[library_idx][peak_idx][0];
                        unannotated_peak.second = specs[library_idx][peak_idx][1];
                        unannotated_peaks.push_back(unannotated_peak);
                    }
                }

                vector<string> original_prefix_array;
                vector<string> original_suffix_array;
                vector<string> random_prefix_array;
                vector<string> random_suffix_array;
                generate_prefix_suffix_peptide(annotation_orig_stripped, original_prefix_array, original_suffix_array);     //D original vs random annotation
                generate_prefix_suffix_peptide(random_annotation, random_prefix_array, random_suffix_array);

                AAJumps aajumps(1);

                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    if( ions_to_extract[i/(peptide_length)].find('b') != -1 ||
                        ions_to_extract[i/(peptide_length)].find('a') != -1){       //D for prefix

                        int charge = 2;
                        vector<float> masses;


                        string orig_prefix = original_prefix_array[i % (peptide_length)];
                        string random_prefix = random_prefix_array[i % (peptide_length)];

                        aajumps.getPRMMasses(create_annotation_ends(orig_prefix).c_str(), masses);
                        float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        aajumps.getPRMMasses(create_annotation_ends(random_prefix).c_str(), masses);
                        float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        float mass_difference = original_mass - random_mass;


                        if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }


                    //D for suffix
                    if(ions_to_extract[i/(peptide_length)].find('y') != -1){        //D a, b and y?
                        int charge = 2;
                        vector<float> masses;


                        string orig_suffix = original_suffix_array[i % (peptide_length)];
                        string random_suffix = random_suffix_array[i % (peptide_length)];

                        aajumps.getPRMMasses(create_annotation_ends(orig_suffix).c_str(), masses);
                        float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        aajumps.getPRMMasses(create_annotation_ends(random_suffix).c_str(), masses);
                        float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        float mass_difference = original_mass - random_mass;


                        if(ions_to_extract[i/(peptide_length)].find("++") != -1){       //D charge 2?

                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(peptide_length)].find('P') != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first;      //D what???
                    }
                }


                Spectrum* new_synth_spec = new Spectrum();      //D create a new decoy spectrum

                //Adding back in the unannotated peaks
                //ion_mass_intensity_pair_vector.clear();
                ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());

                //Adding these peaks to a spectrum
                new_synth_spec->resize(0);
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
                int new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_peaklist_size++;        //D count the number of created peaks
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }


                new_synth_spec->resize(new_peaklist_size);
                int new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){     //D only consider peaks with intensity value greater than 0
                        (*new_synth_spec)[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        (*new_synth_spec)[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }


                //adding projection to this
                psmPtr decoy_psm(new PeptideSpectrumMatch);
                decoy_psm->m_annotation = create_annotation_ends(random_annotation);
                decoy_psm->m_charge = specs[library_idx].psmList.front()->m_charge;
                //D careful to consider. The following statement was commented
                //decoy_psm->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];

                new_synth_spec->psmList.push_back(decoy_psm);
                new_synth_spec->parentCharge = specs[library_idx].parentCharge;
                new_synth_spec->parentMass = specs[library_idx].parentMass;
                new_synth_spec->parentMZ = specs[library_idx].parentMZ;

                unsigned int shared_peaks = 0;
                float score1, score2;
                float cosine_sim = specs[library_idx].scoreMatch(*new_synth_spec, 0.02, shared_peaks, score1, score2, false, false, true);        //D float peakTol = 0.02;
                cout << "cosine_sim = " << cosine_sim << endl;
                if (cosine_sim <= 0.7)       //D float peakTol = 0.02;
                {
                    decoy_library.specs.push_back(*new_synth_spec);
                    //D decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];
                    decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = new_synth_spec;     //D m_spectrum -> a spectrum pointer
                    //delete new_synth_spec;
                    DEBUG_MSG("The decoy spectrum's size = " << decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum->size() << endl);
                    break;
                }
                else
                    delete new_synth_spec;
            }       //D end for retries

            if (!valid_decoy_found){
                cout<<"No Valid Decoy Found"<<endl;
                continue;
            }
        }

        for(int i = 0; i < decoy_library.specs.size(); i++){
            DEBUG_MSG("hot check! " << decoy_library.specs[i].psmList.front()->m_spectrum->size());
        }
        return decoy_library;
    }
    D*/


    SpectralLibrary SpectralLibrary::create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        unsigned int seed = 0;
        srand(seed);

        vector<string> library_peptides;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            annotation = create_annotation_ends(annotation);        //D add .* and *.
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            stripped_annotation += ('0' + specs[library_idx].parentCharge);     //D sequence + parentCharge
            library_peptides.push_back(stripped_annotation);        //D push into library_peptides not the original specs set
        }

        sort(library_peptides.begin(), library_peptides.end());

        SpectralLibrary decoy_library;
        int c2_index = 0;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            if (specs[library_idx].parentCharge == 2)
                c2_index++;
            cout << "c2_index = " << c2_index << endl;

            string annotation = specs[library_idx].psmList.front()->m_annotation;

            //Randomize the annotation
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            string annotation_orig_stripped = stripped_annotation;

            string random_annotation;

            bool valid_decoy_found = false;
            for(int random_retries = 0 ; random_retries < 100; random_retries++){
                random_annotation = create_decoy_peptide(annotation, specs[library_idx].parentCharge);
                //D cout << random_annotation << endl;

                string search_random_annotation = random_annotation;
                search_random_annotation += ('0'  + specs[library_idx].parentCharge);

                if(binary_search(library_peptides.begin(), library_peptides.end(), search_random_annotation)){
                    continue;   //try again, already in library
                }

                cout<<library_idx<<"\t"<<decoy_library.size()<<"\t"<<annotation_orig_stripped<<"\t"<<random_annotation << endl;      //D index << original annotation << decoy annotation

                int peptide_length = getpeptideLength(annotation_orig_stripped);

                //Extracting the ions
                vector<pair <float, float> > ion_mass_intensity_pair_vector;
                specs[library_idx].psmList.front()->annotate(annotation,allIons,model,0,0,0.1);        //D meaningless since annotate is deprecated?
                extractIons(specs[library_idx].psmList.front(), peptide_length, model, ions_to_extract, ion_mass_intensity_pair_vector, 0, 0);

                //D Collecting the peaks that are not annotated in specs[library_idx]
                vector<pair <float, float> > unannotated_peaks;
                for(int peak_idx = 0; peak_idx < specs[library_idx].size(); peak_idx++){
                    bool annotated = false;
                    for(int annotated_idx = 0; annotated_idx < ion_mass_intensity_pair_vector.size(); annotated_idx++){
                        if(specs[library_idx][peak_idx][0] == ion_mass_intensity_pair_vector[annotated_idx].first &&
                            specs[library_idx][peak_idx][1] == ion_mass_intensity_pair_vector[annotated_idx].second){
                            annotated = true;
                            break;
                        }       //D the peak_idx is annotated
                    }

                    if(!annotated){
                        pair<float, float> unannotated_peak;
                        unannotated_peak.first = specs[library_idx][peak_idx][0];
                        unannotated_peak.second = specs[library_idx][peak_idx][1];
                        unannotated_peaks.push_back(unannotated_peak);
                    }
                }

                vector<string> original_prefix_array;
                vector<string> original_suffix_array;
                vector<string> random_prefix_array;
                vector<string> random_suffix_array;
                generate_prefix_suffix_peptide(annotation_orig_stripped, original_prefix_array, original_suffix_array);     //D original vs random annotation
                generate_prefix_suffix_peptide(random_annotation, random_prefix_array, random_suffix_array);

                AAJumps aajumps(1);

                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){     //D just consider a, b and y
                    if( ions_to_extract[i/(peptide_length)].find('b') != -1 ||
                        ions_to_extract[i/(peptide_length)].find('a') != -1){       //D for prefix

                        int charge = 2;
                        vector<float> masses;

                        string orig_prefix = original_prefix_array[i % (peptide_length)];
                        string random_prefix = random_prefix_array[i % (peptide_length)];

                        aajumps.getPRMMasses(create_annotation_ends(orig_prefix).c_str(), masses);
                        float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        aajumps.getPRMMasses(create_annotation_ends(random_prefix).c_str(), masses);
                        float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        float mass_difference = original_mass - random_mass;

                        if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }


                    //D for suffix
                    if(ions_to_extract[i/(peptide_length)].find('y') != -1){        //D a, b and y?
                        int charge = 2;
                        vector<float> masses;

                        string orig_suffix = original_suffix_array[i % (peptide_length)];
                        string random_suffix = random_suffix_array[i % (peptide_length)];

                        aajumps.getPRMMasses(create_annotation_ends(orig_suffix).c_str(), masses);
                        float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        aajumps.getPRMMasses(create_annotation_ends(random_suffix).c_str(), masses);
                        float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                        float mass_difference = original_mass - random_mass;

                        if(ions_to_extract[i/(peptide_length)].find("++") != -1){       //D charge 2?

                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(peptide_length)].find('P') != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first;      //D what???
                    }
                }

                Spectrum* new_synth_spec = new Spectrum();      //D create a new decoy spectrum
                //adding projection to this
                psmPtr decoy_psm(new PeptideSpectrumMatch);
                decoy_psm->m_annotation = create_annotation_ends(random_annotation);
                decoy_psm->m_charge = specs[library_idx].psmList.front()->m_charge;
                //D careful to consider. The following statement was commented
                //decoy_psm->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];

                new_synth_spec->psmList.push_back(decoy_psm);
                new_synth_spec->parentCharge = specs[library_idx].parentCharge;
                new_synth_spec->parentMass = specs[library_idx].parentMass;
                new_synth_spec->parentMZ = specs[library_idx].parentMZ;

                //Adding back in the unannotated peaks
                //ion_mass_intensity_pair_vector.clear();
                //D ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());

                //Adding these peaks to a spectrum
                new_synth_spec->resize(0);
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);

                int new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_peaklist_size++;        //D count the number of created peaks
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }

                new_synth_spec->resize(new_peaklist_size);
                int new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){     //D only consider peaks with intensity value greater than 0
                        (*new_synth_spec)[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        (*new_synth_spec)[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }

                new_synth_spec->filterLowMassPeaks(200.0);
                //D specs[library_idx].filterLowMassPeaks(160.0);
                list<int> peaksToRemove;
                if (annotation_orig_stripped.find("(M,15.994915)") != string::npos)
                {
                    cout << endl << "delpeak" << annotation_orig_stripped << endl;
                    float del_location = new_synth_spec->parentMZ - 63.99/new_synth_spec->parentCharge;
                    for (int i = 0; i < (*new_synth_spec).size(); i++)
                        if ((del_location-0.1 <= (*new_synth_spec)[i][0]) and ((*new_synth_spec)[i][0] <= del_location+0.1))        //D float comparTol = 0.1;
                        {
                            peaksToRemove.push_back(i);
                            cout << i << ' ' << endl;
                        }
                }
                new_synth_spec->removePeaks(peaksToRemove);

                unsigned int shared_peaks = 0;
                float score1, score2;
                float cosine_sim = specs[library_idx].scoreMatch(*new_synth_spec, 0.02, shared_peaks, score1, score2, false, false, true);        //D float peakTol = 0.02;
                cout << "cosine_sim = " << cosine_sim << endl;

                ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);

                new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_peaklist_size++;        //D count the number of created peaks
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }

                new_synth_spec->resize(new_peaklist_size);      //D copy back
                new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){     //D only consider peaks with intensity value greater than 0
                        (*new_synth_spec)[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        (*new_synth_spec)[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }
                /*D
                for (int i = 0; i < (*new_synth_spec).size(); i++)
                {
                    cout << (*new_synth_spec)[i][0] << "," << (*new_synth_spec)[i][1] << "; ";
                }
                cout << endl << endl;

                for (int i = 0; i < specs[library_idx].size(); i++)
                {
                    cout << specs[library_idx][i][0] << "," << specs[library_idx][i][1] << "; ";
                }
                cout << endl << endl;
                D*/

                if (cosine_sim <= 0.4)       //D float peakTol = 0.02;
                {
                    decoy_library.specs.push_back(*new_synth_spec);
                    //D decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];
                    decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = new_synth_spec;     //D m_spectrum -> a spectrum pointer
                    //delete new_synth_spec;
                    valid_decoy_found = true;
                    DEBUG_MSG("The decoy spectrum's size = " << decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum->size() << endl);
                    break;
                }
                else
                    delete new_synth_spec;
            }       //D end for retries

            if (!valid_decoy_found){
                cout<<"No Valid Decoy Found"<< " " << annotation_orig_stripped << " " << specs[library_idx].parentCharge << endl;
                continue;
            }
        }

        for(int i = 0; i < decoy_library.specs.size(); i++){
            DEBUG_MSG("hot check! " << decoy_library.specs[i].psmList.front()->m_spectrum->size());
        }
        return decoy_library;
    }


    string SpectralLibrary::create_decoy_peptide(string peptide, int charge){
        vector<string> deliminated_aminoacids;      //D each string contains a non modification AA or a parenthesis of a modification AA
        string stripped_peptide = remove_annotation_ends(peptide);
        int seed = hashpeptide(peptide, charge) + rand()%1000000000;        //D sum of all characters on the sequence
        //cout<<peptide<<"\tSEED\t"<<seed<<endl;
        srand(seed);

        //cout<<stripped_peptide<<endl;

        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];       //D contains only an amino acid
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }

            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }
        }       //D end for pepidx

        vector<string> randomized_decoy_array;
        vector<string> remaining_aa_to_assign;
        vector<int> empty_indices;

        //Transferring over immobile amino acids
        for(int i = 0 ; i < deliminated_aminoacids.size(); i++){        //D keep R, K and P in the original positions
            if( deliminated_aminoacids[i].find("R") != string::npos ||
                deliminated_aminoacids[i].find("K") != string::npos ||
                deliminated_aminoacids[i].find("P") != string::npos){

                randomized_decoy_array.push_back(deliminated_aminoacids[i]);
                continue;
            }
            else{
                string temp = "";
                randomized_decoy_array.push_back(temp);
            }
        }


        //Now we clear out R, K, P from original list
        //D
        for(int i = 0; i < deliminated_aminoacids.size(); i++){     //D push strings of non R, K and P characters
            if(!( deliminated_aminoacids[i].find("R") != string::npos ||
                deliminated_aminoacids[i].find("K") != string::npos ||
                deliminated_aminoacids[i].find("P") != string::npos)){
                remaining_aa_to_assign.push_back(deliminated_aminoacids[i]);
                empty_indices.push_back(i);
            }
        }

        while(empty_indices.size() > 0){
            int rand_idx = rand()%remaining_aa_to_assign.size();
            randomized_decoy_array[empty_indices[0]] = remaining_aa_to_assign[rand_idx];
            empty_indices.erase(empty_indices.begin());     //D delete at the index 0
            remaining_aa_to_assign.erase(remaining_aa_to_assign.begin() + rand_idx);        //D delete at the random position selected of the non R, K and P ones
        }

        string randomized = "";
        for(int i = 0; i < randomized_decoy_array.size(); i++){
            randomized +=randomized_decoy_array[i];
        }

        return randomized;
    }

    void SpectralLibrary::generate_prefix_suffix_peptide(string peptide, vector<string> &prefix, vector<string> &suffix){
        vector<string> deliminated_aminoacids;
        string stripped_peptide = remove_annotation_ends(peptide);

        //cout<<stripped_peptide<<endl;

        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];       //D temp_str is a character only?
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);     //D vector of strings
                pepidx = pepidx_right;

                continue;
            }

            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }

                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }

                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;

                continue;
            }
        }

        prefix.clear();
        suffix.clear();
        for(int i = 0; i < deliminated_aminoacids.size(); i++){
            string prefix_annotation = "";
            string suffix_annotation = "";
            for(int j = 0; j <= i; j++){
                prefix_annotation += deliminated_aminoacids[j];
                suffix_annotation += deliminated_aminoacids[deliminated_aminoacids.size() - i + j - 1];     //D always contains the last AA
            }
            prefix.push_back(prefix_annotation);
            suffix.push_back(suffix_annotation);
        }
    }


    int SpectralLibrary::add_update_spectrum_to_Library(Spectrum & addition_spectrum){
        if(addition_spectrum.psmList.size() == 0){
            DEBUG_MSG("psm list is empty");
            return -1;
        }
        if(addition_spectrum.psmList.front()->m_organism.size() == 0){
            DEBUG_VAR(addition_spectrum.psmList.front()->m_organism);
            DEBUG_MSG("organism name is empty");
            return -1;
        }

        bool updated = false;
        for(int i = 0; i < this->specs.size(); i++){
            if(this->specs[i].fileName == addition_spectrum.fileName && this->specs[i].scan == addition_spectrum.scan){
                if(this->specs[i].psmList.front()->m_organism.size() > 0){
                    int most_recent_idx_lib = this->specs[i].psmList.front()->m_organism.size()-1;
                    int most_recent_idx_add = addition_spectrum.psmList.front()->m_organism.size()-1;
                    if( this->specs[i].psmList.front()->m_organism[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_organism[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_compound_name[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_compound_name[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_smiles[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_smiles[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_InChI[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_InChI[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_InChI_Aux[most_recent_idx_lib] == addition_spectrum.psmList.front()->m_InChI_Aux[most_recent_idx_add] &&
                        this->specs[i].psmList.front()->m_notes == addition_spectrum.psmList.front()->m_notes){

                        //return i+1;
                        return -2;
                    }
                    else{
                        this->specs[i].psmList.front()->m_submission_metadata.push_back(addition_spectrum.psmList.front()->m_submission_metadata[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_organism.push_back(addition_spectrum.psmList.front()->m_organism[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_compound_name.push_back(addition_spectrum.psmList.front()->m_compound_name[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_smiles.push_back(addition_spectrum.psmList.front()->m_smiles[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_InChI.push_back(addition_spectrum.psmList.front()->m_InChI[most_recent_idx_add]);
                        this->specs[i].psmList.front()->m_InChI_Aux.push_back(addition_spectrum.psmList.front()->m_InChI_Aux[most_recent_idx_add]);
                        updated = true;
                        return i+1;
                    }
                }
                else{
                    return i+1;
                }
            }
        }

        if (!updated){
            this->specs.push_back(addition_spectrum);
            return specs.size();
        }
        return 0;
    }



     void SpectralLibraryGFSearch::train_distribution(SpectralLibrary &library, MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        DEBUG_MSG("START Training");
        DEBUG_MSG("library.size() = " << library.size());

        float annotation_tolerance = 0.1; //D Eventually will be a param value

        //Initializing the bands
        int number_bands = 10;      //D what is bands? - number of lines on the histogram!
        /* for absolute-intensity histogram
        int number_bands = 8;
        float bounds[] = {7500, 12500, 17500, 22500, 27500, 45000, 90000, 1000000000};
        */

        vector<vector< float > > band_deltas;       //D what is band_deltas?
        vector< float > ion_deletion;
        vector< float > ion_insertion;
        vector< float > ion_nondeletion;        //D non
        vector< float > ion_noninsertion;
        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            vector<float> band_delta;
            band_deltas.push_back(band_delta);
        }

        //Clustering all the same annotation spectra        //D why need to cluster?
        map<string, vector<Spectrum *> > same_annotation_clusters;
        for(int spectrum_idx = 0; spectrum_idx < library.size(); spectrum_idx++){
            DEBUG_MSG("library[spectrum_idx].psmList.size() = " << library[spectrum_idx].psmList.size());
            DEBUG_MSG(library[spectrum_idx].psmList.front()->m_annotation);     //D consider the front annotation
            DEBUG_MSG(library[spectrum_idx].psmList.back()->m_annotation);
            DEBUG_MSG("file_scan: " << library[spectrum_idx].fileName << " " << library[spectrum_idx].scan);
            DEBUG_MSG("filesca_n: " << library[spectrum_idx].psmList.back()->m_spectrumFile << " " << library[spectrum_idx].psmList.back()->m_scanNum);
            DEBUG_MSG("filescan: " << library[spectrum_idx].psmList.front()->m_spectrumFile << " " << library[spectrum_idx].psmList.front()->m_scanNum);

            stringstream ss (stringstream::in | stringstream::out);
            ss<<library[spectrum_idx].parentCharge;
            string annotation = library[spectrum_idx].psmList.front()->m_annotation;        //D .front() only?
            //D string annotation = library[spectrum_idx].psmList.back()->m_annotation;
            annotation += ss.str();     //D sequence + parentcharge
            same_annotation_clusters[annotation].push_back(&library[spectrum_idx]);     //D set of spectra of the same annotation and parent charge

            DEBUG_MSG(same_annotation_clusters[annotation].size());     //D print the number of spectra of the same annotation
        }

        map<string, vector<Spectrum *> >::iterator it;
        for ( it=same_annotation_clusters.begin() ; it != same_annotation_clusters.end(); it++ ){
            DEBUG_MSG("TRAININGON\t"<<(*it).first<<"\t"<<(*it).second.size());
            string cluster_annotation = (*it).first;        //D sequence + parentCharge
            int cluster_charge = (*it).second[0]->parentCharge;     //D the parentCharge is also a substring of cluster_annotation

            //TR if(((*it).second.size()) < 10) continue;        //D Why? only consider a set of greater than 10 spectra?

            AAJumps jumps(1);       //D what is the procedure?
            psmPtr training_library_consensus_psm;

            DEBUG_MSG("this->size() = " << this->size() << endl);
            DEBUG_MSG("whether this is m_distilled_spectral_lib" << endl);
            bool found_library_spectrum = false;
            for(int library_idx = 0; library_idx < this->size(); library_idx++){
                if(this->specs[library_idx].psmList.front()->m_annotation == (*it).second[0]->psmList.front()->m_annotation &&
                    this->specs[library_idx].parentCharge == cluster_charge){
                    DEBUG_MSG("FOUND LIBRARY");
                    DEBUG_MSG("this->specs[library_idx].psmList.size() = " << this->specs[library_idx].psmList.size());
                    training_library_consensus_psm = this->specs[library_idx].psmList.front();      //D the refined library spectrum is the consensus spectrum for the spectrum group
                    training_library_consensus_psm->m_spectrum = &this->specs[library_idx];
                    found_library_spectrum = true;
                    break;      //D
                }
            }

            if(!found_library_spectrum) continue;       //D consider annotation clusters which appear on m_distilled_spectral_lib only

            string library_annotation = training_library_consensus_psm->m_annotation;
            training_library_consensus_psm->annotate(training_library_consensus_psm->m_annotation, allIons, model, 0,0, jumps);
            vector<float> library_consensus_ions;       //D a vector of floats

            extractIons(training_library_consensus_psm, create_deliminated_aminoacids(library_annotation).size() ,model,ions_to_extract,library_consensus_ions, 0, 0);
            sqrt_vector(library_consensus_ions);        //D sqrt all elements, but why?
            normalize_extracted_ions(library_consensus_ions);       //D then normalize the vector


            cout << endl << (*it).first << ":";
            for(int ion_idx = 0; ion_idx < library_consensus_ions.size(); ion_idx++){
                cout << library_consensus_ions[ion_idx] << " ";
            }
            cout << endl;


            for(int cluster_spec_idx = 0; cluster_spec_idx < (*it).second.size(); cluster_spec_idx++){
                vector<float> library_training_ions;        //D ERRORRRRR

                float abundance = (*it).second[cluster_spec_idx]->getTotalIonCurrent();     //D getTotalIonCurrent

                if(abundance < 12000.f) continue;
                //TR if(abundance > 12000.f) continue;       //D only consider abundance > 12000

                psmPtr training_library_psm = (*it).second[cluster_spec_idx]->psmList.back();
                training_library_psm->m_spectrum = (*it).second[cluster_spec_idx];
                training_library_psm->annotate(training_library_psm->m_annotation, allIons, model, 0,0, jumps);
                DEBUG_MSG("training_library_psm->m_annotation = " << training_library_psm->m_annotation);
                cout << "training_library_psm->m_annotation = " << training_library_psm->m_annotation << endl;

                //D training_library_consensus_psm vs training_library_psm?
                extractIons(training_library_psm, create_deliminated_aminoacids(library_annotation).size() ,model,ions_to_extract,library_training_ions, 0, 0);
                vector<float> ori_library_training_ions = library_training_ions;
                sqrt_vector(library_training_ions);     //D normalize the sqrts of intensities?
                normalize_extracted_ions(library_training_ions);

                //D calculate cosine SIMILARITIES now!!
                float dot_product = spectrum_similarity(training_library_psm, training_library_consensus_psm, create_deliminated_aminoacids(library_annotation).size(), model, ions_to_extract, allIons);
                float explained_intensity;
                float library_side_dot_product = spectrum_similarity_sqrt_librarypeaks(training_library_consensus_psm, training_library_psm, create_deliminated_aminoacids(library_annotation).size(), model, ions_to_extract, allIons, explained_intensity);
                float full_spectrum_dot = full_spectrum_similarity(*training_library_psm->m_spectrum, *training_library_consensus_psm->m_spectrum);     //D transform spectra to vectors of max-mass elements, then calculate cosine the similarity
                DEBUG_MSG("SIM\t"<<dot_product<<"\tfull spec\t"<<full_spectrum_dot<<"\tLibrarySideSim\t"<<library_side_dot_product<<"\t"<<training_library_psm->m_spectrumFile<<"\t"<<training_library_psm->m_spectrum->scan<<"\t"<<(*training_library_psm->m_spectrum).getTotalIonCurrent()<<"\t"<<cluster_annotation<<"\t"<<training_library_consensus_psm->m_spectrum->parentMZ<<"\t"<<training_library_psm->m_spectrum->parentMZ);

                /*
                cout << endl << cluster_spec_idx << " " << training_library_psm->m_spectrumFile << " " << training_library_psm->m_scanNum << ":";
                for(int ion_idx = 0; ion_idx < library_training_ions.size(); ion_idx++){
                    cout << library_training_ions[ion_idx] << " ";
                }
                cout << endl;
                */

                for(int ion_idx = 0; ion_idx < library_training_ions.size(); ion_idx++){

                    float consensus_ion = library_consensus_ions[ion_idx];      //D what is the different between library_consensus and library_training??
                    float library_ion = library_training_ions[ion_idx];
                    float ori_library_ion = ori_library_training_ions[ion_idx];
                    float max_library_ion = max_ion_intensity(library_consensus_ions);      //D get maximum value of consensus ions

                    if(consensus_ion > 0.0001 && library_ion < 0.0001){
                        //D library_ion < 0.0001 is considered as Deletion Probability???
                        //D float percent_intensity = min(library_ion/max_library_ion, 1.f);        //D Why need min? the value of the ratio is always less than 1.0
                        float percent_intensity = min(consensus_ion/max_library_ion, 1.f);
                        ion_deletion.push_back(percent_intensity);
                        continue;
                    }
                    if(consensus_ion < 0.0001 && library_ion < 0.0001){
                        //D Insertion, not concerned
                        continue;
                    }
                    if(consensus_ion < 0.0001 && library_ion < 0.0001){
                        //Nothing
                        continue;
                    }

                    //D how about the case consensus_ion < 0.0001 && library_ion > 0.0001?


                    float ratio = library_ion/consensus_ion;
                    //DEBUG_MSG("library\t"<<library_ion<<"\t"<<consensus_ion);

                    float percent_intensity = min(consensus_ion/max_library_ion, 1.f);      //D why min is needed here?
                    ion_nondeletion.push_back(percent_intensity);       //D both nondeletion and noninsertion?
                    ion_noninsertion.push_back(percent_intensity);


                    float delta_ion_intensity = (log(library_ion/consensus_ion))/log(2);        //D R_i/L_i? log2 of the ratio

                    //Both are present, lets find a bucket for it
                    int band_idx = min(consensus_ion/max_library_ion, 0.999999f)/(1.0/number_bands);        //D divided by number_bands

                    /*  for absolute-intensity histogram
                    int band_idx;
                    for (band_idx = 0; band_idx < number_bands; band_idx++)
                        if (ori_library_ion <= bounds[band_idx])
                            break;
                    */

                    if(band_idx >= number_bands){
                        cout<<"BAND IDX ERROR"<<"\t"<<consensus_ion<<"\t"<<max_library_ion<<"\t"<<band_idx<<endl;
                        continue;
                    }

                    if(delta_ion_intensity > 1000) continue;        //D only continue if delta_ion_intensity <= 1000? sure thing with the log function?

                    band_deltas[band_idx].push_back(delta_ion_intensity);       //D push the log ratio into its appropriate band

                }
            }

            //Skip Quadratic
            continue;

        }


        DEBUG_MSG("CREATING HISTOGRAMS");
        vector<float> band_averages;        //D averages
        vector<float> band_variances;       //D variances
        //Training Cosine Distributions
        histograms.clear();
        int buckets = 200;      //D bucket size = 200
        float start = -4;
        float end = 4;
        for(int band_idx = 0; band_idx < number_bands; band_idx++){     //D what does bands mean, default value is equal to ten?
            vector<pair<float, float> > histogram = create_histogram(buckets, start, end, band_deltas[band_idx], true);     //D pairs of <index, count>, note: count can be normalized
            histograms.push_back(histogram);
        }

        for(int i = 0; i < buckets; i++){       //D vector<vector<pair<float, float> > > histograms;
            cout<<histograms[0][i].first<<"\t";
            for(int band_idx = 0; band_idx < number_bands; band_idx++){
                cout<<histograms[band_idx][i].second<<"\t";
            }
            cout<<endl;
        }

        //Calculating Average and Variance for each band
        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            double sum = 0.f;
            for(int index = 0; index < band_deltas[band_idx].size(); index++){
                sum += band_deltas[band_idx][index];
                if(band_deltas[band_idx][index] > 100){
                    cout<<band_deltas[band_idx][index]<<endl;
                }
            }
            double mean = sum/band_deltas[band_idx].size();
            band_averages.push_back(mean);      //D band_averages
            cout<<"Band\t"<<band_idx<<"\tAverage\t"<<sum/band_deltas[band_idx].size()<<"\t"<<band_deltas[band_idx].size()<<"\t"<<sum<<endl;
        }


        for(int band_idx = 0; band_idx < number_bands; band_idx++){
            double square_sum = 0.f;
            for(int index = 0; index < band_deltas[band_idx].size(); index++){
                square_sum += band_deltas[band_idx][index] * band_deltas[band_idx][index];
            }
            double variance = square_sum/band_deltas[band_idx].size() - band_averages[band_idx] * band_averages[band_idx];
            band_variances.push_back(variance);
            cout<<"Band\t"<<band_idx<<"\tVariance\t"<<variance<<endl;
        }

        //Training insertion/deletion distributions
        int insertion_deletion_buckets = 11;        //D
        float insertion_delection_start = -0.15f;
        float insertion_delection_end = 0.95;       //D 1.0 for each bucket

        //D non normalized histograms
        //D ion_deletion = {library_ion/max_library_ion}, ion_nondeletion = {consensus_ion/max_library_ion}
        vector < pair < float, float > > ion_deletion_histogram = create_histogram(insertion_deletion_buckets, insertion_delection_start, insertion_delection_end, ion_deletion, false);
        vector < pair < float, float > > ion_nondeletion_histogram = create_histogram(insertion_deletion_buckets, insertion_delection_start, insertion_delection_end, ion_nondeletion, false);


        ion_deletion_histogram.erase(ion_deletion_histogram.begin());       //D from 11 buckets to 10 buckets
        ion_nondeletion_histogram.erase(ion_nondeletion_histogram.begin());

        deletion_probability.clear();       //D vector< pair < float, float > > deletion_probability;

        for(int i = 0; i < ion_deletion_histogram.size(); i++){
            //DEBUG_MSG("DELETION\t"<<ion_deletion_histogram[i].second<<"\t"<<ion_nondeletion_histogram[i].second);

            pair< float, float > deletion_pair;
            deletion_pair.first = ion_deletion_histogram[i].first;
            deletion_pair.second = ion_deletion_histogram[i].second/(ion_nondeletion_histogram[i].second + ion_deletion_histogram[i].second);     //D alright. calculate deletion probabilities
            deletion_probability.push_back(deletion_pair);

            cout<<"INSERTION_DELETION\t"<<deletion_pair.first<<"\t"<<deletion_pair.second<<endl;
        }

        DEBUG_MSG("DONE Training");

    }


    void SpectralLibraryGFSearch::generate_distributions(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        DEBUG_MSG("START SLGF Generation");
        SpectralLibrary temp;
        train_library_dynamicprogramming(model, ions_to_extract, allIons, temp);        //D MS2ScoringModel!!!

        return;

    }

    //D the dynamic procedure in the paper?
    void SpectralLibraryGFSearch::train_library_dynamicprogramming(MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpectralLibrary &library){        //D can remove "library"
        vector<float> band_variance;        //D band likes the order of magnitude? 10000, 100000...
        vector<float> band_mean;


        //Calculating the variance and mean from the band histograms
        //D histograms have been already loaded before (histogram_input_prefix)
        for(int band_idx = 0; band_idx < histograms.size(); band_idx++){
            float sum_mass = 0.f;
            float sum_weighted_mass = 0.f;
            for(int idx = 0; idx < histograms[band_idx].size(); idx++){     //D this size is ~200, note that histograms are already normalized
                sum_mass += histograms[band_idx][idx].second;       //D weights
                sum_weighted_mass += histograms[band_idx][idx].second * histograms[band_idx][idx].first;
            }
            float mean = sum_weighted_mass/sum_mass;
            band_mean.push_back(mean);      //D mean for each band or histogram
            DEBUG_MSG("BAND\t"<<band_idx<<"\t"<<sum_mass<<"\t"<<sum_weighted_mass);     //D sum_mass is always equal to 0 because of normalizing?

            float sum_variance = 0.f;       //D weighted sum_variance
            for(int idx = 0; idx < histograms[band_idx].size(); idx++){
                sum_variance += histograms[band_idx][idx].second * (histograms[band_idx][idx].first - mean)*(histograms[band_idx][idx].first - mean);
            }
            float variance = sum_variance/sum_mass;
            band_variance.push_back(variance);
            DEBUG_MSG("VAR\t"<<band_idx<<"\t"<<variance);
        }




        DEBUG_MSG("this->specs.size() = " << this->specs.size() << endl);
        for(int spec_idx = 0; spec_idx < this->specs.size(); spec_idx++){       //D iterate all spectra in the library?
        //for(int spec_idx = 0; spec_idx < 100; spec_idx++){
            //int thread_id = omp_get_thread_num();

            //if(thread_id == 0){
            //    cout<<"Spec IDX\t"<<spec_idx*6<<" of "<<this->specs.size()<<endl;
            //}


            vector<pair<float, float> > insertion_probability;      //D this vector for what?
            dynamicprogramming_spectrum(this->specs[spec_idx],
                                        model,
                                        ions_to_extract,
                                        allIons,
                                        band_mean,
                                        band_variance,
                                        insertion_probability,
                                        deletion_probability,
                                        this->specs[spec_idx].psmList.front()->SLGF_distribution,
                                        histograms, 4);
        }
    }

    void SpectralLibraryGFSearch::save_distributions(string output_file_path){
        string output_histogram_path = output_file_path + ".histogram";
        string output_deletion_path = output_file_path + ".deletionprob";
        save_cosine_distribution(histograms, output_histogram_path);
        save_deletion_distribution(deletion_probability, output_deletion_path);
    }

    void SpectralLibraryGFSearch::load_distributions(string input_path){
        string input_histogram_path = input_path + ".histogram";
        string input_deletion_path = input_path + ".deletionprob";
        load_cosine_distribution(histograms, input_histogram_path);
        load_deletion_distribution(deletion_probability, input_deletion_path);
    }

    int SpectralLibrary::search_target_decoy_specset_SLGF(SpectralLibrary &decoy,
                                                                         SpecSet searchable_spectra,
                                                                         float parentmz_tolerance,
                                                                         vector<Spectrum *> target_library_ptr,
                                                                         vector<Spectrum *> decoy_library_ptr,
                                                                         int scoring_method,
                                                                         MS2ScoringModel &model,
                                                                         vector<string> &ionsToExtract,
                                                                         string allIons,
                                                                         int do_score_threshold,
                                                                         float score_threshold,
                                                                         PeptideSpectrumMatchSet &output_psm_set){

        vector<int> accepted_fragmentation;
        accepted_fragmentation.push_back(Spectrum::FragType_CID);

        //Preextracting ions for target and decoy libraries
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(target_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(decoy_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        int searched_count = 0;

        #pragma omp parallel for num_threads(8) schedule(guided)
        for(int query_idx = 0; query_idx < searchable_spectra.size(); query_idx++){
            //if(searchable_spectra[query_idx].scan != 9730) continue;
            //if(searchable_spectra[query_idx].scan != 3506) continue;

            //Filtering in acceptable fragmentation types
            bool valid_fragmentation = false;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(searchable_spectra[query_idx].msFragType == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }

            if(!valid_fragmentation) continue;



            //cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\t";
            //cout<<"mslevel\t"<<searchable_spectra[query_idx].msLevel
            //cout<<searchable_spectra[query_idx].parentMass<<"\t"<<searchable_spectra[query_idx].parentMZ<<"\t"<<searchable_spectra[query_idx].parentCharge<<"\t";
            psmPtr targetdecoy_psm(new PeptideSpectrumMatch);
            int target_decoy_search = this->search_target_decoy_SLGF(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psm,
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons);


            if(target_decoy_search == 0){
                targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
                //cout<<"Scan\t"<<targetdecoy_psm->m_scanNum<<"\t"<<"Library IDX\t"<<targetdecoy_psm->m_dbIndex<<"\t"<<targetdecoy_psm->m_annotation<<endl;
                float match_score = targetdecoy_psm->m_score;
                //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
                if( ((do_score_threshold == 0)) || match_score > score_threshold){
                    //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                    #pragma omp critical
                    {
                        searched_count++;
                        cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                        output_psm_set.push_back(targetdecoy_psm);
                    }
                }

            }
            //cout<<endl;

            continue;
        }

        return 0;
    }


    //DW
    float peakTol = 0.02;
    float comparTol = 0.1;
    int tagLen = 3;
    int gap = 0;
    int tagNum = 20;        //D get top 20 tags from each search spectrum
    int TAG_COMMON = 1;
    bool TAG_FILTERING = true;
    pair<float, float> MOD_RANGE(-200.0, 100.0);
    //D pair<float, float> MOD_RANGE(0.0, 0.0);

    /*D
    vector< list <pair<string, float> > > extract_tag_lib(vector<Spectrum> *spec_lib, const char* his_file)
    {
        cout << " extract_tag_lib begin" << endl;
        vector< list< pair<string, float> > > result;       //D each vector for a library spectrum

        AAJumps jumps(1);
        //D
        FILE *fo = fopen(his_file, "w");
        for (int i = 0; i < spec_lib->size(); i++)
        {
            list<Tag> tags;
            //D ExtractTags((*spec_lib)[i], tags, peakTol, tagLen, gap, tagNum);
            vector<unsigned int> *f_heap_times = &(ExtractDenovoTags((*spec_lib)[i], tags, peakTol, tagLen, tagNum, false));
            //D cout << tags.size() << endl;

            string annotation = (*spec_lib)[i].psmList.front()->m_annotation;
            for (int j = 0; j < annotation.size(); j++)
                if (annotation[j] == 'L')
                    annotation[j] = 'I';        //D the denovo tag extracting doesn't support 'L'

            list< pair<string, float> > temp;       //D <tag sequence, tag prefix>
            for (list<Tag>::iterator it = tags.begin(); it != tags.end(); it++)
            {
                string strSequence;
                for (int c = 0; c < it->sequence.size(); c++)
                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
                //D cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << endl;

                if (annotation.find(strSequence) != string::npos)       //D check if the tag sequence actually appears on the correct annotation
                {
                    pair<string, float> p_temp;
                    p_temp.first = strSequence;
                    p_temp.second = it->flankingPrefix;
                    temp.push_back(p_temp);     //D havent't considered the score of the tag yet
                }
            }
            result.push_back(temp);

            if (i%1000 == 0)
                cout << "i = " << i << endl;
            //D
            for (int j = 0; j < f_heap_times->size(); j++)
                fprintf(fo, "%ld ", (*f_heap_times)[j]);
            fprintf(fo, "\n");
        }
        fclose(fo);

        return result;
    }
    D*/


    vector< list <pair<string, float> > > extract_tag_lib(vector<Spectrum> *spec_lib, const char* his_file)
    {
        cout << " extract_tag_lib begin" << endl;
        vector< list< pair<string, float> > > result;       //D each vector for a library spectrum

        AAJumps jumps(1);
        //D
        FILE *fo = fopen(his_file, "w");
        for (int i = 0; i < spec_lib->size(); i++)
        {
            list<Tag> tags;
            //D cout << endl << "spec " << i << " scan " << (*spec_lib)[i].scan << " ";
            //D cout << ExtractCorrectTags((*spec_lib)[i], tags, comparTol, tagLen, false) << endl;
            ExtractCorrectTags((*spec_lib)[i], tags, comparTol, tagLen, false, false);     //D no filtering for extracting correct tags
            //D cout << tags.size() << endl;

            list< pair<string, float> > temp;       //D <tag sequence, tag prefix>
            for (list<Tag>::iterator it = tags.begin(); it != tags.end(); it++)
            {
                string strSequence;
                for (int c = 0; c < it->sequence.size(); c++)
                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
                //D cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << endl;

                pair<string, float> p_temp;
                p_temp.first = strSequence;
                p_temp.second = it->flankingPrefix;
                temp.push_back(p_temp);     //D no consider the score of the tag

                //D cout << p_temp.first << " " << p_temp.second << "; ";
            }
            result.push_back(temp);

            if (i%1000 == 0)
                cout << "i = " << i << endl;
        }
        fclose(fo);

        return result;
    }


    /*DW SpectralLibrary::extract_search_tags ORIGINAL
    vector< list <pair<string, float> > > SpectralLibrary::extract_search_tags(SpecSet search_lib, const char* his_file)
    {
        cout << "extract_search_tags begin" << endl;
        vector< list< pair<string, float> > > result;

        AAJumps jumps(1);
        //D
        FILE *fo = fopen(his_file, "w");

        for (int i = 0; i < search_lib.size(); i++)
        {
            list<Tag> tags;
            //D ExtractTags(search_lib[i], tags, peakTol, tagLen, gap, tagNum);
            vector<unsigned int> *f_heap_times = &(ExtractDenovoTags(search_lib[i], tags, peakTol, tagLen, tagNum, true));      //D do filtering for search spectra
            //D cout << tags.size() << endl;

            list< pair<string, float> > temp;
            for (list<Tag>::iterator it = tags.begin(); it != tags.end(); it++)
            {
                string strSequence;
                for (int c = 0; c < it->sequence.size(); c++)
                    strSequence = strSequence + jumps.aaLetters[it->sequence[c]];
                pair<string, float> p_temp;
                p_temp.first = strSequence;
                p_temp.second = it->flankingPrefix;
                temp.push_back(p_temp);
                //D cout << strSequence.c_str() << " " << strSequence.length() << " " << it->score << endl;
            }
            result.push_back(temp);

            if (i%1000 == 0)
                cout << "i = " << i << endl;
            //D
            for (int j = 0; j < f_heap_times->size(); j++)
                fprintf(fo, "%ld ", (*f_heap_times)[j]);
            fprintf(fo, "\n");
        }
        fclose(fo);

        return result;
    }
    D*/


    //D SpectralLibrary::extract_search_tags modified for combining PRM and normal specs together
    vector< list <pair<string, float> > > SpectralLibrary::extract_search_tags(SpecSet search_lib, const char* his_file)
    {
        cout << "extract_search_tags begin" << endl;
        vector< list< pair<string, float> > > result;

        AAJumps jumps(1);



        FILE *fi = fopen("denovo_tags", "r");
        char line[256];


        map<int, list<pair<string, float> > > scan_tag_map;
        while (fgets(line, sizeof(line), fi)) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
            //S = line;
            //cout << line << endl;
            int scanID = atoi(line);
            fgets(line, sizeof(line), fi);
            //cout << line << endl;
            int N_tag = atoi(line);
            cout << scanID << " " << N_tag << endl;

            for (int i = 0; i < N_tag; i++)
            {
                fgets(line, sizeof(line), fi);
                string S = line;        // YAI 1769.869019 11.384000
                string tag_seq = S.substr(0, 3);    // YAI
                string remain = S.substr(4, S.length()-4);        // 1769.869019 11.384000
                int j;
                for (j = 0; j < remain.length(); j++)
                    if (remain[j] == ' ')
                        break;

                float mass = atof(remain.substr(0, j).c_str());
                cout << remain << " " << j << " " << mass << endl;

                if (scan_tag_map.count(scanID) == 0)
                    scan_tag_map[scanID].clear();

                pair<string, float> p_temp;
                p_temp.first = tag_seq;
                p_temp.second = mass;

                scan_tag_map[scanID].push_back(p_temp);
                cout << p_temp.first << " " << p_temp.second << endl;
            }

        }
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */

        fclose(fi);




        for (int i = 0; i < search_lib.size(); i++)
        {
            list< pair<string, float> > temp = scan_tag_map[search_lib[i].scan];
            result.push_back(temp);

            if (i%1000 == 0)
                cout << "i = " << i << endl;
        }

        return result;
    }


    int common_tag_num(list <pair<string, float> > query_tags, list <pair<string, float> > reference_tags)
    {
        int N = 0;
        for (list< pair<string, float> >::iterator it1 = query_tags.begin(); it1 != query_tags.end(); it1++)
            for (list< pair<string, float> >::iterator it2 = reference_tags.begin(); it2 != reference_tags.end(); it2++)
                if (((it1->first).compare(it2->first) == 0) && (fabs(it1->second-it2->second) <= comparTol))
                {
                    N++;
                }

        return N;
    }


    bool is_tag_common(list <pair<string, float> > query_tags, list <pair<string, float> > reference_tags)
    {
        int N = 0;
        for (list< pair<string, float> >::iterator it1 = query_tags.begin(); it1 != query_tags.end(); it1++)
            for (list< pair<string, float> >::iterator it2 = reference_tags.begin(); it2 != reference_tags.end(); it2++)
                if (((it1->first).compare(it2->first) == 0) && (fabs(it1->second-it2->second) <= comparTol))
                {
                    N++;
                    if (N >= TAG_COMMON)
                        return true;
                }

        return false;
    }


    int SpectralLibrary::search_target_decoy_SLGFNew2(vector< list <pair<string, float> > > target_tag_lib, vector< list <pair<string, float> > > decoy_tag_lib, list <pair<string, float> > query_tags,
                                             SpectralLibrary &decoy,        //D search_target_decoy_SLGF_NEW
                                             Spectrum query_spec,
                                             vector<psmPtr> & output_psms,
                                             int top_psm_number,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons,
                                             int abundance){

    DEBUG_MSG("SEARCHING SCAN "<<query_spec.scan);

    vector<psmPtr> search_results;

    float ISOTOPIC_MZ_ERROR_THRESHOLD = 1.2;

    float query_intensity = query_spec.getTotalIonCurrent();
    psmPtr query_psm(new PeptideSpectrumMatch());
    query_psm->m_spectrum = & query_spec;

    DEBUG_MSG("begin target" << endl);
    int target_pass_num = 0;
    for(int library_idx = 0; library_idx < this->size(); library_idx++){
        float library_mass = specs[library_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = specs[library_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;

        //DW
        //D list<Tag> tags;
        //D ExtractTags(specs[i], tags, peakTol, tagLen, gap, tagNum);
        //D ExtractTags(query_spec, tags, peakTol, tagLen, gap, tagNum);

        //D if (common_tag_num(query_tags, target_tag_lib[library_idx]) < TAG_COMMON)
        if (is_tag_common(query_tags, target_tag_lib[library_idx]) == false)
            continue;
        target_pass_num = target_pass_num + 1;

        float sim = 0.0;
        float percent_intensity = 0.f;

        query_psm->m_annotation = specs[library_idx].psmList.front()->m_annotation;
        DEBUG_MSG("query_psm->m_annotation = " << query_psm->m_annotation << endl);

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);
        DEBUG_MSG("create_deliminated_aminoacids successfully!" << endl);
        DEBUG_MSG("mz_difference = " << mz_difference << endl);

        float explained_intensity = 0.0;
        float rescored_sim = 0.0;

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){        //D SHOULD WE DO FABS(mz_difference) HERE???
            DEBUG_MSG("1 - spectrum_similarity_sqrt_librarypeaks" << endl);
            sim = spectrum_similarity_sqrt_librarypeaks(specs[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);      //D explained_intensity is calculated in spectrum_similarity_sqrt_librarypeaks
            DEBUG_MSG("sim = " << sim << endl);     //D sim is a cosine value
            DEBUG_MSG("specs[library_idx].psmList.front()->SLGF_distribution.size() = " << specs[library_idx].psmList.front()->SLGF_distribution.size() << endl);
            rescored_sim = SLGF_rescore(specs[library_idx].psmList.front()->SLGF_distribution, sim);
            DEBUG_MSG("rescored_sim = " << rescored_sim << endl);       //D rescored_sim is the sum of probabilities a random spectrum having the cosine less than the sim value
            //library_name = specs[library_idx].fileName;
            library_name = specs[library_idx].psmList.front()->m_spectrumFile;
            DEBUG_MSG("library_name = " << library_name << endl);
        }

        //D the if condition cannot happen together with the previous one
        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){      //D SHOULD WE DO FABS(mz_difference) HERE???
        if (target_isocombined.size() > 0){
            DEBUG_MSG("2 - spectrum_similarity_sqrt_librarypeaks_isocombine" << endl);
            DEBUG_MSG("target_isocombined.size() = " << target_isocombined.size() << endl);
            DEBUG_MSG("size() = " << target_isocombined[library_idx].psmList.size() << endl);
            DEBUG_MSG("check target_isocombined[library_idx]" << endl);
            //D library_idx is used for both library and iso_combine spectra?
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(target_isocombined[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);     //D explained_intensity IS NOT CHANGED HERE
            DEBUG_MSG("sim = " << sim << endl);
            DEBUG_MSG("target_isocombined[library_idx].psmList.front()->SLGF_distribution.size() = " << target_isocombined[library_idx].psmList.front()->SLGF_distribution.size() << endl);
            rescored_sim = SLGF_rescore(target_isocombined[library_idx].psmList.front()->SLGF_distribution, sim);
            DEBUG_MSG("rescored_sim = " << rescored_sim << endl);
            //library_name = target_isocombined[library_idx].fileName;
            library_name = target_isocombined[library_idx].psmList.front()->m_spectrumFile;
            DEBUG_MSG("library_name = " << library_name << endl);
        }}

        float final_score = rescored_sim*explained_intensity;       //D SSM_score in the paper


        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;      //D was already assigned to the library annotation
        search_result->m_spectrum = & query_spec;       //D query spectrum
        search_result->m_isDecoy = false;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D a cosine similarity score
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = library_idx + 1;     //D index plus 1
        search_result->m_library_name = get_only_filename(library_name);



        search_results.push_back(search_result);

        DEBUG_MSG("TARGET\t"<<final_score);
    }

    DEBUG_MSG("target_pass_num = " << target_pass_num << endl);
    DEBUG_MSG("end target" << endl);

    DEBUG_MSG("begin decoy" << endl);
    int decoy_pass_num = 0;
    for(int decoy_idx = 0; decoy_idx < decoy.size(); decoy_idx++){
        float library_mass = decoy[decoy_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = decoy[decoy_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;

        //D if (common_tag_num(query_tags, decoy_tag_lib[decoy_idx]) < TAG_COMMON)
        if (is_tag_common(query_tags, decoy_tag_lib[decoy_idx]) == false)
            continue;
        decoy_pass_num = decoy_pass_num + 1;

        float sim = 0.0;
        float percent_intensity = 0.f;


        query_psm->m_annotation = decoy[decoy_idx].psmList.front()->m_annotation;

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);        //D query_psm->m_annotation was already assigned to the decoy spectrum's annotation

        float explained_intensity = 0.0;
        float rescored_sim = 0.0;       //D note: the higher scores the better, all kind of these scores
        rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){
            sim = spectrum_similarity_sqrt_librarypeaks(decoy[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);        //D explained_intensity is the sum of squares of intensities
            rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy[decoy_idx].fileName;
            library_name = decoy[decoy_idx].psmList.front()->m_spectrumFile;
        }

        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){
        if (decoy_isocombined.size() > 0){
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(decoy_isocombined[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);
            rescored_sim = SLGF_rescore(decoy_isocombined[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy_isocombined[decoy_idx].fileName;
            library_name = decoy_isocombined[decoy_idx].psmList.front()->m_spectrumFile;
        }}

        float final_score = rescored_sim*explained_intensity;

        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;
        search_result->m_spectrum = & query_spec;
        search_result->m_isDecoy = true;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D cosine similarity value
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = decoy_idx + 1;
        search_result->m_library_name = get_only_filename(library_name);


        search_results.push_back(search_result);

        DEBUG_MSG("DECOY\t"<<final_score);
    }

    DEBUG_MSG("decoy_pass_num = " << decoy_pass_num << endl);
    DEBUG_MSG("end decoy" << endl);

    //D search_results are mixed by both target and decoy search results
    sort(search_results.begin(), search_results.end(), search_results_comparator_psmPtr);       //D sort descendingly by m_score (more specifically SSM_score)
    if(search_results.size() > 0){
        for(int i = 0; i < top_psm_number; i++){
            if(search_results.size() <= i) break;       //D get at most top_psm_number of psms
            output_psms.push_back(search_results[i]);
        }
        return 0;
    }


    return -1;
}


set<int> get_candidate_spectra(list <pair<string, float> > query_tags, map<string, vector<pair<float, int> > > tag_map)
{
    set<int> spec_idxs;
    for (list <pair<string, float> >::iterator it = query_tags.begin(); it != query_tags.end(); it++)
    {
        string tag_seq = it->first;
        float flankingPrefix = it->second;
        if (tag_map.count(tag_seq) > 0)
        {
            vector<pair<float, int> > temp = tag_map[tag_seq];
            int L = 0;      //D do binary search to find candidate spec Idxs
            int R = temp.size() - 1;
            int M;

            while (L <= R)
            {
                M = (L+R) / 2;
                if (temp[M].first > flankingPrefix + comparTol)
                    R = M - 1;
                else
                    if (temp[M].first < flankingPrefix - comparTol)
                        L = M + 1;
                    else
                        break;
            }

            if (L <= R)     //D there is at least a lib spec matched
            {
                int i = M;
                while ((i >= 0) and (temp[i].first >= flankingPrefix-comparTol))
                {
                    spec_idxs.insert(temp[i].second);       //D found a candidate spectrum ID
                    i--;
                }


                i = M + 1;
                while ((i < temp.size()) and (temp[i].first <= flankingPrefix+comparTol))
                {
                    spec_idxs.insert(temp[i].second);
                    i++;
                }
            }
        }
    }   //D end for iterator it

    return spec_idxs;
}


//D consider modification
set<int> get_candidate_spectra_MOD(list <pair<string, float> > query_tags, map<string, vector<pair<float, int> > > tag_map, Spectrum query_spec, vector<Spectrum> *spec_set, float parent_mass_tolerance)
{
    set<int> spec_idxs;
    for (list <pair<string, float> >::iterator it = query_tags.begin(); it != query_tags.end(); it++)
    {
        string tag_seq = it->first;
        float flankingPrefix = it->second;      //D flanking prefix mass of the query tag
        if (tag_map.count(tag_seq) > 0)
        {
            vector<pair<float, int> > temp = tag_map[tag_seq];
            int L = 0;      //D do binary search to find candidate spec Idxs
            int R = temp.size() - 1;
            int M;

            float upper_mass = flankingPrefix + parent_mass_tolerance + comparTol + MOD_RANGE.second;
            float lower_mass = flankingPrefix - parent_mass_tolerance - comparTol + MOD_RANGE.first;

            while (L <= R)
            {
                M = (L+R) / 2;
                if (temp[M].first > upper_mass)
                    R = M - 1;
                else
                    if (temp[M].first < lower_mass)
                        L = M + 1;
                    else
                        break;
            }

            if (L <= R)     //D there is at least a lib spec matched
            {
                int i = M;
                while ((i >= 0) and (temp[i].first >= lower_mass))
                {
                    //D for the charge issue
                    float mass_difference;
                    for (int query_charge = 2; query_charge <= 4; query_charge++)
                    {
                        float query_mass = query_spec.parentMZ*query_charge - (query_charge-1)*AAJumps::massHion;
                        mass_difference = query_mass - (*spec_set)[temp[i].second].parentMass;
                        //D if (fabs((query_mass-(*spec_set)[temp[i].second].parentMass) - (flankingPrefix-temp[i].first)) <= comparTol)
                        if (fabs(mass_difference - (flankingPrefix-temp[i].first)) <= parent_mass_tolerance + comparTol)
                        {
                            spec_idxs.insert(temp[i].second);       //D found a candidate spectrum ID
                            break;
                        }

                    }
                    i--;
                }


                i = M + 1;
                while ((i < temp.size()) and (temp[i].first <= upper_mass))
                {
                    //D for the charge issue
                    float mass_difference;
                    for (int query_charge = 2; query_charge <= 4; query_charge++)
                    {
                        float query_mass = query_spec.parentMZ*query_charge - (query_charge-1)*AAJumps::massHion;
                        mass_difference = query_mass - (*spec_set)[temp[i].second].parentMass;
                        //D if (fabs((query_mass-(*spec_set)[temp[i].second].parentMass) - (flankingPrefix-temp[i].first)) <= comparTol)
                        if (fabs(mass_difference - (flankingPrefix-temp[i].first)) <= parent_mass_tolerance + comparTol)
                        {
                            spec_idxs.insert(temp[i].second);       //D found a candidate spectrum ID
                            break;
                        }

                    }
                    i++;
                }
            }
        }
    }   //D end for iterator it

    return spec_idxs;
}


int SpectralLibrary::search_target_decoy_SLGFNew3(map<string, vector<pair<float, int> > > target_map, map<string, vector<pair<float, int> > > decoy_map, list <pair<string, float> > query_tags,
                                             SpectralLibrary &decoy,        //D search_target_decoy_SLGF_NEW
                                             Spectrum query_spec,
                                             vector<psmPtr> & output_psms,
                                             int top_psm_number,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons,
                                             int abundance){

    DEBUG_MSG("SEARCHING SCAN "<<query_spec.scan);

//    if ((query_spec.scan != 41558) && (query_spec.scan != 5282))
//        return -1;

    vector<psmPtr> search_results;

    float ISOTOPIC_MZ_ERROR_THRESHOLD = 1.2;

    Spectrum ori_query_spec = query_spec;

    float query_intensity = query_spec.getTotalIonCurrent();
    psmPtr query_psm(new PeptideSpectrumMatch());
    query_psm->m_spectrum = & query_spec;

    DEBUG_MSG("begin target" << endl);
    set<int> target_idxs = get_candidate_spectra_MOD(query_tags, target_map, query_spec, &specs, parentmz_tolerance);
    int target_pass_num = target_idxs.size();

    //D for(int library_idx = 0; library_idx < this->size(); library_idx++){
    for(set<int>::iterator it = target_idxs.begin(); it != target_idxs.end(); it++){
        int library_idx = *it;
        float library_mass = specs[library_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = specs[library_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        //D if(charge != query_charge) continue;

        float sim = 0.0;
        float percent_intensity = 0.f;

        query_psm->m_annotation = specs[library_idx].psmList.front()->m_annotation;
        DEBUG_MSG("query_psm->m_annotation = " << query_psm->m_annotation << endl);
        cout << query_psm->m_spectrum->scan << " query_psm->m_annotation = " << query_psm->m_annotation << endl;


        /*D
        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;
        D*/
        float mass_difference, mz_difference;
        for (query_charge = 2; query_charge <= 4; query_charge++)
        {
            float query_parent_mass = query_mass*query_charge - (query_charge-1)*AAJumps::massHion;
            float library_parent_mass = library_mass*charge - (charge-1)*AAJumps::massHion;
            mass_difference = query_parent_mass - library_parent_mass;
            if (fabs(mass_difference) <= parentmz_tolerance)
                break;
        }

        //D commented for the search with modification
        //D if (fabs(mass_difference) > parentmz_tolerance)
        //D    continue;
        mz_difference = (query_mass - library_mass);

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);
        DEBUG_MSG("create_deliminated_aminoacids successfully!" << endl);
        DEBUG_MSG("mz_difference = " << mz_difference << endl);

        //D float explained_intensity = 0.0;
        float explained_intensity = 1.0;
        float rescored_sim = 0.0;

        string library_name = "";

        DEBUG_MSG("1 - spectrum_similarity_sqrt_librarypeaks" << endl);
        //D sim = spectrum_similarity_sqrt_librarypeaks(specs[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);      //D explained_intensity is calculated in spectrum_similarity_sqrt_librarypeaks

        psmPtr ref_psm = specs[library_idx].psmList.front();

        vector <pair<float, float> > extracted_ion;
        ref_psm->annotate(ref_psm->m_annotation, allIons, model, 0, 0, .50);
        extractIons(ref_psm, deliminated_peptide.size(), model, ionsToExtract, extracted_ion, 0, 0);     //D extractIons(library,peptideLength,model,ionsToExtract,library->m_ion_extraction, 0, 0);

        sort(extracted_ion.begin(), extracted_ion.end(), mass_intensity_pair_mass_comp);

        Spectrum annotated_spec = specs[library_idx];
        annotated_spec.resize(extracted_ion.size());
        for (int i = 0; i < extracted_ion.size(); i++)
        {
            annotated_spec[i][0] = extracted_ion[i].first;
            annotated_spec[i][1] = extracted_ion[i].second;
        }

        unsigned int matchPeaks;
        float score1, score2;
        float temp = 0.5;
        sim = query_spec.scoreMatch3(annotated_spec, temp, matchPeaks, score1, score2, false, true);
        DEBUG_MSG("sim = " << sim << endl);     //D sim is a cosine value
        DEBUG_MSG("specs[library_idx].psmList.front()->SLGF_distribution.size() = " << specs[library_idx].psmList.front()->SLGF_distribution.size() << endl);
        rescored_sim = SLGF_rescore(specs[library_idx].psmList.front()->SLGF_distribution, sim);
        DEBUG_MSG("rescored_sim = " << rescored_sim << endl);       //D rescored_sim is the sum of probabilities a random spectrum having the cosine less than the sim value
        //library_name = specs[library_idx].fileName;
        library_name = specs[library_idx].psmList.front()->m_spectrumFile;
        DEBUG_MSG("library_name = " << library_name << endl);

        Spectrum temp_spec = *query_psm->m_spectrum;
        //D explained_intensity = preprocess_library_ion_extraction(query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons);

        *query_psm->m_spectrum = temp_spec;

        //D explained_intensity = sqrt(score1);
        explained_intensity = score1;

        DEBUG_MSG("explained_intensity = " << explained_intensity << endl);

        float final_score = rescored_sim*explained_intensity;       //D SSM_score in the paper

        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;      //D was already assigned to the library annotation
        search_result->m_spectrum = & query_spec;       //D query spectrum
        search_result->m_isDecoy = false;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D a cosine similarity score
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = library_idx + 1;     //D index plus 1
        search_result->m_library_name = get_only_filename(library_name);

        //D just added
        search_result->m_charge = query_spec.parentCharge;
        search_result->m_mz = query_spec.parentMZ;
        search_result->m_exactmass = query_spec.parentMass;
        //D search_result->m_parentmass_difference = query_spec.parentMass - specs[library_idx].parentMass;
        search_result->m_parentmass_difference = query_spec.parentMZ*specs[library_idx].parentCharge - (specs[library_idx].parentCharge-1)*AAJumps::massHion - specs[library_idx].parentMass;

        //D if ((fabs(search_result->m_parentmass_difference) <= parentmz_tolerance) || (deliminated_peptide.size() > 12))
        search_results.push_back(search_result);

        DEBUG_MSG("TARGET\t"<<final_score);
    }

    DEBUG_MSG("target_pass_num = " << target_pass_num << endl);
    DEBUG_MSG("end target" << endl);

    DEBUG_MSG("begin decoy" << endl);

    query_spec = ori_query_spec;

    set<int> decoy_idxs = get_candidate_spectra_MOD(query_tags, decoy_map, query_spec, &decoy.specs, parentmz_tolerance);
    int decoy_pass_num = decoy_idxs.size();

    //D for(int decoy_idx = 0; decoy_idx < decoy.size(); decoy_idx++){
    for(set<int>::iterator it = decoy_idxs.begin(); it != decoy_idxs.end(); it++){
        int decoy_idx = *it;
        float library_mass = decoy[decoy_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = decoy[decoy_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        //D if(charge != query_charge) continue;

        float sim = 0.0;
        float percent_intensity = 0.f;


        query_psm->m_annotation = decoy[decoy_idx].psmList.front()->m_annotation;
        DEBUG_MSG("query_psm->m_annotation = " << query_psm->m_annotation << endl);
        cout << query_psm->m_spectrum->scan << " query_psm->m_annotation = " << query_psm->m_annotation << endl;

        /*D
        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;
        D*/
        float mass_difference, mz_difference;
        for (query_charge = 2; query_charge <= 4; query_charge++)
        {
            float query_parent_mass = query_mass*query_charge - (query_charge-1)*AAJumps::massHion;
            float library_parent_mass = library_mass*charge - (charge-1)*AAJumps::massHion;
            mass_difference = query_parent_mass - library_parent_mass;
            if (fabs(mass_difference) <= parentmz_tolerance)
                break;
        }

        //D if (fabs(mass_difference) > parentmz_tolerance)
        //D    continue;
        mz_difference = (query_mass - library_mass);

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);        //D query_psm->m_annotation was already assigned to the decoy spectrum's annotation
        int peptideLength = deliminated_peptide.size();

        float explained_intensity = 1.0;
        float rescored_sim = 0.0;       //D note: the higher scores the better, all kind of these scores

        string library_name = "";

        //D sim = spectrum_similarity_sqrt_librarypeaks(decoy[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);        //D explained_intensity is the sum of squares of intensities

        psmPtr ref_psm = decoy[decoy_idx].psmList.front();

        vector <pair<float, float> > extracted_ion;
        ref_psm->annotate(ref_psm->m_annotation, allIons, model, 0, 0, .50);
        extractIons(ref_psm, deliminated_peptide.size(), model, ionsToExtract, extracted_ion, 0, 0);     //D extractIons(library,peptideLength,model,ionsToExtract,library->m_ion_extraction, 0, 0);

        sort(extracted_ion.begin(), extracted_ion.end(), mass_intensity_pair_mass_comp);

        Spectrum annotated_spec = decoy[decoy_idx];
        annotated_spec.resize(extracted_ion.size());
        for (int i = 0; i < extracted_ion.size(); i++)
        {
            annotated_spec[i][0] = extracted_ion[i].first;
            annotated_spec[i][1] = extracted_ion[i].second;
        }


        unsigned int matchPeaks;
        float score1, score2;
        float temp = 0.5;
        sim = query_spec.scoreMatch3(annotated_spec, temp, matchPeaks, score1, score2, false, true);
        rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);
        //library_name = decoy[decoy_idx].fileName;
        library_name = decoy[decoy_idx].psmList.front()->m_spectrumFile;

        Spectrum temp_spec = *query_psm->m_spectrum;
        //D explained_intensity = preprocess_library_ion_extraction(query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons);
        *query_psm->m_spectrum = temp_spec;

        //D explained_intensity = sqrt(score1);
        explained_intensity = score1;

        DEBUG_MSG("explained_intensity = " << explained_intensity << endl);

        float final_score = rescored_sim*explained_intensity;

        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;
        search_result->m_spectrum = & query_spec;
        search_result->m_isDecoy = true;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D cosine similarity value
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = decoy_idx + 1;
        search_result->m_library_name = get_only_filename(library_name);

        //D just added
        search_result->m_charge = query_spec.parentCharge;
        search_result->m_mz = query_spec.parentMZ;
        search_result->m_exactmass = query_spec.parentMass;
        //D search_result->m_parentmass_difference = query_spec.parentMass - decoy[decoy_idx].parentMass;
        search_result->m_parentmass_difference = query_spec.parentMZ*decoy[decoy_idx].parentCharge - (decoy[decoy_idx].parentCharge-1)*AAJumps::massHion - decoy[decoy_idx].parentMass;

        //D if ((fabs(search_result->m_parentmass_difference) <= parentmz_tolerance) || (deliminated_peptide.size() > 12))
        search_results.push_back(search_result);

        DEBUG_MSG("DECOY\t"<<final_score);
    }

    DEBUG_MSG("decoy_pass_num = " << decoy_pass_num << endl);
    DEBUG_MSG("end decoy" << endl);

    query_spec = ori_query_spec;        //D restore the query spec

    //D search_results are mixed by both target and decoy search results
    sort(search_results.begin(), search_results.end(), search_results_comparator_psmPtr);       //D sort descendingly by m_score (more specifically SSM_score)
    if(search_results.size() > 0){
        float top_mod_score = -1.0;
        float top_unmod_score = -1.0;

        for (vector<psmPtr>::iterator it = search_results.begin(); it != search_results.end(); it++)
        {
            psmPtr search_result = *it;
            if ((top_unmod_score < 0) && (fabs(search_result->m_parentmass_difference) <= parentmz_tolerance))
                top_unmod_score = search_result->m_score;
            if ((top_mod_score < 0) && (fabs(search_result->m_parentmass_difference) > parentmz_tolerance))
                top_mod_score = search_result->m_score;
        }

        if ((top_mod_score >= 0) && (top_unmod_score >= 0))
            if (top_mod_score < 1.1*top_unmod_score)        //D not better than 10%
            {
                vector<psmPtr> temp;
                for (vector<psmPtr>::iterator it = search_results.begin(); it != search_results.end(); it++)
                {
                    psmPtr search_result = *it;
                    if(fabs(search_result->m_parentmass_difference) <= parentmz_tolerance)
                        temp.push_back(search_result);
                }
                search_results = temp;
            }

        for(int i = 0; i < top_psm_number; i++){
            if(search_results.size() <= i) break;       //D get at most top_psm_number of psms
            output_psms.push_back(search_results[i]);
        }
        return 0;
    }


    return -1;
}


map<string, vector< pair<float, int> > > create_map(vector< list <pair<string, float> > > tag_lib)
{
    map<string, vector< pair<float, int> > > tag_map;
    for (int i = 0; i < tag_lib.size(); i++)
    {
        for (list< pair<string, float> >::iterator it = tag_lib[i].begin(); it != tag_lib[i].end(); it++)
        {
            string tag_seq = it->first;
            float flankingPrefix = it->second;
            if (tag_map.count(tag_seq) == 0)
                tag_map[tag_seq].resize(0);
            tag_map[tag_seq].push_back(make_pair(flankingPrefix, i));       //D (the prefix of this tag, the index of spec containing this tag)
        }
    }

    for (map<string, vector< pair<float, int> > >::iterator it = tag_map.begin(); it != tag_map.end(); it++)
    {
        string tag_seq = it->first;
        sort(tag_map[tag_seq].begin(), tag_map[tag_seq].end());
        for (int i = 0; i < tag_map[tag_seq].size(); i++)
            cout << std::fixed << std::setprecision(3) << tag_map[tag_seq][i].first << " " << tag_map[tag_seq][i].second << "; ";
        cout << endl;
    }

    return tag_map;
}


    int SpectralLibrary::search_target_decoy_specset_SLGF(SpectralLibrary &decoy,       //D decoy spectrum set
                                                                         SpecSet searchable_spectra,
                                                                         float parentmz_tolerance,
                                                                         vector<Spectrum *> target_library_ptr,
                                                                         vector<Spectrum *> decoy_library_ptr,
                                                                         SpectralLibrary &target_isocombined,
                                                                         SpectralLibrary &decoy_isocombined,
                                                                         vector<Spectrum *> target_isocombined_library_ptr,
                                                                         vector<Spectrum *> decoy_isocombined_library_ptr,
                                                                         int scoring_method,
                                                                         MS2ScoringModel &model,
                                                                         vector<string> &ionsToExtract,
                                                                         string allIons,
                                                                         int do_score_threshold,
                                                                         float score_threshold,
                                                                         PeptideSpectrumMatchSet &output_psm_set,
                                                                         int output_psm_count,
                                                                         int abundance,     //D abundance = 1 low, = 2 high
                                                                         int start_seach_idx,
                                                                         int end_search_idx){



        vector<int> accepted_fragmentation;     //D CID = 0, HCD = 2, None = 100
        accepted_fragmentation.push_back(Spectrum::FragType_CID);
        accepted_fragmentation.push_back(Spectrum::FragType_HCD);      //D consider to remove
        accepted_fragmentation.push_back(Spectrum::FragType_NONE);      //D consider to remove
        DEBUG_MSG("accepted_fragmentation.size() = " << accepted_fragmentation.size() << endl);

        //D Preextracting ions for target and decoy libraries
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
           //target_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            //decoy_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        //D Preextracting ions for target and decoy libraries iso combined
        for(int lib_idx = 0; lib_idx < target_isocombined_library_ptr.size(); lib_idx++){
            //target_isocombined_library_ptr[lib_idx]->setResolution(1.0005, true);
        }

        for(int lib_idx = 0; lib_idx < decoy_isocombined_library_ptr.size(); lib_idx++){
            //decoy_isocombined_library_ptr[lib_idx]->setResolution(1.0005, true);
        }




        //D Preextracting ions for target and decoy libraries
        DEBUG_MSG("target_library_ptr.size() = " << target_library_ptr.size() << endl);
        for(int lib_idx = 0; lib_idx < target_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction(target_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        DEBUG_MSG("decoy_library_ptr.size() = " << decoy_library_ptr.size() << endl);
        for(int lib_idx = 0; lib_idx < decoy_library_ptr.size(); lib_idx++){
            //D DEBUG_MSG("1- " << create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size() << endl);
            //D DEBUG_MSG("2- " << create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size() << endl);
            //D DEBUG_MSG(decoy_library_ptr[lib_idx]->psmList.front()->m_spectrum->size());

            preprocess_library_ion_extraction(decoy_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        //Preextracting ions for target and decoy libraries iso combined
        if (true)
        {
        DEBUG_MSG("target_isocombined_library_ptr.size() = " << target_isocombined_library_ptr.size() << endl);
        for(int lib_idx = 0; lib_idx < target_isocombined_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction_isocombine(target_isocombined_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(target_isocombined_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        }

        if (true)
        {
        DEBUG_MSG("decoy_isocombined_library_ptr.size() = " << decoy_isocombined_library_ptr.size() << endl);
        for(int lib_idx = 0; lib_idx < decoy_isocombined_library_ptr.size(); lib_idx++){
            preprocess_library_ion_extraction_isocombine(decoy_isocombined_library_ptr[lib_idx]->psmList.front(),
                          create_deliminated_aminoacids(decoy_isocombined_library_ptr[lib_idx]->psmList.front()->m_annotation).size(),
                          model,
                          ionsToExtract,
                          allIons);
        }

        }


        DEBUG_MSG("start_seach_idx = " << start_seach_idx << endl);
        DEBUG_MSG("end_search_idx = " << end_search_idx << endl);
        DEBUG_MSG("searchable_spectra.size() = " << searchable_spectra.size() << endl);
        cout << "start_seach_idx = " << start_seach_idx << endl;
        cout << "end_search_idx = " << end_search_idx << endl;
        cout << "searchable_spectra.size() = " << searchable_spectra.size() << endl;

        time_t now = time(0);
        char* dt = ctime(&now);
        cout << "time before filtering = " << dt << endl;

        //D tag filtering
        vector< list <pair<string, float> > > target_tag_lib;
        vector< list <pair<string, float> > > decoy_tag_lib;
        vector< list <pair<string, float> > > search_tag_lib;
        map<string, vector<pair<float, int> > > lib_tag_map;        //D <tag_sequence, <prefix, specID> >
        map<string, vector<pair<float, int> > > decoy_tag_map;

        if (TAG_FILTERING)
        {
            cout << "begin extract lib" << endl;
            const char* target_file = "target_histogram";
            target_tag_lib = extract_tag_lib(&specs, target_file);
            cout << "end extract lib" << endl;
            now = time(0);
            dt = ctime(&now);
            cout << "end tag lib extracting time = " << dt << endl;

            cout << "begin extract decoy lib" << endl;
            const char* decoy_file = "decoy_histogram";
            decoy_tag_lib = extract_tag_lib(&decoy.specs, decoy_file);
            cout << "end extract decoy lib" << endl;
            now = time(0);
            dt = ctime(&now);
            cout << "end decoy lib extracting time = " << dt << endl;

            cout << "begin extract search tags" << endl;
            const char* search_file = "search_histogram";
            search_tag_lib = extract_search_tags(searchable_spectra, search_file);
            cout << "end extract search tags" << endl;
            now = time(0);
            dt = ctime(&now);
            cout << "end search tag extracting time = " << dt << endl;

            lib_tag_map = create_map(target_tag_lib);
            decoy_tag_map = create_map(decoy_tag_lib);
            //D map<string, vector<pair<float, int> > > search_tags = create_map(search_tag_lib);
            now = time(0);
            dt = ctime(&now);
            cout << "end create map time = " << now << endl;
        }

        //Dcout << "SEARCH NOW (double) clock() = " << (double) clock() << endl;
        int searched_count = 0;     //D the number of search spectra which have at least a psm
        #pragma omp parallel for num_threads(1) schedule(guided)
        for(int query_idx = start_seach_idx; query_idx < end_search_idx; query_idx++){
            //Filtering on abundance
            DEBUG_MSG("query_idx = " << query_idx << endl);
            int TIC = searchable_spectra[query_idx].getTotalIonCurrent();
            DEBUG_MSG("TIC = " << TIC << endl);
            /*D
            if(abundance == 1){//Low Abundance
                if(TIC > 12000){
                    continue;
                }
            }
            if(abundance == 2){//High Abundance
                if(TIC <= 12000){
                    continue;
                }
            }
            D*/

            //Filtering in acceptable fragmentation types
            bool valid_fragmentation = false;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                //DEBUG_MSG("valid fragmentation! " << searchable_spectra[query_idx].msFragType << "\t" << accepted_fragmentation[fragmentation_idx] << endl);
                if(searchable_spectra[query_idx].msFragType == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }

            if(!valid_fragmentation) continue;


            //cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\t";
            //cout<<"mslevel\t"<<searchable_spectra[query_idx].msLevel
            //cout<<searchable_spectra[query_idx].parentMass<<"\t"<<searchable_spectra[query_idx].parentMZ<<"\t"<<searchable_spectra[query_idx].parentCharge<<"\t";

            /*
            psmPtr targetdecoy_psm(new PeptideSpectrumMatch);



            int target_decoy_search = this->search_target_decoy_SLGF(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psm,
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons);




            if(target_decoy_search == 0){
                targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
                //cout<<"Scan\t"<<targetdecoy_psm->m_scanNum<<"\t"<<"Library IDX\t"<<targetdecoy_psm->m_dbIndex<<"\t"<<targetdecoy_psm->m_annotation<<endl;
                float match_score = targetdecoy_psm->m_score;
                //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
                if( ((do_score_threshold == 0)) || match_score > score_threshold){
                    //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                    #pragma omp critical
                    {
                        searched_count++;
                        cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                        output_psm_set.push_back(targetdecoy_psm);
                    }
                }

            }*/


           vector<psmPtr> targetdecoy_psms;

           //D two options, either doing tag filtering or not
           int target_decoy_search;
           if (TAG_FILTERING == false)
                target_decoy_search = this->search_target_decoy_SLGFNew(decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psms,       //D output psms
                                                                output_psm_count,       //D the number of top psms from search_target_decoy_SLGFNew
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons,
                                                                abundance);
           //D no tag filtering
           else
                if (TAG_COMMON == 1)        //D can do union of the sets of candidate specIDs
                    target_decoy_search = this->search_target_decoy_SLGFNew3(lib_tag_map, decoy_tag_map, search_tag_lib[query_idx],
                                                                decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psms,       //D output psms
                                                                output_psm_count,       //D the number of top psms from search_target_decoy_SLGFNew
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons,
                                                                abundance);
                else
                    target_decoy_search = this->search_target_decoy_SLGFNew2(target_tag_lib, decoy_tag_lib, search_tag_lib[query_idx],
                                                                decoy,
                                                                searchable_spectra[query_idx],
                                                                targetdecoy_psms,       //D output psms
                                                                output_psm_count,       //D the number of top psms from search_target_decoy_SLGFNew
                                                                parentmz_tolerance,
                                                                target_library_ptr,
                                                                decoy_library_ptr,
                                                                target_isocombined,
                                                                decoy_isocombined,
                                                                target_isocombined_library_ptr,
                                                                decoy_isocombined_library_ptr,
                                                                scoring_method,
                                                                model,
                                                                ionsToExtract,
                                                                allIons,
                                                                abundance);


            if(target_decoy_search == 0){       //D found at least one psm from search_target_decoy_SLGFNew
                for(int i = 0; i < targetdecoy_psms.size(); i++){
                    targetdecoy_psms[i]->m_scanNum = searchable_spectra[query_idx].scan;
                    float match_score = targetdecoy_psms[i]->m_score;       //D use m_score for filtering
                    if( ((do_score_threshold == 0)) || match_score > score_threshold){      //D accept the psm if the SSM_score is strictly greater than the score threshold or no threshold is set
                        #pragma omp critical
                        {
                            searched_count++;
                            cout<<"search_target_decoy_specset_SLGF - Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<searched_count<<"\tof\t"<<searchable_spectra.size()<<endl;
                            output_psm_set.push_back(targetdecoy_psms[i]);
                        }
                    }
                }
            }

            continue;
        }

        now = time(0);
        dt = ctime(&now);
        cout << "end search time = " << dt << endl;
        cout << "Done: int SpectralLibrary::search_target_decoy_specset_SLGF" << endl;

        return 0;
    }


    int SpectralLibrary::search_target_decoy_SLGF(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        float query_intensity = query_spec.getTotalIonCurrent();

        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        //cout<<decoy_start_search_idx<<"\t"<<decoy_end_search_idx<<endl;
        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;

            float sim;// = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;



            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(decoy_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(decoy_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(decoy_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);

                full_cos = sim*explained_intensity;
                orig_dot = sim;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                //sim = rescored_sim;
                full_cos = sim;


                percent_intensity = explained_intensity;

            }

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple_decoy[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple_decoy[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple_decoy[i]);;
            search_results_decoy.push_back(psm);
        }




        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);


        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            float sim;// = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;


            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = target_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(target_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(target_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(target_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);

                //DEBUG_MSG("LIBRARY\t"<<target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim);

                //cout<<rescored_sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim*explained_intensity<<endl;

                orig_dot = sim;
                full_cos = sim*explained_intensity;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                //sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;


            }
            //DEBUG_MSG(target_library_ptr[library_idx]->psmList.front()->m_dbIndex<<"\t"<<fabs(query_spec.parentMZ - library_mass)<<"\t"<<(string)target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<orig_dot<<"\t"<<percent_intensity<<"\t"<<sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best target
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple[i]);;
            search_results.push_back(psm);

            //DEBUG_MSG(psm->m_dbIndex<<"\t"<<psm->m_score<<"\t"<<psm->m_annotation);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }

        //DEBUG_MSG("TARGET\t"<<target_top_scoring<<"\t"<<search_results[0]->m_strict_envelope_score<<"\t"<<search_results[0]->m_unstrict_envelope_score<<"\t"<<search_results[0]->m_annotation);
        //DEBUG_MSG("DECOY\t"<<decoy_top_scoring<<"\t"<<search_results_decoy[0]->m_strict_envelope_score<<"\t"<<search_results_decoy[0]->m_unstrict_envelope_score<<"\t"<<search_results_decoy[0]->m_annotation);

        if(target_top_scoring >= decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results[0]->m_strict_envelope_score;
            float explained_int = search_results[0]->m_unstrict_envelope_score;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<"\t"<<orig_cos<<endl;


            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }


            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = false;
            output_psm->m_fdr = search_results[0]->m_fdr;
            output_psm->m_startMass = deltaD;
            output_psm->m_charge = query_spec.parentCharge;
            output_psm->m_mz = query_spec.parentMZ;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results_decoy[0]->m_strict_envelope_score;
            float explained_int = search_results_decoy[0]->m_unstrict_envelope_score;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = true;
            output_psm->m_fdr = search_results_decoy[0]->m_fdr;
            output_psm->m_startMass = deltaD;
            output_psm->m_charge = query_spec.parentCharge;
            output_psm->m_mz = query_spec.parentMZ;

            return 0;
        }

        return -1;
    }

    //Searching with two models
    int SpectralLibrary::search_target_decoy_SLGF(SpectralLibrary &decoy,
                                             Spectrum query_spec,
                                             psmPtr output_psm,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;

        float query_intensity = query_spec.getTotalIonCurrent();

        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);

        float isotopic_mz_error_threshold = 1.2;

        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        //cout<<decoy_start_search_idx<<"\t"<<decoy_end_search_idx<<endl;
        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){
            //DEBUG_MSG(decoy_library_ptr[library_idx]->parentMZ<<"\t"<<decoy_isocombined_library_ptr[library_idx]->parentMZ);

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;
            float mz_error = fabs(query_spec.parentMZ - library_mass);

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;

            float sim;// = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;

            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(decoy_library_ptr[library_idx]->psmList.front()->m_annotation);



                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(decoy_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons,explained_intensity);


                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(decoy_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                if(mz_error > isotopic_mz_error_threshold){
                    sim = spectrum_similarity_sqrt_librarypeaks_isocombine(decoy_isocombined_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);
                    rescored_sim = SLGF_rescore(decoy_isocombined_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                }



                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);
                full_cos = sim*explained_intensity;
                orig_dot = sim;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;

            }

            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<7>(similarity_tuple) = decoy_library_ptr[library_idx]->psmList.front();
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple_decoy.push_back(similarity_tuple);
        }



        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple_decoy[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple_decoy[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple_decoy[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple_decoy[i]);;
            search_results_decoy.push_back(psm);
        }




        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;

        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);



        for(int library_idx = target_start_search_idx; library_idx <= target_end_search_idx; library_idx++){

            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;
            float mz_error = fabs(query_spec.parentMZ - library_mass);

            if(fabs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            float sim;// = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            float dot_bias;// = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);
            float orig_dot = 0.f;
            float percent_intensity = 0.f;
            float full_cos = 0.f;

            //For rescoring the similarity as a pvalue
            if(scoring_method == SpectralLibrary::MatchScoreType_DotProduct_SLGF){
                psmPtr query_psm(new PeptideSpectrumMatch());
                query_psm->m_spectrum = & query_spec;
                query_psm->m_annotation = target_library_ptr[library_idx]->psmList.front()->m_annotation;
                vector<string> deliminated_aminoacids = create_deliminated_aminoacids(target_library_ptr[library_idx]->psmList.front()->m_annotation);

                float explained_intensity;
                sim = spectrum_similarity_sqrt_librarypeaks(target_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);

                //sim = sim * 1.01;

                float rescored_sim = SLGF_rescore(target_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                if(mz_error > isotopic_mz_error_threshold){
                    sim = spectrum_similarity_sqrt_librarypeaks_isocombine(target_isocombined_library_ptr[library_idx]->psmList.front(), query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons, explained_intensity);
                    rescored_sim = SLGF_rescore(target_isocombined_library_ptr[library_idx]->psmList.front()->SLGF_distribution, sim);
                }

                //explained_intensity = percent_explained_intensity(query_psm, deliminated_aminoacids.size(), model, ionsToExtract, allIons);
                //cout<<rescored_sim<<"\t"<<explained_intensity<<"\t"<<rescored_sim*explained_intensity<<endl;

                orig_dot = sim;
                full_cos = sim*explained_intensity;
                //sim = rescored_sim*explained_intensity*explained_intensity;
                sim = rescored_sim*explained_intensity;
                sim = rescored_sim;
                full_cos = sim;

                percent_intensity = explained_intensity;


            }
            //DEBUG_MSG(target_library_ptr[library_idx]->psmList.front()->m_dbIndex<<"\t"<<fabs(query_spec.parentMZ - library_mass)<<"\t"<<(string)target_library_ptr[library_idx]->psmList.front()->m_annotation<<"\t"<<orig_dot<<"\t"<<percent_intensity<<"\t"<<sim);

            score_results_tuple similarity_tuple;

            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = full_cos;
            tr1::get<2>(similarity_tuple) = target_library_ptr[library_idx]->psmList.front()->m_dbIndex;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            tr1::get<5>(similarity_tuple) = orig_dot;
            tr1::get<6>(similarity_tuple) = percent_intensity;
            tr1::get<8>(similarity_tuple) = sim;
            scores_tuple.push_back(similarity_tuple);
        }


        //Finding the best target
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_dbIndex = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            psm->m_strict_envelope_score = tr1::get<5>(scores_tuple[i]);
            psm->m_unstrict_envelope_score = tr1::get<6>(scores_tuple[i]);
            psm->m_pValue = query_intensity;
            psm->m_fdr = tr1::get<8>(scores_tuple[i]);;
            search_results.push_back(psm);

            //DEBUG_MSG(psm->m_dbIndex<<"\t"<<psm->m_score<<"\t"<<psm->m_annotation);
        }

        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;

        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }

        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }

        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }

        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }

        //DEBUG_MSG("TARGET\t"<<target_top_scoring<<"\t"<<search_results[0]->m_strict_envelope_score<<"\t"<<search_results[0]->m_unstrict_envelope_score<<"\t"<<search_results[0]->m_annotation);
        //DEBUG_MSG("DECOY\t"<<decoy_top_scoring<<"\t"<<search_results_decoy[0]->m_strict_envelope_score<<"\t"<<search_results_decoy[0]->m_unstrict_envelope_score<<"\t"<<search_results_decoy[0]->m_annotation);

        if(target_top_scoring >= decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results[0]->m_strict_envelope_score;
            float explained_int = search_results[0]->m_unstrict_envelope_score;

            //cout<<"Dot Bias\t"<<dot_bias<<"\t"<<search_results[0]->m_annotation<<"\tdot\t"<<dot_product<<"\tdeltaD\t"<<deltaD<<"\t"<<orig_cos<<endl;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }


            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_dbIndex       = search_results[0]->m_dbIndex;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = false;
            output_psm->m_fdr = search_results[0]->m_fdr;
            output_psm->m_startMass = deltaD;

            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            float orig_cos = search_results_decoy[0]->m_strict_envelope_score;
            float explained_int = search_results_decoy[0]->m_unstrict_envelope_score;

            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_SLGF:
                    match_score = dot_product;
                    break;
                default:
                    break;
            }

            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_dbIndex       = search_results_decoy[0]->m_dbIndex;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_strict_envelope_score =  orig_cos;
            output_psm->m_unstrict_envelope_score =  explained_int;
            output_psm->m_pValue = query_intensity;
            output_psm->m_isDecoy = true;
            output_psm->m_fdr = search_results_decoy[0]->m_fdr;
            output_psm->m_startMass = deltaD;

            return 0;
        }

        return -1;
    }


    int SpectralLibrary::search_target_decoy_SLGFNew(SpectralLibrary &decoy,        //D search_target_decoy_SLGF_NEW
                                             Spectrum query_spec,
                                             vector<psmPtr> & output_psms,
                                             int top_psm_number,
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr,
                                             vector<Spectrum *> decoy_library_ptr,
                                             SpectralLibrary &target_isocombined,
                                             SpectralLibrary &decoy_isocombined,
                                             vector<Spectrum *> target_isocombined_library_ptr,
                                             vector<Spectrum *> decoy_isocombined_library_ptr,
                                             int scoring_method,
                                             MS2ScoringModel &model,
                                             vector<string> &ionsToExtract,
                                             string allIons,
                                             int abundance){

    DEBUG_MSG("SEARCHING SCAN "<<query_spec.scan);

    vector<psmPtr> search_results;

    float ISOTOPIC_MZ_ERROR_THRESHOLD = 1.2;

    float query_intensity = query_spec.getTotalIonCurrent();
    psmPtr query_psm(new PeptideSpectrumMatch());
    query_psm->m_spectrum = & query_spec;

    DEBUG_MSG("begin target" << endl);
    for(int library_idx = 0; library_idx < this->size(); library_idx++){
        float library_mass = specs[library_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = specs[library_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;


        float sim = 0.0;
        float percent_intensity = 0.f;

        query_psm->m_annotation = specs[library_idx].psmList.front()->m_annotation;
        DEBUG_MSG("query_psm->m_annotation = " << query_psm->m_annotation << endl);

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);
        DEBUG_MSG("create_deliminated_aminoacids successfully!" << endl);
        DEBUG_MSG("mz_difference = " << mz_difference << endl);

        float explained_intensity = 0.0;
        float rescored_sim = 0.0;

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){        //D SHOULD WE DO FABS(mz_difference) HERE???
            DEBUG_MSG("1 - spectrum_similarity_sqrt_librarypeaks" << endl);
            sim = spectrum_similarity_sqrt_librarypeaks(specs[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);      //D explained_intensity is calculated in spectrum_similarity_sqrt_librarypeaks
            DEBUG_MSG("sim = " << sim << endl);     //D sim is a cosine value
            DEBUG_MSG("specs[library_idx].psmList.front()->SLGF_distribution.size() = " << specs[library_idx].psmList.front()->SLGF_distribution.size() << endl);
            rescored_sim = SLGF_rescore(specs[library_idx].psmList.front()->SLGF_distribution, sim);
            DEBUG_MSG("rescored_sim = " << rescored_sim << endl);       //D rescored_sim is the sum of probabilities a random spectrum having the cosine less than the sim value
            //library_name = specs[library_idx].fileName;
            library_name = specs[library_idx].psmList.front()->m_spectrumFile;
            DEBUG_MSG("library_name = " << library_name << endl);
        }

        //D the if condition cannot happen together with the previous one
        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){      //D SHOULD WE DO FABS(mz_difference) HERE???
        if (target_isocombined.size() > 0){
            DEBUG_MSG("2 - spectrum_similarity_sqrt_librarypeaks_isocombine" << endl);
            DEBUG_MSG("target_isocombined.size() = " << target_isocombined.size() << endl);
            DEBUG_MSG("size() = " << target_isocombined[library_idx].psmList.size() << endl);
            DEBUG_MSG("check target_isocombined[library_idx]" << endl);
            //D library_idx is used for both library and iso_combine spectra?
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(target_isocombined[library_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);     //D explained_intensity IS NOT CHANGED HERE
            DEBUG_MSG("sim = " << sim << endl);
            DEBUG_MSG("target_isocombined[library_idx].psmList.front()->SLGF_distribution.size() = " << target_isocombined[library_idx].psmList.front()->SLGF_distribution.size() << endl);
            rescored_sim = SLGF_rescore(target_isocombined[library_idx].psmList.front()->SLGF_distribution, sim);
            DEBUG_MSG("rescored_sim = " << rescored_sim << endl);
            //library_name = target_isocombined[library_idx].fileName;
            library_name = target_isocombined[library_idx].psmList.front()->m_spectrumFile;
            DEBUG_MSG("library_name = " << library_name << endl);
        }}

        float final_score = rescored_sim*explained_intensity;       //D SSM_score in the paper


        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;      //D was already assigned to the library annotation
        search_result->m_spectrum = & query_spec;       //D query spectrum
        search_result->m_isDecoy = false;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D a cosine similarity score
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = library_idx + 1;     //D index plus 1
        search_result->m_library_name = get_only_filename(library_name);



        search_results.push_back(search_result);

        DEBUG_MSG("TARGET\t"<<final_score);
    }
    DEBUG_MSG("end target" << endl);

    DEBUG_MSG("begin decoy" << endl);
    for(int decoy_idx = 0; decoy_idx < decoy.size(); decoy_idx++){
        float library_mass = decoy[decoy_idx].parentMZ;
        float query_mass = query_spec.parentMZ;

        int charge = decoy[decoy_idx].parentCharge;
        int query_charge = query_spec.parentCharge;

        if(charge != query_charge) continue;


        float sim = 0.0;
        float percent_intensity = 0.f;


        query_psm->m_annotation = decoy[decoy_idx].psmList.front()->m_annotation;

        float mass_difference = (query_mass - library_mass) * charge;
        float mz_difference = (query_mass - library_mass);
        if(fabs(mz_difference) > parentmz_tolerance) continue;

        vector<string> deliminated_peptide = create_deliminated_aminoacids(query_psm->m_annotation);        //D query_psm->m_annotation was already assigned to the decoy spectrum's annotation

        float explained_intensity = 0.0;
        float rescored_sim = 0.0;       //D note: the higher scores the better, all kind of these scores
        rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);

        string library_name = "";

        if(abundance == 1 || (abundance == 2 && mz_difference < ISOTOPIC_MZ_ERROR_THRESHOLD) || abundance == 0){
            sim = spectrum_similarity_sqrt_librarypeaks(decoy[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons,explained_intensity);        //D explained_intensity is the sum of squares of intensities
            rescored_sim = SLGF_rescore(decoy[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy[decoy_idx].fileName;
            library_name = decoy[decoy_idx].psmList.front()->m_spectrumFile;
        }

        if(abundance == 2 && mz_difference > ISOTOPIC_MZ_ERROR_THRESHOLD){
        if (decoy_isocombined.size() > 0){
            sim = spectrum_similarity_sqrt_librarypeaks_isocombine(decoy_isocombined[decoy_idx].psmList.front(), query_psm, deliminated_peptide.size(), model, ionsToExtract, allIons, explained_intensity);
            rescored_sim = SLGF_rescore(decoy_isocombined[decoy_idx].psmList.front()->SLGF_distribution, sim);
            //library_name = decoy_isocombined[decoy_idx].fileName;
            library_name = decoy_isocombined[decoy_idx].psmList.front()->m_spectrumFile;
        }}

        float final_score = rescored_sim*explained_intensity;

        psmPtr search_result(new PeptideSpectrumMatch());
        search_result->m_annotation = query_psm->m_annotation;
        search_result->m_spectrum = & query_spec;
        search_result->m_isDecoy = true;
        search_result->m_score = final_score;       //D SSM_score in the paper
        search_result->m_pValue = sim;      //D cosine similarity value
        search_result->m_spectrumFile = query_spec.fileName;
        search_result->m_dbIndex = decoy_idx + 1;
        search_result->m_library_name = get_only_filename(library_name);


        search_results.push_back(search_result);

        DEBUG_MSG("DECOY\t"<<final_score);
    }
    DEBUG_MSG("end decoy" << endl);

    //D search_results are mixed by both target and decoy search results
    sort(search_results.begin(), search_results.end(), search_results_comparator_psmPtr);       //D sort descendingly by m_score (more specifically SSM_score)
    if(search_results.size() > 0){
        for(int i = 0; i < top_psm_number; i++){
            if(search_results.size() <= i) break;       //D get at most top_psm_number of psms
            output_psms.push_back(search_results[i]);
        }
        return 0;
    }


    return -1;
}

}
