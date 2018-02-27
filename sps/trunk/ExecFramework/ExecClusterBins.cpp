//
// Created by Benjamin Pullman on 7/11/16.
//

/*
 * ExecKLFilter.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: isna
 */

#include "ExecClusterBins.h"
#include "ClusterSet.h"
#include <tr1/unordered_set>
#include <tr1/unordered_map>

namespace specnets {

    enum MatchType {
        None = 0,
        Subset = 1,
        Superset = 2,
        Exact = 3
    };

    //    struct Match
    //    {
    //        MatchType match_type;
    //        float cosine;
    //        float proj_cosine;
    //    };

    struct ClusterMember {
        Spectrum* spec;
        float similarity;
        unsigned int count_phase2;
        unsigned int count_phase3;
        unsigned int count_phase4;

    };

    struct DisCluster {
        vector<ClusterMember> core;
        vector<ClusterMember> sat;
        vector<ClusterMember> wsat;
        vector<ClusterMember> mix;
        vector<pair<unsigned int, unsigned int> > sub_clusters;
    };

    string normalize_peptide(string &peptide) {
        string normalized;
        normalized.reserve(peptide.size());
        string::iterator c;
        for (c = peptide.begin(); c != peptide.end(); c++) {
            if (isalpha(*c)) {
                if (*c == 'L') {
                    *c = 'I';
                }
                normalized += *c;
            }
        }
        return normalized;
    }

    void set_purity(vector<ClusterMember> &cluster_member,
            tr1::unordered_map<string, tr1::unordered_map<unsigned int, string> > &identifications,
            tr1::unordered_map<string, unsigned int> &peptide_counts,
            tr1::unordered_map<string, string> &original_peptide,
            vector<pair<string, string> > &mix_pairs,
            unsigned int &set_id_spectra,
            unsigned int &set_rep_spectra,
            string &representative
            ) {


        vector<ClusterMember>::iterator s;
        for (s = cluster_member.begin(); s != cluster_member.end(); ++s) {
            string previous_mix = "";
            vector<string> peptides;
            vector<string>::iterator peptide;
            stringSplit(identifications[s->spec->fileName][s->spec->scan], peptides, "!");
            for (peptide = peptides.begin(); peptide != peptides.end(); peptide++) {
                string original = *peptide;
                if (original != "PEPTIDE") {
                    string normalize = normalize_peptide(original);
                    if (normalize != previous_mix) {
                        if (previous_mix != "" && representative == "") {
                            mix_pairs.push_back(make_pair(previous_mix, normalize));
                        }
                        previous_mix = normalize;
                        if (peptide_counts.count(normalize) == 0) {
                            original_peptide[normalize] = original;
                            peptide_counts[normalize] = 1;
                        } else {
                            peptide_counts[normalize]++;
                        }
                    }
                }
            }
        }

        unsigned int representative_count = 0;

        if (representative == "") {
            tr1::unordered_map<string, unsigned int>::iterator peptide;
            for (peptide = peptide_counts.begin(); peptide != peptide_counts.end(); ++peptide) {
                if (peptide->second > representative_count) {
                    representative = peptide->first;
                    representative_count = peptide->second;
                }
            }
        }

        for (s = cluster_member.begin(); s != cluster_member.end(); ++s) {
            bool hasID = false;
            bool match = false;
            vector<string> peptides;
            stringSplit(identifications[s->spec->fileName][s->spec->scan], peptides, "!");
            vector<string>::iterator peptide;
            for (peptide = peptides.begin(); peptide != peptides.end(); peptide++) {
                string original = *peptide;
                if (original != "PEPTIDE") {
                    hasID = true;
                    string normalize = normalize_peptide(original);
                    vector<pair<string, string> >::iterator mixture;
                    if (normalize == representative) {
                        match = true;
                    } else {
                        for (mixture = mix_pairs.begin(); mixture != mix_pairs.end(); mixture++) {
                            if (mixture->first == representative) {
                                match = match || (normalize == mixture->second);
                            } else if (mixture->second == representative) {
                                match = match || (normalize == mixture->first);
                            }
                            if (match) {
                                break;
                            }
                        }
                    }
                }
            }
            if (hasID) {
                set_id_spectra++;
            }
            if (match) {
                set_rep_spectra++;
            }
        }

    }

    void purity(DisCluster *cluster,
            tr1::unordered_map<string, tr1::unordered_map<unsigned int, string> > &identifications,
            string &representative,
            string &all_peptides,
            float &all_purity,
            float &core_purity,
            float &sat_purity,
            float &weak_purity,
            float &mix_purity) {

        tr1::unordered_map<string, unsigned int> peptide_counts;
        tr1::unordered_map<string, string> original_peptide;
        vector<pair<string, string> > mix_pairs;
        string normalized_representative = "";


        unsigned int core_id_spectra = 0;
        unsigned int core_rep_spectra = 0;
        unsigned int sat_id_spectra = 0;
        unsigned int sat_rep_spectra = 0;
        unsigned int wsat_id_spectra = 0;
        unsigned int wsat_rep_spectra = 0;
        unsigned int mix_id_spectra = 0;
        unsigned int mix_rep_spectra = 0;

        set_purity(cluster->core,
                identifications,
                peptide_counts,
                original_peptide,
                mix_pairs,
                core_id_spectra,
                core_rep_spectra,
                normalized_representative
                );

        set_purity(cluster->sat,
                identifications,
                peptide_counts,
                original_peptide,
                mix_pairs,
                sat_id_spectra,
                sat_rep_spectra,
                normalized_representative
                );

        set_purity(cluster->wsat,
                identifications,
                peptide_counts,
                original_peptide,
                mix_pairs,
                wsat_id_spectra,
                wsat_rep_spectra,
                normalized_representative
                );


        set_purity(cluster->mix,
                identifications,
                peptide_counts,
                original_peptide,
                mix_pairs,
                mix_id_spectra,
                mix_rep_spectra,
                normalized_representative
                );

        if (core_id_spectra == 0) {
            core_purity = 1;
        } else {
            core_purity = (float) core_rep_spectra / core_id_spectra;
        }

        if (sat_id_spectra == 0) {
            sat_purity = 1;
        } else {
            sat_purity = (float) (core_rep_spectra + sat_rep_spectra) / (core_id_spectra + sat_id_spectra);
        }

        if (wsat_id_spectra == 0) {
            weak_purity = 1;
        } else {
            weak_purity = (float) wsat_rep_spectra / wsat_id_spectra;
        }

        if (mix_id_spectra == 0) {
            mix_purity = 1;
        } else {
            mix_purity = (float) mix_rep_spectra / mix_id_spectra;
        }

        unsigned int total_id_spectra = core_id_spectra + sat_id_spectra + wsat_id_spectra + mix_id_spectra;
        unsigned int total_rep_spectra = core_rep_spectra + sat_rep_spectra + wsat_rep_spectra + mix_rep_spectra;

        if (core_rep_spectra <= .1 * cluster->core.size()) {
            if (core_rep_spectra + sat_rep_spectra <= .1 * (cluster->core.size() + cluster->sat.size())) {
                if (total_id_spectra < .1 * (cluster->core.size() + cluster->sat.size() + cluster->wsat.size() + cluster->mix.size())) {
                    representative = "";
                    normalized_representative = "";
                } else {
                    representative = original_peptide[normalized_representative];
                }
            } else {
                representative = original_peptide[normalized_representative];
            }
        } else {
            representative = original_peptide[normalized_representative];
        }


        if (total_id_spectra == 0) {
            all_purity = 1;
        } else {
            all_purity = (float) total_rep_spectra / total_id_spectra;
        }

        tr1::unordered_set<string> mix_potentials;
        tr1::unordered_set<string>::iterator core_mix;

        vector<pair<string, string> >::iterator mixture;

        for (mixture = mix_pairs.begin(); mixture != mix_pairs.end(); mixture++) {
            if (mixture->first == normalized_representative) {
                mix_potentials.insert(mixture->second);
            } else if (mixture->second == normalized_representative) {
                mix_potentials.insert(mixture->first);
            }
        }

        for (core_mix = mix_potentials.begin(); core_mix != mix_potentials.end(); core_mix++) {
            representative = representative + "!" + original_peptide[*core_mix];
        }

        tr1::unordered_map<string, unsigned int>::iterator peptide;

        for (peptide = peptide_counts.begin(); peptide != peptide_counts.end(); ++peptide) {
            all_peptides += original_peptide[peptide->first];
            all_peptides += " ";
        }



    }

    short output_clusters(const char* cluster_filename,
            const char* spectra_filename,
            const char* cluster_name,
            const char* cluster_file_base,
            vector<DisCluster> &clusters,
            unsigned int c,
            tr1::unordered_map<string, tr1::unordered_map<unsigned int, string> > &identifications,
            tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> > &globalIndex,
            tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> > &fileIndex,
            tr1::unordered_set<unsigned int> &considered_clust) {
        FILE* spectra_out = fopen(spectra_filename, "a");
        FILE* cluster_out = fopen(cluster_filename, "a");

        if (!spectra_out) {
            cerr << "ERROR: cannot open " << spectra_filename << "\n";
            return -1;
        }
        if (!cluster_out) {
            cerr << "ERROR: cannot open " << cluster_filename << "\n";
            return -1;
        }
        //        DEBUG_MSG(filename);
        //        DEBUG_MSG(clusters.size());

        vector<DisCluster>::iterator cluster;
        vector<ClusterMember>::iterator s;

        unsigned int cn = 0;
        float precursor_mass = 0;

        c = 0;

        for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster, ++cn) {
            if (considered_clust.count(cn) == 0) {
                string representative = "";
                string all_peptides = "";
                float all_purity = 0;
                float core_purity = 0;
                float sat_purity = 0;
                float weak_purity = 0;
                float mix_purity = 0;
                unsigned int size = 0;

                unsigned int parentMass;
                unsigned int charge;

                for (s = cluster->core.begin(); s != cluster->core.end(); ++s, ++size) {
                    charge = s->spec->parentCharge;
                    fprintf(spectra_out, "C");
                    precursor_mass = (float) ((s->spec->parentMass) + (s->spec->parentCharge - 1) * specnets::AAJumps::massHion) / s->spec->parentCharge;
                    fprintf(spectra_out, "\t%s-%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\n", cluster_name, c, s->count_phase4, s->spec->fileName.c_str(), cluster_file_base, s->spec->scan, globalIndex[s->spec->fileName][s->spec->scan] + 1, fileIndex[s->spec->fileName][s->spec->scan] + 1, s->spec->parentCharge, precursor_mass, s->similarity, identifications[s->spec->fileName][s->spec->scan].c_str());
                }
                for (s = cluster->sat.begin(); s != cluster->sat.end(); ++s, ++size) {
                    fprintf(spectra_out, "S");
                    precursor_mass = (float) ((s->spec->parentMass) + (s->spec->parentCharge - 1) * specnets::AAJumps::massHion) / s->spec->parentCharge;
                    fprintf(spectra_out, "\t%s-%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\n", cluster_name, c, s->count_phase4, s->spec->fileName.c_str(), cluster_file_base, s->spec->scan, globalIndex[s->spec->fileName][s->spec->scan] + 1, fileIndex[s->spec->fileName][s->spec->scan] + 1, s->spec->parentCharge, precursor_mass, s->similarity, identifications[s->spec->fileName][s->spec->scan].c_str());
                }
                for (s = cluster->wsat.begin(); s != cluster->wsat.end(); ++s, ++size) {
                    fprintf(spectra_out, "W");
                    precursor_mass = (float) ((s->spec->parentMass) + (s->spec->parentCharge - 1) * specnets::AAJumps::massHion) / s->spec->parentCharge;
                    fprintf(spectra_out, "\t%s-%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\n", cluster_name, c, s->count_phase4, s->spec->fileName.c_str(), cluster_file_base, s->spec->scan, globalIndex[s->spec->fileName][s->spec->scan] + 1, fileIndex[s->spec->fileName][s->spec->scan] + 1, s->spec->parentCharge, precursor_mass, s->similarity, identifications[s->spec->fileName][s->spec->scan].c_str());
                }
                for (s = cluster->mix.begin(); s != cluster->mix.end(); ++s, ++size) {
                    fprintf(spectra_out, "M");
                    precursor_mass = (float) ((s->spec->parentMass) + (s->spec->parentCharge - 1) * specnets::AAJumps::massHion) / s->spec->parentCharge;
                    fprintf(spectra_out, "\t%s-%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\n", cluster_name, c, s->count_phase4, s->spec->fileName.c_str(), cluster_file_base, s->spec->scan, globalIndex[s->spec->fileName][s->spec->scan] + 1, fileIndex[s->spec->fileName][s->spec->scan] + 1, s->spec->parentCharge, precursor_mass, s->similarity, identifications[s->spec->fileName][s->spec->scan].c_str());
                }
                purity(&(*cluster), identifications, representative, all_peptides, all_purity, core_purity, sat_purity, weak_purity, mix_purity);
                //                if (size >= 2) {
                fprintf(cluster_out, "%s-%d\t%d\t%d\t%s\t%s\t%f\t%f\t%f\t%f\t%f\n", cluster_name, c, size, charge, representative.c_str(), all_peptides.c_str(), all_purity, core_purity, sat_purity, weak_purity, mix_purity);
                //                }
                c++;
            }
        }

        fclose(spectra_out);
        fclose(cluster_out);

    }

    bool pairGT(pair<float, unsigned int> a, pair<float, unsigned int> b) {
        return a.first > b.first;
    }

    bool sizeGT(DisCluster a, DisCluster b) {
        return a.core.size() > b.core.size();
    }



    // -------------------------------------------------------------------------

    ExecClusterBins::ExecClusterBins(void) :
    ExecBase(), m_inputSpectra(0x0), m_filteredSpectra(0x0), m_dissimilaritySpectra(0x0), m_representativeSpectra(0x0), m_inputComponents(0x0), m_identifications(0x0), m_index(0x0), m_fileIndex(0x0), m_globalIndex(0x0), ownInput(true), ownOutput(true) {
        m_inputSpectra = new SpecSet();
        m_filteredSpectra = new SpecSet();
        m_dissimilaritySpectra = new SpecSet();
        m_representativeSpectra = new SpecSet();
        m_inputComponents = new vector<SpecSet> ();
        m_index = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_fileIndex = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_globalIndex = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_identifications = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, string> >;
        m_name = "ExecClusterBins";
        m_type = "ExecClusterBins";
    }

    // -------------------------------------------------------------------------

    ExecClusterBins::ExecClusterBins(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_filteredSpectra(0x0), m_dissimilaritySpectra(0x0), m_representativeSpectra(0x0), m_inputComponents(0x0), m_identifications(0x0), m_index(0x0), m_fileIndex(0x0), m_globalIndex(0x0), ownInput(true), ownOutput(true) {
        m_inputSpectra = new SpecSet();
        m_filteredSpectra = new SpecSet();
        m_dissimilaritySpectra = new SpecSet();
        m_representativeSpectra = new SpecSet();
        m_inputComponents = new vector<SpecSet> ();
        m_index = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_fileIndex = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_globalIndex = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> >;
        m_identifications = new tr1::unordered_map<string, tr1::unordered_map<unsigned int, string> >;
        m_name = "ExecClusterBins";
        m_type = "ExecClusterBins";
    }

    // -------------------------------------------------------------------------
    //    ExecKLFilter::ExecKLFilter(const ParameterList & inputParams, vector<SpecSet> * inputSpectra) :
    //            ExecBase(inputParams), m_inputSpectra(inputSpectra), m_filteredSpectra(0x0), ownInput(false), ownOutput(true)
    //    {
    //        m_name = "ExecKLFilter";
    //        m_type = "ExecKLFilter";
    //    }

    // -------------------------------------------------------------------------

    ExecClusterBins::~ExecClusterBins(void) {
        //        if( ownInput ){
        //            delete m_inputSpectra;
        //        }
    }

    // -------------------------------------------------------------------------

    ExecBase * ExecClusterBins::clone(const ParameterList & inputParams) const {
        return new ExecClusterBins(inputParams);
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::invoke(void) {

        int top_peaks_n = m_params.getValueInt("TOP_DISSIMILARITY_PEAKS", 5);
        float tolerance = m_params.getValueFloat("TOLERANCE", 0.05);
        float peak_tolerance = m_params.getValueFloat("PEAK_TOLERANCE", 0.1);
        unsigned int req_peak_match = m_params.getValueInt("MIN_COS_PEAK_MATCH", 5);
        float core_threshold = m_params.getValueFloat("COS_CORE_THRESHOLD", 0.7);
        float sat_threshold = m_params.getValueFloat("COS_SAT_THRESHOLD", 0.5);
        unsigned int inner_cluster_min_size = m_params.getValueInt("MIN_INNER_CLUSTER_SIZE", 2);


        vector<SpecSet>::iterator component;
        vector<DisCluster> all_clusters;
        float score1;
        float score2;


        vector<string> path_names;
        stringSplit(m_params.getValue("TAG_CLUSTERS"), path_names, "/");
        string cluster_basename = path_names.back();
        cluster_basename.erase(cluster_basename.end() - 4, cluster_basename.end());

        string cluster_file_base = cluster_basename + ".mgf";

        string dissimilarity_name = m_params.getValue("DISSIMILARITY_FOLDER") + "/" + cluster_file_base;
        string filtered_name = m_params.getValue("FILTERED_FOLDER") + "/" + cluster_file_base;
        string representatives_name = m_params.getValue("REPRESENTATIVES_FOLDER") + "/" + cluster_file_base;

        unsigned int component_count = 0;

        tr1::unordered_set<unsigned int> considered_clust;


        for (component = m_inputComponents->begin(); component != m_inputComponents->end(); component++) {
            vector<DisCluster> tag_clusters;
            DisCluster tag_cluster;

            //initialize distance matrix

            tr1::unordered_map<unsigned int, tr1::unordered_map<unsigned int, float> > cosine;
            tr1::unordered_map<unsigned int, tr1::unordered_map<unsigned int, float> > proj_cosine;
            vector<vector<MatchType> > matrix(component->size(), vector<MatchType>(component->size()));

            for (unsigned int i = 0; i < matrix.size(); ++i) {
                for (unsigned int j = 0; j < matrix.size(); ++j) {
                    matrix[i][j] = None;
                }
            }

            vector<specnets::Spectrum>::iterator spec1;
            vector<specnets::Spectrum>::iterator spec2;

            Spectrum top_peaks;

            unsigned int i = 0;
            unsigned int j = 0;

            component_count++;

            for (spec1 = component->begin(); spec1 != component->end(); ++spec1) {

                ClusterMember start;

                start.spec = &*spec1;
                start.similarity = 1;
                start.count_phase2 = 0;
                start.count_phase3 = 0;
                start.count_phase4 = 0;

                tag_cluster.core.push_back(start);

                float isotope_increment = (float) 1;

                i = (*m_index)[spec1->fileName][spec1->scan];
                (*m_filteredSpectra)[(*m_globalIndex)[spec1->fileName][spec1->scan]] = *spec1;
                (*m_filteredSpectra)[(*m_globalIndex)[spec1->fileName][spec1->scan]].scan = (*m_globalIndex)[spec1->fileName][spec1->scan] + 1;
                //select top 5 for dissimilarity

                //                DEBUG_MSG(spec1->scan);

                top_peaks = *spec1;
                top_peaks.rankFilterPeaks(1, 30);
                top_peaks.selectTopK(5, 0x0);

                (*m_dissimilaritySpectra)[(*m_globalIndex)[spec1->fileName][spec1->scan]] = top_peaks;
                (*m_dissimilaritySpectra)[(*m_globalIndex)[spec1->fileName][spec1->scan]].scan = (*m_globalIndex)[spec1->fileName][spec1->scan] + 1;



                for (spec2 = component->begin(); spec2 != component->end(); ++spec2, ++j) {


                    j = (*m_index)[spec2->fileName][spec2->scan];

                    //                    DEBUG_MSG(spec2->scan);


                    //check pm isotopes

                    float mass_difference = abs(spec1->parentMass - spec2->parentMass);
                    float charge_tolerance = tolerance * spec1->parentCharge;

                    if ((mass_difference <= charge_tolerance)
                            || (mass_difference <= charge_tolerance + isotope_increment && mass_difference >= isotope_increment - charge_tolerance)
                            || (mass_difference <= charge_tolerance + 2 * isotope_increment && mass_difference >= 2 * isotope_increment - charge_tolerance)) {

                        //calculate one way dissimilarity

                        //                        DEBUG_MSG("PM match, calculating dissimilarity");

                        int dissimilarity = top_peaks_n;

                        //                        DEBUG_MSG(top_peaks.size());

                        for (size_t p = 0; p < top_peaks.size(); ++p) {

                            if (spec2->findPeaks(top_peaks[p][0], peak_tolerance) != -1)
                                dissimilarity--;

                        }


                        //                        DEBUG_MSG(dissimilarity);



                        //                        DEBUG_MSG(matrix[j][i].match_type);

                        unsigned int matched_peaks = 0;

                        if (dissimilarity <= 1) {
                            if (matrix[j][i] == None) {
                                float cos = spec1->scoreMatch2(*spec2,
                                        peak_tolerance,
                                        matched_peaks,
                                        score1,
                                        score2,
                                        true,
                                        false,
                                        false);
                                float pcos = spec1->scoreMatch2(*spec2,
                                        peak_tolerance,
                                        matched_peaks,
                                        score1,
                                        score2,
                                        true,
                                        true,
                                        false);
                                //                                    DEBUG_MSG(pcos);
                                if (pcos > sat_threshold && matched_peaks >= req_peak_match) {
                                    if (score1 > score2) {
                                        matrix[i][j] = Subset;
                                        matrix[j][i] = Superset;
                                    } else {
                                        matrix[i][j] = Superset;
                                        matrix[j][i] = Subset;
                                    }
                                    cosine[i][j] = cos;
                                    cosine[j][i] = cos;
                                    proj_cosine[i][j] = pcos;
                                    proj_cosine[j][i] = pcos;
                                }
                            } else {
                                if (cosine[i][j] > sat_threshold) {
//                                    proj_cosine[i][j] = cosine[i][j];
//                                    proj_cosine[j][i] = cosine[j][i];
                                    matrix[i][j] = Exact;
                                    matrix[j][i] = Exact;
                                }
                            }

                        }

                        //                        DEBUG_MSG("Updated matrix");


                    } else {
                        //                        DEBUG_MSG("No PM match, moving on...");
                    }

                }

            }


            tag_clusters.push_back(tag_cluster);


            DEBUG_MSG("OUTPUT PHASE 1")

            output_clusters(m_params.getValue("OUTPUT_CLUSTERS_PHASE1").c_str(), m_params.getValue("OUTPUT_SPECTRA_PHASE1").c_str(), cluster_basename.c_str(), cluster_file_base.c_str(), tag_clusters, component_count, *m_identifications, *m_globalIndex, *m_fileIndex, considered_clust);

            //            DEBUG_MSG("Built matrix");


            vector<pair<float, unsigned int> > best_edges(component->size());

            for (spec1 = component->begin(); spec1 != component->end(); ++spec1) {

                i = (*m_index)[spec1->fileName][spec1->scan];

                float best_edge = 0;

                for (spec2 = component->begin(); spec2 != component->end(); ++spec2) {

                    j = (*m_index)[spec2->fileName][spec2->scan];

                    //                    DEBUG_MSG(matrix[i][j].match_type);

                    if (matrix[i][j] == Exact && cosine[i][j] > best_edge) {

                        best_edge = cosine[i][j];

                    }

                }

                best_edges[i] = make_pair(best_edge, i);
                //                DEBUG_MSG(i);


            }

            sort(best_edges.begin(), best_edges.end(), pairGT);

            // begin clustering - from best edge to worst

            //            DEBUG_MSG((best_edges[0].second)->scan);


            vector<DisCluster> clusters;
            vector<pair<float, unsigned int> >::iterator next_spec;
            vector<Spectrum>::iterator spec;
            vector<unsigned int>::iterator specp;
            vector<ClusterMember>::iterator core;
            vector<pair<float, unsigned int> >::iterator potential_core;
            tr1::unordered_set<unsigned int> considered_spec;
            tr1::unordered_set<unsigned int> considered_spec2;
            vector<specnets::Spectrum>::iterator inner_spec;

            //            for(next_spec = best_edges.begin(); next_spec != best_edges.end(); next_spec++) {
            //
            //                DEBUG_MSG(next_spec->first);
            //                DEBUG_MSG(next_spec->second->scan);
            //
            //            }


            for (next_spec = best_edges.begin(); next_spec != best_edges.end(); ++next_spec) {

                //                DEBUG_MSG(next_spec->first);
                i = next_spec->second;
                //                DEBUG_MSG(&*(next_spec->second)->scan);

                if (considered_spec.count(i) == 0) {

                    unsigned int count = 0;

                    //                    DEBUG_MSG("Considered cluster");

                    DisCluster cluster;
                    ClusterMember start;
                    start.spec = &((*component)[i]);
                    //                    DEBUG_MSG(start.spec->scan);
                    start.similarity = 1;
                    start.count_phase2 = count;
                    start.count_phase3 = count;
                    start.count_phase4 = count;
                    count++;
                    cluster.core.push_back(start);


                    vector<pair<float, unsigned int> > core_edges;

                    for (spec = component->begin(); spec != component->end(); ++spec) {

                        j = (*m_index)[spec->fileName][spec->scan];

                        if (considered_spec.count(j) == 0
                                && matrix[i][j] == Exact
                                && cosine[i][j] >= core_threshold
                                ) {

                            core_edges.push_back(make_pair(cosine[i][j], j));

                        }

                    }

                    //                    DEBUG_MSG(core_edges.size());

                    sort(core_edges.begin(), core_edges.end(), pairGT);

                    for (potential_core = core_edges.begin(); potential_core != core_edges.end(); potential_core++) {

                        unsigned int k = potential_core->second;

                        float clique_cos = 0;
                        unsigned int clique_size = 0;
                        unsigned int edges_size = 0;

                        bool first = true;

                        for (core = cluster.core.begin(); core != cluster.core.end(); ++core) {

                            unsigned int l = (*m_index)[core->spec->fileName][core->spec->scan];

                            if (matrix[k][l] == Exact
                                    && cosine[k][l] >= core_threshold
                                    ) {
                                if (first) {
                                    clique_cos = cosine[k][l];
                                    first = false;
                                }
                                clique_size++;
                            }

                            if (matrix[k][l] != None) {
                                edges_size++;
                            }

                        }

                        if (clique_size > .8 * cluster.core.size()) {

                            ClusterMember mem;
                            mem.spec = &((*component)[k]);
                            mem.similarity = clique_cos;
                            mem.count_phase2 = count;
                            mem.count_phase3 = count;
                            mem.count_phase4 = count;
                            count++;
                            cluster.core.push_back(mem);

                        }

                    }

                    if (cluster.core.size() >= inner_cluster_min_size) {


                        for (core = cluster.core.begin(); core != cluster.core.end(); ++core) {


                            unsigned int k = (*m_index)[core->spec->fileName][core->spec->scan];

                            considered_spec.insert(k);
                            considered_spec2.insert(k);

                            //                            for (inner_spec = component->begin(); inner_spec != component->end(); inner_spec++){
                            //
                            //                                unsigned int l = (*m_index)[inner_spec->fileName][inner_spec->scan];
                            //
                            //                                if (considered_spec.count(l) == 0
                            //                                    && (matrix[k][l] == Exact
                            //                                    && cosine[k][l] >= sat_threshold)
                            //                                    || ((matrix[k][l] == Subset ||
                            //                                        matrix[k][l] == Superset)
                            //                                    && proj_cosine[k][l] >= sat_threshold)) {
                            //
                            //                                    considered_spec.insert(l);
                            //
                            //                                }
                            //
                            //                            }

                        }



                        clusters.push_back(cluster);
                        //                    DEBUG_MSG("Added cluster");
                    }


                }


            }


            //
            vector<DisCluster >::iterator cluster;

            //            DEBUG_MSG("Adding sat");

            //
            //            for(spec = component->begin(); spec != component->end(); ++spec) {
            //
            //                j = (*m_index)[spec->fileName][spec->scan];
            //
            //                if (considered_spec2.count(j) == 0) {
            //
            //                    float best_cluster_cos = 0;
            //                    DisCluster* best_cluster = NULL;
            //
            //                    for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster) {
            //
            //                        float best_edge = 0;
            //                        unsigned int exact_match_edges = 0;
            //
            //                        for (core = cluster->core.begin(); core != cluster->core.end(); ++core) {
            //
            //                            i = (*m_index)[core->spec->fileName][core->spec->scan];
            //
            //                            if (matrix[i][j] == Exact
            //                                && cosine[i][j] >= sat_threshold)
            //                            {
            //                                best_edge = max(best_edge,cosine[i][j]);
            //                                exact_match_edges++;
            //                            }
            //
            //                        }
            //
            //                        if (best_edge >= sat_threshold
            //                                && best_edge > best_cluster_cos
            //                                && exact_match_edges > .75 * cluster->core.size()) {
            //                            best_cluster = &(*cluster);
            //                            best_cluster_cos = best_edge;
            //                        }
            //
            //                    }
            //
            //
            //                    if (best_cluster != NULL) {
            ////                        DEBUG_MSG("Found eligibile sat cluster");
            //                        ClusterMember sat;
            //                        sat.spec = &(*spec);
            //                        sat.similarity = best_cluster_cos;
            //                        sat.count = count;
            //                        count++;
            //                        best_cluster->sat.push_back(sat);
            //                        considered_spec2.insert(j);
            //                    }
            //                }
            //
            //            }

            //            DEBUG_MSG("Adding wsat");
            //
            //            for(spec = component->begin(); spec != component->end(); ++spec) {
            //
            //                j = (*m_index)[spec->fileName][spec->scan];
            //
            //                if (considered_spec2.count(j) == 0) {
            //
            //                    float best_cluster_cos = 0;
            //                    DisCluster* best_cluster = NULL;
            //                    MatchType min_match = Subset;
            //
            //                    for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster) {
            //
            //                        float best_edge = 0;
            //
            //
            //                        for (core = cluster->core.begin(); core != cluster->core.end(); ++core) {
            //
            //                               i = (*m_index)[core->spec->fileName][core->spec->scan];
            //
            //                               if (matrix[i][j] == Exact || matrix[i][j] == min_match)
            //                               {
            //                                   best_edge = max(best_edge,proj_cosine[i][j]);
            //                                   if (matrix[i][j] == Exact) {
            //                                       min_match = Exact;
            //                                   }
            //                               }
            //
            //                        }
            //
            //                        if (best_edge >= sat_threshold
            //                                   && best_edge > best_cluster_cos)
            //                           {
            //                               best_cluster = &*cluster;
            //                               best_cluster_cos = best_edge;
            //                           }
            //
            //                    }
            //
            //
            //                    if (best_cluster != NULL) {
            ////                        DEBUG_MSG("Found eligibile sat cluster");
            //                        ClusterMember sat;
            //                        sat.spec = &(*spec);
            //                        sat.similarity = best_cluster_cos;
            //                        sat.count = count;
            //                        count++;
            //                        best_cluster->wsat.push_back(sat);
            //                        considered_spec2.insert(j);
            //                    }
            //                }
            //
            //            }

            //            DEBUG_MSG("Adding mix");
            //
            //            for(spec = component->begin(); spec != component->end(); ++spec) {
            //
            //                j = (*m_index)[spec->fileName][spec->scan];
            //
            //                if (considered_spec2.count(j) == 0) {
            //
            //                    float best_cluster_cos = 0;
            //                    DisCluster* best_cluster = NULL;
            //                    MatchType min_match = Superset;
            //
            //                    for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster) {
            //
            //                        float best_edge = 0;
            //
            //
            //                        for (core = cluster->core.begin(); core != cluster->core.end(); ++core) {
            //
            //                            i = (*m_index)[core->spec->fileName][core->spec->scan];
            //
            //                            if (matrix[i][j] == Exact || matrix[i][j] == min_match)
            //                            {
            //                                best_edge = max(best_edge,proj_cosine[i][j]);
            //                                if (matrix[i][j] == Exact) {
            //                                    min_match = Exact;
            //                                }
            //                            }
            //
            //                        }
            //
            //                        if (best_edge >= sat_threshold
            //                            && best_edge > best_cluster_cos)
            //                        {
            //                            best_cluster = &*cluster;
            //                            best_cluster_cos = best_edge;
            //                        }
            //
            //                    }
            //
            //
            //                    if (best_cluster != NULL) {
            ////                        DEBUG_MSG("Found eligibile sat cluster");
            //                        ClusterMember sat;
            //                        sat.spec = &(*spec);
            //                        sat.similarity = best_cluster_cos;
            //                        sat.count = count;
            //                        count++;
            //                        best_cluster->mix.push_back(sat);
            //                        considered_spec2.insert(j);
            //                    } else {
            //                        DEBUG_MSG("Need new cluster");
            //                        DisCluster cluster_new;
            //                        ClusterMember start;
            //                        start.spec = &(*spec);
            //                        start.similarity = 1;
            //                        start.count = count;
            //                        cluster_new.core.push_back(start);
            //                        clusters.push_back(cluster_new);
            //                    }
            //                }
            //
            //            }

            //            DEBUG_MSG("Finished adding");

            for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster) {
                all_clusters.push_back(*cluster);
            }

        }

        m_filteredSpectra->SaveSpecSet_mgf(filtered_name.c_str(), true);
        m_filteredSpectra->resize(0);
        m_dissimilaritySpectra->SaveSpecSet_mgf(dissimilarity_name.c_str(), true);
        m_dissimilaritySpectra->resize(0);




        sort(all_clusters.begin(), all_clusters.end(), sizeGT);


        vector<ClusterMember>::iterator mem;

        tr1::unordered_set<unsigned int> considered_clust_clique;
        tr1::unordered_set<unsigned int> merged_clust;
        vector<pair<float, unsigned int> >::iterator next_spec;
        vector<pair<float, unsigned int> >::iterator potential_core;

        DEBUG_MSG("OUTPUT PHASE 2")
        output_clusters(m_params.getValue("OUTPUT_CLUSTERS_PHASE2").c_str(), m_params.getValue("OUTPUT_SPECTRA_PHASE2").c_str(), cluster_basename.c_str(), cluster_file_base.c_str(), all_clusters, component_count, *m_identifications, *m_globalIndex, *m_fileIndex, considered_clust);

        vector<pair<unsigned int, unsigned int> > spec_by_size;

        tr1::unordered_map<unsigned int, tr1::unordered_map<unsigned int, float> > cosine;
        tr1::unordered_map<unsigned int, tr1::unordered_map<unsigned int, float> > proj_cosine;

        tr1::unordered_set<unsigned int> not_core;


        vector<vector<MatchType> > matrix(all_clusters.size(), vector<MatchType>(all_clusters.size()));

        for (unsigned int i = 0; i < matrix.size(); ++i) {
            for (unsigned int j = 0; j < matrix.size(); ++j) {
                matrix[i][j] = None;
            }
        }

        Spectrum top_peaks;



        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {


            float isotope_increment = (float) 1;
            //                DEBUG_MSG(spec1->scan);
            Spectrum spec1 = *all_clusters[idx1].core[0].spec;
            unsigned int i = idx1;
            top_peaks = spec1;
            top_peaks.rankFilterPeaks(1, 30);
            top_peaks.selectTopK(5, 0x0);

            for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {

                Spectrum spec2 = *all_clusters[idx2].core[0].spec;

                unsigned int j = idx2;

                //                    DEBUG_MSG(spec2->scan);

                if (idx1 != idx2) {

                    //check pm isotopes

                    float mass_difference = abs(spec1.parentMass - spec2.parentMass);
                    float charge_tolerance = tolerance * spec1.parentCharge;


                    if ((spec1.parentCharge == spec2.parentCharge) &&
                            ((mass_difference <= charge_tolerance)
                            || (mass_difference <= charge_tolerance + isotope_increment && mass_difference >= isotope_increment - charge_tolerance)
                            || (mass_difference <= charge_tolerance + 2 * isotope_increment && mass_difference >= 2 * isotope_increment - charge_tolerance))) {

                        //calculate one way dissimilarity

                        //                        DEBUG_MSG("PM match, calculating dissimilarity");

                        int dissimilarity = top_peaks_n;

                        //                        DEBUG_MSG(top_peaks.size());

                        for (size_t p = 0; p < top_peaks.size(); ++p) {

                            if (spec2.findPeaks(top_peaks[p][0], peak_tolerance) != -1)
                                dissimilarity--;

                        }


                        //                        DEBUG_MSG(dissimilarity);



                        //                        DEBUG_MSG(matrix[j][i].match_type);

                        unsigned int matched_peaks = 0;

                        if (dissimilarity <= 1) {
                            if (matrix[j][i] == None) {
                                float cos = spec1.scoreMatch2(spec2,
                                        peak_tolerance,
                                        matched_peaks,
                                        score1,
                                        score2,
                                        true,
                                        false,
                                        false);
                                float pcos = spec1.scoreMatch2(spec2,
                                        peak_tolerance,
                                        matched_peaks,
                                        score1,
                                        score2,
                                        true,
                                        true,
                                        false);
                                //                                    DEBUG_MSG(pcos);
                                if (pcos > sat_threshold && matched_peaks >= req_peak_match) {
                                    matrix[i][j] = Subset;
                                    matrix[j][i] = Superset;
                                    cosine[i][j] = cos;
                                    cosine[j][i] = cos;
                                    proj_cosine[i][j] = pcos;
                                    proj_cosine[j][i] = pcos;
                                }
                            } else {
                                if (cosine[i][j] > sat_threshold) {
//                                    proj_cosine[i][j] = cosine[i][j];
//                                    proj_cosine[j][i] = cosine[j][i];
                                    matrix[i][j] = Exact;
                                    matrix[j][i] = Exact;
                                }
                            }

                        }

                        //                        DEBUG_MSG("Updated matrix");


                    } else {
                        //                        DEBUG_MSG("No PM match, moving on...");
                    }

                }
            }

        }

        //            vector<pair<float, unsigned int> > best_edges(all_clusters.size());
        //
        //            for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
        //
        //                float best_edge = 0;
        //
        //                    for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {
        //
        ////                    DEBUG_MSG(matrix[i][j].match_type);
        //
        //                    if (matrix[idx1][idx2] == Exact && cosine[idx1][idx2] > best_edge) {
        //
        //                        best_edge = cosine[idx1][idx2];
        //
        //                    }
        //
        //                }
        //
        //                best_edges[idx1] = make_pair(best_edge,idx1);
        ////                DEBUG_MSG(i);
        //
        //
        //            }
        //
//                    sort(best_edges.begin(), best_edges.end(), pairGT);
        //
        //

        FILE* pairs = fopen(m_params.getValue("PAIRWISE").c_str(), "a");


        for (unsigned int i = 0; i < matrix.size(); ++i) {
            for (unsigned int j = 0; j < matrix.size(); ++j) {
                if (matrix[i][j] != None) {
                    fprintf(pairs, "%s\t%d\t%s\t%d\t%d\t%f\t%f\n",
                            all_clusters[i].core[0].spec->fileName.c_str(),
                            all_clusters[i].core[0].spec->scan,
                            all_clusters[j].core[0].spec->fileName.c_str(),
                            all_clusters[j].core[0].spec->scan,
                            matrix[i][j],
                            cosine[i][j],
                            proj_cosine[i][j]
                            );
                };
            }
        }

        fclose(pairs);


        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
            all_clusters[idx1].sub_clusters.push_back(make_pair(all_clusters[idx1].core.size(),idx1));
        }


        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {

            //                DEBUG_MSG(next_spec->first);
            //                DEBUG_MSG(&*(next_spec->second)->scan);
            unsigned int i = idx1;
            if (considered_clust_clique.count(i) == 0) {

                //                    DEBUG_MSG("Considered cluster");


                vector<pair<float, unsigned int> > core_edges;

                for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {

                    unsigned int j = idx2;


                    if ((considered_clust_clique.count(j) == 0)
                            && (matrix[i][j] == Exact)
                            && (cosine[i][j] >= core_threshold)
                            ) {

                        core_edges.push_back(make_pair(cosine[i][j], j));

                    }

                }

                //                    DEBUG_MSG(core_edges.size());

//                sort(core_edges.begin(), core_edges.end(), pairGT);


                for (potential_core = core_edges.begin(); potential_core != core_edges.end(); potential_core++) {

                    unsigned int k = potential_core->second;

                    float clique_cos = 0;
                    unsigned int clique_size = 0;
                    unsigned int edges_size = 0;


                    for (unsigned int core = 0; core <  all_clusters[i].sub_clusters.size(); core++) {

                        unsigned int l =  all_clusters[i].sub_clusters[core].second;

                        if (matrix[k][l] == Exact
                                && cosine[k][l] >= core_threshold
                                ) {
                            clique_cos = max(clique_cos, cosine[k][l]);
                            clique_size++;

                        }

                        if (matrix[k][l] != None) {
                            edges_size++;
                        }

                    }

                    if (clique_size >= .8 * all_clusters[i].sub_clusters.size()) {

                        unsigned int cluster_size = all_clusters[i].sat.size() + all_clusters[i].core.size();

                        all_clusters[i].sub_clusters.push_back(make_pair(all_clusters[k].core.size(),k));
                        for (mem = all_clusters[k].core.begin(); mem != all_clusters[k].core.end(); ++mem) {
                            mem->count_phase3 = cluster_size + mem->count_phase3;
                            mem->count_phase4 = cluster_size + mem->count_phase4;
                            all_clusters[i].core.push_back((*mem));
                        }
                        for (mem = all_clusters[k].sat.begin(); mem != all_clusters[k].sat.end(); ++mem) {
                            mem->count_phase3 = cluster_size + mem->count_phase3;
                            mem->count_phase4 = cluster_size + mem->count_phase4;
                            all_clusters[i].sat.push_back((*mem));
                        }
                        for (mem = all_clusters[k].wsat.begin(); mem != all_clusters[k].wsat.end(); ++mem) {
                            mem->count_phase3 = cluster_size + mem->count_phase3;
                            mem->count_phase4 = cluster_size + mem->count_phase4;
                            all_clusters[i].wsat.push_back((*mem));
                        }
                        for (mem = all_clusters[k].mix.begin(); mem != all_clusters[k].mix.end(); ++mem) {
                            mem->count_phase3 = cluster_size + mem->count_phase3;
                            mem->count_phase4 = cluster_size + mem->count_phase4;
                            all_clusters[i].mix.push_back((*mem));
                        }
                        considered_clust_clique.insert(k);
                        considered_clust.insert(k);

                    }

                }

                //                    if (cluster.core.size() >= inner_cluster_min_size) {
                //
                //
                //                        for (core = cluster.core.begin(); core != cluster.core.end(); ++core) {
                //
                //
                //                            unsigned int k = (*m_index)[core->spec->fileName][core->spec->scan];
                //
                //                            considered_spec.insert(k);
                //                            considered_spec2.insert(k);
                //
                ////                            for (inner_spec = component->begin(); inner_spec != component->end(); inner_spec++){
                ////
                ////                                unsigned int l = (*m_index)[inner_spec->fileName][inner_spec->scan];
                ////
                ////                                if (considered_spec.count(l) == 0
                ////                                    && (matrix[k][l] == Exact
                ////                                    && cosine[k][l] >= sat_threshold)
                ////                                    || ((matrix[k][l] == Subset ||
                ////                                        matrix[k][l] == Superset)
                ////                                    && proj_cosine[k][l] >= sat_threshold)) {
                ////
                ////                                    considered_spec.insert(l);
                ////
                ////                                }
                ////
                ////                            }
                //
                //                        }
                //
                //                        clusters.push_back(cluster);
                //    //                    DEBUG_MSG("Added cluster");
                //                    }


            }


        }

        DEBUG_MSG("OUTPUT PHASE 3")


        output_clusters(m_params.getValue("OUTPUT_CLUSTERS_PHASE3").c_str(), m_params.getValue("OUTPUT_SPECTRA_PHASE3").c_str(), cluster_basename.c_str(), cluster_file_base.c_str(), all_clusters, component_count, *m_identifications, *m_globalIndex, *m_fileIndex, considered_clust);


        vector<DisCluster >::iterator cluster;

        //            DEBUG_MSG("Adding sat");

        spec_by_size.resize(all_clusters.size());

        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
            spec_by_size[idx1] = make_pair(all_clusters[idx1].core.size(),idx1);
        }

        sort(spec_by_size.begin(), spec_by_size.end(), pairGT);

        for (unsigned int idx1 = 0; idx1 < spec_by_size.size(); ++idx1) {
            unsigned int j = spec_by_size[idx1].second;

            if (considered_clust.count(j) == 0) {

                float best_cluster_cos = 0;
                int best_cluster_idx = -1;


                for (unsigned int idx2 = 0; idx2 < spec_by_size.size(); ++idx2) {

                    unsigned int m = spec_by_size[idx2].second;


                    if (j != m && considered_clust.count(m) == 0) {

                        float best_edge = 0;
                        unsigned int exact_match_edges = 0;
                        unsigned int non_core_edges = 0;

                        bool shallow = false;


                        for (unsigned int core1 = 0; core1 < all_clusters[m].sub_clusters.size(); core1++) {

                            float local_best_edge = 0;

                            unsigned int total_edges = 0;


                            for (unsigned int core2 = 0; core2 < all_clusters[j].sub_clusters.size(); core2++) {

                                unsigned int i = all_clusters[m].sub_clusters[core1].second;
                                unsigned int k = all_clusters[j].sub_clusters[core2].second;

//                            if (matrix[i][k] != Exact) {
//                                exact_match_edges++;
//                            }

                                if (not_core.count(i) != 0) {
                                    non_core_edges += 1;
                                }

                                if (matrix[i][k] == Exact && not_core.count(k) == 0 && not_core.count(i) == 0) {
                                    exact_match_edges++;
                                    total_edges += all_clusters[j].sub_clusters[core2].first;
                                    if (cosine[i][k] > local_best_edge) {
                                        local_best_edge = cosine[i][k];
                                    }
                                }

                                if (matrix[i][k] == Exact && not_core.count(i) != 0 && not_core.count(k) == 0) {
                                    shallow = true;
                                }
                            }

                            best_edge += local_best_edge * total_edges;

                        }

                        if (non_core_edges == 0) {
                            shallow = true;
                        }

                        if (best_edge >= sat_threshold
                                && best_edge > best_cluster_cos) {
                            best_cluster_idx = m;
                            best_cluster_cos = best_edge;
                        }

                    }
                }
                if (best_cluster_idx != -1) {//
                    unsigned int bigger = 0;
                    unsigned int smaller = 0;
                    if (all_clusters[best_cluster_idx].core.size() + all_clusters[best_cluster_idx].sat.size() > all_clusters[j].core.size() + all_clusters[j].sat.size()) {
                        bigger = best_cluster_idx;
                        smaller = j;
                    } else {
                        bigger = j;
                        smaller = best_cluster_idx;
                    }
                    unsigned int cluster_size = all_clusters[bigger].sat.size() + all_clusters[bigger].core.size();
                    for (unsigned int ic = 0; ic < all_clusters[smaller].sub_clusters.size(); ic++) {
                        all_clusters[bigger].sub_clusters.push_back(all_clusters[smaller].sub_clusters[ic]);
                        not_core.insert(all_clusters[smaller].sub_clusters[ic].second);
                    }
                    for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
                        mem->count_phase4 = cluster_size + mem->count_phase4;
                        all_clusters[bigger].sat.push_back((*mem));
                    }
                    for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
                        mem->count_phase4 = cluster_size + mem->count_phase4;
                        all_clusters[bigger].wsat.push_back((*mem));
                    }
                    for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
                        mem->count_phase4 = cluster_size + mem->count_phase4;
                        all_clusters[bigger].wsat.push_back((*mem));
                    }
                    for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
                        mem->count_phase4 = cluster_size + mem->count_phase4;
                        all_clusters[bigger].mix.push_back((*mem));
                    }
                    considered_clust.insert(smaller);
                }
            }
        }


        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
            spec_by_size[idx1] = make_pair(all_clusters[idx1].core.size() + all_clusters[idx1].sat.size(),idx1);
        }

        sort(spec_by_size.begin(), spec_by_size.end(), pairGT);

        //            DEBUG_MSG("Adding wsat");

        for (unsigned int idx1 = 0; idx1 < spec_by_size.size(); ++idx1) {
            unsigned int j = spec_by_size[idx1].second;

                if (considered_clust.count(j) == 0) {

                    float best_cluster_cos = 0;
                    int best_cluster_idx = -1;
                    MatchType min_match_w = Subset;
                    MatchType min_match_m = Superset;


                for (unsigned int idx2 = 0; idx2 < spec_by_size.size(); ++idx2) {

                    unsigned int m = spec_by_size[idx2].second;


                    if (j != m && considered_clust.count(m) == 0) {

                            unsigned int exact_match_edges = 0;
                            unsigned int non_core_edges = 0;
                            bool shallow = false;

                            float best_edge = 0;

                            for (unsigned int core1 = 0; core1 < all_clusters[m].sub_clusters.size(); core1++) {

                                float local_best_edge = 0;
                                unsigned int total_edges = 0;

                                for (unsigned int core2 = 0; core2 < all_clusters[j].sub_clusters.size(); core2++) {

                                    unsigned int i = all_clusters[m].sub_clusters[core1].second;
                                    unsigned int k = all_clusters[j].sub_clusters[core2].second;

                                    if (matrix[i][k] != None) {
                                        exact_match_edges++;
                                        total_edges += all_clusters[j].sub_clusters[core2].first;
                                    }

                                    if ((matrix[i][k] == Exact || matrix[i][k] == min_match_w || matrix[i][k] == min_match_m) && not_core.count(i) == 0 && not_core.count(k) == 0) {
                                        if (proj_cosine[i][k] > local_best_edge) {
                                            local_best_edge = proj_cosine[i][k];
                                        }
                                        if (matrix[i][k] == Exact) {
                                            min_match_w = Exact;
                                            min_match_m = Exact;
                                        }
                                    }

                                    if (not_core.count(i) != 0) {
                                        non_core_edges += 1;
                                    }

                                    if (matrix[i][k] == Exact && not_core.count(i) != 0 && not_core.count(k) == 0) {
                                        shallow = true;
                                    }


                                }

                                best_edge += local_best_edge * total_edges;;

                            }

                            if (non_core_edges == 0) {
                                shallow = true;
                            }

                            if (best_edge >= sat_threshold
                                && best_edge > best_cluster_cos) {
                                best_cluster_idx = m;
                                best_cluster_cos = best_edge;
                            }

                        }

                    }

                    if (best_cluster_idx != -1) {//
                        unsigned int bigger = 0;
                        unsigned int smaller = 0;
                        if (all_clusters[best_cluster_idx].core.size() + all_clusters[best_cluster_idx].sat.size() < all_clusters[j].core.size() + all_clusters[j].sat.size()) {
                           bigger = j;
                           smaller = best_cluster_idx;
                        } else {
                            bigger = best_cluster_idx;
                            smaller = j;
                        }
                        unsigned int cluster_size = all_clusters[bigger].mix.size() + all_clusters[bigger].wsat.size() + all_clusters[bigger].sat.size() + all_clusters[bigger].core.size();
                        for (unsigned int ic = 0; ic < all_clusters[smaller].sub_clusters.size(); ic++) {
                            all_clusters[bigger].sub_clusters.push_back(all_clusters[smaller].sub_clusters[ic]);
                            not_core.insert(all_clusters[smaller].sub_clusters[ic].second);
                        }

                        if (matrix[bigger][smaller] == Subset) {

                            for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].mix.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].mix.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].mix.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].mix.push_back((*mem));
                            }
                        } else {
                            for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].wsat.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].wsat.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].wsat.push_back((*mem));
                            }
                            for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
                                mem->count_phase4 = cluster_size + mem->count_phase4;
                                all_clusters[bigger].mix.push_back((*mem));
                            }
                        }

                        considered_clust.insert(smaller);
                    }
                }



        }

//        //            DEBUG_MSG("Adding mix");
//        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
//
//            if (considered_clust.count(idx1) == 0) {
//
//                j = idx1;
//
//                float best_cluster_cos = 0;
//                int best_cluster_idx = -1;
//                MatchType min_match = Superset;
//
//                for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {
//
//
//
//                    if (idx1 != idx2 && considered_clust.count(idx2) == 0) {
//
//                         i = idx2;
//
//                          float best_edge = 0;
//
//                          unsigned int exact_match_edges = 0;
//
//                        for (unsigned int core = 0; core < all_clusters[j].sub_clusters.size(); core++) {
//
//                            unsigned int k = all_clusters[j].sub_clusters[core];
//
//
//                                if (matrix[i][k] != None) {
//                                    exact_match_edges++;
//                                }
//
//                                if ((matrix[i][k] == Exact || matrix[i][k] == min_match) && not_core.count(k) == 0) {
//                                    best_edge = proj_cosine[i][k];
//                                    if (matrix[i][k] == Exact) {
//                                        min_match = Exact;
//                                    }
//                                }
//                        }
//
//                        if (best_edge >= sat_threshold
//                                && best_edge > best_cluster_cos
//                                && exact_match_edges > .2 * all_clusters[j].sub_clusters.size()) {
//                            best_cluster_idx = idx2;
//                            best_cluster_cos = best_edge;
//                        }
//
//                    }
//                }
//
//
//                if (best_cluster_idx != -1) {//
//                    unsigned int bigger = 0;
//                    unsigned int smaller = 0;
//                    if (all_clusters[best_cluster_idx].core.size() + all_clusters[best_cluster_idx].core.size() > all_clusters[idx1].core.size() + all_clusters[idx1].core.size()) {
//                        bigger = best_cluster_idx;
//                        smaller = idx1;
//                    } else {
//                        bigger = idx1;
//                        smaller = best_cluster_idx;
//                    }
//                    for (unsigned int ic = 0; ic < all_clusters[smaller].sub_clusters.size(); ic++) {
//                        all_clusters[bigger].sub_clusters.push_back(all_clusters[smaller].sub_clusters[ic]);
//                    }
//                    not_core.insert(smaller);
//                    for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
//                        all_clusters[bigger].mix.push_back((*mem));
//                    }
//                    for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
//                        all_clusters[bigger].mix.push_back((*mem));
//                    }
//                    for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
//                        all_clusters[bigger].mix.push_back((*mem));
//                    }
//                    for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
//                        all_clusters[bigger].mix.push_back((*mem));
//                    }
//                    considered_clust.insert(smaller);
//                }
//            }
//
//        }

        //
        //        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
        //            if (considered_clust.count(idx1) == 0) {
        //
        //                considered_clust.insert(idx1);
        //
        //                float isotope_increment = (float) 1;
        //
        //                vector<ClusterMember>::iterator core;
        //
        //                for (core = all_clusters[idx1].core.begin(); core != all_clusters[idx1].core.end(); ++core) {
        //                    Spectrum spec1 = *(*core).spec;
        //
        //                    Spectrum top_peaks = spec1;
        //                    top_peaks.rankFilterPeaks(1,3);
        //                    top_peaks.selectTopK(5, 0x0);
        //
        //                    for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {
        //
        //                        if (considered_clust.count(idx2) == 0) {
        //                            Spectrum spec2 = *all_clusters[idx2].core[0].spec;
        //
        //                            float mass_difference = abs(spec1.parentMass - spec2.parentMass);
        //
        //                            if (spec1.parentCharge == spec2.parentCharge && ((mass_difference < tolerance)
        //                                || (mass_difference < tolerance + isotope_increment && mass_difference > isotope_increment - tolerance)
        //                                || (mass_difference < tolerance + 2*isotope_increment && mass_difference > 2*isotope_increment - tolerance))) {
        //
        //
        //                                int dissimilarity = top_peaks_n;
        //
        //        //                        DEBUG_MSG(top_peaks.size());
        //
        //                                for (size_t p = 0; p < top_peaks.size(); ++p){
        //                                    if (spec2.findPeaks(top_peaks[p][0],tolerance) != -1)
        //                                        dissimilarity--;
        //
        //                                }
        //
        //
        //        //                        DEBUG_MSG(dissimilarity);
        //
        //
        //
        //        //                        DEBUG_MSG(matrix[j][i].match_type);
        //
        //                                unsigned int matched_peaks = 0;
        //
        //                                if (dissimilarity <= 1) {
        //
        //                                        float pcos = spec1.scoreMatch2(spec2,
        //                                                                        tolerance,
        //                                                                        matched_peaks,
        //                                                                        score1,
        //                                                                        score2,
        //                                                                        true,
        //                                                                        true,
        //                                                                        false);
        //
        //
        //    //                                    DEBUG_MSG(pcos);
        //                                        if (pcos > sat_threshold && matched_peaks >= req_peak_match && min(score1,score2) > .5) {
        //
        //                                            for (mem = all_clusters[idx2].core.begin(); mem != all_clusters[idx2].core.end(); ++mem) {
        //                                                all_clusters[idx1].mix.push_back((*mem));
        //                                            }
        //                                            for (mem = all_clusters[idx2].sat.begin(); mem != all_clusters[idx2].sat.end(); ++mem) {
        //                                                all_clusters[idx1].mix.push_back((*mem));
        //                                            }
        //                                            for (mem = all_clusters[idx2].wsat.begin(); mem != all_clusters[idx2].wsat.end(); ++mem) {
        //                                                all_clusters[idx1].mix.push_back((*mem));
        //                                            }
        //                                            for (mem = all_clusters[idx2].mix.begin(); mem != all_clusters[idx2].mix.end(); ++mem) {
        //                                                all_clusters[idx1].mix.push_back((*mem));
        //                                            }
        //
        //                                            considered_clust.insert(idx2);
        //                                            merged_clust.insert(idx2);
        //
        //                                        }
        //                                }
        //
        //                            }
        //
        //
        //                        }
        //                    }
        //                }
        //
        //            }
        //
        //        }
        //

//        m_representativeSpectra->resize(1);
//        (*m_representativeSpectra)[0] = *all_clusters[0].core[0].spec;
//

        DEBUG_MSG("OUTPUT PHASE 4")

        output_clusters(m_params.getValue("OUTPUT_CLUSTERS_PHASE4").c_str(), m_params.getValue("OUTPUT_SPECTRA_PHASE4").c_str(), cluster_basename.c_str(), cluster_file_base.c_str(), all_clusters, component_count, *m_identifications, *m_globalIndex, *m_fileIndex, considered_clust);

//
//        cosine.clear();
//        proj_cosine.clear();
//
//        for (unsigned int i = 0; i < matrix.size(); ++i) {
//            for (unsigned int j = 0; j < matrix.size(); ++j) {
//                matrix[i][j] = None;
//            }
//        }
//
//        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
//
//
//                float isotope_increment = (float) 1;
//                //                DEBUG_MSG(spec1->scan);
//                Spectrum spec1 = *all_clusters[idx1].core[0].spec;
//                i = idx1;
//                top_peaks = spec1;
//                top_peaks.rankFilterPeaks(1, 30);
//                top_peaks.selectTopK(5, 0x0);
//
//                for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {
//
//
//                        Spectrum spec2 = *all_clusters[idx2].core[0].spec;
//
//                        j = idx2;
//
//                        //                    DEBUG_MSG(spec2->scan);
//
//                        if (idx1 != idx2) {
//
//                            //check pm isotopes
//
//                            float mass_difference = abs(spec1.parentMass - spec2.parentMass);
//                            float charge_tolerance = tolerance * spec1.parentCharge;
//
//
//                            if ((spec1.parentCharge == spec2.parentCharge) &&
//                                    ((mass_difference <= charge_tolerance)
//                                    || (mass_difference <= charge_tolerance + isotope_increment && mass_difference >= isotope_increment - charge_tolerance)
//                                    || (mass_difference <= charge_tolerance + 2 * isotope_increment && mass_difference >= 2 * isotope_increment - charge_tolerance))) {
//
//                                //calculate one way dissimilarity
//
//                                //                        DEBUG_MSG("PM match, calculating dissimilarity");
//
//                                int dissimilarity = top_peaks_n;
//
//                                //                        DEBUG_MSG(top_peaks.size());
//
//                                for (size_t p = 0; p < top_peaks.size(); ++p) {
//
//                                    if (spec2.findPeaks(top_peaks[p][0], peak_tolerance) != -1)
//                                        dissimilarity--;
//
//                                }
//
//
//                                //                        DEBUG_MSG(dissimilarity);
//
//
//
//                                //                        DEBUG_MSG(matrix[j][i].match_type);
//
//                                unsigned int matched_peaks = 0;
//
//                                if (dissimilarity <= 3) {
//                                    if (matrix[j][i] == None) {
//                                        float cos = spec1.scoreMatch2(spec2,
//                                                peak_tolerance,
//                                                matched_peaks,
//                                                score1,
//                                                score2,
//                                                true,
//                                                false,
//                                                false);
//                                        float pcos = spec1.scoreMatch2(spec2,
//                                                peak_tolerance,
//                                                matched_peaks,
//                                                score1,
//                                                score2,
//                                                true,
//                                                true,
//                                                false);
//                                        //                                    DEBUG_MSG(pcos);
//                                        if (pcos > sat_threshold && matched_peaks >= req_peak_match && min(score1, score2) > .4) {
//                                            matrix[i][j] = Subset;
//                                            matrix[j][i] = Superset;
//                                            cosine[i][j] = cos;
//                                            cosine[j][i] = cos;
//                                            proj_cosine[i][j] = pcos;
//                                            proj_cosine[j][i] = pcos;
//                                        }
//                                    } else {
//                                        if (cosine[i][j] > sat_threshold) {
////                                            proj_cosine[i][j] = cosine[i][j];
////                                            proj_cosine[j][i] = cosine[j][i];
//                                            matrix[i][j] = Exact;
//                                            matrix[j][i] = Exact;
//                                        }
//                                    }
//
//                                }
//
//                                //                        DEBUG_MSG("Updated matrix");
//
//
//                            } else {
//                                //                        DEBUG_MSG("No PM match, moving on...");
//                            }
//                        }
//
//                    }
//
//        }
//
//        for (unsigned int idx1 = 0; idx1 < all_clusters.size(); ++idx1) {
//                if (considered_clust.count(idx1) == 0) {
//
//                    j = idx1;
//
//                    float best_cluster_cos = 0;
//                    int best_cluster_idx = -1;
//                    MatchType min_match_w = Subset;
//                    MatchType min_match_m = Superset;
//
//
//                    for (unsigned int idx2 = 0; idx2 < all_clusters.size(); ++idx2) {
//
//                        if (j != idx2 && considered_clust.count(idx2) == 0) {
//
//                            i = idx2;
//
//                            unsigned int exact_match_edges = 0;
//
//                            float best_edge = 0;
//
//                            for (unsigned int core1 = 0; core1 < all_clusters[idx2].sub_clusters.size(); core1++) {
//
//                                float local_best_edge = 0;
//                                unsigned int total_edges = 0;
//
//                                for (unsigned int core2 = 0; core2 < all_clusters[j].sub_clusters.size(); core2++) {
//
//                                    unsigned int i = all_clusters[idx2].sub_clusters[core1].second;
//                                    unsigned int k = all_clusters[j].sub_clusters[core2].second;
//
//                                    if (matrix[i][k] != None) {
//                                        exact_match_edges++;
//                                        total_edges += all_clusters[j].sub_clusters[core2].first;
//                                    }
//
//                                    if ((matrix[i][k] == Exact || matrix[i][k] == min_match_w || matrix[i][k] == min_match_m) && not_core.count(k) == 0 && not_core.count(i) == 0) {
//                                        if (proj_cosine[i][k] > local_best_edge) {
//                                            local_best_edge = proj_cosine[i][k];
//                                        }
//                                        if (matrix[i][k] == Exact) {
//                                            min_match_w = Exact;
//                                            min_match_m = Exact;
//                                        }
//                                    }
//
//
//                                }
//
//                                best_edge += local_best_edge * total_edges;;
//
//                            }
//
//                            if (best_edge >= sat_threshold
//                                && best_edge > best_cluster_cos) {
//                                best_cluster_idx = idx2;
//                                best_cluster_cos = best_edge;
//                            }
//
//                        }
//
//                    }
//
//                    if (best_cluster_idx != -1) {//
//                        unsigned int bigger = 0;
//                        unsigned int smaller = 0;
//                        if (all_clusters[best_cluster_idx].core.size() + all_clusters[best_cluster_idx].sat.size() < all_clusters[j].core.size() + all_clusters[j].sat.size()) {
//                           bigger = j;
//                           smaller = best_cluster_idx;
//                        } else {
//                            bigger = best_cluster_idx;
//                            smaller = j;
//                        }
//
//                        for (unsigned int ic = 0; ic < all_clusters[smaller].sub_clusters.size(); ic++) {
//                            all_clusters[bigger].sub_clusters.push_back(all_clusters[smaller].sub_clusters[ic]);
//                            not_core.insert(all_clusters[smaller].sub_clusters[ic].second);
//                        }
//
//                        if (matrix[bigger][smaller] == Subset) {
//
//                            for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
//                                all_clusters[bigger].mix.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
//                                all_clusters[bigger].mix.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
//                                all_clusters[bigger].mix.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
//                                all_clusters[bigger].mix.push_back((*mem));
//                            }
//                        } else {
//                            for (mem = all_clusters[smaller].core.begin(); mem != all_clusters[smaller].core.end(); ++mem) {
//                                all_clusters[bigger].wsat.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].sat.begin(); mem != all_clusters[smaller].sat.end(); ++mem) {
//                                all_clusters[bigger].wsat.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].wsat.begin(); mem != all_clusters[smaller].wsat.end(); ++mem) {
//                                all_clusters[bigger].wsat.push_back((*mem));
//                            }
//                            for (mem = all_clusters[smaller].mix.begin(); mem != all_clusters[smaller].mix.end(); ++mem) {
//                                all_clusters[bigger].mix.push_back((*mem));
//                            }
//                        }
//
//                        considered_clust.insert(smaller);
//                    }
//                }
//
//            }

        m_representativeSpectra->resize(all_clusters.size() - considered_clust.size());

        unsigned int counter = 0;

        for (unsigned int idx = 0; idx < all_clusters.size(); ++idx) {
            if (considered_clust.count(idx) == 0) {
                (*m_representativeSpectra)[counter] = *all_clusters[idx].core[0].spec;
                (*m_representativeSpectra)[counter].scan = counter + 1;
                counter++;
            }
        }

        m_representativeSpectra->SaveSpecSet_mgf(representatives_name.c_str(), true);







        //        DEBUG_MSG(component_count);

        return true;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::loadInputData(void) {

        tr1::unordered_map<string, tr1::unordered_map<unsigned int, unsigned int> > tag_index;

        unsigned int cc_number;

        float low_mass = 200;
        float snr = 3;

        if (m_params.exists("TAG_CLUSTERS")) {
            string line;
            unsigned int scan;
            string fileName;
            string peptide;
            unsigned int tag_component;
            ifstream tag_clusters_file(m_params.getValue("TAG_CLUSTERS").c_str());
            if (tag_clusters_file.is_open()) {
                getline(tag_clusters_file, line);
                cc_number = atoi(line.c_str());
                m_inputComponents->resize(cc_number);

                //                DEBUG_MSG(m_inputComponents->size())

                while (getline(tag_clusters_file, line, '\t')) {
                    tag_component = atoi(line.c_str());
                    getline(tag_clusters_file, fileName, '\t');

                    getline(tag_clusters_file, line, '\t');
                    scan = atoi(line.c_str());

                    getline(tag_clusters_file, peptide);

                    tag_index[fileName][scan] = tag_component;

                    //                    DEBUG_MSG(fileName);
                    //                    DEBUG_MSG(scan);

                    (*m_identifications)[fileName][scan] = peptide;
                }
            }
            tag_clusters_file.close();
        }

        vector<Spectrum>::iterator spec;

        vector<SpecSet>::iterator component;


        for (component = m_inputComponents->begin(); component != m_inputComponents->end(); ++component) {
            component->resize(0);
        }

        if (m_params.exists("INPUT_SPECS")) {
            Spectrum filtered_peaks;
            vector<string> filenames;
            vector<string>::iterator filenames_iter;

            stringSplit(m_params.getValue("INPUT_SPECS"), filenames, ";");

            for (filenames_iter = filenames.begin(); filenames_iter != filenames.end(); ++filenames_iter) {
                m_inputSpectra->Load((*filenames_iter).c_str());
                unsigned int fileIndex = 0;
                for (spec = m_inputSpectra->begin(); spec != m_inputSpectra->end(); ++spec) {
                    //                    DEBUG_MSG(spec->fileName);
                    //                    DEBUG_MSG(spec->scan);
                    if (tag_index.count(spec->fileName) != 0) {
                        tr1::unordered_map<unsigned int, unsigned int> file_dict = tag_index[spec->fileName];
                        if (file_dict.count(spec->scan) != 0) {
                            spec->filterLowMassPeaks(low_mass);
                            spec->removeChargeReducedPrecursors(10, false, 0);
                            spec->filterLowSNR(snr);
                            spec->rankFilterPeaks(10, 100);
                            filtered_peaks = *spec;
                            filtered_peaks.rankFilterPeaks(1, 3);
                            if (filtered_peaks.size() >= 10) {
                                (*m_inputComponents)[file_dict[spec->scan]].push_back(*spec);
                                (*m_fileIndex)[spec->fileName][spec->scan] = fileIndex;
                            }
                        }
                    }
                    fileIndex++;
                }
            }

            //            vector<SpecSet>::iterator component;

            unsigned int read_spec = 0;

            unsigned int globalIndex = 0;

            for (component = m_inputComponents->begin(); component != m_inputComponents->end(); ++component) {
                unsigned int index = 0;
                for (spec = component->begin(); spec != component->end(); ++spec, ++index) {
                    //                    filtered_peaks = *spec;
                    //                    filtered_peaks.rankFilterPeaks(1,5);
                    for (unsigned int j = 0; j < spec->size(); j++)
                        (*spec)[j][1] = sqrt((*spec)[j][1]);
                    spec->normalize2();
                    (*m_index)[spec->fileName][spec->scan] = index;
                    (*m_globalIndex)[spec->fileName][spec->scan] = globalIndex;
                    globalIndex++;
                    read_spec++;
                }
            }

            m_filteredSpectra->resize(globalIndex);
            m_dissimilaritySpectra->resize(globalIndex);

//            DEBUG_MSG(read_spec);
        }

        return true;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::saveOutputData(void) {
        return true;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::saveInputData(std::vector<std::string> & filenames) {
        return true;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::loadOutputData(void) {
        return true;
    }

    // -------------------------------------------------------------------------

    vector<ExecBase*> const & ExecClusterBins::split(int numSplit) {
        return m_subModules;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::merge(void) {
        return true;
    }

    // -------------------------------------------------------------------------

    bool ExecClusterBins::validateParams(std::string & error) {
        return true;
    }


} /* namespace specnets */
