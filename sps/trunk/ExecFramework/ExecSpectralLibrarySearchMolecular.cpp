//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"
#include "SpectrumPairSet.h"

// Header Include
#include "ExecSpectralLibrarySearchMolecular.h"
#include "ExecFdrPeptide.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectralLibrarySearchMolecular::ExecSpectralLibrarySearchMolecular(void){
    m_name = "ExecSpectralLibrarySearchMolecular";
    m_type = "ExecSpectralLibrarySearchMolecular";
  }


  ExecSpectralLibrarySearchMolecular::ExecSpectralLibrarySearchMolecular(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySearchMolecular";
    m_type = "ExecSpectralLibrarySearchMolecular";
    //inputParams.print(cout);

    this->m_params = inputParams;

    //Mass Tolerance when searching
    m_search_parentmass_tolerance = inputParams.getValueFloat("search_parentmass_tolerance",
                                                              inputParams.getValueFloat("tolerance.PM_tolerance", 2));

    m_search_peak_tolerance = inputParams.getValueFloat("search_peak_tolerance",
                                                              inputParams.getValueFloat("tolerance.Ion_tolerance", 0.5));


    //Specifying output psm file
    m_output_psm_filename = inputParams.getValue("RESULTS_DIR", "");

    string output_prefix = inputParams.getValue("RESULTS_DIR_PREFIX", "");
    if(output_prefix.length() > 0){
        m_output_psm_filename = get_only_path(output_prefix) + "/" + get_only_filename(m_output_psm_filename);
    }


    m_score_threshold = inputParams.getValueFloat("score_threshold", inputParams.getValueFloat("SCORE_THRESHOLD", 0.50));

    m_minimum_matched_peaks = inputParams.getValueFloat("MIN_MATCHED_PEAKS_SEARCH", inputParams.getValueFloat("MIN_MATCHED_PEAKS", 0));

    m_top_hits_reported = inputParams.getValueInt("TOP_K_RESULTS", 1);

    //Node Idx and count
    m_node_idx = inputParams.getValueInt("NODEIDX", 0);
    m_node_count = inputParams.getValueInt("NODECOUNT", 1);

    m_analog_search = inputParams.getValueInt("ANALOG_SEARCH",0);

    if(m_analog_search > 0){
        m_analog_search = inputParams.getValueInt("MAX_SHIFT_MASS",0);
    }

    m_library_search_quality = inputParams.getValueInt("SEARCH_LIBQUALITY", 100);

    if(inputParams.getValueInt("FORCE_EXACT_MATCH", 0) == 0){
        m_force_exact_match = false;
    }
    else{
        m_force_exact_match = true;
    }




    if(inputParams.getValueInt("CCMS", 0) == 0){
        //Search Spectra
        stringSplit(inputParams.getValue("searchspectra"), m_input_query_spectra);

        //Annotated MGF files most likely comes from some other spectral library
        //And the annotations are embedded in the file, so don't bother with other
        //annotation file
        string input_annotated_mgf_file = inputParams.getValue("EXISTING_LIBRARY_MGF");
        //Splitting inputs
        stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    }

    if(inputParams.getValueInt("CCMS", 0) == 1){
        //Reading Mapping
        //Reading in file mapping
        int mapping_count = 0;
        map<string,string> file_name_mapping;
        m_mgf_file_names_string = "";
        while(1){
            char buf[100];
            sprintf(buf, "upload_file_mapping%i", mapping_count);

            mapping_count++;

            if(!inputParams.exists(buf))
                break;

            std::string mapping = inputParams.getValue(buf);
            std::string mangled_name = mapping.substr(0, mapping.find("|"));
            std::string original_name = mapping.substr(mapping.find("|")+1);


            if(inputParams.getValueInt("PROTEOSAFE_SPECNETS_SEARCH", 0) == 0){
                if(mangled_name.find("spec-") != string::npos || mangled_name.find("speccontinuous-") != string::npos){
                    file_name_mapping[original_name] = mangled_name;

                    std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
                    m_input_query_spectra.push_back(path_to_spectra + "/" + mangled_name);
                }
            }
            if(mangled_name.find("lib-") != string::npos){
                file_name_mapping[original_name] = mangled_name;
                std::string path_to_lib = inputParams.getValue("EXISTING_LIBRARY_MGF");
                m_mgf_file_names.push_back(path_to_lib + "/" + mangled_name);
                m_mgf_file_names_string += path_to_lib + "/" + mangled_name + " ";
            }
        }

        if(inputParams.getValueInt("PROTEOSAFE_SPECNETS_SEARCH", 0) == 1){
            stringSplit(inputParams.getValue("searchspectra"), m_input_query_spectra);
        }
    }
  }


  ExecSpectralLibrarySearchMolecular::~ExecSpectralLibrarySearchMolecular(void){
  }


  ExecBase * ExecSpectralLibrarySearchMolecular::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySearchMolecular(inputParams);
  }

  bool ExecSpectralLibrarySearchMolecular::invoke(void){
    vector<int> charge_filter;

    ParameterList temp_params;
    temp_params.setValue("FILTER_PRECURSOR_WINDOW", m_params.getValue("FILTER_PRECURSOR_WINDOW", "0"));
    temp_params.setValue("FILTER_STDDEV_PEAK_INT", m_params.getValue("FILTER_STDDEV_PEAK_INT", "0"));
    temp_params.setValue("MIN_PEAK_INT", m_params.getValue("MIN_PEAK_INT", "0"));
    temp_params.setValue("WINDOW_FILTER", m_params.getValue("WINDOW_FILTER", "0"));
    temp_params.setValue("FILTER_SNR_PEAK_INT", m_params.getValue("FILTER_SNR_PEAK_INT", "0.0"));

    DEBUG_MSG("QTOF Filtering Input");
    qtof_filtering_spectrum(searchable_spectra, temp_params);

    for(int i = 0; i < m_libraries.size(); i++){
        //Filtering spectra because of qtof data is messy
        if(m_params.getValueInt("FILTER_LIBRARY", 0) == 1){
            qtof_filtering_spectrum(m_libraries[i], temp_params);
            DEBUG_MSG("QTOF Filtering Library");
        }
        DEBUG_MSG("Done Creating Library Size: " <<m_libraries[i].size());
    }



    //Doing the search
    DEBUG_MSG("SEARCHING.............................");
    vector<psmPtr> all_target_search_results;


    vector<int> accepted_fragmentation;
    accepted_fragmentation.push_back(Spectrum::FragType_CID);

    //Calculating ranges
    int spectra_count = searchable_spectra.size();
    int chunk_size = spectra_count/m_node_count;
    int start_idx = chunk_size * m_node_idx;
    int end_idx = chunk_size * (m_node_idx+1);
    if(m_node_idx == (m_node_count-1)) end_idx = spectra_count;

    DEBUG_MSG("STARTING\t"<<start_idx<<"\t"<<end_idx<<"\t"<<m_node_idx<<"\t"<<m_node_count<<"\t"<<searchable_spectra.size());


    for(int libraries_idx = 0 ; libraries_idx < m_libraries.size(); libraries_idx++){
        DEBUG_MSG("SEARCHING LIB: "<<libraries_idx);
        for(int query_idx = start_idx; query_idx < end_idx; query_idx++){
            //DEBUG_VAR(searchable_spectra[query_idx].parentMass);
            //if (query_idx % 100 == 0 ) DEBUG_MSG("SEARCHING QUERY IDX "<<query_idx<<" of "<<end_idx<< " in library idx " <<libraries_idx<< " of "<<m_libraries.size()<<" with "<<m_all_search_results.size()<<" results so far");
            DEBUG_MSG("SEARCHING QUERY IDX "<<query_idx<<" of "<<end_idx<< " in library idx " <<libraries_idx<< " of "<<m_libraries.size()<<" with "<<m_all_search_results.size()<<" results so far");

            //Skipping empty spectra
            if(searchable_spectra[query_idx].size() == 0)
                continue;

            int target_search_ret = m_libraries[libraries_idx].search(searchable_spectra[query_idx],
                                        m_all_search_results,
                                        m_search_parentmass_tolerance,
                                        m_search_peak_tolerance,
                                        m_top_hits_reported,
                                        m_analog_search,
                                        m_score_threshold,
                                        m_library_search_quality,
                                        m_minimum_matched_peaks,
                                        m_force_exact_match);

        }
    }


    DEBUG_MSG("FILTERING SEARCH RESULTS");
    PeptideSpectrumMatchSet temp_results;
    for(int i = 0; i < m_all_search_results.size(); i++){
        if(m_all_search_results[i]->m_score < m_score_threshold || m_all_search_results[i]->m_shared_peaks < m_minimum_matched_peaks && m_minimum_matched_peaks > 0){
            //m_all_search_results.removePsmSetItem(m_all_search_results[i]);
            //i--;
        }
        else{
            temp_results.push_back(m_all_search_results[i]);
        }
    }
    m_all_search_results = temp_results;

    //Filtering the results again for Unique IDs for each ID
    map<string, psmPtr> top_scoring_psms;
    for(int i = 0; i < m_all_search_results.size(); i++){
        stringstream ss1 (stringstream::in | stringstream::out);
        ss1 << m_all_search_results[i]->m_spectrumFile;
        ss1 << "____";
        ss1 << m_all_search_results[i]->m_scanNum;
        string spectrum_key = ss1.str();

        DEBUG_MSG(spectrum_key);

        map<string, psmPtr>::iterator it = top_scoring_psms.find(spectrum_key);
        if(it == top_scoring_psms.end()){
            top_scoring_psms[spectrum_key] = m_all_search_results[i];
        }
        else{
            if(it->second->m_score < m_all_search_results[i]->m_score){
                top_scoring_psms[spectrum_key] = m_all_search_results[i];
            }
        }
    }

    PeptideSpectrumMatchSet temp_results2;

    for (map<string, psmPtr>::iterator it=top_scoring_psms.begin(); it!=top_scoring_psms.end(); ++it){
        temp_results2.push_back(it->second);
    }
    m_all_search_results = temp_results2;

    //removing sentinal not peptide annotations
    for(int i = 0; i < m_all_search_results.size(); i++){
        if(m_all_search_results[i]->m_annotation == "*.NOTPEPTIDE.*" || m_all_search_results[i]->m_annotation == "NOTPEPTIDE" || m_all_search_results[i]->m_annotation.length() < 2){
            m_all_search_results[i]->m_annotation = "*..*";
        }
    }

    return true;
  }

  bool ExecSpectralLibrarySearchMolecular::loadInputData(void){



    DEBUG_MSG("Loading Libraries");
    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        SpectralLibrary new_loaded_lib;

        DEBUG_MSG("LOADING TARGET MGF\n");
        string mgf_file_name = m_mgf_file_names[file_idx];
        new_loaded_lib.LoadSpecSet_mgf(mgf_file_name.c_str());
        DEBUG_MSG("LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file. Size: "<<new_loaded_lib.size());;

        for(int lib_idx = 0; lib_idx < new_loaded_lib.size(); lib_idx++){
            if(new_loaded_lib[lib_idx].psmList.size() == 0){
                //Creating new blank PSM
                psmPtr psm(new PeptideSpectrumMatch);
                psm->m_annotation = "*.NOTPEPTIDE.*";
                psm->m_spectrumFile = mgf_file_name;
                psm->m_dbIndex = lib_idx + 1;
                new_loaded_lib[lib_idx].psmList.push_back(psm);
            }
        }

        m_libraries.push_back(new_loaded_lib);
    }


    DEBUG_MSG("Loading Search Spectra");

    //Loading searchable spectral
    //Looping through mzxml/pklbin files to load
    for(int file_idx = 0; file_idx < m_input_query_spectra.size(); file_idx++){
        string spectra_file_name = m_input_query_spectra[file_idx];
        string extension = get_extension(spectra_file_name);

        SpecSet temp_specs;
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0 || strcmp(extension.c_str(), "mzML") == 0){
            string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
            temp_specs.Load(pklbin_version.c_str());
            //spectra_file_name = pklbin_version;
            DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            temp_specs.loadPklBin(spectra_file_name.c_str());
            DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a pklbin file");
        }
        else if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
            temp_specs.LoadSpecSet_mgf(spectra_file_name.c_str());
            DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a mgf file");
        }
        else{
            DEBUG_MSG("CANNOTLOAD Search\t"<<spectra_file_name);
        }

        //Writing the file name into the spectra
        for(int spec_idx = 0; spec_idx < temp_specs.size(); spec_idx++){
            temp_specs[spec_idx].fileName = spectra_file_name;
        }

        searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());
    }




    return true;
  }


  bool ExecSpectralLibrarySearchMolecular::saveOutputData(void){
    DEBUG_MSG("OUTPUT PSM\t"<<m_output_psm_filename);

    if(m_output_psm_filename.length() > 0){
        m_all_search_results.saveToFile(m_output_psm_filename.c_str(), true);
    }

    //Saving out spectrum pairset
    string pair_setoutput_name = this->m_params.getValue("PAIRSET_OUTPUT_FILE", "");
    if(pair_setoutput_name.length() > 0){
        SpectrumPairSet pairSet;
        for(int i = 0; i < m_all_search_results.size(); i++){

            int library_idx = m_all_search_results[i]->m_dbIndex;
            int query_idx = m_all_search_results[i]->m_scanNum;
            float cosine_score = m_all_search_results[i]->m_score;

            int library_offset = 10000000;
            //if(searchable_spectra.size() > 0){
            //    library_offset = searchable_spectra.size();
            //}

            SpectrumPair newPair;

            newPair.spec2 = library_idx-1 + library_offset;
            newPair.spec1 = query_idx -1;
            newPair.score1 = cosine_score;
            pairSet.push_back(newPair);
        }

        DEBUG_MSG("SORTING PAIRSET");
        //Sorting
        pairSet.sort_pairs_by_index();
        //pairSet.sort_pairs();
        DEBUG_MSG("SAVING PAIRSET");
        //Saving
        pairSet.saveToBinaryFile(pair_setoutput_name);
    }

    DEBUG_MSG("DONE SAVING PAIRSET");

    return true;
  }


  bool ExecSpectralLibrarySearchMolecular::saveInputData(std::vector<std::string> & filenames){
    DEBUG_MSG("SAVEINPUT");
    string dataDir = m_params.getValue("GRID_DATA_DIR_PARAMS", ".");
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);
    DEBUG_VAR(dataDir);
    DEBUG_VAR(paramFilename);
    filenames.push_back(paramFilename); // Parameter file MUST be first in vector

    return true;
  }

  bool ExecSpectralLibrarySearchMolecular::loadOutputData(void){
    DEBUG_MSG("LOADING OUTPUT DATA\t"<<getName()<<"\t"<<m_output_psm_filename);
    m_all_search_results.loadFromFile(m_output_psm_filename.c_str());
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySearchMolecular::split(int numSplit){
    DEBUG_MSG("SPLIT\t"<<numSplit<<"==================================================");

    //Lets assume split num is just a guideline and not to be strictly followed.
    int num_files = m_input_query_spectra.size();

    if(num_files > 1)
        numSplit = min(num_files,numSplit);


    for(int i = 0; i < numSplit; i++){
        ParameterList childParams(m_params);


        string dataDir = m_params.getValue("INTERMEDIATE_RESULTS_DIR", "intermediateresults");
        string baseDirectory = dataDir + "/";
        string baseFilename = baseDirectory + makeName(getName(), i);
        string psmFilename = baseFilename + ".psms";
        DEBUG_VAR(psmFilename);
        childParams.setValue("RESULTS_DIR", psmFilename);

        childParams.removeParam("GRID_EXECUTION");
        childParams.removeParam("CCMS");

        stringstream ss1 (stringstream::in | stringstream::out);
        stringstream ss2 (stringstream::in | stringstream::out);

        if(num_files == 1){
            //Lets split up the search of this one file
            ss1 << i;
            childParams.setValue("NODEIDX", ss1.str());

            ss2<<numSplit;
            childParams.setValue("NODECOUNT", ss2.str());

            childParams.setValue("searchspectra", m_input_query_spectra[0]);
            childParams.setValue("EXISTING_LIBRARY_MGF", m_mgf_file_names_string);
        }
        if(num_files > 1){
            //NODE 0 of 1 always
            ss1 << 0;
            childParams.setValue("NODEIDX", ss1.str());

            ss2<<1;
            childParams.setValue("NODECOUNT", ss2.str());

            //Set the file name to be searched appropriately
            int files_per_node = (num_files / numSplit) + 1;
            if(num_files == numSplit)
                files_per_node = 1;

            int starting_file_idx = i * files_per_node;
            int ending_file_idx = (i+1) * files_per_node;

            if(ending_file_idx >= num_files)
                ending_file_idx = num_files;

            string search_spectra_output = "";

            for(int file_idx = starting_file_idx; file_idx < ending_file_idx; file_idx++){
                search_spectra_output += m_input_query_spectra[file_idx] + " ";
            }

            DEBUG_VAR(search_spectra_output);

            childParams.setValue("searchspectra", search_spectra_output);
            childParams.setValue("EXISTING_LIBRARY_MGF", m_mgf_file_names_string);
        }


        ExecBase * theClone = new ExecSpectralLibrarySearchMolecular(childParams);

        theClone->setName(makeName(m_name, i));
        m_subModules.push_back(theClone);
    }

    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecSpectralLibrarySearchMolecular::merge(void){
    DEBUG_MSG("MERGING\t"<<m_subModules.size());
    for(int i = 0; i < m_subModules.size(); i++){
        for(int j = 0; j < ((ExecSpectralLibrarySearchMolecular*)m_subModules[i])->m_all_search_results.size(); j++){
            m_all_search_results.push_back(((ExecSpectralLibrarySearchMolecular*)m_subModules[i])->m_all_search_results[j]);
        }
    }
    DEBUG_MSG("MERGED SIZE\t"<<m_all_search_results.size());
    /*
    DEBUG_MSG("TOP K PER ORGANISM FILTERING");
    //Finding top K per organism
    map<string, PeptideSpectrumMatchSet> organism_psm_map;
    for(int i = 0 ; i < m_all_search_results.size(); i++){
        stringstream ofs;
        string organism;
        if(m_all_search_results[i]->m_organism.size() == 0){
            organism = "NONE";
        }
        else{
            organism = m_all_search_results[i]->m_organism[0];
        }

        ofs << m_all_search_results[i]->m_spectrumFile << m_all_search_results[i]->m_scanNum<<organism;
        string spectrum_unique_string = ofs.str();
        //DEBUG_MSG(spectrum_unique_string);

        map<string, PeptideSpectrumMatchSet>::iterator it;
        it = organism_psm_map.find(spectrum_unique_string);
        if(it == organism_psm_map.end()){
            //DEBUG_MSG("NOT FOUND");
            PeptideSpectrumMatchSet psmSet;
            organism_psm_map[spectrum_unique_string] = psmSet;
            organism_psm_map[spectrum_unique_string].push_back(m_all_search_results[i]);
        }
        else{
            //DEBUG_MSG("FOUND");
            organism_psm_map[spectrum_unique_string].push_back(m_all_search_results[i]);

        }
    }

    //Iterating through
    PeptideSpectrumMatchSet output_psms;
    for (map<string, PeptideSpectrumMatchSet>::iterator it=organism_psm_map.begin(); it!=organism_psm_map.end(); ++it){
        DEBUG_MSG(it->first);
        DEBUG_MSG(it->second.size());
        sort(it->second.m_psmSet.begin(), it->second.m_psmSet.end(), search_results_comparator_psmPtr);
        for(int i = 0; i < (min( (unsigned int)m_top_hits_reported, it->second.size())); i++){
            output_psms.push_back(it->second[i]);
        }
    }*/
    //DEBUG_VAR(output_psms.size());
    //DEBUG_VAR(m_all_search_results.size());
    //m_all_search_results = output_psms;

    return true;
  }


  bool ExecSpectralLibrarySearchMolecular::validateParams(std::string & error){
    return true;
  }


}
