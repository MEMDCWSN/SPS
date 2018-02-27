//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibrarySearch.h"
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

  ExecSpectralLibrarySearch::ExecSpectralLibrarySearch(void){
    m_name = "ExecSpectralLibrarySearch";
    m_type = "ExecSpectralLibrarySearch";
  }


  ExecSpectralLibrarySearch::ExecSpectralLibrarySearch(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySearch";
    m_type = "ExecSpectralLibrarySearch";
    //inputParams.print(cout);
    
    this->m_params = inputParams;
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("EXISTING_LIBRARY_MGF");
    string input_annotated_mgf_decoy_file = inputParams.getValue("EXISTING_LIBRARY_MGF_DECOY");

    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
    
    
    //Mass Tolerance when searching
    m_search_parentmass_tolerance = inputParams.getValueFloat("search_parentmass_tolerance", 
                                                              inputParams.getValueFloat("tolerance.PM_tolerance", 2));
                                                              
    m_search_peak_tolerance = inputParams.getValueFloat("search_peak_tolerance", 
                                                              inputParams.getValueFloat("tolerance.Ion_tolerance", 0.5));
    
    //Specfiying amino acid masses file
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile", "");
    
    //Specifying output psm file
    m_output_psm_filename = inputParams.getValue("RESULTS_DIR", "");
    
    string output_prefix = inputParams.getValue("RESULTS_DIR_PREFIX", "");
    if(output_prefix.length() > 0){
        m_output_psm_filename = get_only_path(output_prefix) + "/" + get_only_filename(m_output_psm_filename);
    }
    
    //Specify whether to search the decoy
    m_do_decoy_search = inputParams.getValueInt("search_decoy", 0);
    
    m_score_threshold = inputParams.getValueFloat("score_threshold", inputParams.getValueFloat("SCORE_THRESHOLD", 0.50));
    
    m_do_score_threshold = inputParams.getValueInt("use_score_threshold", 1);
    
    m_score_calculation_method = inputParams.getValueInt("search_scoring_method", 0);
    
    m_top_hits_reported = inputParams.getValueInt("TOP_K_RESULTS", 1);
    
    //Node Idx and count
    m_node_idx = inputParams.getValueInt("NODEIDX", 0);
    m_node_count = inputParams.getValueInt("NODECOUNT", 1);
    
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    stringSplit(input_annotated_mgf_decoy_file, m_mgf_decoy_file_names);
   
    if(inputParams.getValueInt("CCMS", 0) == 0){
        //Search Spectra
        m_input_projection_target_string = inputParams.getValue("searchspectra");
        stringSplit(m_input_projection_target_string, m_input_projection_target_spectra);    
    }
    
    if(inputParams.getValueInt("CCMS", 0) == 1){
        //Reading Mapping
        //Reading in file mapping
        int mapping_count = 0;
        map<string,string> file_name_mapping;
        while(1){
            char buf[100];
            sprintf(buf, "upload_file_mapping%i", mapping_count);
            
            mapping_count++;
            
            if(!inputParams.exists(buf))
                break;
            
            std::string mapping = inputParams.getValue(buf);
            std::string mangled_name = mapping.substr(0, mapping.find("|"));
            std::string original_name = mapping.substr(mapping.find("|")+1);
            
            if(mangled_name.find("spec-") == string::npos) continue;
            
            file_name_mapping[original_name] = mangled_name;
           
            std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
            m_input_projection_target_spectra.push_back(path_to_spectra + "/" + mangled_name);
        }
    }
  }


  ExecSpectralLibrarySearch::~ExecSpectralLibrarySearch(void){
  }


  ExecBase * ExecSpectralLibrarySearch::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySearch(inputParams);
  }

  bool ExecSpectralLibrarySearch::invoke(void){
    cout<<"Invoke"<<endl;
    vector<int> charge_filter;
    
    m_library.createlibrary(-1, 
                          -1, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          false);
                          
    cout<<"Done Creating Library"<<endl;
    if(m_do_decoy_search){
        //Creating a decoy spectral library
        if(m_mgf_decoy_file_names.size() == 0){
            cout<<"Generating Decoy"<<endl;
            m_decoy_library = m_library.create_decoy_spectral_library(model, ionsToExtract, allIons);
        }
        else{
            m_decoy_library.createlibrary(-1, 
                            -1, 
                            model, 
                            ionsToExtract, 
                            allIons, 
                            aminoacidexclusions, 
                            charge_filter,
                            false);
        }
        cout<<"Done Creating Decoy Library"<<endl;
    }
    
    
    //Doing the search
    cout<<"============================================="<<endl;
    cout<<"SEARCH"<<endl;
    vector<psmPtr> all_target_search_results;
    vector<psmPtr> all_decoy_search_results;
    
    cout<<searchable_spectra.size()<<endl;
    
    //Filtering spectra because of qtof data is messy
    
    ParameterList temp_params;
    temp_params.setValue("FILTER_PRECURSOR_WINDOW", "1");
    temp_params.setValue("FILTER_STDDEV_PEAK_INT", "2");
    temp_params.setValue("MIN_PEAK_INT", "50");
    qtof_filtering_spectrum(searchable_spectra, temp_params);
    qtof_filtering_spectrum(m_library, temp_params);
    
    
    
    //Creating Sorted library spectra pointers
    vector<Spectrum *> target_library_ptr;
    vector<Spectrum *> decoy_library_ptr;
    
    sorted_vector_library(target_library_ptr, m_library);
    sorted_vector_library(decoy_library_ptr, m_decoy_library);

    
    vector<int> accepted_fragmentation;
    accepted_fragmentation.push_back(Spectrum::FragType_CID);
    
    //Calculating ranges
    int spectra_count = searchable_spectra.size();
    int chunk_size = spectra_count/m_node_count;
    int start_idx = chunk_size * m_node_idx;
    int end_idx = chunk_size * (m_node_idx+1);
    if(m_node_idx == (m_node_count-1)) end_idx = spectra_count;
    
    DEBUG_MSG("STARTING\t"<<start_idx<<"\t"<<end_idx<<"\t"<<m_node_idx<<"\t"<<m_node_count<<"\t"<<searchable_spectra.size());
    
    
    
        
        
    //PeptideSpectrumMatchSet all_search_results;
    for(int query_idx = start_idx; query_idx < end_idx; query_idx++){
        //Filtering in acceptable fragmentation types
        DEBUG_MSG("SEARCHING QUERY SCAN "<<searchable_spectra[query_idx].scan<<"\t"<<searchable_spectra[query_idx].fileName);
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
        if(m_do_decoy_search == 1){
            psmPtr targetdecoy_psm(new PeptideSpectrumMatch);
            int target_decoy_search = m_library.search_target_decoy(m_decoy_library, 
                                                                searchable_spectra[query_idx], 
                                                                targetdecoy_psm, 
                                                                m_search_parentmass_tolerance, 
                                                                target_library_ptr, 
                                                                decoy_library_ptr, 
                                                                m_score_calculation_method);
            m_all_search_results.push_back(targetdecoy_psm);
        }
                                                               
        
        if(m_do_decoy_search == 0){
            
            int target_search_ret = m_library.search(searchable_spectra[query_idx], 
                                           m_all_search_results, 
                                           m_search_parentmass_tolerance, 
                                           m_search_peak_tolerance,
                                           m_top_hits_reported);
            
            
                                           
            /*targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
            targetdecoy_psm->m_spectrumFile = searchable_spectra[query_idx].fileName;
            float match_score = targetdecoy_psm->m_score;
            //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
            if( ((m_do_score_threshold == 0)) || match_score > m_score_threshold){
                //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                all_search_results.push_back(targetdecoy_psm);
            }
            */
        }
        //cout<<endl;
        
        continue;
    }
    
    if(m_do_decoy_search == 0){
        if(m_do_score_threshold){
            for(int i = 0; i < m_all_search_results.size(); i++){
                if(m_all_search_results[i]->m_score < m_score_threshold){
                    m_all_search_results.removePsmSetItem(m_all_search_results[i]);
                    i--;
                }
            }
        }
    }
    
    //removing sentinal not peptide annotations
    for(int i = 0; i < m_all_search_results.size(); i++){
        if(m_all_search_results[i]->m_annotation == "*.NOTPEPTIDE.*"){
            m_all_search_results[i]->m_annotation = "*..*";
        }
    }

    
    
    if(m_do_decoy_search){
        vector<string> fdr_outputs;
        for(float fdr = 0.f; fdr < 0.05; fdr += 0.001){
            ParameterList fdr_params;
            char buf[1000];
            sprintf(buf,"%f", fdr);
            
            fdr_params.setValue("PEPTIDE_FDR_CUTOFF", buf);
            fdr_params.setValue("TDA_TYPE", "concatenated");
            PeptideSpectrumMatchSet fdr_peptides;
            ExecFdrPeptide calculateFDR(fdr_params, &m_all_search_results, &fdr_peptides);
            if (!calculateFDR.invoke())
            {
                //DEBUG_MSG("Unable to generate fdr results!);
                return false;
            }
            
            stringstream fdr_output (stringstream::in | stringstream::out);
            
            fdr_output<<"FDR\t"<<fdr<<"\t"<<fdr_peptides.size()<<endl;
            fdr_outputs.push_back(fdr_output.str());
        }
        
        for(int i = 0; i < fdr_outputs.size();i++){
            cout<<fdr_outputs[i];
        }
    }
    
    
    
    
    
    return true;
  }

  bool ExecSpectralLibrarySearch::loadInputData(void){      
    
    load_aminoacid_masses();
      
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b-H2O-NH3");
    ionsToExtract.push_back("b-H2O-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("b++-H2O-H2O");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-iso");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y-H2O-NH3");
    ionsToExtract.push_back("y-H2O-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-H2O");
    ionsToExtract.push_back("y++-H2O-H2O");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-iso");
    ionsToExtract.push_back("P++");
    ionsToExtract.push_back("P++-H2O");
    ionsToExtract.push_back("P++-NH3");
    ionsToExtract.push_back("P++-H2O-H2O"); 
    
    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        cout<<"LOADING TARGET MGF\n"<<endl;
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    if(m_do_decoy_search){
        //Loading annotated MGF
        for(int file_idx = 0; file_idx < m_mgf_decoy_file_names.size(); file_idx++){
            string mgf_file_name = m_mgf_decoy_file_names[file_idx];
            cout<<"LOADING DECOY MGF\t"<<mgf_file_name<<endl;
            m_decoy_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
            cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf decoy file"<<endl;
        }
    }

    
    cout<<"Loading Search Spectra"<<endl;
    
    //Loading searchable spectral
    //Looping through mzxml/pklbin files to load
    for(int file_idx = 0; file_idx < m_input_projection_target_spectra.size(); file_idx++){
        string spectra_file_name = m_input_projection_target_spectra[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            SpecSet temp_specs;
            //cout<<"LOADED Search\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            
            string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
            temp_specs.Load(pklbin_version.c_str());
            DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
            //LoadMzxml(spectra_file_name.c_str(), temp_specs, 0, 2);
            //Writing the file name into the spectra
            for(int spec_idx = 0; spec_idx < temp_specs.size(); spec_idx++){
                temp_specs[spec_idx].fileName = spectra_file_name;
            }
            
            searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());
            //projection_spectra.LoadSpecSet_mzxml(spectra_file_name.c_str(), 2);
            //cout<<"MZXML not Supported"<<endl;
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            SpecSet temp_specs;
            temp_specs.loadPklBin(spectra_file_name.c_str());
            cout<<"LOADED Search\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
            //Writing the file name into the spectra
            for(int spec_idx = 0; spec_idx < temp_specs.size(); spec_idx++){
                temp_specs[spec_idx].fileName = spectra_file_name;
            }
            
            searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());

        }
        else if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
            cout<<"LOADED Search\t"<<spectra_file_name<<"\t as a mgf file"<<endl;
            SpecSet temp_specs;
            temp_specs.LoadSpecSet_mgf(spectra_file_name.c_str());
            for(int spec_idx = 0; spec_idx < temp_specs.size(); spec_idx++){
                temp_specs[spec_idx].fileName = spectra_file_name;
            }
            searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());
        }
        else{
            cerr<<"CANNOTLOAD Search\t"<<spectra_file_name<<endl;
        }
    }
    
    
    
    
    return true;
  }


  bool ExecSpectralLibrarySearch::saveOutputData(void){
    DEBUG_MSG("OUTPUT PSM\t"<<m_output_psm_filename);
    if(m_output_psm_filename.length() > 0){
        m_all_search_results.saveToFile(m_output_psm_filename.c_str(), true);
    }
    return true;
  }


  bool ExecSpectralLibrarySearch::saveInputData(std::vector<std::string> & filenames){
    DEBUG_MSG("SAVEINPUT");
    string dataDir = m_params.getValue("GRID_DATA_DIR", ".");
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);
    DEBUG_VAR(dataDir);
    DEBUG_VAR(paramFilename);
    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    
    return true;
  }

  bool ExecSpectralLibrarySearch::loadOutputData(void){
    DEBUG_MSG("LOADING OUTPUT DATA\t"<<getName()<<"\t"<<m_output_psm_filename);
    m_all_search_results.loadFromFile(m_output_psm_filename.c_str());
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySearch::split(int numSplit){
    DEBUG_MSG("SPLIT\t"<<numSplit<<"==================================================");
    
    //Lets assume split num is just a guideline and not to be strictly followed. 
    int num_files = m_input_projection_target_spectra.size();
    
    if(num_files > 1)
        numSplit = num_files;
    
    
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
            
            childParams.setValue("searchspectra", m_input_projection_target_spectra[0]);
        }
        if(num_files > 1){
            //NODE 0 of 1 always
            ss1 << 0;
            childParams.setValue("NODEIDX", ss1.str());
            
            ss2<<1;
            childParams.setValue("NODECOUNT", ss2.str());
            
            //Set the file name to be searched appropriately
            childParams.setValue("searchspectra", m_input_projection_target_spectra[i]);
        }
        
        
        ExecBase * theClone = new ExecSpectralLibrarySearch(childParams);
        
        theClone->setName(makeName(m_name, i));
        m_subModules.push_back(theClone);
    }
    
    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecSpectralLibrarySearch::merge(void){
    DEBUG_MSG("MERGING\t"<<m_subModules.size());
    for(int i = 0; i < m_subModules.size(); i++){
        for(int j = 0; j < ((ExecSpectralLibrarySearch*)m_subModules[i])->m_all_search_results.size(); j++){
            m_all_search_results.push_back(((ExecSpectralLibrarySearch*)m_subModules[i])->m_all_search_results[j]);
        }
    }
    DEBUG_MSG("MERGED SIZE\t"<<m_all_search_results.size());
    
    return true;
  }


  bool ExecSpectralLibrarySearch::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibrarySearch::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
