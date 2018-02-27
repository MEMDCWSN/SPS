//Module Includes
#include "SpectralLibrary.h"
#include "projectionutils.h"
#include "utils.h"

// Header Include
#include "ExecProjectionStatistics.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecProjectionStatistics::ExecProjectionStatistics(void){
    m_name = "ExecProjectionStatistics";
    m_type = "ExecProjectionStatistics";
  }


  ExecProjectionStatistics::ExecProjectionStatistics(const ParameterList & inputParams){
    m_name = "ExecProjectionStatistics";
    m_type = "ExecProjectionStatistics";
    inputParams.print(cout);
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");
    
    //Grabbing the file name of the output mgf file
    m_output_mgf_file_name = inputParams.getValue("output_mgf");
    
    //Grabbing thresholding filters
    m_pvalue_filter = inputParams.getValueFloat("pvalue_cutoff");
    m_envelope_filter = inputParams.getValueFloat("envelope_cutoff");
    
    
    //Pulling out parameters for creating the library
    //MZXML and annotation files come as a pair
    string input_spectrum_names = inputParams.getValue("rawnames");
    string input_psn_names = inputParams.getValue("psmnames");
    string input_target_spectrum_names = inputParams.getValue("targetrawnames");
    string input_target_psn_names = inputParams.getValue("targetpsmnames");
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf");
    string input_target_annotated_mgf_file = inputParams.getValue("targetannotatedmgf");
    
    //Splitting inputs
    stringSplit(input_spectrum_names, m_spectrum_file_names);
    stringSplit(input_psn_names, m_psm_file_names);
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    
    stringSplit(input_target_spectrum_names, m_target_spectrum_file_names);
    stringSplit(input_target_psn_names, m_target_psm_file_names);
    stringSplit(input_target_annotated_mgf_file, m_target_mgf_file_names);
    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //Background sampling interval
    backgroundsampleinterval = inputParams.getValueInt("backgroundsampleinterval");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
  }


  ExecProjectionStatistics::~ExecProjectionStatistics(void){
  }


  ExecBase * ExecProjectionStatistics::clone(const ParameterList & inputParams) const{
    return new ExecProjectionStatistics(inputParams);
  }

  bool ExecProjectionStatistics::invoke(void){
    cout<<"Invoke"<<endl;
    vector<int> charge_filter;
    charge_filter.push_back(2);
    m_specs.createlibrary(m_envelope_filter, 
                          m_pvalue_filter, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter, 
                          true);
    cout<<"Done Creating Library"<<endl;
                          
    m_target_specs.createlibrary(m_envelope_filter, 
                          m_pvalue_filter, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          true);
                          
    cout<<"Done Creating Targets"<<endl;
    
    vector<vector<float > > average_spectra = m_specs.find_global_average_spectrum(model, ionsToExtract, allIons);
    
    calculate_background(average_spectra);
    
    calculate_substitutions(average_spectra);
    
    return true;
  }

  bool ExecProjectionStatistics::loadInputData(void){      
    //Checking if the number of mzxml files lines up with the number of annotation files  
    if(m_spectrum_file_names.size() != m_psm_file_names.size()) return false;
    
    
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    //ionsToExtract.resize(13);
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-iso");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-iso");
    ionsToExtract.push_back("P++");
    ionsToExtract.push_back("P++-H2O");
    ionsToExtract.push_back("P++-NH3");
    
    load_aminoacid_masses();
    
    load_libs();
    
    return true;
  }


  bool ExecProjectionStatistics::saveOutputData(void){
        return true;
  }


  bool ExecProjectionStatistics::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecProjectionStatistics::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecProjectionStatistics::split(int numSplit){
    return m_subModules;
  }

  bool ExecProjectionStatistics::merge(void){
    return true;
  }


  bool ExecProjectionStatistics::validateParams(std::string & error){
    return true;
  }
  
  void ExecProjectionStatistics::calculate_background(vector<vector<float> > average_spectra){
    unsigned int seed = 0;
    srand(seed);
    vector<string> output_lines;
    //Creating some background statistics
    #pragma omp parallel for num_threads(6)
    for(int library_idx = 0; library_idx < m_specs.size(); library_idx++){
        #pragma omp master
        cout<<"library: "<<library_idx*6<<" of "<<m_specs.size()<<endl;
        
        int random;
        #pragma omp critical
        { random = rand();}
        if( (random % backgroundsampleinterval) != 0) continue; 

        for(int search_idx = 0; search_idx < m_specs.size(); search_idx++){
            if(library_idx == search_idx) continue;
            string library_annotation = m_specs[library_idx].psmList.front()->m_annotation;
            string search_annotation = m_specs[search_idx].psmList.front()->m_annotation;
            if(library_annotation.length() != search_annotation.length()) continue;
            
            
            float sim = spectrum_similarity(  m_specs[library_idx].psmList.front(), 
                                        m_specs[search_idx].psmList.front(), 
                                        library_annotation.length() - 4, 
                                        model, 
                                        ionsToExtract, 
                                        allIons);
                                        
            float sim_minus_background =  spectrum_similarity(  m_specs[library_idx].psmList.front(), 
                                                                m_specs[search_idx].psmList.front(), 
                                                                library_annotation.length() - 4, 
                                                                model, 
                                                                ionsToExtract, 
                                                                allIons,
                                                                average_spectra);
            stringstream ss (stringstream::in | stringstream::out);
            ss<<"BACKGROUNDCOSINE:\t"<<library_annotation<<"\t"<<search_annotation<<"\t"<<sim<<"\t"<<sim_minus_background<<endl;
            
            #pragma omp critical
            {
                output_lines.push_back(ss.str());
            }
        }
    }
    for(int outputline_idx = 0; outputline_idx < output_lines.size(); outputline_idx++){
        cout<<output_lines[outputline_idx];
    }
    
  }
      
  void ExecProjectionStatistics::calculate_substitutions(vector<vector<float> > average_spectra){
    vector<string> output_lines;
    //Creating some background statistics
    #pragma omp parallel for num_threads(6)
    for(int library_idx = 0; library_idx < m_specs.size(); library_idx++){
        #pragma omp master
        cout<<"Substitution Progress: "<<library_idx*6<<" of "<<m_specs.size()<<endl;
        for(int search_idx = 0; search_idx < m_specs.size(); search_idx++){
            if(library_idx == search_idx) continue;
            string library_annotation = m_specs[library_idx].psmList.front()->m_annotation;
            string search_annotation = m_specs[search_idx].psmList.front()->m_annotation;
            if(library_annotation.length() != search_annotation.length()) continue;

            int string_difference = getStringDifference(library_annotation, search_annotation);
            if(string_difference != 1) continue;
            
            //Filtering on where the difference is, 
            //Filtering out the beginning of the spectrum because CID fragmentation causes 
            //low intesnsity at low mz ranges
            int difference_idx = getDifferenceIndex(library_annotation, search_annotation);
            int first_third_idx = (library_annotation.length() - 4)/4 + 2;
            int second_third_idx = (library_annotation.length() - 4)/4*3 + 2;
            if(difference_idx < first_third_idx || difference_idx > second_third_idx) continue;
            
            float sim =  spectrum_similarity(  m_specs[library_idx].psmList.front(), 
                                        m_specs[search_idx].psmList.front(), 
                                        library_annotation.length() - 4, 
                                        model, 
                                        ionsToExtract, 
                                        allIons);
                                        
            float sim_minus_background =  spectrum_similarity(  m_specs[library_idx].psmList.front(), 
                                                                m_specs[search_idx].psmList.front(), 
                                                                library_annotation.length() - 4, 
                                                                model, 
                                                                ionsToExtract, 
                                                                allIons,
                                                                average_spectra);
                                                                
            stringstream ss (stringstream::in | stringstream::out);
            ss<<"SUBSTITUTIONCOSINE:\t"<<library_annotation<<"\t"<<search_annotation<<"\t";
            ss<<difference_idx<<"\t"<<library_annotation[difference_idx]<<"\t"<<search_annotation[difference_idx]<<"\t";
            ss<<sim<<"\t"<<sim_minus_background<<endl;
            
            #pragma omp critical
            {
                output_lines.push_back(ss.str());
            }
        }
    }
    for(int outputline_idx = 0; outputline_idx < output_lines.size(); outputline_idx++){
        cout<<output_lines[outputline_idx];
    }
    
  }
  
  void ExecProjectionStatistics::load_libs(){
    //Looping through mzxml/annotation files to load
    for(int file_idx = 0; file_idx < m_spectrum_file_names.size(); file_idx++){
        //Determining file type by file extension, only supporting plkbin and mzxml
        string spectra_file_name = m_spectrum_file_names[file_idx];
        string annotation_file_name = m_psm_file_names[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            //m_specs.LoadSpecSet_mzxml_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            cout<<"MZXML not Supported"<<endl;
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            m_specs.LoadSpecSet_pklbin_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            cout<<"LOADED\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
        }
        else{
            cerr<<"CANNOTLOAD\t"<<spectra_file_name<<endl;
        }
    }

    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_specs.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    //============Loading Targets===================================
    
    //Looping through mzxml/annotation files to load
    for(int file_idx = 0; file_idx < m_target_spectrum_file_names.size(); file_idx++){
        //Determining file type by file extension, only supporting plkbin and mzxml
        string spectra_file_name = m_target_spectrum_file_names[file_idx];
        string annotation_file_name = m_target_psm_file_names[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            //m_target_specs.LoadSpecSet_mzxml_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            cout<<"MZXML not Supported"<<endl;
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            m_target_specs.LoadSpecSet_pklbin_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            cout<<"LOADED\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
        }
        else{
            cerr<<"CANNOTLOAD\t"<<spectra_file_name<<endl;
        }
    }

    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_target_mgf_file_names.size(); file_idx++){
        string mgf_file_name = m_target_mgf_file_names[file_idx];
        m_target_specs.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
  }
  
  void ExecProjectionStatistics::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
}
