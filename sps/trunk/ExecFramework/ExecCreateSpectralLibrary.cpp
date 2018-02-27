//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecCreateSpectralLibrary.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include "utils.h"

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecCreateSpectralLibrary::ExecCreateSpectralLibrary(void){
    m_name = "ExecCreateSpectralLibrary";
    m_type = "ExecCreateSpectralLibrary";
  }


  ExecCreateSpectralLibrary::ExecCreateSpectralLibrary(const ParameterList & inputParams){
    m_name = "ExecCreateSpectralLibrary";
    m_type = "ExecCreateSpectralLibrary";
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
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf");
    
    //Splitting inputs
    stringSplit(input_spectrum_names, m_spectrum_file_names);
    stringSplit(input_psn_names, m_psm_file_names);
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
  }


  ExecCreateSpectralLibrary::~ExecCreateSpectralLibrary(void){
  }


  ExecBase * ExecCreateSpectralLibrary::clone(const ParameterList & inputParams) const{
    return new ExecCreateSpectralLibrary(inputParams);
  }

  bool ExecCreateSpectralLibrary::invoke(void){
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
    
    
    //Setting the parent mass
    for(int spec_idx = 0; spec_idx < m_specs.size(); spec_idx++){
        if(m_specs[spec_idx].parentMass == 0.f){
            if(m_specs[spec_idx].parentCharge != 0){
                 m_specs[spec_idx].parentMass = (m_specs[spec_idx].parentMZ * m_specs[spec_idx].parentCharge) - (AAJumps::massHion * (m_specs[spec_idx].parentCharge - 1));
            }
        }
    }
                          
    cout<<"Done Creating"<<endl;
    
    return true;
  }

  bool ExecCreateSpectralLibrary::loadInputData(void){      
    //Checking if the number of mzxml files lines up with the number of annotation files  
    if(m_spectrum_file_names.size() != m_psm_file_names.size()) return false;
    
    
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    ionsToExtract.resize(13);
    ionsToExtract[0] = "a";
    ionsToExtract[1] = "b";
    ionsToExtract[2] = "b-NH3";
    ionsToExtract[3] = "b-H2O";
    ionsToExtract[4] = "b++";
    ionsToExtract[5] = "b++-iso";
    ionsToExtract[6] = "y";
    ionsToExtract[7] = "y-NH3";
    ionsToExtract[8] = "y-H2O";
    ionsToExtract[9] = "y++";
    ionsToExtract[10] = "y++-iso";
    ionsToExtract[11] = "P++-H2O";
    ionsToExtract[12] = "P++-NH3";
    
    //Looping through mzxml/annotation files to load
    for(int file_idx = 0; file_idx < m_spectrum_file_names.size(); file_idx++){
        //Determining file type by file extension, only supporting plkbin and mzxml
        string spectra_file_name = m_spectrum_file_names[file_idx];
        string annotation_file_name = m_psm_file_names[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            //m_specs.LoadSpecSet_mzxml_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            cout<<"MZXML not Supported"<<endl;
            
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            cout<<"LOADED\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
            m_specs.LoadSpecSet_pklbin_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            
        }
        else if(strcmp(extension.c_str(), "mgf") == 0){
            cout<<"LOADED\t"<<spectra_file_name<<"\t as a mgf file"<<endl;
            m_specs.LoadSpecSet_mgf_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
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
    
    
    
    return true;
  }


  bool ExecCreateSpectralLibrary::saveOutputData(void){
    m_specs.SaveSpecSet_mgf(m_output_mgf_file_name.c_str());
    return true;
  }


  bool ExecCreateSpectralLibrary::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecCreateSpectralLibrary::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecCreateSpectralLibrary::split(int numSplit){
    return m_subModules;
  }

  bool ExecCreateSpectralLibrary::merge(void){
    return true;
  }


  bool ExecCreateSpectralLibrary::validateParams(std::string & error){
    return true;
  }
  
  void ExecCreateSpectralLibrary::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
