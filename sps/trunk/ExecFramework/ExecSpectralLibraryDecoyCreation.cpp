//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibraryDecoyCreation.h"
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

  ExecSpectralLibraryDecoyCreation::ExecSpectralLibraryDecoyCreation(void){
    m_name = "ExecSpectralLibraryDecoyCreation";
    m_type = "ExecSpectralLibraryDecoyCreation";
  }


  ExecSpectralLibraryDecoyCreation::ExecSpectralLibraryDecoyCreation(const ParameterList & inputParams){
    m_name = "ExecSpectralLibraryDecoyCreation";
    m_type = "ExecSpectralLibraryDecoyCreation";
    inputParams.print(cout);
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file", inputParams.getValue("MODEL_FILE"));    
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf", inputParams.getValue("EXISTING_LIBRARY_MGF"));
    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion", "");
    
    
    //Specfiying amino acid masses file
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile", inputParams.getValue("AMINOACIDMASS_FILE"));
    
    //Output file for decoy library
    m_output_decoy_mgf_file_name = inputParams.getValue("outputdecoylibrary", inputParams.getValue("OUTPUT_MGF"));
    
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
  }


  ExecSpectralLibraryDecoyCreation::~ExecSpectralLibraryDecoyCreation(void){
  }


  ExecBase * ExecSpectralLibraryDecoyCreation::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibraryDecoyCreation(inputParams);
  }

  bool ExecSpectralLibraryDecoyCreation::invoke(void){
    cout<<"Invoke"<<endl;
    vector<int> charge_filter;
    //charge_filter.push_back(1);
    //charge_filter.push_back(2);
    //charge_filter.push_back(3);
    //charge_filter.push_back(4);
    //charge_filter.push_back(5);
    //charge_filter.push_back(6);
    
    m_library.createlibrary(-1, 
                          -1, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          false);
    
    cout<<"Done Creating Library"<<endl;
                          

    cout<<"Generating Decoy"<<endl;
    m_decoy_library = m_library.create_decoy_spectral_library(model, ionsToExtract, allIons);
   
    cout<<"Done Creating Decoy Library"<<endl;
    
    if(m_output_decoy_mgf_file_name.length() > 0)
        m_decoy_library.SaveSpecSet_mgf(m_output_decoy_mgf_file_name.c_str());
    else
        cout<<"Invalid output name"<<endl;
    
    return true;
  }

  bool ExecSpectralLibraryDecoyCreation::loadInputData(void){      
    
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
        cout<<"Loading MGF\n"<<endl;
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    return true;
  }


  bool ExecSpectralLibraryDecoyCreation::saveOutputData(void){
        return true;
  }


  bool ExecSpectralLibraryDecoyCreation::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibraryDecoyCreation::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibraryDecoyCreation::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibraryDecoyCreation::merge(void){
    return true;
  }


  bool ExecSpectralLibraryDecoyCreation::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibraryDecoyCreation::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
