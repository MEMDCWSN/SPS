//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectralLibraryBatchRemove.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
//#include "mzxml.h"
#include <fstream>
#include "projectionutils.h"

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectralLibraryBatchRemove::ExecSpectralLibraryBatchRemove(void){
    m_name = "ExecSpectralLibraryBatchRemove";
    m_type = "ExecSpectralLibraryBatchRemove";
  }


  ExecSpectralLibraryBatchRemove::ExecSpectralLibraryBatchRemove(const ParameterList & inputParams){
    m_name = "ExecSpectralLibraryBatchRemove";
    m_type = "ExecSpectralLibraryBatchRemove";
    
    this->m_params = inputParams;
    
    output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");
    library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");
    
    //String Splitting on spaces. 
    stringSplit(inputParams.getValue("SPECTRUMID_TO_DELETE", ""), spectrum_IDs_to_remove);
    
  }


  ExecSpectralLibraryBatchRemove::~ExecSpectralLibraryBatchRemove(void){
  }


  ExecBase * ExecSpectralLibraryBatchRemove::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibraryBatchRemove(inputParams);
  }

  bool ExecSpectralLibraryBatchRemove::invoke(void){
        for(int i = 0; i < m_library.size(); i++){
            //Detecting if spectrum id matches ones to delete
            bool found_flag = false;
            
            for(int delete_id_idx = 0; delete_id_idx < spectrum_IDs_to_remove.size(); delete_id_idx++){
                if(m_library[i].psmList.front()->m_spectrumID == spectrum_IDs_to_remove[delete_id_idx]){
                    cout<<"FOUND: " << spectrum_IDs_to_remove[delete_id_idx]<<endl;
                    cout<<m_library[i].size()<<endl;
                    found_flag = true;
                    break;
                }
            }
            
            if(found_flag == false){
                updated_library.push_back(m_library[i]);
            }
        }
        


        return true;
  }

  bool ExecSpectralLibraryBatchRemove::loadInputData(void){
        string spectra_file_name = library_path;
        string extension = get_extension(spectra_file_name);

        SpecSet existing_lib;

        if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
            DEBUG_MSG("LOADING Lib\t"<<spectra_file_name<<"\t as a mgf file");
            m_library.LoadSpecSet_mgf(spectra_file_name.c_str());
            DEBUG_MSG("DONE");
            
        }
        else{
            DEBUG_MSG("CANNOTLOAD Search\t"<<spectra_file_name);
        }
      
        return true;
  }


  bool ExecSpectralLibraryBatchRemove::saveOutputData(void){
    //Saving out
    cout<<"Saving file"<<endl;
    updated_library.SaveSpecSet(output_mgf_name.c_str());
    
    return true;
  }


  bool ExecSpectralLibraryBatchRemove::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibraryBatchRemove::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibraryBatchRemove::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibraryBatchRemove::merge(void){
    return true;
  }


  bool ExecSpectralLibraryBatchRemove::validateParams(std::string & error){
    return true;
  }

}
