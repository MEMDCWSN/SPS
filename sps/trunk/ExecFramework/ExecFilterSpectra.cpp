//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecFilterSpectra.h"

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

  ExecFilterSpectra::ExecFilterSpectra(void){
    m_name = "ExecFilterSpectra";
    m_type = "ExecFilterSpectra";
  }


  ExecFilterSpectra::ExecFilterSpectra(const ParameterList & inputParams){
    m_name = "ExecFilterSpectra";
    m_type = "ExecFilterSpectra";
    
    this->m_params = inputParams;
    
    output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");
    library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");
  }


  ExecFilterSpectra::~ExecFilterSpectra(void){
  }


  ExecBase * ExecFilterSpectra::clone(const ParameterList & inputParams) const{
    return new ExecFilterSpectra(inputParams);
  }

  bool ExecFilterSpectra::invoke(void){
        int exceeds_peak_max_count = 0;
        int max_peak_count = 500;
        
        for(int i = 0; i < m_library.size(); i++){
            cout<<"PEAK COUNT: "<<m_library[i].size()<<"\t"<<m_library[i].psmList.front()->m_spectrumID<<endl;
            if(m_library[i].size() > max_peak_count){
                exceeds_peak_max_count++;
            }
            
            m_library[i].selectTopK(max_peak_count);
            //cout<<i<<"\t"<<m_library[i].size()<<endl;
            //m_library[i].rankFilterPeaks(10,10);
            cout<<"done with filter"<<endl;
        }
        cout<<"Exceeding max peak: " <<exceeds_peak_max_count<<endl;


        return true;
  }

  bool ExecFilterSpectra::loadInputData(void){
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


  bool ExecFilterSpectra::saveOutputData(void){
    //Saving out
    cout<<"Saving file"<<endl;
    m_library.SaveSpecSet(output_mgf_name.c_str());
    
    return true;
  }


  bool ExecFilterSpectra::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecFilterSpectra::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecFilterSpectra::split(int numSplit){
    return m_subModules;
  }

  bool ExecFilterSpectra::merge(void){
    return true;
  }


  bool ExecFilterSpectra::validateParams(std::string & error){
    return true;
  }

}
