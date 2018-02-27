//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectralLibraryPeakStatistics.h"

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

  ExecSpectralLibraryPeakStatistics::ExecSpectralLibraryPeakStatistics(void){
    m_name = "ExecSpectralLibraryPeakStatistics";
    m_type = "ExecSpectralLibraryPeakStatistics";
  }


  ExecSpectralLibraryPeakStatistics::ExecSpectralLibraryPeakStatistics(const ParameterList & inputParams){
    m_name = "ExecSpectralLibraryPeakStatistics";
    m_type = "ExecSpectralLibraryPeakStatistics";
    
    this->m_params = inputParams;
    string input_annotated_mgf_file = inputParams.getValue("EXISTING_LIBRARY_MGF");
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
  }


  ExecSpectralLibraryPeakStatistics::~ExecSpectralLibraryPeakStatistics(void){
  }


  ExecBase * ExecSpectralLibraryPeakStatistics::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibraryPeakStatistics(inputParams);
  }

  bool ExecSpectralLibraryPeakStatistics::invoke(void){
        map<float, vector<int> > SNR_to_peaks_above;
        map<string, vector<int> > SpectrumID_to_peaks_above;
        int list_length = 0;
        for(float SNR_threshold = 1.0; SNR_threshold < 20.f; SNR_threshold += 1.f){
            vector<int> temp_peaks_above_list;
            list_length = 0;
            
            int too_few_peaks_spectrum = 0;
            for(int library_idx = 0; library_idx < m_libraries.size(); library_idx++){
                for(int i = 0; i < m_libraries[library_idx].size(); i++){
                    float spectrum_noise_level = m_libraries[library_idx][i].getNoiseLevel();
                    if(spectrum_noise_level == 0.f){
                        too_few_peaks_spectrum++;
                        continue;
                    }
                    
                    //DEBUG_VAR(m_libraries[library_idx][i].getNoiseLevel());
                    
                    Spectrum temp = m_libraries[library_idx][i];
                    int before_filter = temp.size();
                    temp.filterLowIntensity(SNR_threshold * spectrum_noise_level);
                    //DEBUG_MSG(before_filter<<"\t"<<temp.size());
                    temp_peaks_above_list.push_back(temp.size());
                    
                    string spectrumid = m_libraries[library_idx][i].psmList.front()->m_spectrumID;
                    if(SpectrumID_to_peaks_above.find(spectrumid) == SpectrumID_to_peaks_above.end()){
                        vector<int> temp_vec;
                        SpectrumID_to_peaks_above[spectrumid] = temp_vec;
                    }
                    SpectrumID_to_peaks_above[spectrumid].push_back(temp.size());
                    
                    list_length++;
                }
            }
            
            SNR_to_peaks_above[SNR_threshold] = temp_peaks_above_list;
            
            DEBUG_MSG("Spectra with under 10 peaks: " << too_few_peaks_spectrum);
        }
        
        //Outputting on a spectrum ID Basis
        cout<<"\t";
        for (std::map<float,vector<int> >::iterator it=SNR_to_peaks_above.begin(); it!=SNR_to_peaks_above.end(); ++it){
            cout << it->first << "\t";
        }
        cout<<endl;
        
        for (std::map<string,vector<int> >::iterator it=SpectrumID_to_peaks_above.begin(); it!=SpectrumID_to_peaks_above.end(); ++it){
            cout<<it->first<<"\t";
            for(int i = 0; i < it->second.size(); i++){
                cout<<it->second[i]<<"\t";
            }
            cout<<endl;
        }
        
        
        /*
        for (std::map<float,vector<int> >::iterator it=SNR_to_peaks_above.begin(); it!=SNR_to_peaks_above.end(); ++it){
            cout << it->first << "\t";
        }
        cout<<endl;
        
        for(int i = 0; i < list_length; i++){
            for (std::map<float,vector<int> >::iterator it=SNR_to_peaks_above.begin(); it!=SNR_to_peaks_above.end(); ++it){
                cout<<it->second[i]<<"\t";
            }
            cout<<endl;
        }*/
        
        return true;
  }

  bool ExecSpectralLibraryPeakStatistics::loadInputData(void){
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
        return true;
  }


  bool ExecSpectralLibraryPeakStatistics::saveOutputData(void){
    return true;
  }


  bool ExecSpectralLibraryPeakStatistics::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibraryPeakStatistics::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibraryPeakStatistics::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibraryPeakStatistics::merge(void){
    return true;
  }


  bool ExecSpectralLibraryPeakStatistics::validateParams(std::string & error){
    return true;
  }

}
