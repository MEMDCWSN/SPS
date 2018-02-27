//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecLibraryView.h"

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

  ExecLibraryView::ExecLibraryView(void){
    m_name = "ExecLibraryView";
    m_type = "ExecLibraryView";
  }


  ExecLibraryView::ExecLibraryView(const ParameterList & inputParams){
    m_name = "ExecLibraryView";
    m_type = "ExecLibraryView";
    inputParams.print(cout);

    input_spectra_data_file = inputParams.getValue("INPUTSPECTRA_DATAFILE", "");

    //checking if CCMS
    if(inputParams.getValueInt("CCMS", 0) == 1){
        //input_spectra_data_file = input_spectra_data_file + "/toolParams-00000.txt";
        input_spectra_data_file = input_spectra_data_file;
    }

    stringSplit(inputParams.getValue("EXISTING_LIBRARY_MGF", ""), existing_library_name);

    results_dir = inputParams.getValue("RESULTS_DIR", ".");

    //Reading in file mapping
    int mapping_count = 0;
    while(1){
      char buf[100];
      sprintf(buf, "upload_file_mapping%i", mapping_count);
      if(!inputParams.exists(buf))
        break;
      std::string mapping = inputParams.getValue(buf);
      std::string mangled_name = mapping.substr(0, mapping.find("|"));
      std::string original_name = mapping.substr(mapping.find("|")+1);

      original_name = get_only_filename(original_name);

      file_name_mapping[original_name] = mangled_name;

      DEBUG_MSG("MAPPING\t"<<original_name<<"\t"<<mangled_name);
      
      if(mangled_name.find("spec-") != string::npos){
          std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
          existing_library_name.push_back(path_to_spectra + "/" + mangled_name);
      }

      mapping_count++;
    }
    
    for (std::map<string,string>::iterator it=file_name_mapping.begin(); it!=file_name_mapping.end(); ++it){
        std::cout << it->first << " => " << it->second << '\n';
    }
    
    
    for(int file_idx = 0; file_idx < existing_library_name.size(); file_idx++){
        string spectra_file_name = existing_library_name[file_idx];
        string extension = get_extension(spectra_file_name);
        
        SpecSet temp_specs;
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
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
            if(temp_specs[spec_idx].fileName.length() == 0)
                temp_specs[spec_idx].fileName = spectra_file_name;
        }
        
        m_libraries.insert(m_libraries.begin(), temp_specs.begin(), temp_specs.end());
    }
    
  }


  ExecLibraryView::~ExecLibraryView(void){
  }


  ExecBase * ExecLibraryView::clone(const ParameterList & inputParams) const{
    return new ExecLibraryView(inputParams);
  }

  bool ExecLibraryView::invoke(void){
    DEBUG_MSG("OUTPUTTING ALL ANNOTATIONS");


    PeptideSpectrumMatchSet output_psms;
    for(int library_idx = 0; library_idx < m_libraries.size(); library_idx++){
        if(m_libraries[library_idx].psmList.size() == 0){
            psmPtr new_psm (new PeptideSpectrumMatch);
            new_psm->m_annotation = "*..*";
            new_psm->m_scanNum = m_libraries[library_idx].scan;
            new_psm->m_dbIndex = library_idx + 1;
            new_psm->m_mz = m_libraries[library_idx].parentMZ;
            new_psm->m_protein = m_libraries[library_idx].instrument_name;
            output_psms.push_back(new_psm);
        }
        else{
            psmPtr temp = m_libraries[library_idx].psmList.front();
            temp->m_scanNum = m_libraries[library_idx].scan;
            temp->m_dbIndex = library_idx + 1;
            if(temp->m_annotation.find("NOTPEPTIDE") != -1){
                temp->m_annotation = "*..*";
            }
            temp->m_protein = m_libraries[library_idx].instrument_name;
            temp->m_mz = m_libraries[library_idx].parentMZ;
            temp->m_notes = temp->m_notes;
            temp->m_library_name = temp->m_spectrumFile;
            DEBUG_VAR(m_libraries[library_idx].fileName);
            if(m_libraries[library_idx].fileName.length() > 0)
                temp->m_spectrumFile = m_libraries[library_idx].fileName;
            temp->m_notes = temp->m_notes;
            output_psms.push_back(temp);
            
        }
    }

    
    DEBUG_MSG("OPENING FILE STREAMS\t"<<results_dir);
    
    output_psms.saveToFile(results_dir.c_str());

    DEBUG_MSG("DONE");            
    return true;
  }

  bool ExecLibraryView::loadInputData(void){
    return true;
  }


  bool ExecLibraryView::saveOutputData(void){
        return true;
  }


  bool ExecLibraryView::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecLibraryView::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecLibraryView::split(int numSplit){
    return m_subModules;
  }

  bool ExecLibraryView::merge(void){
    return true;
  }


  bool ExecLibraryView::validateParams(std::string & error){
    return true;
  }

}
