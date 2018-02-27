//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecCreateSpectralLibraryMolecular.h"

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

  ExecCreateSpectralLibraryMolecular::ExecCreateSpectralLibraryMolecular(void){
    m_name = "ExecCreateSpectralLibraryMolecular";
    m_type = "ExecCreateSpectralLibraryMolecular";
  }


  ExecCreateSpectralLibraryMolecular::ExecCreateSpectralLibraryMolecular(const ParameterList & inputParams){
    m_name = "ExecCreateSpectralLibraryMolecular";
    m_type = "ExecCreateSpectralLibraryMolecular";

    inputParams.print(cout);

    this->m_params = inputParams;

    output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");

    stringSplit(inputParams.getValue("INPUT_PSMS_FILES", "."), input_psms_filenames);


    string library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");

    min_mass_diff = inputParams.getValueFloat("MIN_SHIFT_MASS", -1.f);

    vector<string> possibleDirFiles = directoryContents(library_path,
                                                          "",
                                                          "",
                                                          true);
    for(int i = 0; i < possibleDirFiles.size(); i++){
        DEBUG_VAR(possibleDirFiles[i]);
        existing_library_name.push_back(library_path + "/" + possibleDirFiles[i]);
    }


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
      file_reverse_mapping[mangled_name] = original_name;

      DEBUG_MSG("MAPPING\t"<<original_name<<"\t"<<mangled_name);

      if(mangled_name.find("spec-") != string::npos){
          std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
          input_spectra.push_back(path_to_spectra + "/" + mangled_name);
      }

      mapping_count++;
    }

    for (std::map<string,string>::iterator it=file_name_mapping.begin(); it!=file_name_mapping.end(); ++it){
        std::cout << it->first << " => " << it->second << '\n';
    }

    //Renaming output the same as before
    //if(existing_library_name.size() == 1)
    //    output_mgf_name = output_mgf_name + "/" + get_only_filename(existing_library_name[0]);

  }


  ExecCreateSpectralLibraryMolecular::~ExecCreateSpectralLibraryMolecular(void){
  }


  ExecBase * ExecCreateSpectralLibraryMolecular::clone(const ParameterList & inputParams) const{
    return new ExecCreateSpectralLibraryMolecular(inputParams);
  }

  bool ExecCreateSpectralLibraryMolecular::invoke(void){






    return true;
  }

  bool ExecCreateSpectralLibraryMolecular::loadInputData(void){

    for(int file_idx = 0; file_idx < input_spectra.size(); file_idx++){
        string spectra_file_name = input_spectra[file_idx];
        string extension = get_extension(spectra_file_name);

        SpecSet temp_specs;
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0 || strcmp(extension.c_str(), "mzML") == 0){
            string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
            int ret_val = temp_specs.Load(pklbin_version.c_str());
            DEBUG_MSG(ret_val);
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

        source_spectra.insert(source_spectra.begin(), temp_specs.begin(), temp_specs.end());
    }


    //Loading PSMS data
    for(int psm_idx = 0; psm_idx < input_psms_filenames.size(); psm_idx++){
        string psm_filename = input_psms_filenames[psm_idx];
        PeptideSpectrumMatchSet temp_psms;

        temp_psms.loadFromFile(psm_filename.c_str());

        for(int i = 0 ; i < temp_psms.size(); i++){
            all_psms.push_back(temp_psms[i]);
        }
    }

    DEBUG_MSG("Loaded PSMs, Size: " << all_psms.size());

    for(int psm_idx = 0; psm_idx < all_psms.size(); psm_idx++){
        DEBUG_VAR(all_psms[psm_idx]->m_spectrumFile);
    }

    DEBUG_VAR(source_spectra.size());
    for(int i = 0; i < source_spectra.size(); i++){
        //Unmangling the name
        source_spectra[i].fileName = file_reverse_mapping[get_only_filename(source_spectra[i].fileName)];
    }

    all_psms.addSpectra(&source_spectra);

    for(int i = 0; i < source_spectra.size(); i++){
        if(source_spectra[i].psmList.size() > 0){
            if(abs(source_spectra[i].psmList.front()->m_parentmass_difference) < min_mass_diff)
                continue;

            //DEBUG_MSG(source_spectra[i].psmList->front()->m_dbIndex<<" "<<source_spectra[i].psmList->front()->m_score<<" "<<source_spectra[i].psmList->front()->m_parentmass_difference);
            stringstream ss (stringstream::in | stringstream::out);
            ss<<file_reverse_mapping[source_spectra[i].psmList.front()->m_library_name]<<":"<<source_spectra[i].psmList.front()->m_dbIndex;
            ss<<":Score "<<source_spectra[i].psmList.front()->m_score;
            ss<<":MassDiff "<<source_spectra[i].psmList.front()->m_parentmass_difference;
            ss<<":"<<source_spectra[i].psmList.front()->m_notes;


            source_spectra[i].psmList.front()->m_notes = ss.str();
            m_new_library.push_back(source_spectra[i]);
        }
    }

    return true;
  }


  bool ExecCreateSpectralLibraryMolecular::saveOutputData(void){

    if(output_mgf_name.length() > 0){
        m_new_library.SaveSpecSet_mgf(output_mgf_name.c_str());
    }
    return true;
  }


  bool ExecCreateSpectralLibraryMolecular::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecCreateSpectralLibraryMolecular::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecCreateSpectralLibraryMolecular::split(int numSplit){
    return m_subModules;
  }

  bool ExecCreateSpectralLibraryMolecular::merge(void){
    return true;
  }


  bool ExecCreateSpectralLibraryMolecular::validateParams(std::string & error){
    return true;
  }

}
