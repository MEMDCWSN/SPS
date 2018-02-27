//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraExtractionSingle.h"

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

  ExecSpectraExtractionSingle::ExecSpectraExtractionSingle(void){
    m_name = "ExecSpectraExtractionSingle";
    m_type = "ExecSpectraExtractionSingle";
  }


  ExecSpectraExtractionSingle::ExecSpectraExtractionSingle(const ParameterList & inputParams){
    m_name = "ExecSpectraExtractionSingle";
    m_type = "ExecSpectraExtractionSingle";
    
    
    this->m_params = inputParams;
    
    //output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");
    
    new_lib_results_dir = inputParams.getValue("NEWLIBRARYRESULTS_DIR", ".");
    
    //stringSplit(inputParams.getValue("EXISTING_LIBRARY_MGF", ""), existing_library_name);

    results_dir = inputParams.getValue("RESULTS_DIR", ".");
    
    //string library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");
    
    library_grid_root_dir = inputParams.getValue("LIBUSER_ROOT_DATA_PATH", "");
    library_user_root_dir = inputParams.getValue("LIBUSER_ROOT_VM_DISPLAY_PATH", "");
    library_user = inputParams.getValue("LIBUSER_NAME", "");    //Should check if library path selected matches this
    
    submit_user = inputParams.getValue("user","INVALIDUSER");
    
    
    //vector<string> possibleDirFiles = directoryContents(library_path,
    //                                                      "",
    //                                                      "",
    //                                                      true);
    //for(int i = 0; i < possibleDirFiles.size(); i++){
    //    DEBUG_VAR(possibleDirFiles[i]);
        //existing_library_name.push_back(library_path + "/" + possibleDirFiles[i]);
    //}
    

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

      string original_name_filename = get_only_filename(original_name);

      
      full_file_name_mapping[original_name] = mangled_name;
      file_name_mapping[original_name_filename] = mangled_name;
      file_reverse_mapping[mangled_name] = original_name_filename;

      DEBUG_MSG("MAPPING\t"<<original_name_filename<<"\t"<<mangled_name);
      
      if(mangled_name.find("spec-") != string::npos){
          std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
          input_spectra.push_back(path_to_spectra + "/" + mangled_name);
      }

      mapping_count++;
    }
    
    for (std::map<string,string>::iterator it=full_file_name_mapping.begin(); it!=full_file_name_mapping.end(); ++it){
        if(it->second[0] == 'l' && it->second[1] == 'i' && it->second[2] == 'b' && it->second[3] == '-'){
            //We have a library, find original name
            string full_lib_path = library_grid_root_dir  + it->first;
            existing_library_name.push_back(full_lib_path);
            DEBUG_VAR(full_lib_path);
            DEBUG_VAR(library_user);
            DEBUG_VAR(it->first);
            DEBUG_VAR(it->first.find(library_user));
            
            
            if(full_lib_path.find("*") != string::npos || full_lib_path.find("..") != string::npos || full_lib_path.find("~") != string::npos){
                DEBUG_MSG("Invalid letters in path, exiting");
                exit(1);
            }
            
            //only can write to the owners file
            if(it->first.find(submit_user) != 0 && it->first.find(library_user) != 0){
                DEBUG_MSG("Invalid User");
                exit(1);
            }
            else{
                DEBUG_MSG("Writing to User Space: " + submit_user + " or " + library_user);
            }
        }
    }

    
    //Renaming output the same as before
    //if(existing_library_name.size() == 1)
    //    output_mgf_name = output_mgf_name + "/" + get_only_filename(existing_library_name[0]);
    
  }


  ExecSpectraExtractionSingle::~ExecSpectraExtractionSingle(void){
  }


  ExecBase * ExecSpectraExtractionSingle::clone(const ParameterList & inputParams) const{
    return new ExecSpectraExtractionSingle(inputParams);
  }

  bool ExecSpectraExtractionSingle::invoke(void){
    
    std::string spectrum_sequence;
    std::string compound_name;
    float parentMZ;
    float exact_mass;
    std::string genus;
    std::string species;
    std::string strain;
    std::string instrument;
    std::string ionsource;
    std::string extract_scan_str;
    std::string smiles;
    std::string InChI;
    std::string InChI_Aux;
    std::string notes;
    int precursor_charge;
    std::string ion_mode;
    std::string PI;
    std::string data_collector;
    std::string CASNumber;
    std::string PublicationID;
    std::string acquisition_method;
    
    std::string submission_user;
    std::string submission_ID;
    std::string submission_date;
    
    std::string ncbi_number;
    
    std::string adduct;
    std::string user_species;
    std::string user_genus;
    
    std::string spectrum_novelty;
    int library_quality;
    
    std::string output_new_spectrum_file;
    
    submission_user = m_params.getValue("user");
    submission_ID = m_params.getValue("task");
    submission_date = m_params.getValue("submission_date", get_time_string().c_str());
    
    spectrumID = m_params.getValue("SPECTRUMID", "NO_ID");
    spectrum_sequence = m_params.getValue("ADDSPECTRA_SEQ", "*..*");
    compound_name = m_params.getValue("ADDSPECTRA_COMPOUND_NAME");
    parentMZ = m_params.getValueFloat("ADDSPECTRA_MOLECULEMASS", 0.f);
    genus = m_params.getValue("ADDSPECTRA_GENUS");
    species = m_params.getValue("ADDSPECTRA_SPECIES");
    strain = m_params.getValue("ADDSPECTRA_STRAIN");
    instrument = m_params.getValue("ADDSPECTRA_INSTRUMENT");
    ionsource = m_params.getValue("ADDSPECTRA_IONSOURCE");
    extract_scan_str = m_params.getValue("ADDSPECTRA_EXTRACTSCAN");
    smiles = m_params.getValue("ADDSPECTRA_SMILES");
    InChI = m_params.getValue("ADDSPECTRA_INCHI");
    InChI_Aux = m_params.getValue("ADDSPECTRA_INCHIAUX");
    precursor_charge = m_params.getValueInt("ADDSPECTRA_CHARGE", 0);
    ion_mode = m_params.getValue("ADDSPECTRA_IONMODE", "Positive");
    PI = m_params.getValue("ADDSPECTRA_PI", "UNKNOWN");
    CASNumber = m_params.getValue("ADDSPECTRA_CASNUMBER");
    PublicationID = m_params.getValue("ADDSPECTRA_PUB");
    acquisition_method = m_params.getValue("ADDSPECTRA_ACQUISITION");
    ncbi_number = m_params.getValue("ADDSPECTRA_NCBI_STRAIN", "");
    exact_mass = m_params.getValueFloat("ADDSPECTRA_EXACTMASS", 0.f);
    data_collector = m_params.getValue("ADDSPECTRA_DATACOLLECTOR", "");
    adduct = m_params.getValue("ADDSPECTRA_ADDUCT", "");
    user_species = m_params.getValue("ADDSPECTRA_USERSPECIES", "");
    user_genus = m_params.getValue("ADDSPECTRA_USERGENUS", "");
    spectrum_novelty = m_params.getValue("ADDSPECTRA_INTEREST", "");
    library_quality = m_params.getValueInt("ADDSPECTRA_LIBQUALITY", 0);
    
    output_new_spectrum_file = m_params.getValue("ADDSPECTRA_OUTPUTNEWSPECTRAFILE", "");
    
    //Extracting Spectra
    int extraction_scan = -1;
    if(extract_scan_str != "*"){
        extraction_scan = atoi(extract_scan_str.c_str());
    }
    
    if(extraction_scan != -1){
        output_mgf_name = output_mgf_name + "/" + file_reverse_mapping[get_only_filename(input_spectra[0])] + "_" + extract_scan_str + ".mgf";
        if(input_spectra.size() > 1){
            DEBUG_MSG("Too Many Input Files");
            exit(1);
        }
        DEBUG_VAR(output_mgf_name);
    }
    else{
        output_mgf_name = output_mgf_name + "/" + file_reverse_mapping[get_only_filename(input_spectra[0])] + "_allunannotated_" + genus + "_" + species + "_" + strain + "_" + instrument + ".mgf";
        DEBUG_VAR(output_mgf_name);
    }
    
    
    DEBUG_MSG("OUTPUTTING OLD ANNOTATIONS");
    /*for(int library_idx = 0; library_idx < m_library.size(); library_idx++){
        if(m_library[library_idx].psmList.size() == 0){
            psmPtr new_psm (new PeptideSpectrumMatch);
            new_psm->m_annotation = "*..*";
            new_psm->m_dbIndex = library_idx + 1;
            new_psm->m_scanNum = m_library[library_idx].scan;
            new_psm->m_mz = m_library[library_idx].parentMZ;
            new_psm->m_protein = m_library[library_idx].instrument_name;
            new_psm->m_library_name = get_only_filename(output_mgf_name);
            output_psms.push_back(new_psm);
        }
        else{
            psmPtr temp = m_library[library_idx].psmList.front();
            temp->m_dbIndex = library_idx + 1;
            temp->m_scanNum = m_library[library_idx].scan;
            if(temp->m_annotation.find("NOTPEPTIDE") != -1){
                temp->m_annotation = "*..*";
            }
            temp->m_protein = m_library[library_idx].instrument_name;
            temp->m_mz = m_library[library_idx].parentMZ;
            temp->m_notes = temp->m_notes;
            
            //DEBUG_VAR(temp->m_exactmass);
            //DEBUG_VAR(temp->m_spectrumFile);
            temp->m_library_name = get_only_filename(temp->m_spectrumFile);
            //temp->m_library_name = get_only_filename(output_mgf_name);
            
            
            if(m_library[library_idx].fileName.length() > 0)
                temp->m_spectrumFile = m_library[library_idx].fileName;
            
            
            output_psms.push_back(temp);
        }
    }*/
    
    int found_source_spectrum = -1;
    for(int i = 0; i < source_spectra.size(); i++){
        if(extraction_scan != -1 && extraction_scan != source_spectra[i].scan){
            continue;
        }
        
        found_source_spectrum = 1;
        
        Spectrum extracted_spec;
        extracted_spec = source_spectra[i];
        extracted_spec.msLevel = 2;
        extracted_spec.instrument_name = ionsource + "-" + instrument;
        extracted_spec.fileName = file_reverse_mapping[get_only_filename(source_spectra[i].fileName)];
        
        
        //Correcting for precursor MZ
        if(parentMZ > 1.f){
            extracted_spec.parentMZ = parentMZ;
        }
        
        psmPtr psm(new PeptideSpectrumMatch());
        psm->m_spectrumFile = extracted_spec.fileName;
        psm->m_scanNum = source_spectra[i].scan;
        DEBUG_VAR(source_spectra[i].scan);     
        psm->m_annotation = spectrum_sequence;
        string submission_metadata = submission_user  + ":"  + submission_ID + ":" + submission_date;
        psm->m_submission_metadata = (submission_metadata);
        
        //Determining whether to use user organism
        if(genus == "Other" && user_genus.length() > 0){
            genus = user_genus;
        }
        
        if(species == "Other" && user_species.length() > 0){
            species = user_species;
        }
        
        string full_organism_name = genus + ":" + species;
        if(strain.length() > 0){
            full_organism_name += ":" + strain;
        }
        
        if(ncbi_number.length() > 0){
            full_organism_name = ncbi_number + "|" + full_organism_name;
        }
        
        if(adduct.length() > 0){
            compound_name = compound_name + " [" + adduct + "]";
        }
        
        psm->m_organism = (full_organism_name);
        psm->m_compound_name = (compound_name);
        psm->m_smiles = (smiles);
        psm->m_InChI = (InChI);
        psm->m_InChI_Aux = (InChI_Aux);
        psm->m_notes = PI + ":" + data_collector + ":" + CASNumber + ":" + PublicationID + ":" + acquisition_method;
        
        if ( spectrum_novelty.length() > 0 ) {
            //If novel then scan should be set to 1
            psm->m_notes += ":" + spectrum_novelty;
            //psm->m_scanNum = 1;
            //source_spectra[i].scan = 1;
        }
        
        if(precursor_charge == 0){
            psm->m_charge = extracted_spec.parentCharge;
        }
        else{
            psm->m_charge = precursor_charge;
        }
        psm->m_ionmode = ion_mode;
        psm->m_library_name = get_only_filename(output_mgf_name);
        psm->m_exactmass = exact_mass;
        psm->m_library_quality = library_quality;
        psm->m_spectrumID = spectrumID;
        psm->m_mz = extracted_spec.parentMZ;
        DEBUG_VAR(psm->m_exactmass);
        
        extracted_spec.psmList.push_back(psm);
        int added_index = m_new_library.add_update_spectrum_to_Library(extracted_spec);
        DEBUG_MSG("FOUND SPEC ADDED AT INDEX "<<added_index);
        if(added_index == -2)
            continue;
        
        
    }
    
    if(found_source_spectrum == -1){
        ERROR_MSG("Did not find spectrum in file");
        ERROR_MSG("Please make sure your file is 32-bit centroided uncompressed mzXML and that you are able to plot it on the input form and the scan number is correct.");
        exit(1);
    }
    else{
        DEBUG_MSG("FOUND SPECTRUM IN FILE");
    }
    
    for(int file_idx = 0; file_idx < existing_library_name.size(); file_idx++){
        string spectra_file_name = existing_library_name[file_idx];
        string extension = get_extension(spectra_file_name);
        
        SpecSet existing_lib;
        
        if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
            DEBUG_MSG("LOADING Lib\t"<<spectra_file_name<<"\t as a mgf file");
            existing_lib.LoadSpecSet_mgf(spectra_file_name.c_str());
            DEBUG_MSG("DONE");
            
        }
        else{
            DEBUG_MSG("CANNOTLOAD Search\t"<<spectra_file_name);
        }
        
        //Reading Library Contents
        for(int library_idx = 0; library_idx < existing_lib.size(); library_idx++){
            psmPtr temp = existing_lib[library_idx].psmList.front();
            temp->m_dbIndex = library_idx + 1;
            temp->m_scanNum = existing_lib[library_idx].scan;
            if(temp->m_annotation.find("NOTPEPTIDE") != -1){
                temp->m_annotation = "*..*";
            }
            temp->m_protein = existing_lib[library_idx].instrument_name;
            temp->m_mz = existing_lib[library_idx].parentMZ;
            temp->m_notes = temp->m_notes;
            
            //DEBUG_VAR(temp->m_exactmass);
            //DEBUG_VAR(temp->m_spectrumFile);
            temp->m_library_name = string_replace(spectra_file_name, library_grid_root_dir, library_user_root_dir);
            
            if(existing_lib[library_idx].fileName.length() > 0){
                temp->m_spectrumFile = existing_lib[library_idx].fileName;
            }
            
            switch(temp->m_library_quality){
                case 0:
                    temp->m_origAnnotation = "NO QUALITY";
                    break;
                case 1:
                    temp->m_origAnnotation = "Gold";
                    break;
                case 2:
                    temp->m_origAnnotation = "Silver";
                    break;
                case 3:
                    temp->m_origAnnotation = "Bronze";
                    break;
            }
            
            output_psms.push_back(temp);
        }
        
        
        string organism_name_from_filename = get_only_filename(spectra_file_name);
        organism_name_from_filename = strip_extension(organism_name_from_filename);
        DEBUG_VAR(organism_name_from_filename);
        //string output_library_name = spectra_file_name + "_temp";
        string output_library_name = spectra_file_name;
        
        //Setting organism name 
        for(int i = 0; i < m_new_library.size(); i++){
            m_new_library[i].psmList.front()->m_organism = organism_name_from_filename;
            m_new_library[i].psmList.front()->m_dbIndex = existing_lib.size() + 1;
            
            psmPtr psm = m_new_library[i].psmList.front();
            psm->m_protein = m_new_library[i].instrument_name;
            psm->m_library_name = string_replace(output_library_name, library_grid_root_dir, library_user_root_dir);
            DEBUG_VAR(psm->m_library_name);
            
            switch(psm->m_library_quality){
                case 0:
                    psm->m_origAnnotation = "NO QUALITY";
                    break;
                case 1:
                    psm->m_origAnnotation = "Gold";
                    break;
                case 2:
                    psm->m_origAnnotation = "Silver";
                    break;
                case 3:
                    psm->m_origAnnotation = "Bronze";
                    break;
            }
            
            output_psms_new.push_back(psm);
        }

        DEBUG_MSG(existing_lib.size());
        existing_lib.insert(existing_lib.end(), m_new_library.begin(), m_new_library.end());
        DEBUG_MSG(existing_lib.size());
        
        
        if(output_new_spectrum_file.length() > 0){
            m_new_library.SaveSpecSet_mgf( output_new_spectrum_file.c_str());
        }
        //temp_specs.SaveSpecSet_mgf( (spectra_file_name + "_temp").c_str());
        existing_lib.SaveSpecSet_mgf( output_library_name.c_str());
    }
    
    
    
    return true;
  }

  bool ExecSpectraExtractionSingle::loadInputData(void){
//     for(int file_idx = 0; file_idx < existing_library_name.size(); file_idx++){
//         string spectra_file_name = existing_library_name[file_idx];
//         string extension = get_extension(spectra_file_name);
//         
//         SpecSet temp_specs;
//         if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
//             string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
//             temp_specs.Load(pklbin_version.c_str());
//             //spectra_file_name = pklbin_version;
//             DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
//         }
//         else if(strcmp(extension.c_str(), "pklbin") == 0){
//             temp_specs.loadPklBin(spectra_file_name.c_str());
//             DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a pklbin file");
//         }
//         else if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
//             temp_specs.LoadSpecSet_mgf(spectra_file_name.c_str());
//             DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a mgf file");
//         }
//         else{
//             DEBUG_MSG("CANNOTLOAD Search\t"<<spectra_file_name);
//         }
//         
//         //Writing the file name into the spectra
//         for(int spec_idx = 0; spec_idx < temp_specs.size(); spec_idx++){
//             //temp_specs[spec_idx].fileName = spectra_file_name;
//         }
//         
//         m_library.insert(m_library.begin(), temp_specs.begin(), temp_specs.end());
//     }
    
    
    for(int file_idx = 0; file_idx < input_spectra.size(); file_idx++){
        string spectra_file_name = input_spectra[file_idx];
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
        
        source_spectra.insert(source_spectra.begin(), temp_specs.begin(), temp_specs.end());
    }
    
    return true;
  }


  bool ExecSpectraExtractionSingle::saveOutputData(void){
    DEBUG_MSG("OPENING FILE STREAMS\t"<<results_dir);
    
    output_psms.saveToFile(results_dir.c_str());
    output_psms_new.saveToFile(new_lib_results_dir.c_str());

    //DEBUG_MSG("SAVING " << output_mgf_name);    
    
    //m_new_library.SaveSpecSet_mgf(output_mgf_name.c_str());
    

    DEBUG_MSG("DONE");            
      
    return true;
  }


  bool ExecSpectraExtractionSingle::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraExtractionSingle::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraExtractionSingle::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraExtractionSingle::merge(void){
    return true;
  }


  bool ExecSpectraExtractionSingle::validateParams(std::string & error){
    return true;
  }

}
