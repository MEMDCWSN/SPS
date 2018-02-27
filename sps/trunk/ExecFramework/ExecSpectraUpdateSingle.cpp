//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraUpdateSingle.h"

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

  ExecSpectraUpdateSingle::ExecSpectraUpdateSingle(void){
    m_name = "ExecSpectraUpdateSingle";
    m_type = "ExecSpectraUpdateSingle";
  }


  ExecSpectraUpdateSingle::ExecSpectraUpdateSingle(const ParameterList & inputParams){
    m_name = "ExecSpectraUpdateSingle";
    m_type = "ExecSpectraUpdateSingle";
    
    inputParams.print(cout);
    
    this->m_params = inputParams;
    
    //output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");
    
    new_lib_results_dir = inputParams.getValue("NEWLIBRARYRESULTS_DIR", ".");
    
    //stringSplit(inputParams.getValue("EXISTING_LIBRARY_MGF", ""), existing_library_name);

    results_dir = inputParams.getValue("RESULTS_DIR", ".");
    
    //string library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");
    
    library_user_root_dir = inputParams.getValue("LIBUSER_ROOT_DATA_PATH", "");
    library_user = inputParams.getValue("LIBUSER_NAME", "");    //Should check if library path selected matches this
    
    
    
    
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

      mapping_count++;
    }
    
    for (std::map<string,string>::iterator it=full_file_name_mapping.begin(); it!=full_file_name_mapping.end(); ++it){
        if(it->second.find_first_of("lib-") == 0){
            //We have a library, find original name
            string full_lib_path = library_user_root_dir  + it->first;
            existing_library_name.push_back(full_lib_path);
            DEBUG_VAR(full_lib_path);
            DEBUG_VAR(it->first.find_first_of(library_user));
            
            if(full_lib_path.find("*") != string::npos || full_lib_path.find("..") != string::npos || full_lib_path.find("~") != string::npos || it->first.find_first_of(library_user) != 0){
                DEBUG_MSG("Invalid letters in path, exiting");
                exit(1);
            }
        }
    }

    
    for (std::map<string,string>::iterator it=file_name_mapping.begin(); it!=file_name_mapping.end(); ++it){
        std::cout << it->first << " => " << it->second << '\n';
    } 
    
    //Renaming output the same as before
    //if(existing_library_name.size() == 1)
    //    output_mgf_name = output_mgf_name + "/" + get_only_filename(existing_library_name[0]);
    
  }


  ExecSpectraUpdateSingle::~ExecSpectraUpdateSingle(void){
  }


  ExecBase * ExecSpectraUpdateSingle::clone(const ParameterList & inputParams) const{
    return new ExecSpectraUpdateSingle(inputParams);
  }

  bool ExecSpectraUpdateSingle::invoke(void){
    
    std::string spectrum_sequence;
    std::string compound_name;
    float parentMZ;
    float exact_mass;
    std::string genus;
    std::string species;
    std::string strain;
    std::string instrument;
    std::string ionsource;
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
    
    submission_user = m_params.getValue("user");
    submission_ID = m_params.getValue("task");
    submission_date = m_params.getValue("submission_date", get_time_string().c_str());
    
    spectrum_sequence = m_params.getValue("ADDSPECTRA_SEQ", "*..*");
    compound_name = m_params.getValue("ADDSPECTRA_COMPOUND_NAME");
    parentMZ = m_params.getValueFloat("ADDSPECTRA_MOLECULEMASS", 0.f);
    genus = m_params.getValue("ADDSPECTRA_GENUS");
    species = m_params.getValue("ADDSPECTRA_SPECIES");
    strain = m_params.getValue("ADDSPECTRA_STRAIN");
    instrument = m_params.getValue("ADDSPECTRA_INSTRUMENT");
    ionsource = m_params.getValue("ADDSPECTRA_IONSOURCE");
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
    spectrumID = m_params.getValue("SPECTRUMID", "NO_ID");
    
    /*
    for(int i = 0; i < source_spectra.size(); i++){
        if(extraction_scan != -1 && extraction_scan != source_spectra[i].scan){
            continue;
        }
        
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
        
        
    }*/
    
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
            temp->m_library_name = spectra_file_name;
            
            if(existing_lib[library_idx].fileName.length() > 0)
                temp->m_spectrumFile = existing_lib[library_idx].fileName;
            
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
        DEBUG_VAR(output_library_name);
        
        bool library_updated = false;
        
        //Updating Annotation
        for(int library_idx = 0; library_idx < existing_lib.size(); library_idx++){
            DEBUG_MSG(spectrumID<<"\t"<<existing_lib[library_idx].psmList.front()->m_spectrumID);
            if(spectrumID == existing_lib[library_idx].psmList.front()->m_spectrumID){
                //Checking for new updates
                string new_id_instrument = ionsource + "-" + instrument;
                if(new_id_instrument != existing_lib[library_idx].instrument_name){
                    //Something happens
                }
                psmPtr psm(new PeptideSpectrumMatch());
                *psm = *(existing_lib[library_idx].psmList.front());
                psm->m_annotation = spectrum_sequence;
                string submission_metadata = submission_user  + ":"  + submission_ID + ":" + submission_date;
                psm->m_submission_metadata = (submission_metadata);
                
                string compound_with_adduct = compound_name;
                if(adduct.length() > 0){
                    compound_with_adduct = compound_with_adduct + " [" + adduct + "]";
                }
                
                psm->m_compound_name = (compound_with_adduct);
                psm->m_smiles = (smiles);
                psm->m_InChI = (InChI);
                psm->m_InChI_Aux = (InChI_Aux);
                psm->m_notes = PI + ":" + data_collector + ":" + CASNumber + ":" + PublicationID + ":" + acquisition_method;
                
                if(precursor_charge != 0){
                    psm->m_charge = precursor_charge;
                }
                
                psm->m_ionmode = ion_mode;
                psm->m_exactmass = exact_mass;
                
                psm->m_library_quality = library_quality;
                
                int does_differ = diff_molecular_psm(existing_lib[library_idx].psmList.front(), psm);
                DEBUG_VAR(does_differ);
                if(does_differ != 0){
                    DEBUG_MSG("UPDATED ANNOTATION");
                    library_updated = true;
                    existing_lib[library_idx].psmList.front() = psm;
                    
                    //Making output psm
                    psm->m_dbIndex = library_idx + 1;
                    psm->m_protein = new_id_instrument;
                    psm->m_library_name = output_library_name;
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
                
                
                
            }
        }
        
        //WRITING LIB BACK
        if(library_updated){
            existing_lib.SaveSpecSet_mgf( output_library_name.c_str());
        }
    }
    
    
    
    return true;
  }

  bool ExecSpectraUpdateSingle::loadInputData(void){    
    return true;
  }


  bool ExecSpectraUpdateSingle::saveOutputData(void){
    DEBUG_MSG("OPENING FILE STREAMS\t"<<results_dir);
    
    output_psms.saveToFile(results_dir.c_str());
    output_psms_new.saveToFile(new_lib_results_dir.c_str());

    //DEBUG_MSG("SAVING " << output_mgf_name);    
    
    //m_new_library.SaveSpecSet_mgf(output_mgf_name.c_str());
    

    DEBUG_MSG("DONE");            
      
    return true;
  }


  bool ExecSpectraUpdateSingle::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraUpdateSingle::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraUpdateSingle::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraUpdateSingle::merge(void){
    return true;
  }


  bool ExecSpectraUpdateSingle::validateParams(std::string & error){
    return true;
  }

}
