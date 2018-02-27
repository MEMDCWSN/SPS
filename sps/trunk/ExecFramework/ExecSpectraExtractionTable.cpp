//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraExtractionTable.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
//#include "mzxml.h"
#include <fstream>
#include "projectionutils.h"
#include "DelimitedTextReader.h"

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectraExtractionTable::ExecSpectraExtractionTable(void){
    m_name = "ExecSpectraExtractionTable";
    m_type = "ExecSpectraExtractionTable";
  }


  ExecSpectraExtractionTable::ExecSpectraExtractionTable(const ParameterList & inputParams){
    m_name = "ExecSpectraExtractionTable";
    m_type = "ExecSpectraExtractionTable";
    
    //inputParams.print(cout);
    
    this->m_params = inputParams;
    
    //output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");
    
    new_lib_results_dir = inputParams.getValue("NEWLIBRARYRESULTS_DIR", ".");
    
    //stringSplit(inputParams.getValue("EXISTING_LIBRARY_MGF", ""), existing_library_name);

    results_dir = inputParams.getValue("RESULTS_DIR", ".");
    
    //string library_path = inputParams.getValue("EXISTING_LIBRARY_MGF", "");
    
    library_grid_root_dir = inputParams.getValue("LIBUSER_ROOT_DATA_PATH", "");
    library_user_root_dir = inputParams.getValue("LIBUSER_ROOT_VM_DISPLAY_PATH", "");
    library_user = inputParams.getValue("LIBUSER_NAME", "");    //Should check if library path selected matches this
    
    input_table_filename = inputParams.getValue("INPUT_TABLE_PATH",""); 
        
    submit_user = inputParams.getValue("user","INVALIDUSER"); 

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
          path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
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
            
            
            if(full_lib_path.find("*") != string::npos || full_lib_path.find("..") != string::npos || full_lib_path.find("~") != string::npos ){
                ERROR_MSG("Invalid letters in path, exiting");
                exit(1);
            }
            
            //only can write to the owners file
            if(it->first.find(submit_user) != 0){
                ERROR_MSG("Invalid User File");
                exit(1);
            }
            else{
                DEBUG_MSG("Writing to User Space: " + submit_user);
            }
        }
    }

    
    for (std::map<string,string>::iterator it=file_name_mapping.begin(); it!=file_name_mapping.end(); ++it){
        DEBUG_MSG( it->first << " => " << it->second );
    } 
    

    
  }


  ExecSpectraExtractionTable::~ExecSpectraExtractionTable(void){
  }


  ExecBase * ExecSpectraExtractionTable::clone(const ParameterList & inputParams) const{
    return new ExecSpectraExtractionTable(inputParams);
  }

  bool ExecSpectraExtractionTable::invoke(void){
    load_new_spectra();
    
    add_to_library();
    
    return true;
      
    

    
    
    return true;
  }

  bool ExecSpectraExtractionTable::loadInputData(void){
    
    
    
    vector<string> required_headers;
    vector<int> required_headers_int;
    
    DelimitedTextReader::loadDelimitedFile(input_table_filename.c_str(), "\t", "#", annotation_headers_map, annotation_lines, required_headers, required_headers_int);
    
    for(int i = 0; i < annotation_lines.size(); i++){
        string output_string = "";
        for(int j = 0; j < annotation_lines[i].size(); j++){
            output_string += annotation_lines[i][j] + "\t";
        }
        DEBUG_MSG(output_string);
    }
    
    return true;
  }


  bool ExecSpectraExtractionTable::saveOutputData(void){
      
    //return true;
    DEBUG_MSG("WRITING OUTPUT\t"<<results_dir);
    DEBUG_MSG("WRITING OUTPUT\t"<<new_lib_results_dir);
    
    output_psms.saveToFile(results_dir.c_str());
    output_psms_new.saveToFile(new_lib_results_dir.c_str());

    //DEBUG_MSG("SAVING " << output_mgf_name);    
    
    //m_new_library.SaveSpecSet_mgf(output_mgf_name.c_str());
    

    DEBUG_MSG("DONE");            
      
    return true;
  }


  bool ExecSpectraExtractionTable::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraExtractionTable::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraExtractionTable::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraExtractionTable::merge(void){
    return true;
  }


  bool ExecSpectraExtractionTable::validateParams(std::string & error){
    return true;
  }
  
  int ExecSpectraExtractionTable::load_new_spectra(){
      bool error_state = false;
      
      //Validation of annotation_headers_map
        if(annotation_headers_map.find("FILENAME") == annotation_headers_map.end() ||
            annotation_headers_map.find("SEQ") == annotation_headers_map.end() ||
            annotation_headers_map.find("COMPOUND_NAME") == annotation_headers_map.end() ||
            annotation_headers_map.find("MOLECULEMASS") == annotation_headers_map.end() ||
            annotation_headers_map.find("EXACTMASS") == annotation_headers_map.end() ||
            annotation_headers_map.find("GENUS") == annotation_headers_map.end() ||
            annotation_headers_map.find("SPECIES") == annotation_headers_map.end() ||
            annotation_headers_map.find("STRAIN") == annotation_headers_map.end() ||
            annotation_headers_map.find("INSTRUMENT") == annotation_headers_map.end() ||
            annotation_headers_map.find("IONSOURCE") == annotation_headers_map.end() ||
            annotation_headers_map.find("EXTRACTSCAN") == annotation_headers_map.end() ||
            annotation_headers_map.find("SMILES") == annotation_headers_map.end() ||
            annotation_headers_map.find("INCHI") == annotation_headers_map.end() ||
            annotation_headers_map.find("INCHIAUX") == annotation_headers_map.end() ||
            annotation_headers_map.find("CHARGE") == annotation_headers_map.end() ||
            annotation_headers_map.find("IONMODE") == annotation_headers_map.end() ||
            annotation_headers_map.find("PI") == annotation_headers_map.end() ||
            annotation_headers_map.find("DATACOLLECTOR") == annotation_headers_map.end() ||
            annotation_headers_map.find("PUBMED") == annotation_headers_map.end() ||
            annotation_headers_map.find("ACQUISITION") == annotation_headers_map.end() ||
            annotation_headers_map.find("ADDUCT") == annotation_headers_map.end() ||
            annotation_headers_map.find("INTEREST") == annotation_headers_map.end() ||
            annotation_headers_map.find("LIBQUALITY") == annotation_headers_map.end() ||
            annotation_headers_map.find("SPECTRUMID") == annotation_headers_map.end() ||
            annotation_headers_map.find("CASNUMBER") == annotation_headers_map.end()){
            ERROR_MSG("HEADERS MALFORMED.\nPlease check that the headers are identical to the template file and no whitespace was introduced.");
            exit(1);
        }
      
      map<string, SpecSet> specset_cache;
      
      for(int annotation_index = 0; annotation_index < annotation_lines.size(); annotation_index++){
        DEBUG_MSG("Original "<<annotation_lines[annotation_index][annotation_headers_map["FILENAME"]]);
        //DEBUG_VAR(atof(annotation_lines[annotation_index][annotation_headers_map["MOLECULEMASS"]].c_str()));
        //DEBUG_VAR(annotation_headers_map["MOLECULEMASS"]);
         
        
        
        
        std::string original_input_file = get_only_filename(annotation_lines[annotation_index][annotation_headers_map["FILENAME"]]);
        
        std::string spectrum_sequence = annotation_lines[annotation_index][annotation_headers_map["SEQ"]];
        std::string compound_name = annotation_lines[annotation_index][annotation_headers_map["COMPOUND_NAME"]];
        float parentMZ = atof(annotation_lines[annotation_index][annotation_headers_map["MOLECULEMASS"]].c_str());
        float exact_mass = atof(annotation_lines[annotation_index][annotation_headers_map["EXACTMASS"]].c_str());
        std::string genus = annotation_lines[annotation_index][annotation_headers_map["GENUS"]];
        std::string species = annotation_lines[annotation_index][annotation_headers_map["SPECIES"]];
        std::string strain = annotation_lines[annotation_index][annotation_headers_map["STRAIN"]];
        std::string instrument = annotation_lines[annotation_index][annotation_headers_map["INSTRUMENT"]];
        std::string ionsource = annotation_lines[annotation_index][annotation_headers_map["IONSOURCE"]];
        std::string extract_scan_str = annotation_lines[annotation_index][annotation_headers_map["EXTRACTSCAN"]];
        std::string smiles = annotation_lines[annotation_index][annotation_headers_map["SMILES"]];
        std::string InChI = annotation_lines[annotation_index][annotation_headers_map["INCHI"]];
        std::string InChI_Aux = annotation_lines[annotation_index][annotation_headers_map["INCHIAUX"]];
        std::string notes;
        int precursor_charge = atoi(annotation_lines[annotation_index][annotation_headers_map["CHARGE"]].c_str());
        std::string ion_mode = annotation_lines[annotation_index][annotation_headers_map["IONMODE"]];
        std::string PI = annotation_lines[annotation_index][annotation_headers_map["PI"]];
        std::string data_collector = annotation_lines[annotation_index][annotation_headers_map["DATACOLLECTOR"]];
        std::string CASNumber = annotation_lines[annotation_index][annotation_headers_map["CASNUMBER"]];
        std::string PublicationID = annotation_lines[annotation_index][annotation_headers_map["PUBMED"]];
        std::string acquisition_method = annotation_lines[annotation_index][annotation_headers_map["ACQUISITION"]];
        std::string adduct = annotation_lines[annotation_index][annotation_headers_map["ADDUCT"]];
        std::string spectrum_novelty = annotation_lines[annotation_index][annotation_headers_map["INTEREST"]];
        int library_quality = atoi(annotation_lines[annotation_index][annotation_headers_map["LIBQUALITY"]].c_str());
        std::string spectrumID = annotation_lines[annotation_index][annotation_headers_map["SPECTRUMID"]];
        
        
        std::string submission_user;
        std::string submission_ID;
        std::string submission_date;
        
        submission_user = m_params.getValue("user");
        submission_ID = m_params.getValue("task");
        submission_date = m_params.getValue("submission_date", get_time_string().c_str());
        
        //Extracting Spectra
        int extraction_scan = -1;
        if(extract_scan_str != "*"){
            extraction_scan = atoi(extract_scan_str.c_str());
        }
        
        
        //Checking if input file name has a mangled equivalent
        if(file_name_mapping.find(original_input_file) == file_name_mapping.end()){
            ERROR_MSG(original_input_file<<" not uploaded, try again");
            error_state = true;
        }
        
        //Loading Input File
        string spectra_file_name = path_to_spectra + "/" + file_name_mapping[original_input_file];
        DEBUG_MSG("Mangled Name: "<<spectra_file_name);
        
        string extension = get_extension(spectra_file_name);
        
        //Checking if file is in the cache
        map<string,SpecSet>::iterator it = specset_cache.find(spectra_file_name);
        
        
        
        SpecSet source_spectra;
        
        if(it == specset_cache.end()){
            DEBUG_MSG("NOT IN CACHE");
            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
                string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
                source_spectra.Load(pklbin_version.c_str());
                DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
            }
            else if(strcmp(extension.c_str(), "pklbin") == 0){
                source_spectra.loadPklBin(spectra_file_name.c_str());
                DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a pklbin file");
            }
            else if(strcmp(extension.c_str(), "mgf") == 0 || strcmp(extension.c_str(), "MGF") == 0){
                source_spectra.LoadSpecSet_mgf(spectra_file_name.c_str());
                DEBUG_MSG("LOADED Search\t"<<spectra_file_name<<"\t as a mgf file");
            }
            else{
                DEBUG_MSG("CANNOTLOAD Search\t"<<spectra_file_name);
            }
            
            //Writing the file name into the spectra
            for(int spec_idx = 0; spec_idx < source_spectra.size(); spec_idx++){
                source_spectra[spec_idx].fileName = spectra_file_name;
            }
            
            specset_cache[spectra_file_name] = source_spectra;
        }
        else{
            DEBUG_MSG("IN CACHE");
            source_spectra = specset_cache[spectra_file_name];
        }
        
        
        
        
        
        
        //Looking for spectrum
        int found_source_spectrum = -1;
        DEBUG_VAR(extraction_scan);
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
            
            
            string full_organism_name = genus + ":" + species;
            if(strain.length() > 0){
                full_organism_name += ":" + strain;
            }
            
            if(adduct.length() > 0){
                compound_name = compound_name + " [" + adduct + "]";
            }
            
            psm->m_organism = (full_organism_name);
            //psm->m_organism = "GNPS";
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
            DEBUG_VAR(psm->m_organism);
            
            extracted_spec.psmList.clear();
            extracted_spec.psmList.push_back(psm);
            int added_index = m_new_library.add_update_spectrum_to_Library(extracted_spec);
            DEBUG_MSG("FOUND SPEC ADDED AT INDEX "<<added_index);
            if(added_index == -2)
                continue;
        }
        
        if(found_source_spectrum == -1){
            ERROR_MSG("DID NOT FIND SPECTRUM IN FILE\t" << original_input_file << "\tSCAN: " << extraction_scan);
            error_state = true;
        }
        else{
            DEBUG_MSG("FOUND SPECTRUM IN FILE");
        }
    }
    
    if(error_state){
        ERROR_MSG("ERRORS, FAILURE");
        exit(1);
    }
    
    return 0;
  }
  
  void ExecSpectraExtractionTable::add_to_library(){
    if(existing_library_name.size() != 1){
        ERROR_MSG("ERROR, Too many Libraries");
        exit(1);
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
        
        
        std::string output_new_spectrum_file = m_params.getValue("ADDSPECTRA_OUTPUTNEWSPECTRAFILE", "");
        
        if(output_new_spectrum_file.length() > 0){
            m_new_library.SaveSpecSet_mgf( output_new_spectrum_file.c_str());
        }
        
        //temp_specs.SaveSpecSet_mgf( (spectra_file_name + "_temp").c_str());
        existing_lib.SaveSpecSet_mgf( output_library_name.c_str());
    }
  }

}
