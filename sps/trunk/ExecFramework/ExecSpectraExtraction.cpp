//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraExtraction.h"

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

  ExecSpectraExtraction::ExecSpectraExtraction(void){
    m_name = "ExecSpectraExtraction";
    m_type = "ExecSpectraExtraction";
  }


  ExecSpectraExtraction::ExecSpectraExtraction(const ParameterList & inputParams){
    m_name = "ExecSpectraExtraction";
    m_type = "ExecSpectraExtraction";
    inputParams.print(cout);

    input_spectra_data_file = inputParams.getValue("INPUTSPECTRA_DATAFILE", "");

    //checking if CCMS
    if(inputParams.getValueInt("CCMS", 0) == 1){
        //input_spectra_data_file = input_spectra_data_file + "/toolParams-00000.txt";
        input_spectra_data_file = input_spectra_data_file;
    }


    submission_user = inputParams.getValue("submission_user");
    submission_ID = inputParams.getValue("submission_ID");
    submission_date = inputParams.getValue("submission_date", get_time_string().c_str());

    output_mgf_name = inputParams.getValue("OUTPUT_MGF", "");

    existing_library_name = inputParams.getValue("EXISTING_LIBRARY_MGF", "");

    results_dir = inputParams.getValue("RESULTS_DIR", ".");
    new_lib_results_dir = inputParams.getValue("NEWLIBRARYRESULTS_DIR", ".");
    spectra_dir = inputParams.getValue("SPECTRA_DIR", ".");

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

      mapping_count++;
    }
    
    for (std::map<string,string>::iterator it=file_name_mapping.begin(); it!=file_name_mapping.end(); ++it){
        std::cout << it->first << " => " << it->second << '\n';
    }
  }


  ExecSpectraExtraction::~ExecSpectraExtraction(void){
  }


  ExecBase * ExecSpectraExtraction::clone(const ParameterList & inputParams) const{
    return new ExecSpectraExtraction(inputParams);
  }

  bool ExecSpectraExtraction::invoke(void){



    return true;
  }

  bool ExecSpectraExtraction::loadInputData(void){

    DEBUG_MSG("LOADING INPUT");


    ifstream input_spectra_metadata(input_spectra_data_file.c_str(), ios::binary);
    string line;
    SpectralLibrary extracted_specs;

    if(existing_library_name.length() > 0){
        cout<<"Loading Existing library"<<endl;
        extracted_specs.LoadSpecSet_mgf(existing_library_name.c_str());
    }


    cout<<"Opened\t"<<input_spectra_data_file<<endl;

    stringstream new_annotations_string (stringstream::in | stringstream::out);
    stringstream all_annotations_string (stringstream::in | stringstream::out);

    //Results Headers
    new_annotations_string<<"Filename\t"<<"Scan\t"<<"LibraryIdx\t"<<"Peptide\t"<<"Organism\t"<<"Smiles\t"<<"Inchi\t"<<"InchiAux\t"<<"SpectrumQuality\t"<<"CompoundName\t"<<"Charge\t"<<"IonMode\t"<<"Notes"<<endl;
    all_annotations_string<<"Filename\t"<<"Scan\t"<<"LibraryIdx\t"<<"Peptide\t"<<"Organism\t"<<"Smiles\t"<<"Inchi\t"<<"InchiAux\t"<<"SpectrumQuality\t"<<"CompoundName\t"<<"Charge\t"<<"IonMode\t"<<"Notes"<<endl;

    if(!input_spectra_metadata.is_open()){
        DEBUG_MSG("COULD NOT OPEN FILE: "<<input_spectra_data_file);
    }
    else{
        getline (input_spectra_metadata ,line);
        //cout<<line<<endl;
        int line_counts = 0;
        while(input_spectra_metadata.good()){
            getline (input_spectra_metadata ,line);

            //cout<<line<<endl;

            if(line.length() < 1)
                continue;   //empty line
            DEBUG_MSG("PARSING META DATA LINE");
            stringstream metadata_line_stream;
            metadata_line_stream<<line;

            std::string input_mz_file_name;
            std::string spectrum_sequence;
            std::string compound_name;
            float molecule_mass;
            std::string molecule_mass_str;
            std::string protein;
            std::string organism;
            float ITOL;
            std::string ITOL_str;
            std::string ITOLU;
            float TOL;
            std::string TOL_str;
            std::string TOLU;
            std::string source_instrument;
            int extract_scan;
            std::string extract_scan_str;
            int spectrum_quality;
            std::string spectrum_quality_str;
            std::string smiles;
            std::string InChI;
            std::string InChI_Aux;
            std::string notes;
            std::string precursor_charge_str;
            int precursor_charge;
            std::string ion_mode;



            /*metadata_line_stream>>input_mz_file_name;
            metadata_line_stream>>spectrum_sequence;
            metadata_line_stream>>compound_name;
            metadata_line_stream>>molecule_mass;
            metadata_line_stream>>protein;
            metadata_line_stream>>organism;
            metadata_line_stream>>ITOL;
            metadata_line_stream>>ITOLU;
            metadata_line_stream>>TOL;
            metadata_line_stream>>TOLU;
            metadata_line_stream>>source_instrument;
            metadata_line_stream>>extract_scan;
            metadata_line_stream>>spectrum_quality;
            metadata_line_stream>>smiles;
            metadata_line_stream>>InChI;
            metadata_line_stream>>InChI_Aux;*/

            getline(metadata_line_stream, input_mz_file_name, '\t');
            getline(metadata_line_stream, spectrum_sequence, '\t');
            getline(metadata_line_stream, compound_name, '\t');
            getline(metadata_line_stream, molecule_mass_str, '\t');
            getline(metadata_line_stream, protein, '\t');
            getline(metadata_line_stream, organism, '\t');
            getline(metadata_line_stream, ITOL_str, '\t');
            getline(metadata_line_stream, ITOLU, '\t');
            getline(metadata_line_stream, TOL_str, '\t');
            getline(metadata_line_stream, TOLU, '\t');
            getline(metadata_line_stream, source_instrument, '\t');
            getline(metadata_line_stream, extract_scan_str, '\t');
            getline(metadata_line_stream, spectrum_quality_str, '\t');
            getline(metadata_line_stream, smiles, '\t');
            getline(metadata_line_stream, InChI, '\t');
            getline(metadata_line_stream, InChI_Aux, '\t');
            getline(metadata_line_stream, notes, '\t');
            getline(metadata_line_stream, precursor_charge_str, '\t');
            getline(metadata_line_stream, ion_mode, '\t');

            molecule_mass = atof(molecule_mass_str.c_str());
            ITOL = atof(ITOL_str.c_str());
            TOL = atof(TOL_str.c_str());
            spectrum_quality = atoi(spectrum_quality_str.c_str());
            precursor_charge = atoi(ion_mode.c_str());

            if(extract_scan_str.compare("*") == 0){
                extract_scan = ExecSpectraExtraction::ALL;
                compound_name = "ORGANISM:"+organism;
            }
            else{
                extract_scan = atoi(extract_scan_str.c_str());
            }

            ion_mode = remove_char(ion_mode, '\r');
            ion_mode = remove_char(ion_mode, '\n');


            input_mz_file_name = get_only_filename(input_mz_file_name);
            int dotpos = input_mz_file_name.find_last_of('.');
            string extension = input_mz_file_name.substr(dotpos+1);


            SpecSet specs;
            vector<short> mslevel;


            string mangled_name;



            if( file_name_mapping.find(input_mz_file_name) != file_name_mapping.end()){
                mangled_name = spectra_dir + "/" + file_name_mapping[input_mz_file_name];
                DEBUG_MSG("Mangled from "<<input_mz_file_name<<" to "<<mangled_name);
            }
            else{
                mangled_name = input_mz_file_name;
            }

            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
                //LoadMzxml(mangled_name.c_str(), specs, &mslevel, 0);
                string pklbin_version = strip_extension(mangled_name) + ".pklbin";
                specs.Load(pklbin_version.c_str());
                DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<mangled_name);
            }
            if(strcmp(extension.c_str(), "pklbin") == 0 || strcmp(extension.c_str(), "pklbin") == 0){
                specs.Load(mangled_name.c_str());
                DEBUG_MSG("Loaded "<<mangled_name<< "\t"<<specs.size()); 
            }
            else if(strcmp(extension.c_str(), "mgf") == 0){
                specs.LoadSpecSet_mgf(mangled_name.c_str());
            }

            DEBUG_MSG("LOADING file " + input_mz_file_name + " as mangled " + mangled_name);
            line_counts++;

            DEBUG_MSG("LOOKING FOR SCAN: "<<extract_scan<<"\t"<<source_instrument);
            bool found = false;
            for(int i = 0; i < specs.size(); i++){
                if( (specs[i].scan == extract_scan || extract_scan == ExecSpectraExtraction::ALL ) && specs[i].msLevel != 1 ){
                    
                    if(! (extract_scan == ExecSpectraExtraction::ALL))
                        DEBUG_MSG("FOUND SCAN "<<specs[i].scan<<"\tMSLEVEL\t"<<specs[i].msLevel);
                    
                    found = true;
                    Spectrum extracted_spec;
                    bool scan_found = false;
                    int found_index = -1;
                    scan_found = true;
                    found_index = i;
                    extracted_spec = specs[i];

                    extracted_spec.instrument_name = source_instrument;
                    extracted_spec.ITOL = ITOL;
                    extracted_spec.ITOLU = ITOLU;
                    extracted_spec.TOL = TOL;
                    extracted_spec.TOLU = TOLU;
                    extracted_spec.msLevel = 2;
                    extracted_spec.spectrum_quality = spectrum_quality;
                    
                    //Correcting for precursor MZ
                    if(molecule_mass > 1.f){
                        extracted_spec.parentMZ = molecule_mass;
                    }

                    
                    //if(input_mz_file_name.find_last_of("/") != string::npos)
                    extracted_spec.fileName = input_mz_file_name;

                    psmPtr psm(new PeptideSpectrumMatch());

                    psm->m_scanNum = specs[i].scan;
                    psm->m_annotation = spectrum_sequence;
                    string submission_metadata = submission_user  + ":"  + submission_ID + ":" + submission_date;
                    psm->m_submission_metadata = (submission_metadata);
                    
                    psm->m_organism = (organism);
                    psm->m_compound_name = (compound_name);
                    psm->m_smiles = (smiles);
                    psm->m_InChI = (InChI);
                    psm->m_InChI_Aux = (InChI_Aux);
                    
                    psm->m_notes = notes;
                    if(precursor_charge == 0){
                        psm->m_charge = extracted_spec.parentCharge;
                    }
                    else{
                        psm->m_charge = precursor_charge;
                    }
                    psm->m_ionmode = ion_mode;

                    extracted_spec.psmList.push_back(psm);

                    //DEBUG_MSG("PUSHING BACK PSM");
                    int added_index = extracted_specs.add_update_spectrum_to_Library(extracted_spec);
                    if(added_index == -2)
                        continue;
                    
                    
                    
                    new_annotations_string<<input_mz_file_name<<"\t"<<extracted_spec.scan<<"\t"<<added_index<<"\t";;

                    if(extracted_spec.psmList.front()->m_annotation == "*.NOTPEPTIDE.*" || extracted_spec.psmList.front()->m_annotation == "NOTPEPTIDE")
                        new_annotations_string<<"*..*"<<"\t";
                    else
                        new_annotations_string<<extracted_spec.psmList.front()->m_annotation<<"\t";


                    new_annotations_string<<extracted_spec.psmList.front()->m_organism<<"\t";
                    new_annotations_string<<extracted_spec.psmList.front()->m_smiles<<"\t";
                    new_annotations_string<<extracted_spec.psmList.front()->m_InChI<<"\t";
                    new_annotations_string<<extracted_spec.psmList.front()->m_InChI_Aux<<"\t";
                    new_annotations_string<<extracted_spec.spectrum_quality<<"\t";
                    new_annotations_string<<extracted_spec.psmList.front()->m_compound_name<<"\t";
                    new_annotations_string<<extracted_spec.psmList.front()->m_charge<<"\t";
                    if(extracted_spec.psmList.front()->m_ionmode.length() == 0)
                        new_annotations_string<<" "<<"\t";
                    else
                        new_annotations_string<<extracted_spec.psmList.front()->m_ionmode<<"\t";    
                    
                    if(extracted_spec.psmList.front()->m_notes.length() == 0)
                        new_annotations_string<<" "<<endl;
                    else
                        new_annotations_string<<extracted_spec.psmList.front()->m_notes<<endl;

                    if(extract_scan != ExecSpectraExtraction::ALL)
                        break;
                }
            }
            
            if(found == false){
                DEBUG_MSG("SPECTRUM IN ANNOTATION NOTFOUND");
            }


            //if(extracted_spec.psmList.size() == 0)
            //    continue;

            //Checking if Spectra are already present in library

            //extracted_specs.push_back(extracted_spec);

            DEBUG_MSG("SPEC SIZE: " << extracted_specs.size());


        }
        DEBUG_VAR(line_counts);
    }

    DEBUG_MSG("OUTPUTTING ALL ANNOTATIONS");


    for(int library_idx = 0; library_idx < extracted_specs.size(); library_idx++){
      all_annotations_string<<extracted_specs[library_idx].fileName<<"\t"<<extracted_specs[library_idx].scan<<"\t"<<library_idx+1<<"\t";

      if(extracted_specs[library_idx].psmList.front()->m_annotation == "*.NOTPEPTIDE.*" || extracted_specs[library_idx].psmList.front()->m_annotation == "NOTPEPTIDE")
          all_annotations_string<<"*..*"<<"\t";
      else
          all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_annotation<<"\t";


      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_organism<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_smiles<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_InChI<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_InChI_Aux<<"\t";
      all_annotations_string<<extracted_specs[library_idx].spectrum_quality<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_compound_name<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_charge<<"\t";
      
      if(extracted_specs[library_idx].psmList.front()->m_ionmode.length() == 0)
          all_annotations_string<<" "<<"\t";
      else
          all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_ionmode<<"\t";    
       
      if(extracted_specs[library_idx].psmList.front()->m_notes.length() == 0)
          all_annotations_string<<" "<<endl;
      else
          all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_notes<<endl;
      
      
      //all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_ionmode<<"\t";
      //all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_notes<<endl;
    }

    //cout<<"Writing\t"<<(results_dir + "/new_annotations.tsv")<<endl;
    //cout<<"Writing\t"<<(results_dir + "/all_annotations.tsv")<<endl;
    DEBUG_MSG("OPENING FILE STREAMS\t"<<results_dir<<"\t"<<new_lib_results_dir);
    //fstream new_file_stream( (results_dir + "/new_annotations.tsv").c_str(), fstream::out | fstream::binary);
    //fstream all_file_stream( (results_dir + "/all_annotations.tsv").c_str(), fstream::out | fstream::binary);
    fstream all_file_stream( (results_dir).c_str(), fstream::out | fstream::binary);
    fstream new_file_stream( (new_lib_results_dir).c_str(), fstream::out | fstream::binary);

    DEBUG_MSG("WRITING STREAMS");

    all_file_stream<<all_annotations_string.str();
    all_file_stream.close();

    new_file_stream<<new_annotations_string.str();
    new_file_stream.close();

    DEBUG_MSG("SAVING MGF");

    output_mgf_name += "/library.mgf";
    extracted_specs.SaveSpecSet_mgf(output_mgf_name.c_str());

    DEBUG_MSG("DONE");

    return true;
  }


  bool ExecSpectraExtraction::saveOutputData(void){
        return true;
  }


  bool ExecSpectraExtraction::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraExtraction::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraExtraction::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraExtraction::merge(void){
    return true;
  }


  bool ExecSpectraExtraction::validateParams(std::string & error){
    return true;
  }

}
