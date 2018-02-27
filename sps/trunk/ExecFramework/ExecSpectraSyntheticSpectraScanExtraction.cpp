//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraSyntheticSpectraScanExtraction.h"

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

  ExecSpectraSyntheticSpectraScanExtraction::ExecSpectraSyntheticSpectraScanExtraction(void){
    m_name = "ExecSpectraSyntheticSpectraScanExtraction";
    m_type = "ExecSpectraSyntheticSpectraScanExtraction";
  }


  ExecSpectraSyntheticSpectraScanExtraction::ExecSpectraSyntheticSpectraScanExtraction(const ParameterList & inputParams){
    m_name = "ExecSpectraSyntheticSpectraScanExtraction";
    m_type = "ExecSpectraSyntheticSpectraScanExtraction";
    
    this->m_params = inputParams;
    
    input_table_filename = inputParams.getValue("INPUT_TABLE", "");
    output_table_filename = inputParams.getValue("OUTPUT_FILE", "");
  }


  ExecSpectraSyntheticSpectraScanExtraction::~ExecSpectraSyntheticSpectraScanExtraction(void){
  }


  ExecBase * ExecSpectraSyntheticSpectraScanExtraction::clone(const ParameterList & inputParams) const{
    return new ExecSpectraSyntheticSpectraScanExtraction(inputParams);
  }

  bool ExecSpectraSyntheticSpectraScanExtraction::invoke(void){
        if(annotation_headers_map.find("Filename") == annotation_headers_map.end() ||
            annotation_headers_map.find("monoisotopic") == annotation_headers_map.end() ||
            annotation_headers_map.find("min_precursor_int") == annotation_headers_map.end() ||
            annotation_headers_map.find("ppm") == annotation_headers_map.end() ||
            annotation_headers_map.find("compoundname") == annotation_headers_map.end()){
            ERROR_MSG("HEADERS MALFORMED.\nPlease check that the headers are identical to the template file and no whitespace was introduced.");
            exit(1);
        }
        
        fstream output_file_stream( (output_table_filename).c_str(), fstream::out | fstream::binary);
        
        output_file_stream<<"FILENAME"<<"\t";
        output_file_stream<<"SCAN"<<"\t";
        output_file_stream<<"ADDUCT"<<"\t";
        output_file_stream<<"PPMERROR"<<"\t";
        output_file_stream<<"ExperimentalMZ"<<"\t";
        output_file_stream<<"COMPOUNDNAME"<<"\t";
        output_file_stream<<"PRECURSOR_INT"<<"\t";
        output_file_stream<<"PPM"<<"\t";
        output_file_stream<<"Peak_After_SNR_Filter(5.0)"<<"\t";
        output_file_stream<<"TotalOtherHits"<<"\t";
        output_file_stream<<"TheoreticalMZ"<<"\t";
        output_file_stream<<"DeltaMZ"<<"\t";
        //output_file_stream<<"Charge"<<"\t";
        
        for(int i = 0; i < annotation_headers_map.size(); i++){
            output_file_stream<<"_"<<"\t";
        }
        
        //Headers for Batch Template
        output_file_stream<<"FILENAME"<<"\t";
        output_file_stream<<"SEQ"<<"\t";
        output_file_stream<<"COMPOUND_NAME"<<"\t";
        output_file_stream<<"MOLECULEMASS"<<"\t";
        output_file_stream<<"INSTRUMENT"<<"\t";
        output_file_stream<<"IONSOURCE"<<"\t";
        output_file_stream<<"EXTRACTSCAN"<<"\t";
        output_file_stream<<"SMILES"<<"\t";
        output_file_stream<<"INCHI"<<"\t";
        output_file_stream<<"INCHIAUX"<<"\t";
        output_file_stream<<"CHARGE"<<"\t";
        output_file_stream<<"IONMODE"<<"\t";
        output_file_stream<<"PUBMED"<<"\t";
        output_file_stream<<"ACQUISITION"<<"\t";
        output_file_stream<<"EXACTMASS"<<"\t";
        output_file_stream<<"DATACOLLECTOR"<<"\t";
        output_file_stream<<"ADDUCT"<<"\t";
        output_file_stream<<"INTEREST"<<"\t";
        output_file_stream<<"LIBQUALITY"<<"\t";
        output_file_stream<<"GENUS"<<"\t";
        output_file_stream<<"SPECIES"<<"\t";
        output_file_stream<<"STRAIN"<<"\t";
        output_file_stream<<"CASNUMBER"<<"\t";
        output_file_stream<<"PI"<<"\t";
        
        
        output_file_stream<<"\n";
        
        set<string> all_input_files;
        map<float, int> mass_occurences_in_samples;
        map<float, float> mass_occurences_ppm_tolerance;
        
        //Cataloging files and values to look for
        for(int annotation_index = 0; annotation_index < annotation_lines.size(); annotation_index++){
            std::string spectra_file_name = get_only_filename(annotation_lines[annotation_index][annotation_headers_map["Filename"]]);
            all_input_files.insert(spectra_file_name);
            
            float ppm = atof(annotation_lines[annotation_index][annotation_headers_map["ppm"]].c_str());
            
            float hydrogen_mass = 1.007825;
            float monoisotopic_mass = atof(annotation_lines[annotation_index][annotation_headers_map["monoisotopic"]].c_str());
            
            //float m_h_adduct_mz = atof(annotation_lines[annotation_index][annotation_headers_map["mz"]].c_str());
            float m_h_adduct_mz = monoisotopic_mass + hydrogen_mass;
            float sodium_adduct_mz = monoisotopic_mass + 22.989218;
            float potasium_adduct_mz = monoisotopic_mass + 38.963158;
            float double_charged_m_2h = (monoisotopic_mass + 2 * hydrogen_mass)/2;
            
            mass_occurences_in_samples[m_h_adduct_mz] = 0;
            mass_occurences_in_samples[sodium_adduct_mz] = 0;
            mass_occurences_in_samples[potasium_adduct_mz] = 0;
            mass_occurences_in_samples[double_charged_m_2h] = 0;
            
            mass_occurences_ppm_tolerance[m_h_adduct_mz] = ppm;
            mass_occurences_ppm_tolerance[sodium_adduct_mz] = ppm;
            mass_occurences_ppm_tolerance[potasium_adduct_mz] = ppm;
            mass_occurences_ppm_tolerance[double_charged_m_2h] = ppm;
            
        }
        
        //Iterating over all files
        for (set<string>::iterator it=all_input_files.begin(); it!=all_input_files.end(); ++it){
            string spectra_file_name = *it;
            
            SpecSet source_spectra;
            string extension = get_extension(spectra_file_name);
            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0 || strcmp(extension.c_str(), "mzML") == 0){
                string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
                source_spectra.Load(pklbin_version.c_str());
                DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
            }
            else{
                DEBUG_MSG("Cannot load: "<<spectra_file_name);
                continue;
            }
            
            for(int spec_idx = 0; spec_idx < source_spectra.size(); spec_idx++){
                float spec_mz = source_spectra[spec_idx].parentMZ;
                
                //Iterating over all masses
                for (std::map<float,int>::iterator it=mass_occurences_in_samples.begin(); it!=mass_occurences_in_samples.end(); ++it){
                    float search_mass = it->first;
                    float ppm_error = abs(spec_mz - search_mass)/search_mass * 1000000;
                    float ppm_tolerance = mass_occurences_ppm_tolerance[search_mass];
                    
                    if(ppm_error < ppm_tolerance){
                        mass_occurences_in_samples[search_mass]++;
                    }
                }
            }
        }
        
        map<string, float> compound_to_max_abundance;
        map<string, string> compound_to_plot_string;
        
        for(int annotation_index = 0; annotation_index < annotation_lines.size(); annotation_index++){
            DEBUG_VAR(annotation_index);
            DEBUG_VAR(annotation_lines[annotation_index].size());
            std::string spectra_file_name = get_only_filename(annotation_lines[annotation_index][annotation_headers_map["Filename"]]);
            std::string compound_name = (annotation_lines[annotation_index][annotation_headers_map["compoundname"]]);
            std::string pubmed_id = (annotation_lines[annotation_index][annotation_headers_map["PubChem ID"]]);
            float monoisotopic_mass = atof(annotation_lines[annotation_index][annotation_headers_map["monoisotopic"]].c_str());
            float ppm = atof(annotation_lines[annotation_index][annotation_headers_map["ppm"]].c_str());
            float min_precursor_int = atof(annotation_lines[annotation_index][annotation_headers_map["min_precursor_int"]].c_str());
            
            
            DEBUG_MSG(spectra_file_name<<"\t"<<monoisotopic_mass<<"\t"<<compound_name);
            
            SpecSet source_spectra;
            string extension = get_extension(spectra_file_name);
            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0  || strcmp(extension.c_str(), "mzML") == 0){
                string pklbin_version = strip_extension(spectra_file_name) + ".pklbin";
                source_spectra.Load(pklbin_version.c_str());
                DEBUG_MSG("Loaded "<<pklbin_version<<" instead of " <<spectra_file_name);
            }
            else{
                DEBUG_MSG("Cannot load: "<<spectra_file_name);
                continue;
            }
            
            
            string orig_string = "";
            for(int i = 0; i < annotation_lines[annotation_index].size(); i++){
                orig_string += annotation_lines[annotation_index][i] + "\t";
            }
            
            bool found_flag = false;
            bool found_flag_m_h = false;
            //Look through spectra
            for(int spec_idx = 0; spec_idx < source_spectra.size(); spec_idx++){
                float spec_precursor = source_spectra[spec_idx].precursor_intensity;//getTotalIonCurrent
                
                if(spec_precursor < min_precursor_int)
                    continue;
                
                float spec_mz = source_spectra[spec_idx].parentMZ;
                int scan = source_spectra[spec_idx].scan;
                
                float hydrogen_mass = 1.007825;
                
                float m_h_adduct = monoisotopic_mass + hydrogen_mass;
                float sodium_adduct_mz = monoisotopic_mass + 22.989218;
                float potasium_adduct_mz = monoisotopic_mass + 38.963158;                
                
                
                float double_charged_m_2h = (monoisotopic_mass + 2 * hydrogen_mass)/2;
                

                
                source_spectra[spec_idx].filterLowSNR(5.0);
                float peaks_after_snr_filter = source_spectra[spec_idx].size();
                
                map<string, float> adducts_searched;
                adducts_searched["M+H"] = m_h_adduct;
                //adducts_searched["[M+Na]"] = sodium_adduct_mz;
                //adducts_searched["[M+K]"] = potasium_adduct_mz;
                //adducts_searched["[M+2H]"] = double_charged_m_2h;
                
                for (std::map<string,float>::iterator it=adducts_searched.begin(); it!=adducts_searched.end(); ++it){
                    float mass_search = it->second;
                    
                    float ppm_error = abs(spec_mz - mass_search)/mass_search * 1000000;
                    
                    if(ppm > ppm_error){
                        output_file_stream<<spectra_file_name<<"\t";
                        output_file_stream<<scan<<"\t";
                        output_file_stream<<it->first<<"\t";
                        output_file_stream<<ppm_error<<"\t";
                        output_file_stream<<spec_mz<<"\t";
                        output_file_stream<<compound_name<<"\t";
                        output_file_stream<<spec_precursor<<"\t";
                        output_file_stream<<ppm_error<<"\t";
                        output_file_stream<<peaks_after_snr_filter<<"\t";
                        output_file_stream<<mass_occurences_in_samples[mass_search]<<"\t";
                        output_file_stream<<mass_search<<"\t";
                        output_file_stream<<abs(spec_mz - mass_search)<<"\t";
                        output_file_stream<<orig_string<<"\t";
                        
                        
                        //Outputting template information
                        
                        output_file_stream<<spectra_file_name<<"\t";
                        output_file_stream<<"*..*"<<"\t";
                        output_file_stream<<compound_name<<"\t";
                        output_file_stream<<mass_search<<"\t";
                        output_file_stream<<"qTof"<<"\t";
                        output_file_stream<<"LC-ESI"<<"\t";
                        output_file_stream<<scan<<"\t";
                        output_file_stream<<"SMILES"<<"\t";
                        output_file_stream<<"INCHI"<<"\t";
                        output_file_stream<<"INCHIAUX"<<"\t";
                        
                        if(it->first == "[M+2H]"){
                            output_file_stream<<"2"<<"\t";
                        }
                        else{
                            output_file_stream<<"1"<<"\t";
                        }
                        
                        output_file_stream<<"Positive"<<"\t";
                        output_file_stream<<"PUBMED"<<"\t";
                        output_file_stream<<"Commercial"<<"\t";
                        output_file_stream<<monoisotopic_mass<<"\t";
                        output_file_stream<<"DATACOLLECTOR"<<"\t";
                        output_file_stream<<it->first<<"\t"; //ADDUCT
                        output_file_stream<<"N/A"<<"\t";
                        output_file_stream<<"1"<<"\t";
                        output_file_stream<<"N/A"<<"\t";
                        output_file_stream<<"N/A"<<"\t";
                        output_file_stream<<"N/A"<<"\t";
                        output_file_stream<<"CASNUMBER"<<"\t";
                        output_file_stream<<"Dorrestein"<<"\t";
                        
                        
                        
                        output_file_stream<<"\n";
                        
                        if(it->first == "M+H"){
                            found_flag_m_h = true;
                        }
                        
                        
                        //Saving into max map
                        string compound_key = spectra_file_name  + "_" + compound_name + "_" + it->first;
                        stringstream ss (stringstream::in | stringstream::out);
                        ss<< "specplot "<<"--infile "<<spectra_file_name;
                        ss<<" "<<"--spectrumscan "<<scan;
                        ss<<" "<<"--outfile " <<"processing_files/spectrum_images/"<<string_replace_all(string_replace_all(pubmed_id, " ", ""), "/", "")<<"_"<<scan<<".png";
                        string plot_string = string_replace_all(ss.str(), "!", "");
                        plot_string = string_replace_all(plot_string, "\"", "");
                        plot_string = string_replace_all(plot_string, "(", "");
                        plot_string = string_replace_all(plot_string, ")", "");
                        plot_string = string_replace_all(plot_string, ";", "");
                        plot_string = string_replace_all(plot_string, "'", "");
                        
                        
                        if(compound_to_max_abundance.find(compound_key) != compound_to_max_abundance.end()){
                            //The key is in here, so check abundance
                            if(spec_precursor > compound_to_max_abundance[compound_key]){
                                compound_to_max_abundance[compound_key] = spec_precursor;
                                compound_to_plot_string[compound_key] = plot_string;
                            }
                        }
                        else{
                            compound_to_max_abundance[compound_key] = spec_precursor;
                            compound_to_plot_string[compound_key] = plot_string;
                        }
                    }
                }
            }
            
            if(found_flag_m_h == false){
                output_file_stream<<spectra_file_name<<"\t";
                output_file_stream<<"NotFound M+H"<<"\t";
                output_file_stream<<"NotFound"<<"\t";
                output_file_stream<<"NotFound"<<"\t";
                output_file_stream<<"NotFound"<<"\t";
                output_file_stream<<compound_name<<"\t";
                output_file_stream<<"NotFound"<<"\t";
                output_file_stream<<"\n";
            }
        }
        
        //Outputting all the plotting lines
        for (std::map<string, string>::iterator it=compound_to_plot_string.begin(); it!=compound_to_plot_string.end(); ++it){
            cout<<it->second<<endl;
        }
        
        return true;
  }

  bool ExecSpectraSyntheticSpectraScanExtraction::loadInputData(void){
        vector<string> required_headers;
        vector<int> required_headers_int;
        
        DelimitedTextReader::loadDelimitedFile(input_table_filename.c_str(), "\t", "#", annotation_headers_map, annotation_lines, required_headers, required_headers_int);
        
        for(int i = 0; i < annotation_lines.size(); i++){
            string output_string = "";
            for(int j = 0; j < annotation_lines[i].size(); j++){
                output_string += annotation_lines[i][j] + "\t";
            }
        }
    
        return true;
  }


  bool ExecSpectraSyntheticSpectraScanExtraction::saveOutputData(void){
    return true;
  }


  bool ExecSpectraSyntheticSpectraScanExtraction::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraSyntheticSpectraScanExtraction::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraSyntheticSpectraScanExtraction::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraSyntheticSpectraScanExtraction::merge(void){
    return true;
  }


  bool ExecSpectraSyntheticSpectraScanExtraction::validateParams(std::string & error){
    return true;
  }

}
