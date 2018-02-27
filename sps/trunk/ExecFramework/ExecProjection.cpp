//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecProjection.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include "utils.h"

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecProjection::ExecProjection(void){
    m_name = "ExecProjection";
    m_type = "ExecProjection";
  }


  ExecProjection::ExecProjection(const ParameterList & inputParams){
    m_name = "ExecProjection";
    m_type = "ExecProjection";
    inputParams.print(cout);

    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");

    //Grabbing the file name of the output mgf file
    m_output_mgf_file_name = inputParams.getValue("output_mgf");

    //Grabbing thresholding filters
    m_pvalue_filter = inputParams.getValueFloat("pvalue_cutoff");
    m_envelope_filter = inputParams.getValueFloat("envelope_cutoff");


    //Pulling out parameters for creating the library
    //MZXML and annotation files come as a pair
    string input_spectrum_names = inputParams.getValue("rawnames");
    string input_psn_names = inputParams.getValue("psmnames");

    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf");

    //Splitting inputs
    stringSplit(input_spectrum_names, m_spectrum_file_names);
    stringSplit(input_psn_names, m_psm_file_names);
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);

    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");


    //Filename for targets of projections
    m_projection_targets_file_name = inputParams.getValue("targetprojectionfile");

    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
  }


  ExecProjection::~ExecProjection(void){
  }


  ExecBase * ExecProjection::clone(const ParameterList & inputParams) const{
    return new ExecProjection(inputParams);
  }

  bool ExecProjection::invoke(void){
    cout<<"Invoke"<<endl;
    vector<int> charge_filter;
    charge_filter.push_back(2);

    m_specs.createlibrary(m_envelope_filter,
                          m_pvalue_filter,
                          model,
                          ionsToExtract,
                          allIons,
                          aminoacidexclusions,
                          charge_filter,
                          true);

    cout<<"Done Creating Library"<<endl;


    for(int target_idx = 0; target_idx < m_input_projection_target_sequences.size(); target_idx++){
        Spectrum projection_spectrum;
        int projection_status = m_specs.projection(m_input_projection_target_sequences[target_idx],
                                                   model, ionsToExtract, allIons, projection_spectrum);

        if(projection_status == -1) //No projection
            continue;

        projection_spectra.push_back(projection_spectrum);
    }

    cout<<"Projections: "<<projection_spectra.size()<<endl;;


    return true;
  }

  bool ExecProjection::loadInputData(void){
    //Checking if the number of mzxml files lines up with the number of annotation files
    if(m_spectrum_file_names.size() != m_psm_file_names.size()) return false;


    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b-H2O-NH3");
    ionsToExtract.push_back("b-H2O-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("b++-H2O-H2O");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-iso");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y-H2O-NH3");
    ionsToExtract.push_back("y-H2O-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-H2O");
    ionsToExtract.push_back("y++-H2O-H2O");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-iso");
    ionsToExtract.push_back("P++");
    ionsToExtract.push_back("P++-H2O");
    ionsToExtract.push_back("P++-NH3");
    ionsToExtract.push_back("P++-H2O-H2O");



    //Looping through mzxml/annotation files to load
    for(int file_idx = 0; file_idx < m_spectrum_file_names.size(); file_idx++){
        //Determining file type by file extension, only supporting plkbin and mzxml
        string spectra_file_name = m_spectrum_file_names[file_idx];
        string annotation_file_name = m_psm_file_names[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            //m_specs.LoadSpecSet_mzxml_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            //cout<<"LOADED\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            cout<<"MZXML not Supported"<<endl;
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            m_specs.LoadSpecSet_pklbin_with_annotation(spectra_file_name.c_str(), annotation_file_name.c_str());
            cout<<"LOADED\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
        }
        else{
            cerr<<"CANNOTLOAD\t"<<spectra_file_name<<endl;
        }
    }

    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_specs.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }


    //Loading Targets
    ifstream data_input_file_stream(m_projection_targets_file_name.c_str(), ios::binary);
    while(!data_input_file_stream.eof()) {
            string target_annotation;
            data_input_file_stream>>target_annotation;
            if(target_annotation.length() > 2){
                    m_input_projection_target_sequences.push_back(target_annotation);
            }
    }

    return true;
  }


  bool ExecProjection::saveOutputData(void){
        projection_spectra.SaveSpecSet_mgf(m_output_mgf_file_name.c_str());
        return true;
  }


  bool ExecProjection::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecProjection::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecProjection::split(int numSplit){
    return m_subModules;
  }

  bool ExecProjection::merge(void){
    return true;
  }


  bool ExecProjection::validateParams(std::string & error){
    return true;
  }

  void ExecProjection::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }

}
