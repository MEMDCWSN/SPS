//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecCCMSMetabolomicsSpecnetsParams.h"
#include "ExecFdrPeptide.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unistd.h>


using namespace specnets;
using namespace std;


namespace specnets
{

  ExecCCMSMetabolomicsSpecnetsParams::ExecCCMSMetabolomicsSpecnetsParams(void){
    m_name = "ExecCCMSMetabolomicsSpecnetsParams";
    m_type = "ExecCCMSMetabolomicsSpecnetsParams";
  }


  ExecCCMSMetabolomicsSpecnetsParams::ExecCCMSMetabolomicsSpecnetsParams(const ParameterList & inputParams){
    m_name = "ExecCCMSMetabolomicsSpecnetsParams";
    m_type = "ExecCCMSMetabolomicsSpecnetsParams";
    //inputParams.print(cout);
    m_params = inputParams;
    
    m_params.setValue("PAIRS_MATCH_MODE", "cosine");
    m_params.setValue("PAIRS_MIN_COSINE", inputParams.getValue("PAIRS_MIN_COSINE" , "0.5"));
    m_params.setValue("MAX_MOD_MASS", "9999");
    m_params.setValue("MIN_RATIO", "0.4");
    m_params.setValue("MIN_MATCHED_PEAKS", inputParams.getValue("MIN_MATCHED_PEAKS" , "6"));
    m_params.setValue("CLUSTER_MODEL", "LTQ_TRYP");
    m_params.setValue("MIN_SPECTRUM_QUALITY", "0.0");
    m_params.setValue("CORRECT_PM", "no");
    m_params.setValue("GUESS_CHARGE", "no");
    m_params.setValue("GRID_NUMNODES", "-1");
    m_params.setValue("MIN_OVERLAP_AREA", "0.45");
    m_params.setValue("TOLERANCE_PEAK", inputParams.getValue("tolerance.Ion_tolerance" , "0.5"));
    m_params.setValue("TOLERANCE_PM", inputParams.getValue("tolerance.PM_tolerance" , "1.0"));
    m_params.setValue("RESOLUTION", "0.1");
    m_params.setValue("FILTER_TRIGS", "no");
    m_params.setValue("MIN_MOD_MASS", "-100");
    m_params.setValue("MAX_PVALUE", "1");
    m_params.setValue("PARTIAL_OVERLAPS", "0");
    m_params.setValue("FILTER_PRECURSOR_WINDOW", inputParams.getValue("FILTER_PRECURSOR_WINDOW" , "1"));
    m_params.setValue("MIN_PEAK_INT", inputParams.getValue("MIN_PEAK_INT" , "50.0"));
    m_params.setValue("FILTER_STDDEV_PEAK_INT", inputParams.getValue("FILTER_STDDEV_PEAK_INT" , "2.0"));
    
    //If we want to filter precursor window, also filter after clustering
    if(inputParams.getValueInt("WINDOW_FILTER" , 1) > 0){
        m_params.setValue("CLUST_RANK_FILTER", "6");
        m_params.setValue("RANK_FILTER", "6");
        m_params.setValue("RANK_FILTER_RADIUS", "50.0");
    }
    
    //Setting the Threaded Parallelism
    if(inputParams.getValueInt("THREADEDPARALLELISM" , 0) > 0){
        m_params.setValue("GRID_NUMNODES", inputParams.getValue("THREADEDPARALLELISM" , "2"));
        m_params.setValue("GRID_TYPE", "threaded");
    }
    
    if(inputParams.getValue("RUN_MSCLUSTER", "off") == "on")
        m_params.setValue("CLUSTER_MIN_SIZE", inputParams.getValue("CLUSTER_MIN_SIZE" , "1"));
    else
        m_params.setValue("CLUSTER_MIN_SIZE", "0");

    //Reading File Names
    int mapping_count = 0;
    map<string,string> file_name_mapping;
    std::vector<std::string> search_spectra_names;
    while(1){
        char buf[100];
        sprintf(buf, "upload_file_mapping%i", mapping_count);
        
        mapping_count++;
        
        if(!m_params.exists(buf))
            break;
        
        std::string mapping = m_params.getValue(buf);
        std::string mangled_name = mapping.substr(0, mapping.find("|"));
        std::string original_name = mapping.substr(mapping.find("|")+1);
        
        if( mangled_name.find("spec-") == string::npos &&
            mangled_name.find("specone-") == string::npos &&
            mangled_name.find("spectwo-") == string::npos &&
            mangled_name.find("specthree-") == string::npos &&
            mangled_name.find("specfour-") == string::npos &&
            mangled_name.find("specfive-") == string::npos &&
            mangled_name.find("specsix-") == string::npos) continue;
        
        file_name_mapping[original_name] = mangled_name;
        
        std::string mangled_extension = get_extension(mangled_name);
        if (mangled_extension != "mzXML" && mangled_extension != "mzxml" && mangled_extension != "mgf" && mangled_extension != "pklbin" && mangled_extension != "mzML" && mangled_extension != "gz"){
            std::cout<<"ERROR: INVALID INPUT FORMAT: "<<original_name<<std::endl;
            continue;
        }
        
        std::string path_to_spectra = m_params.getValue("SPECTRA_DIR", "");
        search_spectra_names.push_back(path_to_spectra + "/" + mangled_name);
    }
    
    //Constructing input spectra
    string all_input_spectra_specnetsformat  = "";
    for(int i = 0; i < search_spectra_names.size(); i++){
        all_input_spectra_specnetsformat += search_spectra_names[i];
        all_input_spectra_specnetsformat += ";";
    }
    
    m_params.setValue("INPUT_SPECS_MS", all_input_spectra_specnetsformat.c_str());
    

    char buf[1024];
    char* res = getcwd(buf, 1024);
    cout<<"CWD"<<"\t"<<buf<<endl;;
    
    string exe_dir = buf;
    m_params.setValue("EXE_DIR", exe_dir);
    m_params.setValue("REPORT_DIR", "report");
    
    
    
    
  }


  ExecCCMSMetabolomicsSpecnetsParams::~ExecCCMSMetabolomicsSpecnetsParams(void){
  }


  ExecBase * ExecCCMSMetabolomicsSpecnetsParams::clone(const ParameterList & inputParams) const{
    return new ExecCCMSMetabolomicsSpecnetsParams(inputParams);
  }

  bool ExecCCMSMetabolomicsSpecnetsParams::invoke(void){
    string param_outfile = m_params.getValue("RESULTS_DIR");
    
    m_params.setValue("EXE_DIR", get_only_path(m_params.getValue("MSCLUSTER_DIR")));
    
    DEBUG_MSG("PEPNOVO_DIR"<<"\t"<<m_params.getValue("PEPNOVO_DIR"));
    DEBUG_MSG("MSCLUSTER_DIR"<<"\t"<<m_params.getValue("MSCLUSTER_DIR"));
    DEBUG_MSG(get_only_path(m_params.getValue("MSCLUSTER_DIR")));
    
    
    DEBUG_MSG("SAVING\t"<<param_outfile);
    m_params.writeToFile(param_outfile);
    
    return true;
  }

  bool ExecCCMSMetabolomicsSpecnetsParams::loadInputData(void){      
    
   
    
    return true;
  }


  bool ExecCCMSMetabolomicsSpecnetsParams::saveOutputData(void){

    return true;
  }


  bool ExecCCMSMetabolomicsSpecnetsParams::saveInputData(std::vector<std::string> & filenames){

    return true;
  }

  bool ExecCCMSMetabolomicsSpecnetsParams::loadOutputData(void){

    return true;
  }

  std::vector<ExecBase *> const & ExecCCMSMetabolomicsSpecnetsParams::split(int numSplit){

    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecCCMSMetabolomicsSpecnetsParams::merge(void){
    
    return true;
  }


  bool ExecCCMSMetabolomicsSpecnetsParams::validateParams(std::string & error){
    return true;
  }
  
  
}
