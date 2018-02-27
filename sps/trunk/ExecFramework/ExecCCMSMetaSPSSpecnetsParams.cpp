//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecCCMSMetaSPSSpecnetsParams.h"

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

  ExecCCMSMetaSPSSpecnetsParams::ExecCCMSMetaSPSSpecnetsParams(void)
  {
    m_name = "ExecCCMSMetaSPSSpecnetsParams";
    m_type = "ExecCCMSMetaSPSSpecnetsParams";
  }

  ExecCCMSMetaSPSSpecnetsParams::ExecCCMSMetaSPSSpecnetsParams(const ParameterList & inputParams)
  {
    m_name = "ExecCCMSMetaSPSSpecnetsParams";
    m_type = "ExecCCMSMetaSPSSpecnetsParams";
    //inputParams.print(cout);
    m_params = inputParams;

    // pk tol
    // pm tol Da
    // ppm tol Da
    // instrument type
    // merge same prec
    // fasta database
  }

  ExecCCMSMetaSPSSpecnetsParams::~ExecCCMSMetaSPSSpecnetsParams(void)
  {
  }

  ExecBase * ExecCCMSMetaSPSSpecnetsParams::clone(const ParameterList & inputParams) const
  {
    return new ExecCCMSMetaSPSSpecnetsParams(inputParams);
  }

  bool ExecCCMSMetaSPSSpecnetsParams::invoke(void)
  {
    if (!m_params.exists("Default.instrument"))
    {
      ERROR_MSG("Parameter \'Default.instrument\' not specified!!");
      return false;
    }
    string IT_ID = "ESI-ION-TRAP";
    string FT_ID = "FT-HYBRID";
    string Inst_ID = m_params.getValue("Default.instrument");
    if (Inst_ID == IT_ID)
    {
      m_params.setValue("DECONV_MS2", "0");
      m_params.setValue("INSTRUMENT_TYPE", "IT");
    }
    else if (Inst_ID == FT_ID)
    {
      m_params.setValue("DECONV_MS2", "1");
      m_params.setValue("INSTRUMENT_TYPE", "FT");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << Inst_ID << " for parameter \'Default.instrument\'!!");
      return false;
    }

    if (!m_params.exists("metasps.fragmentation"))
    {
      ERROR_MSG("Parameter \'metasps.fragmentation\' not specified!!");
      return false;
    }

    if (m_params.getValue("custom_preprocess.Enable_Clustering", "off") == "on")
    {
      m_params.setValue("CLUSTER_MIN_SIZE", "1");
    }
    else
    {
      m_params.setValue("CLUSTER_MIN_SIZE", "0");
    }

    string Mult_ID = "multiple";
    string CID_ID = "CID";
    string HCD_ID = "HCD";
    string ETD_ID = "ETD";
    string Frag_ID = m_params.getValue("metasps.fragmentation");
    if (Frag_ID == Mult_ID)
    {
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("MERGE_SAME_PREC", "1");
    }
    else if (Frag_ID == CID_ID)
    {
      m_params.setValue("ACTIVATION", "CID");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
    }
    else if (Frag_ID == HCD_ID)
    {
      m_params.setValue("ACTIVATION", "HCD");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
    }
    else if (Frag_ID == ETD_ID)
    {
      m_params.setValue("ACTIVATION", "ETD");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("MERGE_SAME_PREC", "0");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << Frag_ID << " for parameter \'Default.fragmentation\'!!");
      return false;
    }

    const string pmTolParam = "tolerance.PM_tolerance";
    const string pmTolUnitParam = "tolerance_unit.PM_unit";
    const string ionTolParam = "tolerance.Ion_tolerance";
    const string ionTolUnitParam = "tolerance_unit.Ion_unit";

    if (!m_params.exists(pmTolParam))
    {
      ERROR_MSG("Parameter \'" << pmTolParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(pmTolUnitParam))
    {
      ERROR_MSG("Parameter \'" << pmTolUnitParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(ionTolParam))
    {
      ERROR_MSG("Parameter \'" << ionTolParam << "\' not specified!!");
      return false;
    }
    if (!m_params.exists(ionTolUnitParam))
    {
      ERROR_MSG("Parameter \'" << ionTolUnitParam << "\' not specified!!");
      return false;
    }

    float pmTol = m_params.getValueFloat(pmTolParam);
    string pmTolUnits = m_params.getValue(pmTolUnitParam);
    const float ppmFac = 1000000.0;

    std::transform(pmTolUnits.begin(),
                   pmTolUnits.end(),
                   pmTolUnits.begin(),
                   ::tolower);

    if (pmTolUnits == "da")
    {
      m_params.setValue("TOLERANCE_PM", parseFloat(pmTol, 5));
    }
    else if (pmTolUnits == "ppm")
    {
      m_params.setValue("TOLERANCE_PM_PPM", parseFloat(pmTol, 5));
      m_params.setValue("TOLERANCE_PM",
                        parseFloat(((pmTol / ppmFac) * 2000.0), 5));
    }
    else
    {
      ERROR_MSG("Unknown tolerance type \'" << pmTolUnits << "\' !!");
      return false;
    }

    float ionTol = m_params.getValueFloat(ionTolParam);
    string ionTolUnits = m_params.getValue(ionTolUnitParam);

    std::transform(ionTolUnits.begin(),
                   ionTolUnits.end(),
                   ionTolUnits.begin(),
                   ::tolower);

    if (ionTolUnits == "da")
    {
      m_params.setValue("TOLERANCE_PEAK", parseFloat(ionTol, 5));
    }
    else if (ionTolUnits == "ppm")
    {
      m_params.setValue("TOLERANCE_PEAK_PPM", parseFloat(ionTol, 5));
      m_params.setValue("TOLERANCE_PEAK",
                        parseFloat(((ionTol / ppmFac) * 2000.0), 5));
    }
    else
    {
      ERROR_MSG("Unknown tolerance type \'" << ionTolUnits << "\' !!");
      return false;
    }

    if (!m_params.exists("Default.complexity"))
    {
      ERROR_MSG("Parameter \'Default.complexity\' not specified!!");
      return false;
    }

    const string Simple_ID = "purified simple mixture";
    const string Complex_ID = "complex mixture";
    string complexity_ID = m_params.getValue("Default.complexity");

    if (complexity_ID == Simple_ID)
    {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "1");
      m_params.setValue("PARTIAL_OVERLAPS", "1");
      m_params.setValue("MAX_PVALUE", "0.05");
      m_params.setValue("CLUSTALW_MINSCORE", "500");
    }
    else if (complexity_ID == Complex_ID)
    {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "2");
      m_params.setValue("PARTIAL_OVERLAPS", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("CLUSTALW_MINSCORE", "50000");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << complexity_ID << " for parameter \'Default.complexity\'!!");
      return false;
    }

    m_params.setValue("REPORT_SERVER", "http://rodney.ucsd.edu/cgi-bin/");
    m_params.setValue("REPORT_USER", "aguthals");
    m_params.setValue("REPORT_DYNAMIC", "1");
    m_params.setValue("REPORT_DIR", "report");
    m_params.setValue("GRID_NUMCPUS", "1");
    m_params.setValue("GRID_SGE_EXE_DIR", "/opt/ge2011.11/bin/linux-x64");
    m_params.setValue("MIN_OVERLAP_AREA", "0.450000");
    m_params.setValue("FILTER_TRIGS", "no");
    m_params.setValue("MAX_MOD_MASS", "100");
    m_params.setValue("MIN_MOD_MASS", "-100");
    m_params.setValue("MIN_RATIO", "0.35");
    m_params.setValue("MIN_SPECTRUM_QUALITY", "0");
    m_params.setValue("MIN_MATCHED_PEAKS", "6");
    m_params.setValue("RESOLUTION", m_params.getValue("custom_alignment.Resolution"));
    m_params.setValue("TAG_LEN", "5");
    m_params.setValue("MAX_NUM_TAGS", "400");
    m_params.setValue("MAX_NUM_MODS", "2");
    m_params.setValue("MIN_MATCHED_PEAKS_DB", "6");
    m_params.setValue("MAX_PARSIMONY", "1");
    m_params.setValue("MIN_METACONTIG_SIZE", "1");
    m_params.setValue("MIN_METACONTIG_SCORE", "3.0");
    m_params.setValue("SPSPATH_MIN_NUM_PEAKS", "5");
    m_params.setValue("SPSPATH_MIN_NUM_SPECS", "2");
    m_params.setValue("PENALTY_PTM", "-2000");
    m_params.setValue("PENALTY_SAME_VERTEX", "-1000000");
    m_params.setValue("SPEC_TYPE_MSMS", "0");

    if (m_params.getValue("custom_preprocess.Correct_PM", "off") == "on")
    {
      m_params.setValue("CORRECT_PM", "yes");
    }
    else
    {
      m_params.setValue("CORRECT_PM", "no");
    }

    if (m_params.getValue("custom_preprocess.Correct_Charge", "off") == "on")
    {
      m_params.setValue("GUESS_CHARGE", "yes");
      m_params.setValue("FIX_CHARGE_ZEROS","0");
    }
    else
    {
      m_params.setValue("GUESS_CHARGE", "no");
      m_params.setValue("FIX_CHARGE_ZEROS","1");
    }

    addParamFromProteoSAFe("DECONV_MS2",
                           "custom_preprocess.Deconvolute_MS",
                           true);
    addParamFromProteoSAFe("MIN_SPECTRUM_QUALITY",
                           "custom_preprocess.Min_Spec_Quality",
                           false);

    addParamFromProteoSAFe("PARTIAL_OVERLAPS",
                           "custom_alignment.Partial_Overlaps",
                           true);
    addParamFromProteoSAFe("MAX_PVALUE", "custom_alignment.Max_Pvalue", false);
    addParamFromProteoSAFe("MIN_MATCHED_PEAKS",
                           "custom_alignment.Min_Matched_Peaks",
                           false);
    addParamFromProteoSAFe("MIN_RATIO",
                           "custom_alignment.Min_Align_Ratio",
                           false);
    addParamFromProteoSAFe("MIN_OVERLAP_AREA",
                           "custom_alignment.Min_Ovlp_Area",
                           false);
    addParamFromProteoSAFe("MAX_MOD_MASS",
                           "custom_alignment.Max_Mod_Mass",
                           false);
    addParamFromProteoSAFe("MIN_MOD_MASS",
                           "custom_alignment.Min_Mod_Mass",
                           false);

    addParamFromProteoSAFe("SPS_MIN_EDGES_TO_COMPONENT",
                           "custom_assembly.Min_Edges_To_Component",
                           false);

    if (m_params.getValue("custom_assembly.Filter_Trig", "off") == "on")
    {
      m_params.setValue("FILTER_TRIGS", "yes");
    }
    else
    {
      m_params.setValue("FILTER_TRIGS", "no");
    }

    addParamFromProteoSAFe("SPSPATH_MIN_NUM_PEAKS",
                           "custom_assembly.Path_Min_Peaks",
                           false);
    addParamFromProteoSAFe("SPSPATH_MIN_NUM_SPECS",
                           "custom_assembly.Path_Min_Specs",
                           false);
    addParamFromProteoSAFe("PENALTY_PTM", "custom_assembly.Penalty_Ptm", false);
    addParamFromProteoSAFe("MIN_METACONTIG_SCORE",
                           "custom_assembly.Min_Meta_Score",
                           false);

    addParamFromProteoSAFe("MIN_MATCHED_PEAKS_DB",
                           "custom_homology.Min_Matched_Peaks_DB",
                           false);
    addParamFromProteoSAFe("TAG_LEN", "custom_homology.Tag_Len", false);
    addParamFromProteoSAFe("MAX_NUM_TAGS",
                           "custom_homology.Max_Num_Tags",
                           false);
    addParamFromProteoSAFe("MAX_NUM_MODS",
                           "custom_homology.Max_Num_Mods",
                           false);
    addParamFromProteoSAFe("CLUSTALW_MINSCORE",
                           "custom_homology.Clustalw_Min_Score",
                           false);

    //Reading File Names
    int mapping_count = 0;
    map<string, string> file_name_mapping;
    std::vector<std::string> search_spectra_names;
    std::vector<std::string> original_names;
    while (1)
    {
      char buf[100];
      sprintf(buf, "upload_file_mapping%i", mapping_count);

      mapping_count++;

      if (!m_params.exists(buf))
        break;

      std::string mapping = m_params.getValue(buf);
      std::string mangled_name = mapping.substr(0, mapping.find("|"));
      std::string original_name = mapping.substr(mapping.find("|") + 1);

      if (mangled_name.find("spec-") == string::npos)
        continue;

      file_name_mapping[original_name] = mangled_name;

      std::string path_to_spectra = m_params.getValue("SPECTRA_DIR", "");
      search_spectra_names.push_back(path_to_spectra + "/" + mangled_name);
      original_names.push_back(original_name);
    }

    //Constructing input spectra
    string all_input_spectra_specnetsformat = "";
    string set_filename_input = "";
    for (int i = 0; i < search_spectra_names.size(); i++)
    {
      all_input_spectra_specnetsformat += search_spectra_names[i];
      all_input_spectra_specnetsformat += ";";
      set_filename_input += original_names[i];
      set_filename_input += ";";
    }

    m_params.setValue("INPUT_SPECS_MS",
                      all_input_spectra_specnetsformat.c_str());

    m_params.setValue("SET_FILENAMES", set_filename_input.c_str());

    char buf[1024];
    getcwd(buf, 1024);
    cout << "CWD" << "\t" << buf << endl;

    m_params.print(cerr);
    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::loadInputData(void)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::saveOutputData(void)
  {
    string param_outfile = m_params.getValue("RESULTS_DIR");

    DEBUG_MSG("EXE_DIR"<<"\t"<<m_params.getValue("EXE_DIR"));

    DEBUG_MSG("SAVING\t"<<param_outfile);
    m_params.writeToFile(param_outfile);
    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::saveInputData(std::vector<std::string> & filenames)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::loadOutputData(void)
  {

    return true;
  }

  std::vector<ExecBase *> const & ExecCCMSMetaSPSSpecnetsParams::split(int numSplit)
  {

    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::merge(void)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::validateParams(std::string & error)
  {
    return true;
  }

  void ExecCCMSMetaSPSSpecnetsParams::addParamFromProteoSAFe(const string & spsParam,
                                                             const string & proteosafeParam,
                                                             const bool mapYesNoToInt)
  {
    string paramVal = m_params.getValue(proteosafeParam, "default");

    if (paramVal != "default")
    {
      if (mapYesNoToInt && paramVal == "yes")
      {
        m_params.setValue(spsParam, "1");
      }
      else if (mapYesNoToInt && paramVal == "no")
      {
        m_params.setValue(spsParam, "0");
      }
      else if (!mapYesNoToInt)
      {
        m_params.setValue(spsParam, paramVal);
      }
      else
      {
        ERROR_MSG("Invalid parameter combo \'" << proteosafeParam << "=" << paramVal);
      }
    }
  }

}
