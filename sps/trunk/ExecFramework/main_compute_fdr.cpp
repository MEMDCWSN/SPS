//
// main_compute_fdr - Computes PSMs using FDR from a set of input PSMs
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "DelimitedTextReader.h"
#include "ExecFdrPeptide.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

#include <stdlib.h>

using namespace specnets;
using namespace std;
using namespace sps;

// -------------------------------------------------------------------------
inline int convertStringToInt(string & stringValue)
{
  int returnValue = -1;
  returnValue = atoi(stringValue.c_str());
  return returnValue;
}

// -------------------------------------------------------------------------
inline float convertStringToFloat(string & stringValue)
{
  float returnValue = -1.0;
  returnValue = atof(stringValue.c_str());
  return returnValue;
}

// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_compute_fdr <psm_file> <out_psm_file> [options]" << endl;
  cerr << "   -decoy filename  Input decoy psm file for 2 file [target is psm_file]" << endl;
  cerr << "   -length          Minimum length for PSM" << endl;
  cerr << "   -moda            Input IDs are in MODa (multiple results) style file" << endl;
  cerr << "   -msgfdb          Input PSMs are in MSGFDB style file" << endl;
  cerr << "   -msgfplus        Input PSMs are in MSGF+ format file" << endl;
  cerr << "   -msplit          Input IDs are in MSplit format file" << endl;
  cerr << "   -msplit2         Input IDs are in MSplit format file (use 2nd IDs)" << endl;
  cerr << "   -peptide         Peptide level FDR" << endl;
  cerr << "   -qscore          Use Qscore to compute FDR" << endl;
  cerr << "   -pmtol N         Parent mass tolerance for variant group calculation" << endl;
  cerr << "   -proteinfilter   Use protein level FDR to filter PSMs by protein" << endl;
  cerr << "   -removedecoys    Remove decoy PSMs from the output" << endl;
  cerr << "   -threshold N     FDR cutoff threshold (default is 0.1)" << endl;
  cerr << "   -type T          Type of FDR: concatenated, merged, split, split_mods, split_charge, split_tryp, separate)" << endl;
  cerr << "                         (default is concatenated)" << endl;
  cerr << "   -keepall         Keep all PSMs" << endl;
  cerr << "   -usepvalue       Use pvalue for sorting instead of score" << endl;
  cerr << "   -usepvaluerep    Use pvalue for determination of better PSM instead of score" << endl;
  cerr << "   -variant         Use variant groups for clustering instead of scan numbers" << endl;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    displayUsage();
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("decoy", "INPUT_DECOY_PSMS", 1));
  listOptions.push_back(CommandLineParser::Option("length", "MINIMUM_PSM_LENGTH", 1));
  listOptions.push_back(CommandLineParser::Option("threshold", "PEPTIDE_FDR_CUTOFF", 1));
  listOptions.push_back(CommandLineParser::Option("type", "TDA_TYPE", 1));
  listOptions.push_back(CommandLineParser::Option("msgfdb", "LOAD_MSGFDB_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msgfplus", "LOAD_MSGFPLUS_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("moda", "LOAD_MODA_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msplit", "LOAD_MSPLIT_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msplit2", "LOAD_MSPLIT_STYLE_FILE2", 0));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("pmtol", "PARENT_MASS_TOLERANCE", 1));
  listOptions.push_back(CommandLineParser::Option("qscore", "QSCORE_FDR", 0));
  listOptions.push_back(CommandLineParser::Option("removedecoys", "REMOVE_DECOYS_IN_OUTPUT", 0));
  listOptions.push_back(CommandLineParser::Option("usepvalue", "USE_PVALUE_SORT", 0));
  listOptions.push_back(CommandLineParser::Option("usepvaluerep", "USE_PVALUE_REPLACE", 0));
  listOptions.push_back(CommandLineParser::Option("variant", "VARIANT_FDR", 0));
  listOptions.push_back(CommandLineParser::Option("peptide", "PEPTIDE_FDR", 0));
  listOptions.push_back(CommandLineParser::Option("keepall", "KEEP_ALL_PSMS", 0));
  listOptions.push_back(CommandLineParser::Option("proteinfilter", "FILTER_BY_PROTEIN_FDR", 0));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
    cerr << parserError << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  if (commandLineParams.exists("QSCORE_FDR") &&
      commandLineParams.exists("USE_PVALUE_COMPARE")) {
    cerr << "-qscore flag incompatible with -usepvalue" << endl;
    displayUsage();
	  return -1;
  }

  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading target PSM file [" << argv[1] << "]");
  if (commandLineParams.exists("LOAD_MSGFDB_STYLE_FILE")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSGFDB");
    if (!psmSet.loadMSGFDBResultsFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
    // This should be looked at, but MScore and p-value are reversed here?
    for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
      float temp = psmSet[iPsm]->m_pValue;
      psmSet[iPsm]->m_pValue = psmSet[iPsm]->m_score;
      psmSet[iPsm]->m_score = temp;

      vector<float> modifications;
      psmSet[iPsm]->getModifications(modifications);
      psmSet[iPsm]->m_numMods = modifications.size();
        
    }
  } else if (commandLineParams.exists("LOAD_MSGFPLUS_STYLE_FILE")) {
      DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSGF+");
      if (!psmSet.loadMSGFPlusResultsFile(argv[1])) {
        ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
        return -2;
      }
  } else if (commandLineParams.exists("LOAD_MODA_STYLE_FILE")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MODA");
    // Load moda intermediate file
    vector<vector<string> > lines;
    DelimitedTextReader::loadDelimitedFileNoHeader(argv[1],
                              "\t",
                              "",
                              lines);
    DEBUG_VAR(lines.size());
    for (int i = 0; i < lines.size(); i++) {
#if 0
      for (int j = 0; j < lines[i].size(); j++) {
        DEBUG_MSG(lines[i][j]);
      }
      DEBUG_VAR(lines[i].size());
#endif      
      if (lines[i][0].find(">>") != string::npos) {
        psmPtr newPsm(new PeptideSpectrumMatch);
        //DEBUG_VAR(lines[i][0]);
        newPsm->m_spectrumFile = lines[i][0].substr(2);
        //DEBUG_VAR(lines[i][1]);
        newPsm->m_scanNum = convertStringToInt(lines[i][2]);
        newPsm->m_charge = convertStringToInt(lines[i][4]);
        i++;
#if 0
        DEBUG_VAR(lines[i].size());
        for (int j = 0; j < lines.size(); j++) {
          DEBUG_MSG(lines[i][j]);
        }
#endif      
        //DEBUG_VAR(lines[i][2]);
        newPsm->m_score = convertStringToInt(lines[i][2]);
        //DEBUG_VAR(lines[i][3]);
        newPsm->m_pValue = 1.0 - convertStringToFloat(lines[i][3]);
        //DEBUG_VAR(lines[i][4]);
        newPsm->m_origAnnotation = lines[i][4];
        
        string strippedAnno;
        newPsm->stripPrecedingAndFollowing(newPsm->m_origAnnotation, strippedAnno);
        PeptideSpectrumMatch::inspectToSpecNets(strippedAnno,
                                                newPsm->m_annotation);
        
        vector<float> modifications;
        newPsm->getModifications(modifications);
        newPsm->m_numMods = modifications.size();
        
        //DEBUG_VAR(lines[i][5]);
        newPsm->m_protein = lines[i][5];
        //DEBUG_VAR(lines[i][6]);
        string startEnd = lines[i][6];
        //DEBUG_VAR(startEnd);
        int tpos = startEnd.find("~");
        //DEBUG_VAR(tpos);
        string stringStart = startEnd.substr(0,tpos);
        //DEBUG_VAR(stringStart);
        int start = convertStringToInt(stringStart);
        newPsm->m_startIndex = start;
        string stringEnd = startEnd.substr(tpos+1);
        //DEBUG_VAR(stringEnd);
        int end = convertStringToInt(stringEnd);
        newPsm->m_endIndex = end;
        psmSet.push_back(newPsm);
      }
    }                              
  } else if (commandLineParams.exists("LOAD_MSPLIT_STYLE_FILE")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSPLIT");
    if (!psmSet.loadMsplitResultsFile(argv[1])) {
      ERROR_MSG("Loading PSM file [" << argv[1] << "]");
      return -2;
    }
  } else if (commandLineParams.exists("LOAD_MSPLIT_STYLE_FILE2")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSPLIT2");
    if (!psmSet.loadMsplitResultsFile2(argv[1])) {
      ERROR_MSG("Loading PSM file [" << argv[1] << "]");
      return -2;
    }
  } else {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "]");
    if (!psmSet.loadFromFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
  }
  //exit(-1);  
  //psmSet.saveToFile("adebug.txt", true, true);

  DEBUG_VAR(psmSet.size());
  DEBUG_VAR(commandLineParams.exists("QSCORE_FDR"));
  
  // Make sure the isDecoy fields are set on any decoy PSMs
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    if (psmSet[iPsm]->m_protein.find("XXX") != string::npos ||
        psmSet[iPsm]->m_protein.find("REV") != string::npos ||
        psmSet[iPsm]->m_protein.find("DECOY") != string::npos) {
      psmSet[iPsm]->m_isDecoy = true;
    }

    // While we are at it.. lets set the scores to the FDR values
    //    if we are doing a "QScore" FDR
    // Set to (1 - FDR) So that 0 becomes 1 and goes down from there
    //     higher scores are better
    if (commandLineParams.exists("QSCORE_FDR")) {
      psmSet[iPsm]->m_score = 1 - psmSet[iPsm]->m_fdr;
    }
  }
  psmSet.saveToFile("bdebug.txt", true, true);

  PeptideSpectrumMatchSet	psmSetDecoy;
  if (commandLineParams.exists("INPUT_DECOY_PSMS")) {
    string decoyFile = commandLineParams.getValue("INPUT_DECOY_PSMS");
    DEBUG_MSG("Loading decoy PSM file [" << decoyFile << "]");
    if (!psmSetDecoy.loadFromFile(decoyFile.c_str())) {
      ERROR_MSG("Loading decoy PSM file [" << decoyFile << "]");
      return -2;
    }
    
    // Make sure the isDecoy fields are set on the decoy PSMs
    for (int iDecoy = 0; iDecoy < psmSetDecoy.size(); iDecoy++) {
      psmSetDecoy[iDecoy]->m_isDecoy = true;
    }
    psmSetDecoy.saveToFile("debug_psm_decoy.txt", true, true);
  }
  DEBUG_VAR(psmSetDecoy.size());

  
  int minLength = commandLineParams.getValueInt("MINIMUM_PSM_LENGTH", 0);
  DEBUG_VAR(minLength);
    
  PeptideSpectrumMatchSet	psmSetFiltered;
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    string annoClean;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[iPsm]->m_annotation, annoClean);
    //DEBUG_VAR(annoClean);
    if (annoClean.size() >= minLength) {
      psmSetFiltered.push_back(psmSet[iPsm]);
    }
  }
  DEBUG_VAR(psmSetFiltered.size());

  PeptideSpectrumMatchSet	psmSetDecoyFiltered;
  for (int iPsm = 0; iPsm < psmSetDecoy.size(); iPsm++) {
    string annoClean;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetDecoy[iPsm]->m_annotation, annoClean);
    if (annoClean.length() >= minLength) {
      psmSetDecoyFiltered.push_back(psmSetDecoy[iPsm]);
    }
  }
  DEBUG_VAR(psmSetDecoyFiltered.size());


  if (!commandLineParams.exists("TDA_TYPE")) {
    commandLineParams.setValue("TDA_TYPE", "concatenated");
  }

  if (commandLineParams.exists("USE_PVALUE_SORT")) {
    commandLineParams.setValue("PEPTIDE_FDR_USE_PVALUE_SORT", "1");
  }
  
  if (commandLineParams.exists("USE_PVALUE_REPLACE")) {
    commandLineParams.setValue("PEPTIDE_FDR_USE_PVALUE_REPLACE", "1");
  }

  if (!commandLineParams.exists("PEPTIDE_FDR_CUTOFF")) {
    commandLineParams.setValue("PEPTIDE_FDR_CUTOFF", "0.01");
  }
    
  if (commandLineParams.exists("REMOVE_DECOYS_IN_OUTPUT")) {
    commandLineParams.setValue("PEPTIDE_FDR_REMOVE_DECOYS", "1");
  }

  commandLineParams.setValue("FDR_TYPE", "1");
  if (commandLineParams.exists("KEEP_ALL_PSMS")) {
    commandLineParams.setValue("FDR_TYPE", "0");
  } else if (commandLineParams.exists("PEPTIDE_FDR")) {
    commandLineParams.setValue("FDR_TYPE", "2");
  } else if (commandLineParams.exists("VARIANT_FDR")) {
    commandLineParams.setValue("FDR_TYPE", "3");
  }

  string parentMassTol = commandLineParams.getValue("PARENT_MASS_TOLERANCE", "1.0");
  DEBUG_VAR(parentMassTol);
  commandLineParams.setValue("TOLERANCE_PM", parentMassTol);

  PeptideSpectrumMatchSet psmSetFdr;
  if (psmSetDecoyFiltered.size() != 0) {
    DEBUG_TRACE;
    ExecFdrPeptide moduleFdrPeptide(commandLineParams,
                                    &psmSetFiltered,
                                    &psmSetDecoyFiltered,
                                    &psmSetFdr);

    bool returnStatus = moduleFdrPeptide.invoke();
    if (!returnStatus) {
      ERROR_MSG("Problem invoking ExecFdrPeptide!");
      return -1;
    }
  } else {
    DEBUG_TRACE;
    ExecFdrPeptide moduleFdrPeptide(commandLineParams,
                                    &psmSetFiltered,
                                    &psmSetFdr);

    bool returnStatus = moduleFdrPeptide.invoke();
    if (!returnStatus) {
      ERROR_MSG("Problem invoking ExecFdrPeptide!");
      return -1;
    }
  }

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[2];
  DEBUG_MSG("Saving file [" << outputFile << "]");
  psmSetFdr.saveToFile(outputFile.c_str(), true, true);

  return 0;
}

