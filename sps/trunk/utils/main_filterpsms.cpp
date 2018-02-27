//
// util_parsimony - Create parsimonious PSM set
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"
#include "tuple.h"

#include <stdlib.h>

const float OVERLAP_REGION = 0.5;
using namespace specnets;
using namespace std;
using namespace sps;

struct TupeCompare : public std::binary_function<tuple<int,int,int>, tuple<int,int,int>, bool>
{
    bool operator()(const tuple<int,int,int>& lhs, const tuple<int,int,int>& rhs) const
    {
      if (lhs.m0 < rhs.m0) {
        return true;
      }
      if (lhs.m0 == rhs.m0 && lhs.m1 < rhs.m1) {
        return true;
      }
      
      if (lhs.m0 == rhs.m0 && lhs.m1 == rhs.m1 && lhs.m2 < rhs.m2) {
        return true;
      }
      return false;
    }
};

  inline int convertStringToInt(const char * stringValue)
  {
    int returnValue = -1;
    returnValue = atoi(stringValue);
    return returnValue;
  }

  inline float convertStringToFloat(const char * stringValue)
  {
    float returnValue = -1.0;
    returnValue = atof(stringValue);
    return returnValue;
  }

  inline int convertStringToInt(string & stringValue)
  {
    int returnValue = -1;
    returnValue = atoi(stringValue.c_str());
    return returnValue;
  }

  inline float convertStringToFloat(string & stringValue)
  {
    float returnValue = -1.0;
    returnValue = atof(stringValue.c_str());
    return returnValue;
  }



// -------------------------------------------------------------------------
string convertAnnotation(string annotation)
{
  string convertedAnnotation; 
  for (int iChar = 0; iChar < annotation.length(); iChar++) {
    if (annotation[iChar] == 'L') {
      convertedAnnotation += 'I';
    } else if (annotation[iChar] == 'K') {
      convertedAnnotation += 'Q';
    } else {
      convertedAnnotation += annotation[iChar];
    }
  }
  return convertedAnnotation;
}

// -------------------------------------------------------------------------
float computeTotalMass(psmPtr p, AAJumps & jumps)
{
  string cleanAnno;
  PeptideSpectrumMatch::getUnmodifiedPeptide(p->m_annotation, cleanAnno);
  float totalMass = jumps.getPeptideMass(cleanAnno);
	  
  vector<float> modifications;
  vector<unsigned int> positions;
  p->getModificationsAndPositions(modifications, positions);
  
  for (int i = 0; i < modifications.size(); i++) {
    totalMass += modifications[i];
  }
  
  return totalMass;
  
}

// ------------------------------------------------------------------------
void findMatchingLocations(DB_fasta & db, 
                           DB_index & index, 
                           string annotation,
                           list<tuple<int, int, int> > & locationsReturn)
{
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(annotation, cleanAnno);
    string convertedAnno = convertAnnotation(cleanAnno);
    
    // Find all the tag matches
    list<pair<int, int> > * locations;
    index.find(convertedAnno.c_str(), &locations);
    if (locations == 0x0) {
      WARN_MSG("No locations found for [" << convertedAnno << "]");
      return;
    }
    list<pair<int, int> >::iterator itr = locations->begin();
    list<pair<int, int> >::iterator itrEnd = locations->end();
    for (; itr != itrEnd; itr++) {
        int dbIndex = itr->first;
        int start = itr->second;
        int end = start + convertedAnno.length() - 1;
        string sequence = db.getSequence(itr->first);
        if (sequence.length() < end) {
          continue;
        }
        string seqAnno = sequence.substr(start, convertedAnno.length());
        if (seqAnno != convertedAnno) {
          continue;
        }
        
        locationsReturn.push_back(make_tuple<int, int, int>(dbIndex, start, end));
    }
    
    return;
}

// -------------------------------------------------------------------------
bool overlap(int first1, int second1,
             int first2, int second2,
             float percent)
{
  float overlap1 = 0.0;
  float overlap2 = 0.0;
  if (first1 >= first2 && first1 <= second2) {
    overlap1 = (float)(second2 - first1) / (float)(second1 - first1);
    //DEBUG_VAR(overlap1);
  } else if (second1 >= first2 && second1 <= second2) {
    overlap1 = (float)(second1 - first2) / (float)(second1 - first1);
    //DEBUG_VAR(overlap1);
  }
  if (first2 >= first1 && first2 <= second1) {
    overlap2 = (float)(second1 - first2) / (float)(second2 - first2);
    //DEBUG_VAR(overlap2);
  } else if (second2 >= first1 && second2 <= second1) {
    overlap2 = (float)(second2 - first1) / (float)(second2 - first2);
    //DEBUG_VAR(overlap2);
  }
  float overlap = max(overlap1, overlap2);
  //DEBUG_VAR(overlap);
  if (overlap >= percent) {
    return true;
  }
  return false;
}                            

// -------------------------------------------------------------------------
string getSequence(DB_fasta & dbFasta, int index, int start, int end)
{
  string sequence = dbFasta.getSequence(index);
  string seqAnno = sequence.substr(start, end - start + 1);
  return seqAnno;
} 
                            
// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    cerr << "Usage: main_filter_psms psm_file_in psm_file_out [options]" << endl;
    cerr << "     -moda              Input IDs are in MODa (multiple results) style file" << endl;
    cerr << "     -modaf             Input IDs are in MODa (filtered results) style file" << endl;
    cerr << "     -msgfdb            Input PSMs are in MSGFDB style file" << endl;
    cerr << "     -msgfp             Input PSMs are in MSGF-PLUS style file" << endl;
    cerr << "     -charge C          Output only PSMs with charge = C" << endl;
    cerr << "     -tryptic T         Output only PSMs that are tryptic (t=1) or non-tryptic T = 0" << endl;
    cerr << "     -firstscan S       Output only PSMs that have scan number >= S" << endl;
    cerr << "     -lastscan S        Output only PSMs that have scan number >= S" << endl;
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("charge", "CHARGE_VALUE", 1));
  listOptions.push_back(CommandLineParser::Option("tryptic", "IS_TRYPTIC", 1));
  listOptions.push_back(CommandLineParser::Option("firstscan", "FIRST_SCAN", 1));
  listOptions.push_back(CommandLineParser::Option("lastscan", "LAST_SCAN", 1));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("msgfdb", "LOAD_MSGFDB_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msgfp", "LOAD_MSGF_PLUS_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("moda", "LOAD_MODA_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("modaf", "LOAD_MODA_STYLE2_FILE", 0));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_filter_psms psm_file_in psm_file_out [options]" << endl;
    cerr << "     -moda              Input IDs are in MODa (multiple results) style file" << endl;
    cerr << "     -modaf             Input IDs are in MODa (filtered results) style file" << endl;
    cerr << "     -msgfdb            Input PSMs are in MSGFDB style file" << endl;
    cerr << "     -msgfp             Input PSMs are in MSGF-PLUS style file" << endl;
    cerr << "     -charge c          Output only PSMs with charge = c" << endl;
    cerr << "     -tryptic t         Output only PSMs that are tryptic (t=1) or non-tryptic t = 0" << endl;
    cerr << "     -firstscan S       Output only PSMs that have scan number >= S" << endl;
    cerr << "     -lastscan S        Output only PSMs that have scan number >= S" << endl;
    cerr << "Invalid options" << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);


  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading target PSM file [" << argv[1] << "]");
  if (commandLineParams.exists("LOAD_MSGFDB_STYLE_FILE")) {
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
  } else if (commandLineParams.exists("LOAD_MSGF_PLUS_STYLE_FILE")) {
    if (!psmSet.loadMSGFPlusResultsFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
    for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
      vector<float> modifications;
      psmSet[iPsm]->getModifications(modifications);
      psmSet[iPsm]->m_numMods = modifications.size();
        
    }
  } else if (commandLineParams.exists("LOAD_MODA_STYLE_FILE")) {
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
        
        PeptideSpectrumMatch::inspectToSpecNets(newPsm->m_origAnnotation,
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
  } else if (commandLineParams.exists("LOAD_MODA_STYLE2_FILE")) {
    // Load moda intermediate file
    vector<vector<string> > lines;
    vector<string> header;
    vector<string> requiredHeader;  // empty
    vector<int> requiredHeaderIndex;  // empty
    DelimitedTextReader::loadDelimitedFile(argv[1],
                              "\t",
                              "",
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
    DEBUG_VAR(lines.size());
    for (int i = 0; i < lines.size(); i++) {
#if 0
      for (int j = 0; j < lines[i].size(); j++) {
        DEBUG_MSG(lines[i][j]);
      }
      DEBUG_VAR(lines[i].size());
#endif      
      psmPtr newPsm(new PeptideSpectrumMatch);
      //DEBUG_VAR(lines[i][0]);
      newPsm->m_spectrumFile = lines[i][0];
      //DEBUG_VAR(lines[i][1]);
      newPsm->m_dbIndex = convertStringToInt(lines[i][1]);
      
      //DEBUG_VAR(lines[i][2]);
      newPsm->m_scanNum = convertStringToInt(lines[i][2]);
      //DEBUG_VAR(lines[i][3]);
      newPsm->m_mz = convertStringToFloat(lines[i][3]);
      //DEBUG_VAR(lines[i][4]);
      newPsm->m_charge = convertStringToInt(lines[i][4]);
      
      //DEBUG_VAR(lines[i][8]);
      newPsm->m_score = convertStringToInt(lines[i][8]);
      //DEBUG_VAR(lines[i][9]);
      newPsm->m_pValue = 1.0 - convertStringToFloat(lines[i][9]);
      //DEBUG_VAR(lines[i][10]);
      newPsm->m_origAnnotation = lines[i][10];
        
      PeptideSpectrumMatch::inspectToSpecNets(newPsm->m_origAnnotation,
                                              newPsm->m_annotation);
        
      vector<float> modifications;
      newPsm->getModifications(modifications);
      newPsm->m_numMods = modifications.size();
        
      //DEBUG_VAR(lines[i][12]);
      newPsm->m_protein = lines[i][12];
      
      //DEBUG_VAR(lines[i][13]);
      int ppos = lines[i][13].find(".");
      string startEnd = lines[i][13].substr(ppos+1);
      //DEBUG_VAR(startEnd);
      int tpos = startEnd.find("~");
      //DEBUG_VAR(tpos);
      string stringStart = startEnd.substr(0,tpos);
      //DEBUG_VAR(stringStart);
      int start = convertStringToInt(stringStart);
      //DEBUG_VAR(start);
      newPsm->m_startIndex = start;

      ppos = startEnd.find(".");
      string stringEnd = startEnd.substr(tpos+1,ppos);
      //DEBUG_VAR(stringEnd);
      int end = convertStringToInt(stringEnd);
      newPsm->m_endIndex = end;
      //DEBUG_VAR(end);
      
      psmSet.push_back(newPsm);
    }                              
  } else {
    DEBUG_TRACE;
    if (!psmSet.loadFromFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
  }
  //exit(-1);  

  DEBUG_VAR(psmSet.size());

  int firstScan = -1;
  DEBUG_VAR(commandLineParams.getValue("FIRST_SCAN"));
  if (commandLineParams.exists("FIRST_SCAN")) {
    firstScan = commandLineParams.getValueInt("FIRST_SCAN");
  }
  DEBUG_VAR(firstScan);
  int lastScan = -1;
  DEBUG_VAR(commandLineParams.getValue("LAST_SCAN"));
  if (commandLineParams.exists("LAST_SCAN")) {
    lastScan = commandLineParams.getValueInt("LAST_SCAN");
  }
  DEBUG_VAR(lastScan);

  int charge = -1;
  DEBUG_VAR(commandLineParams.getValue("CHARGE_VALUE"));
  if (commandLineParams.exists("CHARGE_VALUE")) {
    charge = commandLineParams.getValueInt("CHARGE_VALUE");
  }
  DEBUG_VAR(charge);

  int isTryptic = -1;
  DEBUG_VAR(commandLineParams.getValue("IS_TRYPTIC"));
  if (commandLineParams.exists("IS_TRYPTIC")) {
    isTryptic = commandLineParams.getValueBool("IS_TRYPTIC");
  }
  DEBUG_VAR(isTryptic);

  PeptideSpectrumMatchSet	psmSetOutput;
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    if (charge != -1 && psmSet[iPsm]->m_charge != charge) {
      continue;
    }
    if (isTryptic != -1 && psmSet[iPsm]->PeptideSpectrumMatch::isTryptic() != isTryptic) {
      continue;
    }
    if (firstScan != -1 && psmSet[iPsm]->m_scanNum < firstScan) {
      continue;
    }
    if (lastScan != -1 && psmSet[iPsm]->m_scanNum > lastScan) {
      continue;
    }
    psmSetOutput.push_back(psmSet[iPsm]);
  }
  DEBUG_TRACE;

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[2];

  DEBUG_VAR(psmSetOutput.size());
  DEBUG_MSG("Saving file [" << outputFile << "]");
  if (!psmSetOutput.saveToFile(outputFile.c_str(), true, false)) {
    ERROR_MSG("Error saving to PSM file [ " << outputFile << "]");
    return -1;
  }
  
  
  return 0;
}

