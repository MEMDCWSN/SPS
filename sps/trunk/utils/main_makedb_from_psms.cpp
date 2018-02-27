//
// main_makedb_from_psms - filters DB by psms proteins
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"
#include "tuple.h"

#include <algorithm>
#include <stdlib.h>

using namespace specnets;
using namespace std;
using namespace sps;


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
string reverseString(string & inputString)
{
  string reversedString = inputString;
  int lenStr = inputString.length();
  for (int j = 0; j < lenStr; j++) {
    reversedString[lenStr - j - 1] = inputString[j];
  }
  return reversedString;
}

// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_compute_variants <psm_file> <database_file> <out_file> [options]" << endl;
//  cerr << "   -msgfdb               Input PSMs are in MSGFDB format file" << endl;
//  cerr << "   -msgfplus             Input PSMs are in MSGF+ format file" << endl;
}

                            
// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 4) {
    displayUsage();
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  //listOptions.push_back(CommandLineParser::Option("addflank", "ADD_FLANKING_AAS", 0));
  //listOptions.push_back(CommandLineParser::Option("aafile", "AMINO_ACID_MASS_FILE", 1));
  
  CommandLineParser clp(argc, argv, 3, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
    cerr << parserError << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading PSM file [" << argv[1] << "]");
  if (!psmSet.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -2;
  }


  DEBUG_MSG("Loading Database file [" << argv[2] << "]");
  DB_fasta dbFasta;
  dbFasta.Load(argv[2]);

  vector<vector<string> > lines;
  DelimitedTextReader::loadDelimitedFileNoHeader(argv[2],
                            "",
                            "",
                            lines);
  DEBUG_VAR(lines.size());


  int minAnnoLength = 10;
  for (int i = 0; i < psmSet.size(); i++) {
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[i]->m_annotation, cleanAnno);
    if (cleanAnno.length() < minAnnoLength) {
      minAnnoLength = cleanAnno.length();
    }
  }
  DEBUG_VAR(minAnnoLength);

  // "Normalize" the AA's in the database
  DB_fasta dbCopy;  // Create the DB first
  dbCopy = dbFasta;   // Then copy.. if you don't.. you seg fault! (dunno why)
  dbCopy.replaceAA('L', 'I');
  dbCopy.replaceAA('K', 'Q');
  DB_index indexAnno(dbCopy, 1024, 1024, minAnnoLength);

  map<int, list<pair<int,int> > > regions;

  set<string> proteinNames;
  set<string> proteinNamesDecoy;
  //--------------------------------------------------------------------
  // Find matches in the database and put thier names in the list
  //--------------------------------------------------------------------
  for (int i = 0; i < psmSet.size(); i++) {
    bool isDecoy = psmSet[i]->m_isDecoy;

    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[i]->m_annotation, cleanAnno);
    string convertedAnno = convertAnnotation(cleanAnno);

    // Going to assume reverse decoys (or we can't do anything)
    if (isDecoy) {
      convertedAnno = reverseString(convertedAnno);
      DEBUG_MSG("Decoy PSM: " << convertedAnno);
    }

    // Find all the tag matches
    list<tuple<int, int, int> > locations;
    findMatchingLocations(dbCopy, indexAnno, convertedAnno, locations);
    list<tuple<int, int, int> >::iterator itr = locations.begin();
    list<tuple<int, int, int> >::iterator itrEnd = locations.end();
    for (; itr != itrEnd; itr++) {
        //DEBUG_MSG(itr->m0 << "  " << itr->m1 << "  " << itr->m2);
        string protein = dbCopy.getID(itr->m0);
        if (isDecoy) {
          proteinNamesDecoy.insert(protein);
          DEBUG_MSG("Inserting Decoy Protein: " << protein);
        } else {
          proteinNames.insert(protein);
          DEBUG_MSG("Inserting Protein: " << protein);
        }
        // If this was originally a decoy also save it as a decoy
        //   This way we can recover the reversed decoy later and add it instead
    }
  }
  DEBUG_VAR(proteinNames.size());
  DEBUG_VAR(proteinNamesDecoy.size());

  map<string, string> mapDB;

  string description;
  string proteinName;
  string protein;
  for (int i = 0; i < dbCopy.size(); i++) {
    string dbProteinName = dbCopy.getID(i);
    string dbProtein = dbCopy.getSequence(i);
    DEBUG_VAR(dbProteinName);
    DEBUG_VAR(dbProtein);

    if (proteinNames.find(dbProteinName) != proteinNames.end()) {
        mapDB[dbProteinName] = dbProtein;
        DEBUG_MSG("   FOUND");
        DEBUG_VAR(dbProteinName);
        DEBUG_VAR(dbProtein);
    }
    // If it was one of the decoys insert the decoy
    if (proteinNamesDecoy.find(dbProteinName) != proteinNamesDecoy.end()) {
      // We'll need to reverse this entry because it was a decoy PSM
      string decoyProtein = reverseString(dbProtein);
      string decoyProteinName = "YYY_" + dbProteinName;
      mapDB[decoyProteinName] = decoyProtein;
      DEBUG_MSG("   DECOY");
      DEBUG_VAR(decoyProteinName);
      DEBUG_VAR(decoyProtein);
    }
  }
  DEBUG_TRACE;
  
  ofstream ofs(argv[3]);
  if (!ofs) {
    ERROR_MSG("Error opening output file [ " << argv[3] << "]");
    return -1;
  }

  map<string, string>::iterator itr = mapDB.begin();
  map<string, string>::iterator itrEnd = mapDB.end();
  for ( ; itr != itrEnd; itr++) {
    ofs << ">" << itr->first << endl;
    ofs << itr->second << endl;
  }
  ofs.close();
  
  return 0;
}
