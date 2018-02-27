//
// main_compute_variants - computes variants and peptide variant regions
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"
#include "tuple.h"

#include <algorithm>
#include <stdlib.h>

const float OVERLAP_REGION = 0.5;
using namespace specnets;
using namespace std;
using namespace sps;

// -------------------------------------------------------------------------
struct sortPsmsByMass
{
    bool operator () (const psmPtr & lhs, const psmPtr & rhs)
    {
          return lhs->m_exactmass <= rhs->m_exactmass;
    }
};

// -------------------------------------------------------------------------
float closestMassInVector(vector<float> & massVec, float massToFind)
{
  float closest = -10000000;  // invalid group
  float closestMass = -1.0; // invalid mass
  std::vector<float>::iterator lower;
  std::vector<float>::iterator upper;

  lower = std::lower_bound(massVec.begin(), massVec.end(), massToFind);
  upper =  std::upper_bound(massVec.begin(), massVec.end(), massToFind);

  if (lower != massVec.end()) {
    closestMass = *lower;
    //DEBUG_VAR(closestMass);
    closest = closestMass;
  }
  if (upper != massVec.end()) {
    //DEBUG_VAR(*upper);
    if (abs(*upper - massToFind) < abs(closestMass - massToFind)) {
      closest = *upper;
    }
  }

  //DEBUG_VAR(closest);
  return closest;
}

// -------------------------------------------------------------------------
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
void displayUsage(void)
{
  cerr << "Usage: main_compute_variants <psm_file> <database_file> <out_file> <out_pvr_file> <out_pep_var> [options]" << endl;
  cerr << "   -addflank             Will add flanking AAs to the PSMs (database must be specified)" << endl;
  cerr << "   -aafile <filename>    Amino acid mass file (if none specified default AA masses will be used)" << endl;
  cerr << "   -decoy                Specifically a decoy file (normally decoys are ignored)" << endl;
  cerr << "   -msgfdb               Input PSMs are in MSGFDB format file" << endl;
  cerr << "   -msgfplus             Input PSMs are in MSGF+ format file" << endl;
  cerr << "   -moda                 Input IDs are in MODa (multiple results) style file" << endl;
  cerr << "   -msplit               Input IDs are in MSplit format file" << endl;
  cerr << "   -msplit2              Input IDs are in MSplit format file (use 2nd IDs)" << endl;
  cerr << "   -generic <col>        Input IDs are in Generic TSV format file" << endl;
  cerr << "                            parameter is column number (0 index) for annotation string" << endl;
  cerr << "                            string should be in standard specnets format, e.g. mods are (A,16)" << endl;
  cerr << "   -protein <col>        Used with -generic to specify the (0 index) column for the protein" << endl;
  cerr << "   -nodupe               Remove duplicate variants in output" << endl;
  cerr << "   -mergedupe            Merge duplicate variants in output" << endl;
  cerr << "   -noreg                Don't compute PVR region info" << endl;
  cerr << "   -novar                Don't compute variant info (assumed to already be present in file)" << endl;
  cerr << "                         database file will not be used but a dummy name must be given" << endl;
  cerr << "   -outdir <directory>   Output directory for output files (if none specified current directory will be used)" << endl;
  cerr << "   -pmtol N              Parent mass tolerance for variant group calculation" << endl;
  cerr << "   -pmtolunits           Parent mass tolerance units (either 'da' or 'ppm')" << endl;
  cerr << "   -showallregions       Show all PVR regions even if they are duplicates" << endl;
  cerr << "   -variants <filename>  Use previously computed variants to filter PSMs" << endl;
  cerr << "   -useregions           Use existing regions" << endl;
  cerr << "   -removedecoys         Remove decoys from the output" << endl;
}

                            
// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 6) {
    displayUsage();
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("addflank", "ADD_FLANKING_AAS", 0));
  listOptions.push_back(CommandLineParser::Option("aafile", "AMINO_ACID_MASS_FILE", 1));
  listOptions.push_back(CommandLineParser::Option("decoy", "COMPUTING_DECOYS", 0));
  listOptions.push_back(CommandLineParser::Option("msgfdb", "LOAD_MSGFDB_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msgfplus", "LOAD_MSGFPLUS_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("moda", "LOAD_MODA_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msplit", "LOAD_MSPLIT_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msplit2", "LOAD_MSPLIT_STYLE_FILE2", 0));
  listOptions.push_back(CommandLineParser::Option("generic", "LOAD_GENERIC_TSV_FILE", 1));
  listOptions.push_back(CommandLineParser::Option("nodupe", "NO_DUPLICATES", 0));
  listOptions.push_back(CommandLineParser::Option("mergedupe", "MERGE_DUPLICATES", 1));
  listOptions.push_back(CommandLineParser::Option("noreg", "DONT_COMPUTE_REGIONS", 0));
  listOptions.push_back(CommandLineParser::Option("novar", "DONT_COMPUTE_VARIANTS", 0));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));  
  listOptions.push_back(CommandLineParser::Option("protein", "GENERIC_TSV_PROTEIN_COLUMN", 1));
  listOptions.push_back(CommandLineParser::Option("pmtol", "PARENT_MASS_TOLERANCE", 1));
  listOptions.push_back(CommandLineParser::Option("pmtolunits", "PARENT_MASS_TOLERANCE_UNITS", 1));
  listOptions.push_back(CommandLineParser::Option("variants", "USE_EXISTING_VARIANTS", 1));
  listOptions.push_back(CommandLineParser::Option("useregions", "USE_EXISTING_REGIONS", 0));
  listOptions.push_back(CommandLineParser::Option("removedecoys", "REMOVE DECOYS", 0));
  listOptions.push_back(CommandLineParser::Option("showallregions", "SHOW_ALL_REGIONS", 0));
  
  CommandLineParser clp(argc, argv, 5, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
    cerr << parserError << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  float parentMassTol = commandLineParams.getValueFloat("PARENT_MASS_TOLERANCE", 1.0);
  string pmtolunits = commandLineParams.getValue("PARENT_MASS_TOLERANCE_UNITS", "da");
  std::transform(pmtolunits.begin(),
                 pmtolunits.end(),
                 pmtolunits.begin(),
                 ::tolower);

  if (pmtolunits == "ppm") {
    parentMassTol = parentMassTol * 1.5 / 1000.0;
  } else if (pmtolunits != "da") {
    cerr << "Parent mass tolerance units must be 'da' or 'ppm'" << endl;
    displayUsage();
    return -1;
  }
  DEBUG_VAR(parentMassTol);
  
  // Scale up parent mass tolerance so we can "round" to it
  if (parentMassTol == 0.0) {
	ERROR_MSG("Parent Mass Tolerance not set!");
	return -1;
  }
  DEBUG_VAR(parentMassTol);

  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading PSM file [" << argv[1] << "]");
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
  } else if (commandLineParams.exists("LOAD_MSPLIT_STYLE_FILE")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSPLIT");
    if (!psmSet.loadMsplitResultsFile(argv[1])) {
      ERROR_MSG("Loading PSM file [" << argv[1] << "]");
      return -2;
    }
  } else if (commandLineParams.exists("LOAD_MSPLIT_STYLE_FILE2")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as MSPLIT2");
    if (!psmSet.loadMsplitResultsFile(argv[1])) {
      ERROR_MSG("Loading PSM file [" << argv[1] << "]");
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
      if (lines[i][0].find(">>") != string::npos) {
        psmPtr newPsm(new PeptideSpectrumMatch);
        //DEBUG_VAR(lines[i][0]);
        newPsm->m_spectrumFile = lines[i][0].substr(2);
        //DEBUG_VAR(lines[i][1]);
        newPsm->m_scanNum = convertStringToInt(lines[i][2]);
        newPsm->m_charge = convertStringToInt(lines[i][4]);
        i++;

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
  } else if (commandLineParams.exists("LOAD_GENERIC_TSV_FILE")) {
    DEBUG_MSG("Loading PSM file [" << argv[1] << "] as GENERIC");
    int annotationIndex = commandLineParams.getValueInt("LOAD_GENERIC_TSV_FILE", 0);
    string filename;
    vector<string> header;
    vector<vector<string> > lines;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;
    DelimitedTextReader::loadDelimitedFile(argv[1],
                              "\t",
                              "",
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);

    DEBUG_VAR(lines.size());
    if (lines.size() == 0) {
      ERROR_MSG("Input file [" << argv[1] << "] is empty");
      return -3;
    }
    if (annotationIndex >= lines[0].size()) {
      ERROR_MSG("Invalid peptide index [" << annotationIndex << "]");
      ERROR_MSG("Input file has [" << lines[0].size() << "] columns");
      return -4;
    }
    
    int proteinIndex = commandLineParams.getValueInt("GENERIC_TSV_PROTEIN_COLUMN", -1);
    if (proteinIndex >= lines[0].size()) {
      ERROR_MSG("Invalid protein index [" << proteinIndex << "]");
      ERROR_MSG("Input file has [" << lines[0].size() << "] columns");
      return -4;
    }
    
    for (int i = 0; i < lines.size(); i++) {
      psmPtr newPsm(new PeptideSpectrumMatch);
      //DEBUG_VAR(lines[i][0]);
      newPsm->m_origAnnotation = lines[i][annotationIndex];
      newPsm->stripPrecedingAndFollowing(newPsm->m_origAnnotation,
                                         newPsm->m_annotation);
      
#if 0
      DEBUG_VAR(lines[i].size());
      for (int j = 0; j < lines.size(); j++) {
        DEBUG_MSG(lines[i][j]);
      }
#endif      
      vector<float> modifications;
      newPsm->getModifications(modifications);
      newPsm->m_numMods = modifications.size();
      
      //DEBUG_VAR(lines[i][5]);
      if (proteinIndex != -1) {
        newPsm->m_protein = lines[i][proteinIndex];
      }
      psmSet.push_back(newPsm);
    }                              
    
  } else {
    DEBUG_TRACE;
    if (!psmSet.loadFromFile(argv[1])) {
      ERROR_MSG("Loading PSM file [" << argv[1] << "]");
      return -2;
    }
  }

  string aaFile = commandLineParams.getValue("AMINO_ACID_MASS_FILE");
  DEBUG_VAR(aaFile);

  
  bool doingDecoys = commandLineParams.exists("COMPUTING_DECOYS");
  DEBUG_VAR(doingDecoys);

  bool noVariants = commandLineParams.exists("DONT_COMPUTE_VARIANTS");
  DEBUG_VAR(noVariants);

  // Load amino acid masses
  DEBUG_MSG("Loading Database ...");
  AAJumps jumps(1);
  if (!aaFile.empty()) {
    jumps.loadJumps(aaFile.c_str());
  }

  // First thing I'm going to do is set decoy flags
  //   DECOYS DO NOT WORK IN THIS UTILITY (for... reasons...)
  // Make sure the isDecoy fields are set on any decoy PSMs
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    if (psmSet[iPsm]->m_protein.find("XXX") != string::npos ||
        psmSet[iPsm]->m_protein.find("REV") != string::npos ||
        psmSet[iPsm]->m_protein.find("DECOY") != string::npos) {
      psmSet[iPsm]->m_isDecoy = true;
    }
  }  
  
  bool useExistingVariants = commandLineParams.exists("USE_EXISTING_VARIANTS");
  bool useExistingRegions = commandLineParams.exists("USE_EXISTING_REGIONS");

  bool removeDecoys = commandLineParams.exists("REMOVE DECOYS");

  map<pair<string, float>, int> variantGroups;

  // Put all variant PSMs into a set with their common un-moded peptide strings
  map<string, set<psmPtr, sortPsmsByMass> > mapVariantToPsmSet;

  //--------------------------------------------------------------------
  // Load existing variants (if using them)
  //--------------------------------------------------------------------
  PeptideSpectrumMatchSet psmExistingVariants;
  map<int, psmPtr> mapExistingVariants;
  vector<float> variantMassVector;
  DEBUG_VAR(useExistingVariants);
  if (useExistingVariants) {
    string filename = commandLineParams.getValue("USE_EXISTING_VARIANTS");
    DEBUG_VAR(filename);
    if (!psmExistingVariants.loadFromFile(filename.c_str())) {
      ERROR_MSG("Loading existing variant PSM file [" << filename << "]");
      return -2;
    }

    // Set up map, vector, etc for later lookups
    for (int i = 0; i < psmExistingVariants.size(); i++) {
      mapExistingVariants[psmExistingVariants[i]->m_variantGroup] = psmExistingVariants[i];
      variantMassVector.push_back(psmExistingVariants[i]->m_exactmass);
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmExistingVariants[i]->m_annotation, cleanAnno);

      mapVariantToPsmSet[cleanAnno].insert(psmExistingVariants[i]);

      // Store the variant
      pair<string, float> variantPair =
        make_pair<string, float>(cleanAnno, psmExistingVariants[i]->m_exactmass);
      variantGroups[variantPair] = psmExistingVariants[i]->m_variantGroup;
      //DEBUG_MSG(cleanAnno << "  " <<  psmExistingVariants[i]->m_exactmass << "  " << psmExistingVariants[i]->m_variantGroup);
    }

    // sort the vector of variant masses so we can find quickly
    sort(variantMassVector.begin(), variantMassVector.end());
  }

#if 0
  for (int i = 0; i < variantMassVector.size(); i++) {
    DEBUG_VAR(variantMassVector[i]);
  }
#endif

  DEBUG_VAR(mapExistingVariants.size());
  DEBUG_VAR(variantMassVector.size());

  PeptideSpectrumMatchSet psmSetStarting;
  if (noVariants) {

    for (int i = 0; i < psmSet.size(); i++) {
      // Set up map, vector, etc for later lookups
      mapExistingVariants[psmSet[i]->m_variantGroup] = psmSet[i];
      variantMassVector.push_back(psmSet[i]->m_exactmass);
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[i]->m_annotation, cleanAnno);
      mapVariantToPsmSet[cleanAnno].insert(psmSet[i]);
      // Store the variant
      pair<string, float> variantPair =
        make_pair<string, float>(cleanAnno, psmSet[i]->m_exactmass);
      variantGroups[variantPair] = psmSet[i]->m_variantGroup;
      //DEBUG_MSG(cleanAnno << "  " <<  psmExistingVariants[i]->m_exactmass << "  " << psmExistingVariants[i]->m_variantGroup);
      psmSetStarting.push_back(psmSet[i]);
    }

  } else if (useExistingVariants) {

    //--------------------------------------------------------------------
    // Use the existing variants
    //--------------------------------------------------------------------
    for (int i = 0; i < psmSet.size(); i++) {
      if (psmSet[i]->m_isDecoy && !doingDecoys) {
        continue;
      }
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[i]->m_annotation, cleanAnno);
      //DEBUG_VAR(cleanAnno);
      float totalMass = computeTotalMass(psmSet[i], jumps);
      //DEBUG_VAR(totalMass);
      psmSet[i]->m_exactmass = totalMass;


      set<psmPtr,sortPsmsByMass>::iterator itrs = mapVariantToPsmSet[cleanAnno].begin();
      set<psmPtr,sortPsmsByMass>::iterator itrsEnd = mapVariantToPsmSet[cleanAnno].end();
      bool found = false;
      for ( ; itrs != itrsEnd; itrs++) {
        //DEBUG_VAR((*itrs)->m_exactmass);
        if (abs((*itrs)->m_exactmass - totalMass) < parentMassTol) {

          //DEBUG_MSG("FOUND");
          psmSet[i]->m_variantGroup = (*itrs)->m_variantGroup;
          psmSet[i]->m_exactmass = (*itrs)->m_exactmass; // set exact mass to the variant one
          if (useExistingRegions) {
            psmPtr exsitingPsm = mapExistingVariants[psmSet[i]->m_variantGroup];
            psmSet[i]->m_peptideRegionGroup = exsitingPsm->m_peptideRegionGroup;
            psmSet[i]->m_peptideRegion = exsitingPsm->m_peptideRegion;
            psmSet[i]->m_peptideRegionLength = exsitingPsm->m_peptideRegionLength;
          }
          found = true;
          break;
        }
      }

      if (found) {
        psmSetStarting.push_back(psmSet[i]);
      }
    }
  } else {
    //--------------------------------------------------------------------
	  // Compute all the variants
    //--------------------------------------------------------------------

    // Put all PSMs into a set with their common un-moded peptide strings
    map<string, set<psmPtr, sortPsmsByMass> > mapPeptideToPsmSet;
    for (int i = 0; i < psmSet.size(); i++) {
      if (psmSet[i]->m_isDecoy && !doingDecoys) {
        psmSetStarting.push_back(psmSet[i]);
        continue;
      }

      //DEBUG_VAR(psmSet[i]->m_annotation);
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSet[i]->m_annotation, cleanAnno);
      //DEBUG_VAR(cleanAnno);
      float totalMass = computeTotalMass(psmSet[i], jumps);
      psmSet[i]->m_exactmass = totalMass;

      mapPeptideToPsmSet[cleanAnno].insert(psmSet[i]);
    }
    DEBUG_VAR(mapPeptideToPsmSet.size());

    int currentVariantGroup = 0;

    map<int, float> mapVariantToCount;
    map<int, float> mapVariantToMass;

    // Compute variant groups by closeness within the peptide set
    map<string, set<psmPtr, sortPsmsByMass> >::iterator itrm = mapPeptideToPsmSet.begin();
    map<string, set<psmPtr, sortPsmsByMass> >::iterator itrmEnd = mapPeptideToPsmSet.end();
    for ( ; itrm != itrmEnd; itrm++) {
      //DEBUG_VAR(itrm->first);
      set<psmPtr,sortPsmsByMass> newPsmSet;
      set<psmPtr,sortPsmsByMass>::iterator itrs = itrm->second.begin();
      set<psmPtr,sortPsmsByMass>::iterator itrsEnd = itrm->second.end();
      float startingMass = 0.0;
      //DEBUG_VAR(startingMass);
      for ( ; itrs != itrsEnd; itrs++) {
        float totalMass = (*itrs)->m_exactmass;
        //DEBUG_VAR(totalMass);
        if (abs(totalMass - startingMass) <= parentMassTol) {
          (*itrs)->m_variantGroup = currentVariantGroup;
          //DEBUG_VAR(currentVariantGroup);
        } else {
          currentVariantGroup++;
          (*itrs)->m_variantGroup = currentVariantGroup;
          //DEBUG_VAR(currentVariantGroup);
          startingMass = totalMass;
          //DEBUG_VAR(startingMass);
        }
        // Want to set the variant masses to average for variant group
        mapVariantToCount[currentVariantGroup]++;
        mapVariantToMass[currentVariantGroup] += totalMass;

        (*itrs)->m_peptideRegionGroup = "-1";
        (*itrs)->m_peptideRegion = "N/A";
        (*itrs)->m_peptideRegionLength = -1;

        newPsmSet.insert(*itrs);
      }

      // Change set to the new variant set
      mapPeptideToPsmSet[itrm->first] = newPsmSet;
    }
    DEBUG_VAR(mapPeptideToPsmSet.size());

    map<int, float>::iterator itrVMap = mapVariantToMass.begin();
    map<int, float>::iterator itrVMapEnd = mapVariantToMass.end();
    for ( ; itrVMap != itrVMapEnd; itrVMap++) {
      if (mapVariantToCount[itrVMap->first] == 0.0) {
        WARN_MSG("Variant #" << itrVMap->first << "has 0 mass!");
        continue;
      }
      //DEBUG_MSG(itrVMap->first << "  " <<
      //          itrVMap->second << "  " <<
      //          mapVariantToCount[itrVMap->first]);
      mapVariantToMass[itrVMap->first] /= mapVariantToCount[itrVMap->first];
      //DEBUG_MSG(itrVMap->first << "  " << mapVariantToMass[itrVMap->first]);
    }
    DEBUG_VAR(mapVariantToMass.size());

    // Copy to starting PSM set
    for (itrm = mapPeptideToPsmSet.begin(); itrm != itrmEnd; itrm++) {
      set<psmPtr, sortPsmsByMass>::iterator itrs = itrm->second.begin();
      set<psmPtr, sortPsmsByMass>::iterator itrsEnd = itrm->second.end();

      // Set the exact mass to the average variant mass for the group
      (*itrs)->m_exactmass = mapVariantToMass[(*itrs)->m_variantGroup];

      pair<string, float> variantPair =
        make_pair<string, float>(itrm->first, (*itrs)->m_exactmass);
      variantGroups[variantPair] = (*itrs)->m_variantGroup;
      for ( ; itrs != itrsEnd; itrs++) {
        psmSetStarting.push_back(*itrs);
      }
    }

  } // else (compute all variants)

  DEBUG_VAR(psmSetStarting.size());
#if 0
  if (!psmSetStarting.saveToFile("adebug.txt", true, true)) {
    ERROR_MSG("Error saving to PSM file [adebug.txt]");
    return -1;
  }
#endif

  // Count number of variants per peptide
  map<string, set<int> > mapPeptideToVariantSet;
  for (int i = 0; i < psmSetStarting.size(); i++) {
    if (psmSetStarting[i]->m_isDecoy && !doingDecoys) {
      continue;
    }

    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetStarting[i]->m_annotation, cleanAnno);
    //DEBUG_VAR(cleanAnno);
    
    mapPeptideToVariantSet[cleanAnno].insert(psmSetStarting[i]->m_variantGroup);
    //DEBUG_VAR(mapPeptideToVariantSet[cleanAnno].size());
  }
  DEBUG_VAR(mapPeptideToVariantSet.size());

  string pepVarFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    pepVarFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    pepVarFile += "/";
  }
  pepVarFile += argv[5];
  ofstream ofsPepVar(pepVarFile.c_str());
  if (!ofsPepVar) {
    ERROR_MSG("Error opening peptide variant file [ " << pepVarFile << "]");
    return -1;
  }
  
  //---------------------------------------------------------------------------
	// Output histogram of peptide variant counts
  //---------------------------------------------------------------------------
  ofsPepVar << "Peptide\tVariant_Count" << endl;
  map<string, set<int> >::iterator itrPep = mapPeptideToVariantSet.begin();
  map<string, set<int> >::iterator itrPepEnd = mapPeptideToVariantSet.end();
  for (; itrPep != itrPepEnd; itrPep++) {
    ofsPepVar << itrPep->first << "\t";
    ofsPepVar << itrPep->second.size() << endl;
  }
  ofsPepVar.close();
  DEBUG_TRACE;
  
  if (commandLineParams.exists("DONT_COMPUTE_REGIONS") || useExistingRegions) {
    string outputFile;
    if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
      outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
      outputFile += "/";
    }
    outputFile += argv[3];
    DEBUG_MSG("Saving file [" << outputFile << "]");
    psmSetStarting.saveToFile(outputFile.c_str(), true, true);
    return(0);
  }
  
  DB_fasta dbFasta;
  DEBUG_MSG("Loading Database ...");
  if (dbFasta.Load(argv[2]) <= 0) {
    ERROR_MSG("Loading FASTA file [" << argv[2] << "]");
    return -3;
  }
  
  int minAnnoLength = 10;
  for (int i = 0; i < psmSetStarting.size(); i++) {
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetStarting[i]->m_annotation, cleanAnno);
    if (cleanAnno.length() < minAnnoLength) {
      minAnnoLength = cleanAnno.length(); 
    }
  }
  DEBUG_VAR(minAnnoLength);

  // "Normalize" the AA's in the database
  DB_fasta dbCopy;  // Create the DB first
  if (!commandLineParams.exists("DONT_COMPUTE_REGIONS")) {
    dbCopy = dbFasta;   // Then copy.. if you don't.. you seg fault! (dunno why)
    dbCopy.replaceAA('L', 'I');
    dbCopy.replaceAA('K', 'Q');
  }
  DB_index indexAnno(dbCopy, 1024, 1024, minAnnoLength);

  map<int, list<pair<int,int> > > regions;
  
  //--------------------------------------------------------------------
  // Find matches in the database and compute their start and end points
  //      then put the regions into map by protein (index)
  //--------------------------------------------------------------------
  for (int i = 0; i < psmSetStarting.size(); i++) {
    if (psmSetStarting[i]->m_isDecoy && !doingDecoys) {
      continue;
    }
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetStarting[i]->m_annotation, cleanAnno);
    string convertedAnno = convertAnnotation(cleanAnno);

    // Find all the tag matches
    list<tuple<int, int, int> > locations;
    findMatchingLocations(dbCopy, indexAnno, convertedAnno, locations);
    list<tuple<int, int, int> >::iterator itr = locations.begin();
    list<tuple<int, int, int> >::iterator itrEnd = locations.end();
    for (; itr != itrEnd; itr++) {
        //DEBUG_MSG(itr->m0 << "  " << itr->m1 << "  " << itr->m2);
        psmSetStarting[i]->m_protein = dbCopy.getID(itr->m0);
        psmSetStarting[i]->m_dbIndex = itr->m0;
        psmSetStarting[i]->m_startIndex = itr->m1;
        psmSetStarting[i]->m_endIndex = itr->m2;
        
        pair<int, int> region;
        region.first = itr->m1;
        region.second = itr->m2;
        regions[itr->m0].push_back(region);
    }
  }
  DEBUG_VAR(regions.size()); 
  
  //--------------------------------------------------------------------
	// Merge all the regions that overlap
  //--------------------------------------------------------------------
  map<int, list<pair<int,int> > >::iterator itrRegions = regions.begin();
  map<int, list<pair<int,int> > >::iterator itrRegionsEnd = regions.end();
  for (; itrRegions != itrRegionsEnd; itrRegions++) {
  
     //DEBUG_VAR(itrRegions->second.size());
     list<pair<int,int> >::iterator itrListEnd = itrRegions->second.end();
     list<pair<int,int> >::iterator itrFirst = itrRegions->second.begin();

     while (itrFirst != itrListEnd) {
       
       bool merged = false;     
       list<pair<int,int> >::iterator itrList = itrFirst;
       
       for (itrList++; itrList != itrListEnd; itrList++) {
       
         //DEBUG_MSG(itrFirst->first << "  " << itrFirst->second << "  " <<
         //                     itrList->first << "  " << itrList->second);
         if (overlap(itrFirst->first, itrFirst->second,
                     itrList->first, itrList->second, OVERLAP_REGION)) {
           //DEBUG_MSG("OVERLAP");
           itrFirst->first = min(itrFirst->first, itrList->first);
           itrFirst->second = max(itrFirst->second, itrList->second);
           //DEBUG_MSG(itrFirst->first << "  "  << itrFirst->second);
           itrList = itrRegions->second.erase(itrList);
           merged = true;
           break;
         }
         
       } // for (; itrList != itrListEnd; itrList++) {
       
       if (!merged) {
         //DEBUG_MSG("NEXT COMPARE");
         itrFirst++;
         if (itrFirst == itrListEnd) {
           break;
         }
       }
       
     } // while (itrFirst != itrRegions->second.end());
     //DEBUG_VAR(itrRegions->second.size());
  }

  DEBUG_VAR(regions.size());
  DEBUG_VAR(variantGroups.size());

#if 0
  int totalRegions = 0;
  itrRegions = regions.begin();
  for (; itrRegions != itrRegionsEnd; itrRegions++) {
     string sequence = dbFasta.getSequence(itrRegions->first);
     totalRegions += itrRegions->second.size();
     DEBUG_VAR(itrRegions->second.size()); 
     list<pair<int,int> >::iterator itrList = itrRegions->second.begin();
     list<pair<int,int> >::iterator itrListEnd = itrRegions->second.end();
     for (; itrList != itrListEnd; itrList++) {
       string seqAnno = sequence.substr(itrList->first, itrList->second - itrList->first + 1);
       DEBUG_MSG(itrList->first << "  " << itrList->second << "  " << seqAnno);
     }
  }
  DEBUG_VAR(totalRegions);
#endif
	
  //--------------------------------------------------------------------
  // Count how many variants are in each region
  //--------------------------------------------------------------------
  int maxRegionCount = 0;
  map<tuple<int,int,int>, int, TupeCompare> regionCounts;
  map<pair<string, float>, int>::iterator itrVar = variantGroups.begin();
  map<pair<string, float>, int>::iterator itrVarEnd = variantGroups.end();
  for ( ;itrVar != itrVarEnd; itrVar++) {
    string variantString = itrVar->first.first;
    //DEBUG_VAR(variantString);
    string variantConverted = convertAnnotation(variantString);
    list<tuple<int, int, int> > locations;
    findMatchingLocations(dbCopy, indexAnno, variantConverted, locations);
    list<tuple<int, int, int> >::iterator itr = locations.begin();
    list<tuple<int, int, int> >::iterator itrEnd = locations.end();
    //DEBUG_VAR(locations.size());
    for (; itr != itrEnd; itr++) {
    
       list<pair<int,int> >::iterator itrList = regions[itr->m0].begin();
       list<pair<int,int> >::iterator itrListEnd = regions[itr->m0].end();
       for (; itrList != itrListEnd; itrList++) {
         if (overlap(itr->m1, itr->m2,
                     itrList->first, itrList->second, OVERLAP_REGION)) {
           //DEBUG_MSG("   " << itr->m0 << "  " << itrList->first);
           tuple<int,int,int> regionId = make_tuple(itr->m0, itrList->first, itrList->second);
           regionCounts[regionId]++;
           if (regionCounts[regionId] > maxRegionCount) {
             maxRegionCount = regionCounts[regionId];
           }
           //DEBUG_VAR(regionCounts.size());
           //DEBUG_VAR(maxRegionCount);
         }
       }
    }
  }  
  DEBUG_VAR(regionCounts.size());
  DEBUG_VAR(maxRegionCount);

  map<tuple<int,int,int>, int, TupeCompare>::iterator itrReg = regionCounts.begin();
  map<tuple<int,int,int>, int, TupeCompare>::iterator itrRegEnd = regionCounts.end();
#if 0
  for ( ;itrReg != itrRegEnd; itrReg++) {
    string sequence = dbFasta.getSequence(itrReg->first.m0);

    string seqAnno = sequence.substr(itrReg->first.m1, 
                                     itrReg->first.m2 - itrReg->first.m1 + 1);
    //DEBUG_MSG(itrReg->first.m0 << "  " << seqAnno << " " << itrReg->second);
  }
#endif

  map<int, string> mapPvrRegions;

  bool addFlankingAas = commandLineParams.exists("ADD_FLANKING_AAS");
  DEBUG_VAR(addFlankingAas);

  map<int,int> regionCounts2;

  //--------------------------------------------------------------------
	// Put all PSMs into all groups (duplicate PSMs)
  //--------------------------------------------------------------------
  PeptideSpectrumMatchSet psmSetOutput;
  map<int, set<int> > mapRegionCountHistogram;  // histogram of region counts
  for (int i = 0; i < psmSetStarting.size(); i++) {
    if (psmSetStarting[i]->m_isDecoy && !doingDecoys) {
      continue;
    }

    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetStarting[i]->m_annotation, cleanAnno);
    string convertedAnno = convertAnnotation(cleanAnno);
    
    map<tuple<int,int,int>, int, TupeCompare>::iterator itrReg = regionCounts.begin();
    map<tuple<int,int,int>, int, TupeCompare>::iterator itrRegEnd = regionCounts.end();
    for (int iRegion = 1; itrReg != itrRegEnd; itrReg++, iRegion++) {
      string regionAnno = getSequence(dbCopy, itrReg->first.m0, 
                                      itrReg->first.m1, itrReg->first.m2);
      string regionAnnoOrig = getSequence(dbFasta, itrReg->first.m0, 
                                      itrReg->first.m1, itrReg->first.m2);
    
      if (regionAnno.find(convertedAnno) != string::npos) {
      
        unsigned int start = itrReg->first.m1 + regionAnno.find(convertedAnno);
        unsigned int end = start + convertedAnno.length() - 1;
      
        psmPtr pDupe(new PeptideSpectrumMatch);
        *pDupe = *psmSetStarting[i];
        
        stringstream ss;
        ss << iRegion;
        pDupe->m_peptideRegionGroup = ss.str();
        // I'm not sure why this is necessary exactly.. but it is
        //   Need to recount the region counts
        regionCounts2[iRegion]++;

        pDupe->m_peptideRegion = regionAnnoOrig;
        pDupe->m_peptideRegionLength = regionAnnoOrig.size();
        mapPvrRegions[iRegion] = regionAnnoOrig;
        pDupe->m_dbIndex = itrReg->first.m0;
        pDupe->m_protein = dbFasta.getID(itrReg->first.m0);
        pDupe->m_startIndex = start;
        pDupe->m_endIndex = end;

        if (!aaFile.empty()) {
          pDupe->changeGapAnnosToSingle();
          pDupe->addFixedCysteineMods();
        }

        if (addFlankingAas) {               
          // While we are at it.. lets put preceding and following AA's on anno
          string dbSequence = dbFasta.getSequence(pDupe->m_dbIndex);

          string prevAA = "_";
          if (start != 0) {
            prevAA = dbSequence[start - 1];
          }
          string postAA = "_";
          if (end < dbSequence.length() - 1) {
            postAA = dbSequence[end + 1];
          }
          //DEBUG_VAR(pDupe->m_annotation);
          pDupe->m_annotation = prevAA + "." + pDupe->m_annotation + "." + postAA;
          //DEBUG_VAR(pDupe->m_annotation);
        }

        psmSetOutput.push_back(pDupe);

        // Make histogram at the same time
        //mapRegionCountHistogram[itrReg->second].insert(iRegion);
      }
    }
  }
  DEBUG_VAR(psmSetOutput.size());


  bool noDupes = commandLineParams.exists("NO_DUPLICATES");
  DEBUG_VAR(noDupes);

  PeptideSpectrumMatchSet psmSetOutput2;
  if (noDupes) {
    map<pair<string, float>, psmPtr> mapVariantDupes;
    for (int i = 0; i < psmSetOutput.size(); i++) {

      //DEBUG_VAR(psmSetOutput[i]->m_annotation);
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetOutput[i]->m_annotation, cleanAnno);
      //DEBUG_VAR(cleanAnno);

      float totalMass = computeTotalMass(psmSetOutput[i], jumps);
      //DEBUG_VAR(totalMass);
      
      int massToCompare = (int)(totalMass / parentMassTol);
      //DEBUG_VAR(massToCompare);

      // Compute the variants
      pair<string, float> variantPair = 
        make_pair<string, float>(cleanAnno, massToCompare);

      if (mapVariantDupes.find(variantPair) == mapVariantDupes.end()) {
        mapVariantDupes[variantPair] = psmSetOutput[i];
      } else if (mapVariantDupes[variantPair]->m_pValue > psmSetOutput[i]->m_pValue) {
        mapVariantDupes[variantPair] = psmSetOutput[i];
      }
    }
    DEBUG_VAR(mapVariantDupes.size());
    map<pair<string, float>, psmPtr>::iterator itrDupe = mapVariantDupes.begin();
    map<pair<string, float>, psmPtr>::iterator itrDupeEnd = mapVariantDupes.end();
    for ( ; itrDupe != itrDupeEnd; itrDupe++) {
      psmSetOutput2.push_back(itrDupe->second);
    }
  } else {
    for (int i = 0; i < psmSetOutput.size(); i++) {
      psmSetOutput2.push_back(psmSetOutput[i]);
    }
  } // if (noDupes) {
  DEBUG_VAR(psmSetOutput2.size());

  string outputFile;
  string histFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
    histFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    histFile += "/";
  }

  outputFile += argv[3];
  DEBUG_MSG("Saving file [" << outputFile << "]");
  if (!psmSetOutput2.saveToFile(outputFile.c_str(), true, true)) {
    ERROR_MSG("Error saving to PSM file [ " << outputFile << "]");
    return -1;
  }
  
  histFile += argv[4];
  ofstream ofsHist(histFile.c_str());
  if (!ofsHist) {
    ERROR_MSG("Error opening histogram file [ " << histFile << "]");
    return -1;
  }

  ofsHist << "RegionNum\tRegionCount\tRegionAnno\tStart\tEnd\tProtein" << endl;
  itrReg = regionCounts.begin();
  itrRegEnd = regionCounts.end();
  for (int iRegion = 1; itrReg != itrRegEnd; itrReg++, iRegion++) {

    string regionAnnoOrig = getSequence(dbFasta, itrReg->first.m0,
                                    itrReg->first.m1, itrReg->first.m2);

    ofsHist << iRegion << "\t";
    ofsHist << regionCounts2[iRegion] << "\t";
    ofsHist << regionAnnoOrig << "\t";
    ofsHist << itrReg->first.m1 << "\t";
    ofsHist << itrReg->first.m2 << "\t";
    ofsHist << dbFasta.getID(itrReg->first.m0) << endl;
  }
  ofsHist.close();
  DEBUG_TRACE;

  //--------------------------------------------------------------------
  // Create the merged result file (if desired)
  //--------------------------------------------------------------------
  bool mergeDupes = commandLineParams.exists("MERGE_DUPLICATES");
  DEBUG_VAR(mergeDupes);

  bool showAllRegions = commandLineParams.exists("SHOW_ALL_REGIONS");
  DEBUG_VAR(showAllRegions);

  if (mergeDupes) {

    map<pair<string, float>, set<string> > mapVariantSetRegions;

    PeptideSpectrumMatchSet psmSetOutput3;
    map<pair<string, float>, psmPtr> mapVariantDupes;
    for (int i = 0; i < psmSetOutput2.size(); i++) {

      //DEBUG_VAR(psmSetOutput[i]->m_annotation);
      string cleanAnno;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSetOutput2[i]->m_annotation, cleanAnno);
      //DEBUG_VAR(cleanAnno);

      float totalMass = computeTotalMass(psmSetOutput2[i], jumps);
      //DEBUG_VAR(totalMass);

      int massToCompare = (int)(totalMass / parentMassTol);
      //DEBUG_VAR(massToCompare);

      // Compute the variants
      pair<string, float> variantPair =
        make_pair<string, float>(cleanAnno, massToCompare);

      if (mapVariantDupes.find(variantPair) == mapVariantDupes.end()) {
        mapVariantDupes[variantPair] = psmSetOutput2[i];
        mapVariantSetRegions[variantPair].insert(psmSetOutput2[i]->m_peptideRegion);
      } else if (mapVariantDupes[variantPair]->m_scanNum == psmSetOutput2[i]->m_scanNum &&
                 mapVariantDupes[variantPair]->m_protein != psmSetOutput2[i]->m_protein) {
        mapVariantDupes[variantPair]->m_peptideRegionGroup += ";";
        mapVariantDupes[variantPair]->m_peptideRegionGroup += psmSetOutput2[i]->m_peptideRegionGroup;
        mapVariantSetRegions[variantPair].insert(psmSetOutput2[i]->m_peptideRegion);
        mapVariantDupes[variantPair]->m_peptideRegion += ";";
        mapVariantDupes[variantPair]->m_peptideRegion += psmSetOutput2[i]->m_peptideRegion;
        mapVariantDupes[variantPair]->m_protein += ";";
        mapVariantDupes[variantPair]->m_protein += psmSetOutput2[i]->m_protein;
      } else {
        //DEBUG_VAR(psmSetOutput2[i]->m_scanNum);
        //DEBUG_VAR(psmSetOutput2[i]->m_protein);
        //DEBUG_VAR(mapVariantDupes[variantPair]->m_peptideRegionGroup);
        //DEBUG_VAR(mapVariantDupes[variantPair]->m_peptideRegion);
        //DEBUG_VAR(mapVariantDupes[variantPair]->m_protein);
        //DEBUG_VAR(mapVariantSetRegions[variantPair].size());
      }
    }

    map<pair<string, float>, psmPtr>::iterator itrDupe = mapVariantDupes.begin();
    map<pair<string, float>, psmPtr>::iterator itrDupeEnd = mapVariantDupes.end();
    for ( ; itrDupe != itrDupeEnd; itrDupe++) {

      if (!showAllRegions) {
        string cleanAnno;
        PeptideSpectrumMatch::getUnmodifiedPeptide(itrDupe->second->m_annotation, cleanAnno);
        float totalMass = computeTotalMass(itrDupe->second, jumps);
        int massToCompare = (int)(totalMass / parentMassTol);
        pair<string, float> variantPair =
          make_pair<string, float>(cleanAnno, massToCompare);

        set<string>::iterator itrSet = mapVariantSetRegions[variantPair].begin();
        set<string>::iterator itrSetEnd = mapVariantSetRegions[variantPair].end();
        itrDupe->second->m_peptideRegion = "";
        for (int i = 0; itrSet != itrSetEnd; itrSet++, i++) {
          itrDupe->second->m_peptideRegion += *itrSet;
          set<string>::iterator itrSetNext = itrSet;
          itrSetNext++;
          if (itrSetNext != itrSetEnd) {
            itrDupe->second->m_peptideRegion += ";";
          }
        }
      }
      psmSetOutput3.push_back(itrDupe->second);
    }
    DEBUG_VAR(psmSetOutput3.size());

    string outputFile;
    if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
      outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
      outputFile += "/";
    }
    outputFile += commandLineParams.getValue("MERGE_DUPLICATES");
    DEBUG_MSG("Saving file [" << outputFile << "]");
    if (!psmSetOutput3.saveToFile(outputFile.c_str(), true, true)) {
      ERROR_MSG("Error saving to PSM file [ " << outputFile << "]");
      return -1;
    }

  } // if (mergeDupes) {

  //---------------------------------------------------------------------------
	// Output histogram of region counts and all regions that make them up
  //---------------------------------------------------------------------------
  
  return 0;
}
