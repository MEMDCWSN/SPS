//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include "tuple.h"

#include <stdlib.h>

#define DEBUG_ 0
#define DEBUG_PERCENT_MATCH 0

#define THRESHOLD_MATCH 0.8

using namespace std;
using namespace sps;
using namespace specnets;

struct Variant {
  string   peptide;
  float    mass;
  int   match1;
  int   match2;
  int   match3;
  float   pm1;
  float   pm2;
  float   pm3;
  float   pvalue1;
  float   pvalue2;
  float   pvalue3;
  int   scan1;
  int   scan2;
  int   scan3;
  string   anno1;
  string   anno2;
  string   anno3;

  Variant() {
    match1 = 0;
    match2 = 0;
    match3 = 0;
    pm1 = 0.0;
    pm2 = 0.0;
    pm3 = 0.0;
    pvalue1 = -1.0;
    pvalue2 = -1.0;
    pvalue3 = -1.0;
    scan1 = -1;
    scan2 = -1;
    scan3 = -1;
    anno1="";
    anno2="";
    anno3="";
  }
};

// -------------------------------------------------------------------------
string convertAnnotation(string annotation)
{
  string convertedAnnotation; 
  for (int iChar = 0; iChar < annotation.length(); iChar++) {
    if (annotation[iChar] == 'I') {
      convertedAnnotation += 'L';
    } else if (annotation[iChar] == 'Q') {
      convertedAnnotation += 'K';
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
  if (DEBUG_) DEBUG_VAR(cleanAnno);
  float totalMass = jumps.getPeptideMass(cleanAnno);
  if (DEBUG_) DEBUG_VAR(totalMass);

  vector<float> modifications;
  vector<unsigned int> positions;
  p->getModificationsAndPositions(modifications, positions);

  for (int i = 0; i < modifications.size(); i++) {
    totalMass += modifications[i];
  }

  if (DEBUG_) DEBUG_VAR(totalMass);
  return totalMass;

}

// -------------------------------------------------------------------------
float getPercentMatch(string & string1, string & string2)
{
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string1.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string2.size());
  int dynaArray[string1.size()+1][string2.size()+1];
  for (int i = 0; i < string1.size() + 1; i++) {
    for (int j = 0; j < string2.size() + 1; j++) {
      dynaArray[i][j] = 0;
    }
  }
  for (int i = 1; i <= string1.size(); i++) {
    for (int j = 1; j <= string2.size(); j++) {
      int up = dynaArray[i-1][j];
      int left = dynaArray[i][j-1];
      int diag = dynaArray[i-1][j-1];
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(i << "  " << j);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(string1[i-1] << "  " << string2[j-1]);
      if (string1[i-1] == string2[j-1]) {
        diag++;
      }
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      int max1 = max(left, up);
      dynaArray[i][j] = max(max1, diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(dynaArray[i][j]);
    }
  }

  if (DEBUG_PERCENT_MATCH) {
    cout << string1 << "  " << string2 << endl;
    for (int i = 0; i <= string1.size(); i++) {
      for (int j = 0; j <= string2.size(); j++) {
        cout << dynaArray[i][j] << "  ";
      }
      cout << endl;
    }
    cout << endl;
  }
  int match = dynaArray[string1.size()][string2.size()];
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(match);
  float percent = (float)match / (float)min(string1.size(), string2.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  percent = min(percent, (float)1.0);
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  return percent;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 3 && argc != 4) {
    cerr << "Usage: main_comparepsms psm_file1 psm_file2 [psm_file2]" << endl;
    return -1;
  }
  
  // Load amino acid masses
  DEBUG_MSG("Loading Database ...");
  AAJumps jumps(1);
/*
  if (!aaFile.empty()) {
    jumps.loadJumps(aaFile.c_str());
  }
*/

  PeptideSpectrumMatchSet	psmSet1;
  if (!psmSet1.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -2;
  }
  if (DEBUG_) DEBUG_VAR(psmSet1.size());
  
  PeptideSpectrumMatchSet	psmSet2;
  if (!psmSet2.loadFromFile(argv[2])) {
    ERROR_MSG("Loading PSM file [" << argv[2] << "]");
    return -3;
  }
  if (DEBUG_) DEBUG_VAR(psmSet2.size());

  PeptideSpectrumMatchSet	psmSet3;
  if (argc == 4 && !psmSet3.loadFromFile(argv[3])) {
    ERROR_MSG("Loading PSM file [" << argv[3] << "]");
    return -4;
  }
  if (DEBUG_) DEBUG_VAR(psmSet3.size());
  
  float parentMassTol = 2.5;
  int pmTol = (int)(parentMassTol * 100.0 + 0.5);
  if (DEBUG_) DEBUG_VAR(pmTol);

  vector<Variant> variants;

  //----------------------------------------------------------
  // Set all the variant matches
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSet1.size(); iPSM++ ) {
    if (psmSet1[iPSM]->m_isDecoy) {
      continue;
    }
    //if (DEBUG_) DEBUG_VAR(psmSet1[iPSM]->m_annotation);
    string convertAnno = convertAnnotation(psmSet1[iPSM]->m_annotation);
    //if (DEBUG_) DEBUG_VAR(convertAnno);
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno, cleanAnno);
    float mass = computeTotalMass(psmSet1[iPSM], jumps);

    Variant newVariant;
    newVariant.peptide = cleanAnno;
    newVariant.anno1 = psmSet1[iPSM]->m_annotation;
    newVariant.match1 = 1;
    newVariant.mass = mass;
    newVariant.pm1 = mass;
    newVariant.pvalue1 = psmSet1[iPSM]->m_pValue;
    newVariant.scan1 = psmSet1[iPSM]->m_scanNum;

    variants.push_back(newVariant);
  }

  for (int iPSM = 0; iPSM < psmSet2.size(); iPSM++ ) {
    if (psmSet2[iPSM]->m_isDecoy) {
      continue;
    }
    //if (DEBUG_) DEBUG_VAR(psmSet2[iPSM]->m_annotation);
    string convertAnno = convertAnnotation(psmSet2[iPSM]->m_annotation);
    //if (DEBUG_) DEBUG_VAR(convertAnno);
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno, cleanAnno);
    float mass = computeTotalMass(psmSet2[iPSM], jumps);
    //if (DEBUG_) DEBUG_VAR(mass);

    // Compare to all previous variants
    bool found = false;
    int j = 0;
    for ( ; j < variants.size() ; j++) {

      // Mass within tolerance (faster to compute than getPercentMatch()
      if (abs(variants[j].mass - mass) > parentMassTol) {
        continue;
      }

      float percentMatch = getPercentMatch(variants[j].peptide, cleanAnno);
      if (percentMatch < THRESHOLD_MATCH) {
        continue;
      }
      found = true;
      break;
    }

    if (found == true) {
      variants[j].anno2 = psmSet2[iPSM]->m_annotation;
      variants[j].match2 = 1;
      variants[j].pm2 = mass;
      variants[j].pvalue2 = psmSet2[iPSM]->m_pValue;
      variants[j].scan2 = psmSet2[iPSM]->m_scanNum;
    } else {
      Variant newVariant;
      newVariant.peptide = cleanAnno;
      newVariant.mass = mass;
      newVariant.anno2 = psmSet2[iPSM]->m_annotation;
      newVariant.match2 = 1;
      newVariant.pm2 = mass;
      newVariant.pvalue2 = psmSet2[iPSM]->m_pValue;
      newVariant.scan2 = psmSet2[iPSM]->m_scanNum;
      variants.push_back(newVariant);
    }
  }

  for (int iPSM = 0; iPSM < psmSet3.size(); iPSM++ ) {
    if (psmSet3[iPSM]->m_isDecoy) {
      continue;
    }
    //if (DEBUG_) DEBUG_VAR(psmSet3[iPSM]->m_annotation);
    string convertAnno = convertAnnotation(psmSet3[iPSM]->m_annotation);
    //if (DEBUG_) DEBUG_VAR(convertAnno);
    string cleanAnno;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno, cleanAnno);
    float mass = computeTotalMass(psmSet3[iPSM], jumps);
    //if (DEBUG_) DEBUG_VAR(mass);

    // Compare to all previous variants
    bool found = false;
    int j = 0;
    for ( ; j < variants.size() ; j++) {

      // Mass within tolerance (faster to compute than getPercentMatch()
      if (abs(variants[j].mass - mass) > parentMassTol) {
        continue;
      }

      float percentMatch = getPercentMatch(variants[j].peptide, cleanAnno);
      if (percentMatch < THRESHOLD_MATCH) {
        continue;
      }
      found = true;
      break;
    }

    if (found == true) {
      Variant newVariant;
      variants[j].anno3 = psmSet3[iPSM]->m_annotation;
      variants[j].match3 = 1;
      variants[j].pm3 = mass;
      variants[j].pvalue3 = psmSet3[iPSM]->m_pValue;
      variants[j].scan3 = psmSet3[iPSM]->m_scanNum;
    } else {
      Variant newVariant;
      newVariant.peptide = cleanAnno;
      newVariant.anno3 = psmSet3[iPSM]->m_annotation;
      newVariant.match3 = 1;
      newVariant.mass = mass;
      newVariant.pm3 = mass;
      newVariant.pvalue3 = psmSet3[iPSM]->m_pValue;
      newVariant.scan3 = psmSet3[iPSM]->m_scanNum;
      variants.push_back(newVariant);
    }
  }

  //----------------------------------------------------------
  // Output results
  //----------------------------------------------------------
  if (argc == 4) {
    cout << "anno\tmass\tmatch1\tmatch2\tmatch3\tpm1\tpm2\tpm3\tpvalue1\tpvalue2\tpvalue3\tscan1\tscan2\tscan3\tanno1\tanno2\tanno3" << endl;
  } else {
    cout << "anno\tmass\tmatch1\tmatch2\tpm1\tpm2\tpvalue1\tpvalue2\tscan1\tscan2\tanno1\tanno2" << endl;
  }
  for (int k = 0 ; k < variants.size() ; k++) {

      cout << variants[k].peptide << "\t";
      cout << variants[k].mass << "\t";

      cout << variants[k].match1 << "\t";
      cout << variants[k].match2 << "\t";
      if (argc == 4) {
        cout << variants[k].match3 << "\t";
      }
      cout << variants[k].pm1 << "\t";
      cout << variants[k].pm2 << "\t";
      if (argc == 4) {
        cout << variants[k].pm3 << "\t";
      }
      cout << variants[k].pvalue1 << "\t";
      cout << variants[k].pvalue2 << "\t";
      if (argc == 4) {
        cout << variants[k].pvalue3 << "\t";
      }
      cout << variants[k].scan1 << "\t";
      cout << variants[k].scan2 << "\t";
      if (argc == 4) {
        cout << variants[k].scan3 << "\t";
      }
      cout << variants[k].anno1 << "\t";
      cout << variants[k].anno2 << "\t";
      if (argc == 4) {
        cout << variants[k].anno3 << "\t";
      }
      cout << endl;
  }

  cout << endl;
  return 0;
}
