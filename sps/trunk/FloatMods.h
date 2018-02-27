/*
 * Mass-spectrometry clustering Project
 * This program is property of University of California at San Diego
 *
 * Supervisor: Pavel Pevzner
 * coded by Dumitru Brinza 
 *
 * ppevzner@cs.ucsd.edu
 * dima@cs.ucsd.edu
 *
 */

#ifndef FLOATMODS_H_
#define FLOATMODS_H_

#include <string>
#include <vector>
#include <sstream>

// FloatMods includes a set of classes and functions required 
// for floating modifications
// 
// example of usage:
//
// 1. TFloatMods ModFloat;
// 2. ModFloat.genome = genome;
// 3. ModFloat.FloatMods(peptide, variation, MaxMods, spectrum*, floated_pep, topK);
//
// 1. -- declare a variable of class TfloatMods
// 2. -- set the database, database should be one string, 
//       proteins in the database are separated by space, unknown are X
// 3. -- call Floating function which returns floated peptide
//       peptide -- sequence of peptide
//       variation -- how much to float to the left/right
//       MaxMods -- peptides with more than MaxMods modifications will not be floated
//       spectrum* -- pointer to the spectrum of a peptide

//---------------------------------------------------------------------------------
struct peak
{
public:
 float mass;
 float intensity;
 peak()
 {
  mass = 0;
  intensity = 0;
 }
 void operator=(peak& data)
 {
	mass = data.mass;
	intensity = data.intensity;
 }
};
//---------------------------------------------------------------------------------
class spectra
{
public:
 float total_intensity;
 unsigned short n_peaks;
 float parent_mass;
 peak *spectrum;
 spectra()
 {
  total_intensity = 0;
  spectrum = NULL;
  n_peaks = 0;
 }
};
//---------------------------------------------------------------------------------
struct TDB_Title{
vector<std::string> title;
vector<int> pos;
inline void set (std::string t, int pz){title.push_back(t); pos.push_back(pz);}
inline std::string get(int p)
{ int i=0; for(;i<pos.size()-1;i++) if(p<pos[i+1]) break; if(p>pos[i+1]) i++; return title[i];}
};
//---------------------------------------------------------------------------------
   struct TPeptide
  {
   vector<char>  amino; // 0 - if it is a mod
   vector<float> mass;
   vector<short> mods;
   
   bool init(std::string peptide); // returns false if peptide does not have mods
   void GenerateMixedMasses(vector<float> & masses); // generates a set of masses, where mods are mixed with aa
   void GenerateCandidateMasses(vector<float> & masses); // generate all possible modifications of existing masses 
   bool ValidPeptide(); // returns true if peptides modifications are in acceptable bounds
   bool ValidateCase3();
   bool ValidPeptidePrefix(int max_mod_id);
   std::string getSequence(short withmods); // returns peptide string without/with mods or with mixing close mods in [], 0,1,2
   std::string MixCloseMods();

  };
//---------------------------------------------------------------------------------
struct TPepScore{
  std::string peptide;
  double intensity;
  TPepScore(std::string seq, double intens){peptide=seq;intensity=intens;}
  
  bool operator<(const TPepScore& data) const
	 {
	  if(intensity > data.intensity) return true;
    return false;
	 }
};

struct TPepScoreList{
vector<TPepScore> list;
int top;

void add(TPepScore e){
list.push_back(e);
if(list.size()>100000)
{
sort(list.begin(), list.end());
vector<TPepScore> temp;
temp.push_back(list[0]);
int k=1;
for(int i=0;k<top && i<list.size()-1;i++) if(list[i].peptide!=list[i+1].peptide) {k++;temp.push_back(list[i+1]);}
list.clear();
for(int i=0;i<top && i<temp.size();i++) list.push_back(temp[i]);
}
}
void getList(vector<std::string> &v)
{
sort(list.begin(), list.end());
int k=0;
if(list.size()>0){v.push_back(list[0].peptide); k=1;}
for(int i=0;k<top && i<list.size()-1;i++) if(list[i].peptide!=list[i+1].peptide) {v.push_back(list[i+1].peptide);k++;}
}
void add(std::string seq, double intens){ add(TPepScore(seq, intens));}
};

//---------------------------------------------------------------------------------
class TFloatMods{

private:

bool        initial_pep_not_valid;
std::string  best_peptide;
double      max_intensity;
spectra *    spect;


vector<std::string> extended_peptides;
vector<std::vector <char> > pref_permutations;
vector<std::vector <char> > suff_permutations;

std::string genome;
std::string genome_orig;
TPepScoreList pep_list;

public:


TFloatMods(std::string &g, std::string &og, int t)
{
 pep_list.top = t;
 genome = g;
 genome_orig = og;
}

void FloatMods(std::string peptide, int variation, int maxMods, spectra * spec, vector<std::string> &list, int & mappos);
double ComputeExplainedIntensity(vector<float> &masses);
void generatePrefixExtensions(TPeptide tpep, int variation, int fix, int x, int y, int mod_id);
void generateSuffixExtensions(TPeptide tpep, int mod_id, int left, int variation, int fix, int y);
void permute(TPeptide tpep,int mod_id, int left, int right, int fix);
inline void EvaluateRounds(TPeptide & tpep);
inline void EvaluateOne(TPeptide & tpep);
void generatePermutations(vector<char> from, vector<char> to, int x, bool flag);
void SrinkPeptide(TPeptide & tpep);
void PrepareShrinkSuffix(TPeptide & tpep, int next_mod_id);
void SrinkPrefix(TPeptide tpep, int perm_id, int mods_id);
void SrinkSuffix(TPeptide tpep, int perm_id, int mods_id);

};

//---------------------------------------------------------------------------------

double ms(std::string s);
double ms(char c);
inline std::string ftos(float i){std::stringstream s;s << i;return s.str();}
inline std::string itos(int i){std::stringstream s;s << i;return s.str();}

//---------------------------------------------------------------------------------

#endif /*FLOATMODS_H_*/

