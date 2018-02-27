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


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "FloatMods.h"


using namespace std;

//---------------------------------------------------------------------

void TFloatMods::FloatMods(string peptide, int variation, int maxMods, spectra * spec, vector<std::string> &list, int & mappos)
 {

   spect = spec;
   TPeptide tpep;
   pep_list.list.clear();
   
   if(!tpep.init(peptide)) {list.push_back(peptide);  return;}

   vector<float> masses;
   tpep.GenerateCandidateMasses(masses);
     
   max_intensity = ComputeExplainedIntensity(masses); 
   best_peptide = peptide;
   initial_pep_not_valid = true;
  
   if(tpep.ValidPeptide()) {initial_pep_not_valid = false;} else max_intensity = 0;
   //add original peptide to the list of top peptides
   pep_list.add(best_peptide, max_intensity);
     
   //cout << best_peptide << "\t" <<  max_intensity  << endl;
     
   //initialize with empty set of shrinking peptides
   extended_peptides.clear();
   // not completely empty
   extended_peptides.push_back(best_peptide);

   //make sure that extending or srinking may happen
   bool hasPositiveMass = false;
   bool hasNegativeMass = false;
     
   int num_poz_mods = 0;
   int num_neg_mods = 0;
     
   for(int j=0;j<tpep.amino.size();j++)
   {
   if(tpep.amino[j]==0) 
     if (tpep.mass[j]>0) {hasPositiveMass = true; num_poz_mods++;} else {hasNegativeMass = true;num_neg_mods++;}
   }
   
   if(num_poz_mods+num_neg_mods>maxMods) {list.push_back(tpep.MixCloseMods());return;}
 
   // extension may happen
   if(hasPositiveMass)
    {
     string pepm = tpep.getSequence(0);
     for(int i=0;i<pepm.length();i++) if(pepm[i]=='Q') pepm[i]='K'; else if(pepm[i]=='I') pepm[i]='L';
     string::size_type loc = genome.find(pepm,0);
     if(loc != string::npos) mappos = loc;
     char pf = 'X', sf = 'X';
     //Map to genome
     
    while( loc != string::npos && (genome[loc-1]!=pf || genome[loc+pepm.length()]!=sf) )
     {
     
      // Extend, permute, Evaluate good configurations, generate set for shrinking
      generatePrefixExtensions(tpep, variation, variation, loc,loc+pepm.length()-1,0);
      pf = genome[loc-1];
      sf = genome[loc+pepm.length()];
     }
     } 
  
     // permute original set of modifications, evaluate good configurations
     // generate set for shrinking
     permute(tpep,0, variation, variation, variation);
  
     // at this point all valid extensions were tested
     bool ext = false;
     if(extended_peptides.size()>0)
     {     
     ext = true;
     vector<string> extPep;
     for(int j=0;j<extended_peptides.size();j++)
     {
     extPep.push_back(extended_peptides[j]);
     }
     
     sort(extPep.begin(), extPep.end());
     
     string prev = extPep[0];
     if(tpep.init(prev)) SrinkPeptide(tpep);
     for(int j=1;j<extPep.size();j++)
     if(prev!=extPep[j])
     {
      prev=extPep[j];
      if(tpep.init(prev)) SrinkPeptide(tpep);
     }
     }
  
    pep_list.getList(list);
    
    for(int i=0;i<list.size();i++)
    {
    if(tpep.init(list[i])) list[i] = tpep.MixCloseMods();
    }
  
}

//---------------------------------------------------------------------------------

double TFloatMods::ComputeExplainedIntensity(vector<float> &masses)
{
      double explained_intensity = 0;
    	 int x=0,y=0;
    	 while(x<spect->n_peaks&&y<masses.size()){
    	 
    	 while(y<masses.size()&&masses[y]<spect->spectrum[x].mass && fabs(-masses[y]+spect->spectrum[x].mass)>0.5)y++;
    	 while(x<spect->n_peaks&&masses[y]>spect->spectrum[x].mass && fabs(masses[y]-spect->spectrum[x].mass)>0.5)x++;
    	 
    	 
    	 if(fabs(masses[y]-spect->spectrum[x].mass)<=0.5){
    	 explained_intensity +=spect->spectrum[x].intensity;
    	 x++; } }
    	 return explained_intensity/spect->total_intensity;
}
//---------------------------------------------------------------------

 // Takes an imput PEPTIDE string and produce a vector of masses and characters 
 // corresponding to each possition, as well as a vector of modifications pointing
 // to their location in first vector. Returns false if there no mods.

bool TPeptide::init(string peptide)
{ 
    amino.clear();
    mass.clear();
    mods.clear();
    
    string::size_type loc1=peptide.find( "[", 0 ),loc2 = 0;
    
    if(loc1==string::npos) return false;
    int k = 0;
     while(loc1!=string::npos)
      {	
       for(int i=loc2;i<loc1;i++){amino.push_back(peptide[i]); mass.push_back(ms(peptide[i]));k++;}
       loc2 = peptide.find( "]", loc1+1 )+1;
       amino.push_back(0); mass.push_back(atof(peptide.substr(loc1+1).c_str()));
       mods.push_back(k);
       k++;
       loc1 = peptide.find( "[", loc2 );
      }
    for(int i=loc2;i<peptide.length();i++){amino.push_back(peptide[i]); mass.push_back(ms(peptide[i]));}
  return true;
}
//---------------------------------------------------------------------
bool TPeptide::ValidateCase3()
{
 float imass;
 char aa  = 0;
 int size = amino.size();

 //cout << "?2";
 //for(int i = 0;i<size;i++) cout << amino_mod[i];
 //cout << endl;


 bool include = false;
 for(int i = 0;i<size;i++)
 {
  imass = 0;
  while(i<size && amino[i]==0)
  {
   imass+=mass[i];
   i++;
  }

   if(imass>200) return false; // check condition #3 = modification shold be smaller than +100   
  
  }  
  
  return true;
}
//----------------------------------------------------------------------------------------
bool TPeptide::ValidPeptidePrefix(int max_mod_id)
{
 float imass;
 float aa  = 0;
 int size = amino.size();

 //cout << "?1";
 //for(int i = 0;i<size;i++) cout << amino_mod[i];
 //cout << endl;


 int k = 0;
 for(int i = 0;(i<size&&k<=max_mod_id);i++)
 {
  imass = 0;
  while((i<size && k<=max_mod_id) && amino[i]==0)
  {
   imass+=mass[i];
   i++;k++;
  }

   if(aa == 0)
   {
    if(imass<-1.0) return false; // check condition #1 = prefix modification shold be positive 
   } 
   else
   if(aa+imass<57 and fabs(aa+imass)>1.0 ) return false; // check condition #2 = aa + modification shold be larger than +57
     
   if(imass>200) return false; // check condition #3 = modification shold be smaller than +100   
   
   if(i<size)aa=mass[i];
  }  
  
  return true;
}
//-------------------------------------------------------------------------------
bool TPeptide::ValidPeptide()
{
 float imass;
 float aa  = 0;
 int size = amino.size();
 
 for(int i = 0;i<size;i++)
 {
  imass = 0;
  while( i<size && amino[i]==0 )
  {
   imass+=mass[i];
   i++;
  }

   if(aa == 0)
   {
    if(imass<-1.0) return false; // check condition #1 = prefix modification shold be positive 
   } 
   else
   if(aa+imass<57 and (fabs(aa + imass)>1.0 or i==size)) return false;// check condition #2 = aa + modification shold be larger than +57
     
   if(imass>200) return false;// check condition #3 = modification shold be smaller than +100   
   
   if(i<size)aa=mass[i];
  }  
  return true;
}
//-------------------------------------------------------------------------------------------------
void TFloatMods::generateSuffixExtensions(TPeptide tpep, int mod_id, int left, int variation, int fix, int y)
{

/*  cout << "<s>" << mod_id << " " << variation << " ";
  for(int i = 0;i<amino_mod.size();i++) cout << amino_mod[i];
  cout << endl;*/

   if(mod_id==tpep.mods.size()){ 
   permute(tpep,0, left, variation, fix);
   // extend to the right 
   //generateSuffixExtensions(amino_mod, hash_mod, 0, fix, fix, y);
   return;
   }

   if(mod_id<tpep.mods.size()) generateSuffixExtensions(tpep, mod_id+1, left, variation, fix, y);



   int k = 0;
   for(;mod_id<tpep.mods.size();mod_id++)
   {
      int pos = tpep.mods[mod_id];
      
       
//     string tme="";
//  for(int i=0;i<amino_mod.size();i++)tme+=amino_mod[i];
//  cout << "\next " << pos << " - "  << hash_mod.size() << "-" << mod_id << " " << tme << endl;
   
    
   // extend to the right
   
    int right = variation+k;
    for(int i=pos+1;(i<tpep.amino.size() && i<=pos+right);i++)if(tpep.amino[i]==0)right++;

    if(pos+right>=tpep.amino.size()-1)
    {
    float mod_mass = tpep.mass[tpep.mods[mod_id]];
    while(mod_mass>-40 && genome[y+1]!=' ')
    {
     float amino_mass = ms(genome_orig[y+1]);
     mod_mass-=amino_mass;
     if(mod_mass>-100)
      {
      k++;
      tpep.amino.push_back(genome_orig[y+1]);
      tpep.mass.push_back(amino_mass);
      //for(int i=0;i<hash_mod.size();i++) hash_mod[i]++;
      tpep.mass[tpep.mods[mod_id]]=mod_mass;
       
//           tme="";
//  for(int i=0;i<amino_mod.size();i++)tme+=amino_mod[i];
//  cout << "\nint " << pos << " - "  << tme << endl;

       
      generateSuffixExtensions(tpep, mod_id+1, left, variation+k, fix, y+1);    
      //permute(amino_mod,hash_mod,0, left, variation+k, variation);
      
      y++;
      } 
      }
     }
    }
}
//-----------------------------------------------------------------------------------------
void TFloatMods::generatePrefixExtensions(TPeptide tpep, int variation, int fix, int x, int y, int mod_id)
{
  /*cout << "<p>" << mod_id << " " << variation << " ";
  for(int i = 0;i<amino_mod.size();i++) cout << amino_mod[i];
  cout << endl;
*/
   if(mod_id==tpep.mods.size()){ 
   //permute(amino_mod,hash_mod,0, variation, fix, fix);
   // extend to the right 
   generateSuffixExtensions(tpep, 0, variation, fix, fix, y);
   return;
   }

   if(mod_id<tpep.mods.size()) generatePrefixExtensions(tpep, variation, fix, x, y, mod_id+1);

   int k = 0;
   for(;mod_id<tpep.mods.size();mod_id++)
   {
    int pos = tpep.mods[mod_id];
    
   // extend to the left
   
    int left = variation+k;
    for(int i=pos-1;(i>=0 && i>=pos-left);i--)if(tpep.amino[i]==0)left++;
       
    if(pos<=left)
    {
    float mod_mass = tpep.mass[tpep.mods[mod_id]];
    while(mod_mass>-40 && genome[x-1]!=' ')
    {
     float amino_mass = ms(genome_orig[x-1]);
     mod_mass-=amino_mass;
     if(mod_mass>-100)
      {
      k++;
      tpep.amino.insert(tpep.amino.begin(),genome_orig[x-1]);
      tpep.mass.insert(tpep.mass.begin(),amino_mass);
      for(int i=0;i<tpep.mods.size();i++) tpep.mods[i]++;
      tpep.mass[tpep.mods[mod_id]]=mod_mass;
       
//    permute(amino_mod,hash_mod,0, variation+k, variation, variation);
    
      // extend to the right 
       generatePrefixExtensions(tpep, variation+k, fix, x-1, y, mod_id+1);
      //generateSuffixExtensions(amino_mod, hash_mod, mod_id+1, variation+k, variation, y);
      
      x--;
      } 
      }
      }
     else return;
    
    }
}
//-----------------------------------------------------------------------------------------
string TPeptide::getSequence(short withmods)
{
  stringstream ss;
  string tmp="";
  float imass = 0;
  for(int i=0;i<amino.size();i++){
   if(amino[i]!=0) 
   {
    if(withmods==2 && imass!=0)
     {
       ss << "[" << imass << "]";
       imass = 0;
     }
    ss << amino[i]; 
   } 
   else 
   {
    imass+=mass[i];
    if(withmods==1) ss << "[" << mass[i] << "]";
   }
  }
  if(withmods==2 && imass!=0) ss << "[" << imass << "]";

  ss >> tmp; 
  return tmp;
}

//-----------------------------------------------------------------------------------------
void TFloatMods::generatePermutations(vector<char> from, vector<char> to, int x, bool flag)
{
if(from.size()>1)
{
 to.push_back(from[x]); from.erase(from.begin()+x);
 for(int i=0;i<from.size();i++)  generatePermutations( from,  to, i, flag);
}
else
{
 to.push_back(from[x]);
 if(!flag) pref_permutations.push_back(to);
 else suff_permutations.push_back(to);
}
}

//-----------------------------------------------------------------------------------------
string TPeptide::MixCloseMods()
{
 float imass;
 string aa  = "";
 int size = amino.size();
 string pep = "";
 for(int i = 0;i<size;i++)
 {
  imass = 0;
  while(i<size && amino[i]==0)
  {
   imass+=mass[i];
   i++;
  }
  if(imass!=0)
  {
   if(imass>0.7) aa +="["+ftos(imass)+"]";
   if(imass<-0.7) 
   {
   if(aa!=""&&fabs(ms(aa)+imass)<0.7) aa="";
   else
   aa +="["+ftos(imass)+"]";
   }
   }
  
  pep += aa;
  if(i<size) {aa = " "; aa[0]=amino[i];} else aa = "";
  }
  pep += aa;
 
 return pep;
}
//-----------------------------------------------------------------------------------------
void TFloatMods::SrinkSuffix(TPeptide tpep, int perm_id, int mods_id)
{
 if(mods_id==suff_permutations[perm_id].size()) 
 {
 if(tpep.ValidPeptide()) EvaluateRounds(tpep); 
 return;
 }
 int mod_id = suff_permutations[perm_id][mods_id];
 float imass = tpep.mass[tpep.mods[mod_id]];
 int i = tpep.amino.size()-1;
 while(imass<200&&i>=0)
 {
  SrinkSuffix(tpep, perm_id, mods_id+1);
  while(i>0 && tpep.amino[i]==0)i--;
  imass+=tpep.mass[i];
  if(imass<200)
  {

  tpep.mass[tpep.mods[mod_id]] = imass;
  tpep.amino.erase(tpep.amino.begin()+i);
  tpep.mass.erase(tpep.mass.begin()+i);
  for(int j=0;j<tpep.mods.size();j++)if(tpep.mods[j]>i)tpep.mods[j]--;

  i--;
  }
 }
}
//-----------------------------------------------------------------------------------------
void TFloatMods::PrepareShrinkSuffix(TPeptide & tpep, int next_mod_id)
{
  int j;
  vector<char> from, to;
  for(j=next_mod_id;j<tpep.mods.size();j++) from.push_back(j);

  //generate a list of orders for mass consuming
  suff_permutations.clear();
  for(j=0;j<from.size();j++) generatePermutations( from,  to, j,1);
  
  //srink their suffix  
  for(j=0;j<suff_permutations.size();j++) SrinkSuffix(tpep,j,0);
}
//-----------------------------------------------------------------------------------------
void TFloatMods::SrinkPrefix(TPeptide tpep, int perm_id, int mods_id)
{

 if(mods_id==pref_permutations[perm_id].size())
 {
  if(tpep.ValidPeptidePrefix(mods_id-1))
  {
  if(mods_id>=tpep.mods.size()) {
  EvaluateRounds(tpep);
  }
    else  
    {
    PrepareShrinkSuffix(tpep, mods_id);
    }
  }
  
  return;
 }
 
 int mod_id = pref_permutations[perm_id][mods_id];
 float imass = tpep.mass[tpep.mods[mod_id]];
 int i=0;
// cout << mod_id << " " << perm_id << " "<< pref_permutations[perm_id].size() << " " << imass << endl;
 while(imass<200&& i<tpep.amino.size())
 {
  SrinkPrefix(tpep, perm_id, mods_id+1);
  while(i<tpep.amino.size() && tpep.amino[i]==0)i++;
  if(i==tpep.amino.size()) break;
  imass+=tpep.mass[i];
  if(imass<200)
  {
  tpep.mass[tpep.mods[mod_id]] = imass;
  tpep.amino.erase(tpep.amino.begin()+i);
  tpep.mass.erase(tpep.mass.begin()+i);
  for(int j=0;j<tpep.mods.size();j++)if(tpep.mods[j]>i)tpep.mods[j]--;

  }
 }
}
//-----------------------------------------------------------------------------------------
void TFloatMods::SrinkPeptide(TPeptide & tpep)
{
  int j;
  
  //finding middle of the peptide
  int k = (tpep.amino.size() - tpep.mods.size())/2;
  int middle =0;
  for(j=0;(j<tpep.amino.size()&&k>=0);j++)if(tpep.amino[j]==0)middle++; else k--;
  vector<char> from, to;
  for(j=0;j<middle;j++) from.push_back(j);
  //generate a list of orders for mass consuming
  pref_permutations.clear();
  for(j=0;j<middle;j++) generatePermutations( from,  to, j,0);
  //srink their prefix  
  if(middle>0)
  {
  for(j=0;j<pref_permutations.size();j++) SrinkPrefix(tpep,j,0);
    //suffix is done inside   
  }
  else PrepareShrinkSuffix(tpep, 0);
}
//-----------------------------------------------------------------------------------------
inline void TFloatMods::EvaluateRounds(TPeptide & tpep)
  {

 /*cout << "+";
 for(int i = 0;i<amino_mod.size();i++) cout << amino_mod[i];
 cout << endl;
 */
  EvaluateOne(tpep);
  for(int i = 0;i<tpep.amino.size();i++) 
  {
   if(tpep.amino[i]==0 && fabs(tpep.mass[i])<30)
   {
    float imass = tpep.mass[i];
    if(((int)imass)!=0) tpep.mass[i]=((int)imass);
    EvaluateOne(tpep);
    if(((int)(imass+0.95))!=0) tpep.mass[i]=((int)(imass+0.95));
    EvaluateOne(tpep);
    tpep.mass[i] = imass;
   }
  }
 } 
//-----------------------------------------------------------------------------------------
inline void TFloatMods::EvaluateOne(TPeptide & tpep)
{
  vector<float> masses;
  tpep.GenerateCandidateMasses(masses);
  
  double loc_int = ComputeExplainedIntensity(masses); 
  
  pep_list.add(tpep.getSequence(2), loc_int);
  
//  cout << " - " << tpep.getSequence(2) << " " << initial_pep_not_valid << " " << max_intensity << " " << loc_int << endl;
  
  if(initial_pep_not_valid){initial_pep_not_valid=false;max_intensity=loc_int;best_peptide=tpep.getSequence(2); }
    
  if(max_intensity + 0.01 < loc_int) { max_intensity=loc_int; best_peptide=tpep.getSequence(2); }
  else
  {
  int mod1=0,mod2=0;
  
  bool seenmods = false;
  for(int i=0;i<tpep.amino.size();i++)
  {
  if(tpep.amino[i]==0) seenmods = true;
  else
  if(tpep.amino[i]!=0)
  {
   if(seenmods) mod2++;
   seenmods = false;
  } 
  }
  if(seenmods) mod2++;
  
  for(int i=2;i<best_peptide.length();i++){if(best_peptide[i]==']'&&(i==best_peptide.length()-1 || best_peptide[i+1]!='[' )) mod1++;}

  //cout << mod1 << " " << mod2 << " ";

  if(fabs(max_intensity-loc_int)<0.01 && /*s2.length()>s1.length()*/ mod1>mod2) {max_intensity=loc_int;best_peptide=tpep.getSequence(2);}
  }
  
 //cout << " -> " << best_peptide << endl;
 
 
}
//-----------------------------------------------------------------------------------------
void TFloatMods::permute(TPeptide tpep, int mod_id, int left, int right, int fix)
{
  
 if(mod_id==tpep.mods.size()) return;
 int pos = tpep.mods[mod_id];
  
 //set left bound
 int left1 = left;
 for(int i=pos-1;(i>=0 && i>=pos-left1);i--)if(tpep.amino[i]==0)left1++;
 int start = pos-left1;
 if(start>0){
   left1 = fix; 
   for(int i=pos-1;(i>=0 && i>=pos-left1);i--)if(tpep.amino[i]==0)left1++;
   start=pos-left1;
   }
 if(start<0) start=0; 
 
 //set right bound
 int right1 = right;
 for(int i=pos+1;(i<tpep.amino.size() && i<=pos+right1);i++)if(tpep.amino[i]==0)right1++;
 int end = pos+right1;
 if(end<tpep.amino.size()-1)
 {
   right1 = fix; 
   for(int i=pos+1;(i<tpep.amino.size() && i<=pos+right1);i++)if(tpep.amino[i]==0)right1++;
   end=pos+right1;
 }
 if(end>tpep.amino.size()-1) end = tpep.amino.size()-1;

 //possition the modification at the begining   of sliding segment
 char temp = tpep.amino[pos];
 float tmass = tpep.mass[pos];
 for(int i=pos-1;i>=start;i--) {tpep.amino[i+1]=tpep.amino[i];tpep.mass[i+1]=tpep.mass[i];}
 tpep.amino[start]=temp;tpep.mass[start]=tmass;
 for(int i=0;i<tpep.mods.size();i++) if(tpep.mods[i]>=start&&tpep.mods[i]<pos)tpep.mods[i]++;
 tpep.mods[mod_id]=start;
  /*
  string tme="";
  for(int i=0;i<amino_mod.size();i++)tme+=amino_mod[i];
   cout << "\n?? " << pos << " - " << left1 << "," << right1 << " = " << start << " " << tme << endl;
 */
 for(int i=start;i<=end;i++)
 {
  //!!! Notice we are Shrinking only last level of extension
  //check shrinked version of peptide
  //cout << "*" << tpep.getSequence(2) << " " << tpep.ValidateCase3() << " " << mod_id << " " << tpep.mods.size()-1 << endl;
  if(tpep.ValidateCase3() && mod_id == tpep.mods.size()-1) extended_peptides.push_back(tpep.getSequence(2));
   
  //Evaluate any acceptable configuration 
 if(tpep.ValidPeptide() && mod_id != tpep.mods.size()) EvaluateRounds(tpep); 
  /*
  tme="";
  for(int h=0;h<amino_mod.size();h++)tme+=amino_mod[h];
  cout << "\n$$ ---" << pos << " - " << left1 << "," << right1 << " = " << start << " " << tme << endl;
*/
  permute(tpep, mod_id+1, left, right, fix);

  while(i<=end)
  {
  i++;
  if(tpep.mods[mod_id]>=tpep.amino.size()-1){i=end+2;break;}
  if(tpep.amino[tpep.mods[mod_id]+1]==0)
    for(int t=0;t<tpep.mods.size();t++) if(tpep.mods[t]==tpep.mods[mod_id]+1) { tpep.mods[t]--; break;}
  tpep.mods[mod_id]++;
  
  char temp = tpep.amino[tpep.mods[mod_id]-1];
  tpep.amino[tpep.mods[mod_id]-1] = tpep.amino[tpep.mods[mod_id]];
  tpep.amino[tpep.mods[mod_id]] = temp;
  
  float tmass = tpep.mass[tpep.mods[mod_id]-1];
  tpep.mass[tpep.mods[mod_id]-1] = tpep.mass[tpep.mods[mod_id]];
  tpep.mass[tpep.mods[mod_id]] = tmass;

  if(tpep.amino[tpep.mods[mod_id]-1]!=0) {i--;break;}
  }
 
  }
}
//--------------------------------------------------------------
void TPeptide::GenerateMixedMasses(vector<float> & masses)
{
  float imass = 0;
  for(int i=0;i<amino.size();i++)
  {
   if(amino[i]==0) imass+=mass[i];
   else
   {
   if(imass!=0) { masses.push_back(imass); }
   imass = mass[i];
   }
  }
  masses.push_back(imass);
}
//------------------------------------------------------------------------
void TPeptide::GenerateCandidateMasses(vector<float> & masses)
{

  vector<float> pure_masses;
  GenerateMixedMasses(pure_masses);
  
  for(int i=0;i<pure_masses.size();i++)
  {
   double b_mass = 0, y_mass = 0;
   for(int j=0;j<pure_masses.size();j++)
 	{
 	 if(j<i) b_mass+=pure_masses[j]; else y_mass+=pure_masses[j];
 	 }
         if(b_mass>0){
 			for(int z=1;z<3;z++)
 			{
 			masses.push_back((b_mass+(double)z*ms("H+"))/(double(z))+1.00);
 			masses.push_back((b_mass+(double)z*ms("H+"))/(double(z)));
 			masses.push_back((b_mass+(double)z*ms("H+")-ms("H2O"))/(double(z)));
 			masses.push_back((b_mass+(double)z*ms("H+")-ms("NH3"))/(double(z)));
 			//masses.push_back((b_mass+(double)z*ms("H+")-ms("CO"))/(double(z)));
 			}
 		     }
 	if(y_mass>0){
 			for(int z=1;z<3;z++)
 			{	
 			masses.push_back((y_mass+(double)z*ms("H+")+ms("H2O"))/(double(z))+1.00);
 			masses.push_back((y_mass+(double)z*ms("H+")+ms("H2O"))/(double(z)));
			masses.push_back((y_mass+(double)z*ms("H+"))/(double(z)));
			masses.push_back((y_mass+(double)z*ms("H+")+ms("H2O")-ms("NH3"))/(double(z)));
 			}
 		    }
   }
 sort(masses.begin(), masses.end());
}
//------------------------------------------------------------------------

double ms(char c)
{
 switch(c)
 {
   case 'A': return 71.037113787;
   case 'R': return 156.101111026;
   case 'D': return 115.026943031;
   case 'N': return 114.042927446;
   case 'C': return 160.030648200;
   case 'E': return 129.042593095;
   case 'Q': return 128.058577510; 
   case 'G': return 57.021463723;
   case 'H': return 137.058911861;
   case 'I': return 113.084063979;
   case 'L': return 113.084063979;
   case 'K': return 128.094963016;
   case 'M': return 131.040484605;
   case 'F': return 147.068413915;
   case 'P': return 97.052763851;
   case 'S': return 87.032028409;
   case 'T': return 101.047678473; 
   case 'W': return 186.079312952;
   case 'Y': return 163.063328537;
   case 'V': return 99.068413915;
 }

return 0;

}
//---------------------
double ms(string s)
{
if(s.length()==1) return ms(s[0]);
if(s=="H+")  return 1.0072763;
if(s=="H2O") return 18.010564686;
if(s=="NH3") return 17.026549101;
if(s=="CO")  return 27.994914622;

return 0;
/*
Mass(C) = 12
Mass(H)  = 1.007825032
Mass(N) = 14.003074005
Mass(O) = 15.994914622
Mass(S) = 31.97207069
*/
}
//--------------------------------------------------------------


