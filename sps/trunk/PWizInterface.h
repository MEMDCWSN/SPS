////////////////////////////////////////////////////////////////////////////////
#ifndef __PWIZ_INTERFACE_H__
#define __PWIZ_INTERFACE_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>

#include "spectrum.h"
#include "SpecSet.h"
#include "aminoacid.h"
#include "PeptideSpectrumMatchSet.h"

//#if defined(__linux__)
//#if !defined(__MINGW32__)
////////////////////////////////////////////////////////////////////////////////
// This flag is defined in the Makefiles, and is set to 0 if the user types make pwiz=no. This defines the beginning of a
// conditional compilation section, based on the developer input.
#if INCLUDE_PWIZ == 1


// Includes from ProteoWizard

#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"
#include "pwiz/data/identdata/DefaultReaderList.hpp"
#include "pwiz/data/identdata/IdentDataFile.hpp"
//#include "pwiz/analysis/peptideid/PeptideID_pepXML.hpp"
#include "pwiz/data/tradata/TraData.hpp"


// Namespaces defined by proteowizard
using namespace pwiz;
using namespace pwiz::msdata;
using namespace pwiz::cv;
using namespace pwiz::data;
using namespace pwiz::identdata;
using namespace pwiz::util;
//using namespace pwiz::peptideid;
using namespace pwiz::tradata;

#endif
////////////////////////////////////////////////////////////////////////////////
// Namespaces defined by SPS
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
 /*! \brief SearchModData helper class

   Holds the data gathered when searching for mods in the ProteoWizard data structures.
   This class holds information for a single modification, which includes the associated mass shift and the residues involved.

   */
class SearchModData {

 public:

    //! \name CONSTRUCTORS
    //@{

    /*! \brief The exemplar constructor.

     Default contructor
     */
  SearchModData() {};

    /*! \brief constructor

     Accepts mass and a residue as a C++ string
     */
  SearchModData(double m, string s) : massShift(m) {residues.push_back(s);};

    /*! \brief constructor

     Accepts mass and a residue as a C char
     */
  SearchModData(double m, char c) : massShift(m) {string aux;aux=c;residues.push_back(aux);};

    //@}


    /*! \brief mod mass shift
     Holds mass shift information (assotiated with the residue(s) )
     */
  double massShift;

  // residues it applyes to
    /*! \brief residues it applyes to
     Contains the string of residues that the mass shift refers to. Usually only one residue will be found in this data structure.
     */
  vector<string> residues;


    /*! \brief Check for existence
     Checks if the residue already exists in the list
    @param res The residue to look for
    @return Boolean value; true if the residue was found
     */
  bool exists(char res)
  {
    for(int j = 0 ; j < residues.size() ; j++) {
      string aux;
      aux = res;
      if(residues[j] == aux)
        return true;
    }
    return false;
  }

    /*! \brief Check for existence
     Checks if the residue already exists in the list, also given a mass shift
    @param mass The mass value to check
    @param res The residue to look for
    @return Boolean value; true if the residue was found
     */
  bool exists(double &mass, char res)
  {
    if( fabs(mass-massShift) > 0.001)
      return false;
    return exists(res);
  };

    /*! \brief Check for existence
     Matches a list of residues agains the internal list. Also compares with a given mass value
    @param mass The mass value to check
    @param res The residue list to check
    @return Boolean value; true if the residue was found
     */
  bool exists(double &mass, vector<char> &res)
  {
    if( fabs(mass-massShift) > 0.001)
      return false;
    for(int i = 0 ; i < res.size() ; i++)
      if(exists(res[i]))
        return true;
    return false;
  };

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief SearchModsData helper class

   Holds the data gathered when searching for mods in the ProteoWizard data structures.

   */
class SearchModsData {

    /*! \brief Gets aa mass value

    @param R The mass value to check
    @return (double) mass value
     */
  double getMass(char R)
  {
    for(int aaIndex = 0; aaIndex < AAcount; aaIndex++)
      if(AAletters[aaIndex] == R)
        return AAmasses[aaIndex];
    return 0.0;
  }

 public:

    /*! \brief Holds fixed mods
     */

  vector<SearchModData>  fixedMods;

    /*! \brief Holds variable mods
     */
  vector<SearchModData>  variable;

    /*! \brief Check for existence
     Checks if the residue already exists in the list, also given a mass shift
    @param mass The mass value to check
    @param res The residue to look for
    @return Boolean value; true if the residue was found
     */
  bool exists(double &mass, char residue)
  {
    for(int i = 0 ; i < fixedMods.size() ; i++)
      if(fixedMods[i].exists(mass, residue))
        return true;
    return false;
  };

    /*! \brief Search for mod
     Checks if a modification exists on the list. First looks in the fixed mods list and then in the variable mods list.
    @param mass The mass value to check
    @param res The residues list to look for
    @return Boolean value; true if the mod was found
     */
  bool fixed(double &mass, vector<char> &residues)
  {
    for(int i = 0 ; i < fixedMods.size() ; i++)
      if(fixedMods[i].exists(mass, residues))
        return true;
    for(int i = 0 ; i < variable.size() ; i++)
      if(variable[i].exists(mass, residues))
        return false;
    return true;
  };


    /*! \brief Remove mods
     */
  void remove(AAJumps &jumpsStd, AAJumps *jumps)
  {
    vector<SearchModData>  newFixedMods, AB;


    for(int i = 0 ; i < fixedMods.size() ; i++) {
      double A = fixedMods[i].massShift;
      for(int j = 0 ; j < fixedMods[i].residues.size() ; j++) {
        char S = fixedMods[i].residues[j][0];
        double Smass = 0.0, T = 0.0;
        if(jumps)
          jumps->getAAref(S, Smass);
        else
          Smass = getMass(S);
        jumpsStd.getAAref(S, T);
        double B = A + T - Smass;

        //cout << "S: " << S << endl;
        //cout << "A: " << A << endl;
        //cout << "T: " << T << endl;
        //cout << "Sm: " << Smass << endl;
        //cout << "B: " << B << endl;
        //cout << "-----------" << endl;

        if( fabs(B-A) > 0.001) {
          // remove the mass
          fixedMods[i].residues.erase(fixedMods[i].residues.begin() + j);
          if(fixedMods[i].residues.size() == 0) {
            fixedMods.erase(fixedMods.begin() + i);
            i--;
          }
          // create a new one, but first check if it already exists
          if(!exists(B,S) && (fabs(B) > 0.001) ) {
            SearchModData m(B,S);
            newFixedMods.push_back(m);
          }
          // end this for cycle
          break;
        }
      }
    }
    AB.reserve( newFixedMods.size() + fixedMods.size() ); // preallocate memory
    AB.insert( AB.end(), newFixedMods.begin(), newFixedMods.end() );
    AB.insert( AB.end(), fixedMods.begin(), fixedMods.end() );
    fixedMods = AB;
  };

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief SearchModsData helper class to hold pwiz IDdata

   */
struct IdData {

    /*! \brief CHARGE
     */
  int m_charge;

    /*! \brief SCORE
     */
  float m_score;

    /*! \brief scan
     */
  int m_scan;

    /*! \brief protein
     */
  string m_protein;

    /*! \brief spectrum file
     */
  string m_spectrumFile;


    /*! \brief Default constructor
     */
  IdData() : m_scan(-1) {};

};
////////////////////////////////////////////////////////////////////////////////
//#if defined(__linux__)
//#if !defined(__MINGW32__)
#if INCLUDE_PWIZ == 1


////////////////////////////////////////////////////////////////////////////////
 /*! \brief PWizInterface class

  ProteoWizard wrapper class
  */

class PWizInterface {

    /*! \brief convertMSData2SpecSet

    Top level method to convert data from the ProteoWizard's data structure to SPS's SpecSet

    @param msd The MSData object, by reference
    @param specSet The SpecSet object, by reference
    @param msl The MS Level, used when readim msxml files, to filter out elements with a lower then required ms level.
     */
  int convertMSData2SpecSet(MSData &msd, SpecSet &specSet, int msl);

    /*! \brief convertIdentData2psmSet

    Top level method to convert ProteoWizard's identity data to SPS's psmSet.

    @param iddList the pwiz's iddList data structure
    @param psmSet the psmSet data structure
    @param searchModsData the structure containing the found mods
    @param outputFixed defines if the fixed mods are to be included
     */
  int convertIdentData2psmSet(vector<IdentDataPtr> &iddList, PeptideSpectrumMatchSet &psmSet, vector<SearchModsData> &searchModsData, bool outputFixed);

    /*! \brief populateFromIdent

    Populate parts of psmSet data strcture from a ProteoWizard's IdentData structure.

    @param iddList the pwiz's iddList data structure
    @param psmSet the psmSet data structure
    @param searchModsData the structure containing the found mods
    @param outputFixed defines if the fixed mods are to be included
     */
  int populateFromIdent(PeptideSpectrumMatchSet &psmSet, IdentData& idd, SearchModsData &searchModsData, bool outputFixed);

    /*! \brief populateFromSpectrumIdentification

    Populate parts of psmSet data structure from a ProteoWizard's SpectrumIdentificationList structure.

    @param psmSet the psmSet data structure
    @param searchModsData the structure containing the found mods
    @param SpectrumIdentificationListPtr the pwiz's ID list data structure
    @param outputFixed defines if the fixed mods are to be included
     */
  int populateFromSpectrumIdentification(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationListPtr &spectrumIdentificationListPtr, bool outputFixed);

    /*! \brief populateSpectrumIdentificationResult

    Populate parts of psmSet data structure from a ProteoWizard's spectrumIdentificationResult structure.

    @param psmSet the psmSet data structure
    @param searchModsData the structure containing the found mods
    @param spectrumIdentificationResult the pwiz's ID list data structure
    @param outputFixed defines if the fixed mods are to be included
     */
  int populateSpectrumIdentificationResult(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationResult & spectrumIdentificationResult, bool outputFixed);

    /*! \brief populateSpectrumIdentificationItem

    Populate parts of psmSet data structure from a ProteoWizard's SpectrumIdentificationItem structure.

    @param psmSet the psmSet data structure
    @param searchModsData the structure containing the found mods
    @param spectrumIdentificationItem the pwiz's ID list data structure
    @param outputFixed defines if the fixed mods are to be included
     */
  int populateSpectrumIdentificationItem(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationItem &spectrumIdentificationItem, IdData &data, bool outputFixed);

    /*! \brief populatePeptideEvidence

    Populate parts of psmSet data structure from a ProteoWizard's PeptideEvidence structure.

    @param psmSet the psmSet data structure
    @param peptideEvidence the structure containing the found mods
    @param spectrumIdentificationItem the pwiz's ID list data structure
    @param outputFixed defines if the fixed mods are to be included
     */
  int populatePeptideEvidence(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, PeptideEvidence &peptideEvidence, IdData &data, bool outputFixed);

    /*! \brief joinSequence

    Joins the sequence parts that compose the final peptide sequence.

    @param vec vector containing the sequence bits
    @param str Final peptide sequence
     */
  void joinSequence(vector<string> &vec, string &str);

    /*! \brief splitSequence

    Splits the peptide sequence into sequence parts.

    @param vec vector containing the sequence bits
    @param str Initial peptide sequence
     */
  void splitSequence(vector<string> &vec, string &str);

    /*! \brief populateModification

    Populates the modifications internal structure from PWiz's data structura that contains info about modifications.

    @param psmSet the psmSet data structure
    @param modification the structure containing the found mods
    @param spectrumIdentificationItem the pwiz's ID list data structure
    @param outputFixed defines if the fixed mods are to be included
     */
  int populateModification(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, pwiz::identdata::Modification & modification, vector<string> &vec, bool outputFixed);

    /*! \brief acquireSearchModification

    Cycles though all the modification sub structures and processes each one individually

    @param search the helper SearchModsData data structure to be filled out
    @param protocols the structure containing the found mods
     */
  int acquireSearchModification(SearchModsData &search, std::vector<SpectrumIdentificationProtocolPtr> &protocols);

  //int readPeptideData(PeptideID_pepXml &peptide, const string &fn);
    /*! \brief getScan

    Finds the scan number int the data string

    @param str Strign to search
    @param scan Where to return the scan #
     */
  bool getScan(string &str, int &scan, char *scanStr);

    /*! \brief findIndex

    Finds the index of a given attribute in the SpectrumPtr data structure

    @param spectrum Attribute list
    @param id Attribute to look for
     */
  int findIndex(SpectrumPtr &spectrum, int id);
  //
  //int getPeptideIDSize(PeptideID &peptide);


  //////////////////////////////////////////////////////////////////////////////
  // Structure dump to file
  //////////////////////////////////////////////////////////////////////////////

  void dumpMSData(MSData &msd, int msl);
  //
  void dumpIdentVector(stringstream &ss, string p, vector<IdentDataPtr> &iddList);
  //
  void dumpIdent(stringstream &ss, string p, IdentDataPtr &idd);
  //
  void dumpIdentifiable(stringstream &ss, string p, Identifiable &identifiable);
  //
  void dumpCV(stringstream &ss, string p, std::vector<CV> &cvs);
  //
  void dumpAnalysisSoftwareVector(stringstream &ss, string p, std::vector<AnalysisSoftwarePtr> &analysisSoftwareList);
  //
  void dumpAnalysisSoftware(stringstream &ss, string p, AnalysisSoftwarePtr &analysisSoftwareList);
  //
  void dumpProvider(stringstream &ss, string p, Provider &provider);
  //
  void dumpContactVector(stringstream &ss, string p, std::vector<pwiz::identdata::ContactPtr> &contact);
  //
  void dumpContact(stringstream &ss, string p, pwiz::identdata::ContactPtr &contact);
  //
  void dumpAnalysisSampleCollection(stringstream &ss, string p,  AnalysisSampleCollection &analysisSampleCollection);
  //
  void dumpSequenceCollection(stringstream &ss, string p,  SequenceCollection &sequenceCollection);
  //
  void dumpAnalysisCollection(stringstream &ss, string p, AnalysisCollection &analysisCollection);
  //
  void dumpAnalysisProtocolCollection(stringstream &ss, string p, AnalysisProtocolCollection &analysisProtocolCollection);
  //
  void dumpDataCollection(stringstream &ss, string p, DataCollection &dataCollection);
  //
  void dumpBibliographicReferenceVector(stringstream &ss, string p, std::vector<BibliographicReferencePtr> &bibliographicReference);
  //
  void dumpBibliographicReference(stringstream &ss, string p, BibliographicReferencePtr &bibliographicReference);
  //
  void dumpDBSequencesVector(stringstream &ss, string p, std::vector<DBSequencePtr> &dbSequences);
  //
  void dumpDBSequence(stringstream &ss, string p, DBSequencePtr &dbSequences);
  //
  void dumpPeptides(stringstream &ss, string p, std::vector< pwiz::identdata::PeptidePtr> &peptides);
  //
  void dumpPeptide(stringstream &ss, string p,  pwiz::identdata::PeptidePtr &peptide);
  //
  void dumpPeptideEvidenceVector(stringstream &ss, string p, std::vector<PeptideEvidencePtr> &peptideEvidence);
  //
  void dumpPeptideEvidence(stringstream &ss, string p, PeptideEvidencePtr &peptideEvidence);
  //
  void dumpSpectrumIdentificationVector(stringstream &ss, string p, std::vector<SpectrumIdentificationPtr> &spectrumIdentification);
  //
  void dumpSpectrumIdentification(stringstream &ss, string p, SpectrumIdentificationPtr &spectrumIdentification);
  //
  void dumpProteinDetectionProtocol(stringstream &ss, string p, std::vector<ProteinDetectionProtocolPtr> &proteinDetectionProtocol);
  //
  void dumpProteinDetectionList(stringstream &ss, string p, ProteinDetectionList &proteinDetectionList);
  //
  void dumpProteinDetection(stringstream &ss, string p, ProteinDetection &proteinDetection);
  //
  void dumpSpectrumIdentificationProtocolVector(stringstream &ss, string p, std::vector<SpectrumIdentificationProtocolPtr> &spectrumIdentificationProtocol);
  //
  void dumpSpectrumIdentificationProtocol(stringstream &ss, string p, SpectrumIdentificationProtocolPtr &spectrumIdentificationProtocol);
  //
  void dumpInputs(stringstream &ss, string p, identdata::Inputs &inputs);
  //
  void dumpSourceFile(stringstream &ss, string p, std::vector<identdata::SourceFilePtr> &sourceFile);
  //
  void dumpSearchDatabaseVector(stringstream &ss, string p, std::vector<SearchDatabasePtr> &searchDatabase);
  //
  void dumpSearchDatabase(stringstream &ss, string p, SearchDatabasePtr &searchDatabase);
  //
  void dumpSpectraDataVector(stringstream &ss, string p, std::vector<SpectraDataPtr> &spectraData);
  //
  void dumpSpectraData(stringstream &ss, string p, SpectraDataPtr &spectraData);
  //
  void dumpAnalysisData(stringstream &ss, string p, AnalysisData &analysisData);
  //
  void dumpContactRolePtr(stringstream &ss, string p, ContactRolePtr &contactRolePtr);
  //
  void dumpCVParamsVector(stringstream &ss, string p, std::vector<CVParam> &cvParams);
  //
  void dumpCVParam(stringstream &ss, string p, CVParam &CvParams);
  //
  void dumpIdentifiableParamContainer(stringstream &ss, string p, IdentifiableParamContainer &identifiableParamContainer);
  //
  void dumpParamContainer(stringstream &ss, string p, ParamContainer &paramContainer);
  //
  void dumpParamGroup(stringstream &ss, string p, std::vector<ParamGroupPtr> &paramGroupPtrs);
  //
  void dumpUserParam(stringstream &ss, string p, std::vector<UserParam> &userParams);
  //
  void dumpModificationsVector(stringstream &ss, string p, std::vector<pwiz::identdata::ModificationPtr> &modifications);
  //
  void dumpModification(stringstream &ss, string p, pwiz::identdata::ModificationPtr &modification);
  //
  void dumpProteinsVector(stringstream &ss, string p, std::vector<ProteinPtr> &proteinPtrs);
  //
  void dumpProtein(stringstream &ss, string p, ProteinPtr &proteinPtrs);
  //
  void dumpRetentionTimesVector(stringstream &ss, string p, std::vector<RetentionTime> &retentionTimes);
  //
  void dumpRetentionTime(stringstream &ss, string p, RetentionTime &retentionTime);
  //
  void dumpEvidence(stringstream &ss, string p, Evidence &evidence);
  //
  void dumpSubstitutionModificationsVector(stringstream &ss, string p,  vector<SubstitutionModificationPtr> &substitutionModification);
  //
  void dumpSubstitutionModification(stringstream &ss, string p,  SubstitutionModificationPtr &substitutionModification);
  //
  void dumpSample(stringstream &ss, string p, pwiz::identdata::SamplePtr &sample);
  //
  void dumpContactRoleVector(stringstream &ss, string p, std::vector<ContactRolePtr> &contactRoles);
  //
  void dumpSoftware(stringstream &ss, string p, pwiz::tradata::SoftwarePtr &software);
  //
  void dumpTranslationTable(stringstream &ss, string p, TranslationTablePtr &translationTablePtr);
  //
  void dumpSpectrumIdentificationList(stringstream &ss, string p, SpectrumIdentificationListPtr &spectrumIdentificationList);
  //
  void dumpEnzymes(stringstream &ss, string p, Enzymes &enzymes);
  //
  void dumpEnzyme(stringstream &ss, string p, EnzymePtr &enzyme);
  //
  void dumpDatabaseTranslation(stringstream &ss, string p, DatabaseTranslationPtr &databaseTranslation);
  //
  void dumpSearchModificationVector(stringstream &ss, string p, std::vector<SearchModificationPtr> &searchModifications);
  //
  void dumpSearchModification(stringstream &ss, string p, SearchModificationPtr &searchModifications);
  //
  void dumpMassTableVector(stringstream &ss, string p, std::vector<MassTablePtr> &massTables);
  //
  void dumpMassTable(stringstream &ss, string p, MassTablePtr &massTable);
  //
  void dumpFilterVector(stringstream &ss, string p, std::vector<FilterPtr> &filters);
  //
  void dumpFilter(stringstream &ss, string p, FilterPtr &filter);
  //
  void dumpResiduesVector(stringstream &ss, string p, std::vector<ResiduePtr> &residues);
  //
  void dumpResidue(stringstream &ss, string p, ResiduePtr &residue);
  //
  void dumpAmbiguousResiduesVector(stringstream &ss, string p, std::vector<AmbiguousResiduePtr> &ambiguousResidues);
  //
  void dumpAmbiguousResidue(stringstream &ss, string p, AmbiguousResiduePtr &ambiguousResidue);
  //
  void dumpTranslationTableVector(stringstream &ss, string p, std::vector<TranslationTablePtr> &translationTables);
  //
  void dumpMeasureVector(stringstream &ss, string p, std::vector<MeasurePtr> &measure);
  //
  void dumpMeasure(stringstream &ss, string p, MeasurePtr &measure);
  //
  void dumpSpectrumIdentificationItemVector(stringstream &ss, string p, std::vector<SpectrumIdentificationItemPtr> &spectrumIdentificationItems);
  //
  void dumpSpectrumIdentificationItem(stringstream &ss, string p, SpectrumIdentificationItemPtr &spectrumIdentificationItem);
  //
  void dumpIonTypeVector(stringstream &ss, string p, std::vector<IonTypePtr> &ionType);
  //
  void dumpIonType(stringstream &ss, string p, IonTypePtr &ionType);
  //
  void dumpFragmentArrayVector(stringstream &ss, string p, std::vector<FragmentArrayPtr> &fragmentArray);
  //
  void dumpFragmentArray(stringstream &ss, string p, FragmentArrayPtr &fragmentArray);
  //
  void dumpSpectrumIdentificationResult(stringstream &ss, string p, SpectrumIdentificationResultPtr &spectrumIdentificationResult);
  //
  void dumpSpectrumIdentificationResultVector(stringstream &ss, string p, std::vector<SpectrumIdentificationResultPtr> &spectrumIdentificationResult);
  //
  //void dumpPeptideID(stringstream &ss, string p, PeptideID &peptide);

  //
  string getCVParamName(int id);

  double findpercursormz(pwiz::msdata::SpectrumPtr &spectrum);

  double findPercursorIntensity(SpectrumPtr &spectrum);

  double findRetentionTime(SpectrumPtr &spectrum);

  //////////////////////////////////////////////////////////////////////////////
  // Data members
  //////////////////////////////////////////////////////////////////////////////

  // iddList
  vector<IdentDataPtr> iddList;

  string inputFilename;


 public:

    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Default contructor
     */
  PWizInterface()  {};
    //@}

    //! \name DESTRUCTOR
    //@{
  ~PWizInterface() {};
    //@}

  /*! \brief acquaireMods

   Gets the list of modifications from the ProteoWzard data structure

   */
  int acquaireMods(vector<SearchModsData> &searchModsData);

  // file open
  /*! \brief openIdent

   Gets information from an identity file

   */
  bool openIdent(string &inName);

  /*! \brief dump

   Dumps the contents of the proteowizard data structures to the screen

   */
  bool dump(void) {stringstream ss;string empty;dumpIdentVector(ss, empty, iddList);cout << ss.str() << endl;};

  /*! \brief loadDataUsingPWiz

   Load data using the ProteoWizard libraries. This method loads specifically spectra data

   */
  int loadDataUsingPWiz(const string &inName, SpecSet &specs, int msl);

  /*! \brief loadDataUsingPWiz

   Load data using the ProteoWizard libraries. This method loads specifically identity data

   */
  int loadDataUsingPWiz(PeptideSpectrumMatchSet &psmSet, vector<SearchModsData> &searchModsData, bool outputFixed);

};
////////////////////////////////////////////////////////////////////////////////
#else
////////////////////////////////////////////////////////////////////////////////
class PWizInterface {


 public:
  // constructors and destructor
  PWizInterface()  {};
  ~PWizInterface() {};

  // file open
  /*! \brief openIdent

   Gets information from an identity file

   */
  bool openIdent(string &inName)  {return false;};

  int acquaireMods(vector<SearchModsData> &searchModsData) {return false;};

  /*! \brief dump

   Dumps the contents of the proteowizard data structures to the screen

   */
  bool dump(void) {};

  /*! \brief loadDataUsingPWiz

   Load data using the ProteoWizard libraries. This method loads specifically spectra data

   */
  int loadDataUsingPWiz(const string &inName, SpecSet &specs, int msl) {return -2;};

  /*! \brief loadDataUsingPWiz

   Load data using the ProteoWizard libraries. This method loads specifically identity data

   */
  int loadDataUsingPWiz(PeptideSpectrumMatchSet &psmSet, vector<SearchModsData> &searchModsData, bool outputFixed)  {return -2;};
};
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////

