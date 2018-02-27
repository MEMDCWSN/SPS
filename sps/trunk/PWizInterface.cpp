////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "PWizInterface.h"

////////////////////////////////////////////////////////////////////////////////

//#if defined(__linux__)
//#if !defined(__MINGW32__)
#if INCLUDE_PWIZ == 1


#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"
#include "pwiz/data/identdata/DefaultReaderList.hpp"
#include "pwiz/data/identdata/IdentDataFile.hpp"
//#include "pwiz/analysis/peptideid/PeptideID_pepXML.hpp"
#include "pwiz/data/tradata/TraData.hpp"

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
using namespace std;
using namespace specnets;

////////////////////////////////////////////////////////////////////////////////
//#if defined(__linux__)
//#if !defined(__MINGW32__)
#if INCLUDE_PWIZ == 1


////////////////////////////////////////////////////////////////////////////////
/*int PWizInterface::getPeptideIDSize(PeptideID &peptide)
{
  int peptideSize = 0;
  PeptideID_pepXml::Iterator it = peptide.begin();
  for( ; it != peptide.end() ; it++)
    peptideSize++;
  return peptideSize;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpPeptideID(stringstream &ss, string p, PeptideID &peptide)
{
  PeptideID_pepXml::Iterator it = peptide.begin();
  int peptideSize = 0;
  for( ; it != peptide.end() ; it++)   {
    ss << p << "nativeID:         " << it->nativeID << endl;
    ss << p << "sequence:         " << it->sequence << endl;
    ss << p << "protein_descr:    " << it->protein_descr << endl;
    ss << p << "mz:               " << it->mz << endl;
    ss << p << "retentionTimeSec: " << it->retentionTimeSec << endl;
    ss << p << "normalizedScore:  " << it->normalizedScore << endl;

    peptideSize++;
  }
  ss << p << "peptide size : " << peptideSize << endl;
}*/
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpMSData(MSData &msd, int msl)
{
  SpectrumList& spectrumList = *msd.run.spectrumListPtr;
  SpectrumPtr spectrum;
  //Read m/z and intensity values from the spectra
  const bool getBinaryData = true;
  size_t numSpectra = spectrumList.size();

  for (int i = 0 ; i < numSpectra ; ++i) {
    // define specnets spectrum object
    specnets::Spectrum spect;

    ///////// Get the peak list ///////
    // Get spectrum binary data
  	spectrum = spectrumList.spectrum(i, getBinaryData);
  	// Define pairs data holder
  	vector<MZIntensityPair> pairs;
  	// Get MZ intensity pairs
  	spectrum->getMZIntensityPairs(pairs);
    // Add the peaks
  	for (vector<MZIntensityPair>::const_iterator it = pairs.begin(), end = pairs.end(); it!=end; ++it) {
      // add a peak
      cout << it->mz << " ; " << it->intensity << endl;
    }


    ///////// Data ///////
    int index;

    // MS Level
    if((index = findIndex(spectrum, MS_ms_level)) >= 0) {
      cout << "MS level: " << spectrum->cvParams[index].value.c_str() << endl;
    }

    // parent M/Z
    if((index = findIndex(spectrum, MS_base_peak_m_z)) >= 0)
      cout << "Parent MZ: " << spectrum->cvParams[index].value.c_str() << endl;
    else
      if((index = findIndex(spectrum, MS_base_peak)) >= 0)
        cout << "Parent MZ: " << spectrum->cvParams[index].value.c_str() << endl;


    // parent charge
    if((index = findIndex(spectrum, MS_charge_state)) >= 0)
      cout << "Parent cbharge: " << spectrum->cvParams[index].value.c_str() << endl;

    // base peak intensity
    if((index = findIndex(spectrum, MS_base_peak_intensity)) >= 0)
      cout << "base peak intensity: " << spectrum->cvParams[index].value.c_str() << endl;

    // scan number
    int scan;
    spect.scan = -1;
    if(getScan(spectrum->id, scan, "scan"))
       cout << "scan: " << scan << endl;

    // File name
    size_t  found = inputFilename.find_last_of("/\\");
    spect.fileName  = inputFilename.substr(found + 1);

    // parent mass
    float parentMass          =  spect.parentMZ * spect.parentCharge - AAJumps::massHion * (spect.parentCharge - 1.0);
      cout << "parent mass: " << parentMass << endl;

    // MS Fragmentation model
    int msFragType            = specnets::Spectrum::FragType_CID;
    if((index = findIndex(spectrum, MS_electron_transfer_dissociation)) >= 0)
      cout << "Fragmentation model: " << specnets::Spectrum::FragType_ETD << endl;
    if((index = findIndex(spectrum, MS_high_energy_collision_induced_dissociation)) >= 0)
      cout << "Fragmentation model: "  << specnets::Spectrum::FragType_HCD << endl;

  }

}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIdentVector(stringstream &ss, string p, vector<IdentDataPtr> &iddList)
{
  ss << p << "iddlist size : " << iddList.size() << endl;

  for (size_t i=0; i < iddList.size(); ++i) {
    IdentData& idd = *iddList[i];
    dumpIdent(ss, p, iddList[i]);
  }
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIdent(stringstream &ss, string p, IdentDataPtr &idd)
{
  p += "Ident.";

  ss << p << "creationDate : " << idd->creationDate << endl;
  dumpCV(ss, p, idd->cvs);
  dumpAnalysisSoftwareVector(ss, p, idd->analysisSoftwareList);
  dumpProvider(ss, p, idd->provider);
  dumpContactVector(ss, p, idd->auditCollection);
  dumpAnalysisSampleCollection(ss, p, idd->analysisSampleCollection);
  dumpSequenceCollection(ss, p, idd->sequenceCollection);
  dumpAnalysisCollection(ss, p, idd->analysisCollection);
  dumpAnalysisProtocolCollection(ss, p, idd->analysisProtocolCollection);
  dumpDataCollection(ss, p, idd->dataCollection);
  dumpBibliographicReferenceVector(ss, p, idd->bibliographicReference);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIdentifiable(stringstream &ss, string p, Identifiable &identifiable)
{
  p += "Identifiable.";
  ss << p << "id :    " << identifiable.id << endl;
  ss << p << "name :  " << identifiable.name << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpCV(stringstream &ss, string p, std::vector<CV> &cvs)
{
  p += "CV.";
  ss << p << ".size() : " << cvs.size() << endl;
  for(int j = 0 ; j < cvs.size() ; j++)
    ss << p << cvs[j].id << " - " << cvs[j].URI << " - " << cvs[j].fullName << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisSoftwareVector(stringstream &ss, string p, std::vector<AnalysisSoftwarePtr> &analysisSoftwareList)
{
  ss << p << "analysisSoftwareList.size() : " << analysisSoftwareList.size() << endl;
  for(int j = 0 ; j < analysisSoftwareList.size() ; j++)
    dumpAnalysisSoftware(ss, p, analysisSoftwareList[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisSoftware(stringstream &ss, string p, AnalysisSoftwarePtr &analysisSoftware)
{
  if(analysisSoftware == 0) return;
  p +=  "AnalysisSoftware.";
  dumpIdentifiable(ss, p, *(analysisSoftware));
  ss << p << "version :         " << analysisSoftware->version << endl;
  ss << p << "URI :             " << analysisSoftware->URI << endl;
  ss << p << "customizations :  " << analysisSoftware->customizations << endl;
  dumpContactRolePtr(ss, p, analysisSoftware->contactRolePtr);
  dumpParamContainer(ss, p, analysisSoftware->softwareName);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpContactRoleVector(stringstream &ss, string p, std::vector<ContactRolePtr> &contactRoles)
{
  ss << p << "ContactRoles.size() : " << contactRoles.size() << endl;
  for(int j = 0 ; j < contactRoles.size() ; j++)
    dumpContactRolePtr(ss, p, contactRoles[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpContactRolePtr(stringstream &ss, string p, ContactRolePtr &contactRolePtr)
{
  if(contactRolePtr == 0) return;
  p +=  "ContactRole.";
  dumpCVParam(ss, p, *contactRolePtr);
  dumpContact(ss, p, contactRolePtr->contactPtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpCVParamsVector(stringstream &ss, string p, std::vector<CVParam> &cvParams)
{
  ss << p << "cvParams.size() : " << cvParams.size() << endl;
  for(int j = 0 ; j < cvParams.size() ; j++)
    dumpCVParam(ss, p, cvParams[j]);
}
////////////////////////////////////////////////////////////////////////////////
string PWizInterface::getCVParamName(int id)
{
  stringstream ret;
  ret << id << " -> ";

  switch(id) {
  case -1:
    ret << "CVID_Unknown";
    break;
  case 1001460:
    ret << "unknown modification";
    break;
  case 1001348:
    ret << "MS_FASTA_format";
    break;
  case 1001073:
    ret << "MS_database_type_amino_acid";
    break;
  case 1001088:
    ret << "MS_protein_description";
    break;
  case 1000589:
    ret << "MS_contact_email";
    break;
  case 1000615:
    ret << "MS_ProteoWizard";
    break;
  case 1001208:
    ret << "MS_Sequest";
    break;
  case 1001456:
    ret << "MS_analysis_software";
    break;
  case 1000776:
    ret << "MS_scan_number_only_nativeID_format";
    break;
  case 1001211:
    ret << "MS_parent_mass_type_mono";
    break;
  case 1001256:
    ret << "MS_fragment_mass_type_mono";
    break;
  case 1001251:
    ret << "MS_Trypsin";
    break;
  case 1001121:
    ret << "MS_number_of_matched_peaks";
    break;
  case 1001362:
    ret << "MS_number_of_unmatched_peaks";
    break;
  case 1001155:
    ret << "MS_Sequest_xcorr";
    break;
  case 1001156:
    ret << "MS_Sequest_deltacn";
    break;
  }
  return ret.str();
}

void PWizInterface::dumpCVParam(stringstream &ss, string p, CVParam &CvParam)
{
  //if(CvParam == 0) return;
  p +=  "CVParam.";
  ss << p << "cvid :   " << getCVParamName(CvParam.cvid) << endl;
  ss << p << "value :  " << CvParam.value << endl;
  ss << p << "units :  " << CvParam.units << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpProvider(stringstream &ss, string p, Provider &provider)
{
  p +=  "Provider.";
  dumpIdentifiable(ss, p, provider);
  dumpContactRolePtr(ss, p, provider.contactRolePtr);
  dumpAnalysisSoftware(ss, p, provider.analysisSoftwarePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpContactVector(stringstream &ss, string p, std::vector<pwiz::identdata::ContactPtr> &contact)
{
  ss << p << "Contact.size() : " << contact.size() << endl;
  for(int j = 0 ; j < contact.size() ; j++)
    dumpContact(ss, p, contact[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpContact(stringstream &ss, string p, pwiz::identdata::ContactPtr &contact)
{
  if(contact == 0) return;
  p += "Contact.";
  dumpIdentifiableParamContainer(ss, p, *contact);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisSampleCollection(stringstream &ss, string p,  AnalysisSampleCollection &analysisSampleCollection)
{
  p += "AnalysisSampleCollection.";
  ss << p << "samples.size() : " << analysisSampleCollection.samples.size() << endl;
  for(int j = 0 ; j < analysisSampleCollection.samples.size() ; j++)
    dumpSample(ss, p,  analysisSampleCollection.samples[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSample(stringstream &ss, string p, pwiz::identdata::SamplePtr &sample)
{
  if(sample == 0) return;
  p += "Sample.";
  dumpIdentifiableParamContainer(ss, p, *sample);
  dumpContactRoleVector(ss, p, sample->contactRole);
  for(int j = 0 ; j < sample->subSamples.size() ; j++)
    dumpSample(ss, p,  sample->subSamples[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSequenceCollection(stringstream &ss, string p,  SequenceCollection &sequenceCollection)
{
  p += "SequenceCollection.";
  dumpDBSequencesVector(ss, p, sequenceCollection.dbSequences);
  dumpPeptides(ss, p, sequenceCollection.peptides);
  dumpPeptideEvidenceVector(ss, p, sequenceCollection.peptideEvidence);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisCollection(stringstream &ss, string p, AnalysisCollection &analysisCollection)
{
  p += "AnalysisCollection.";
  dumpSpectrumIdentificationVector(ss, p, analysisCollection.spectrumIdentification);
  dumpProteinDetection(ss, p, analysisCollection.proteinDetection);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisProtocolCollection(stringstream &ss, string p, AnalysisProtocolCollection &analysisProtocolCollection)
{
  p += "AnalysisProtocolCollection.";
  dumpSpectrumIdentificationProtocolVector(ss, p, analysisProtocolCollection.spectrumIdentificationProtocol);
  dumpProteinDetectionProtocol(ss, p, analysisProtocolCollection.proteinDetectionProtocol);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpDataCollection(stringstream &ss, string p, DataCollection &dataCollection)
{
  p += "DataCollection.";
  dumpInputs(ss, p, dataCollection.inputs);
  dumpAnalysisData(ss, p, dataCollection.analysisData);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpBibliographicReferenceVector(stringstream &ss, string p, std::vector<BibliographicReferencePtr> &bibliographicReferences)
{
  ss << p << "bibliographicReferences.size() : " << bibliographicReferences.size() << endl;
  for(int j = 0 ; j < bibliographicReferences.size() ; j++)
    dumpBibliographicReference(ss, p,  bibliographicReferences[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpBibliographicReference(stringstream &ss, string p, BibliographicReferencePtr &bibliographicReference)
{
  p += "BibliographicReference.";
  dumpIdentifiable(ss, p, *bibliographicReference);
  ss << p << "authors :     " << bibliographicReference->authors << endl;
  ss << p << "publication : " << bibliographicReference->publication << endl;
  ss << p << "publisher :   " << bibliographicReference->publisher << endl;
  ss << p << "editor :      " << bibliographicReference->editor << endl;
  ss << p << "year :        " << bibliographicReference->year << endl;
  ss << p << "volume :      " << bibliographicReference->volume << endl;
  ss << p << "issue :       " << bibliographicReference->issue << endl;
  ss << p << "pages :       " << bibliographicReference->pages << endl;
  ss << p << "title :       " << bibliographicReference->title << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpDBSequencesVector(stringstream &ss, string p, std::vector<DBSequencePtr> &dbSequences)
{
  ss << p << "dbSequences.size() : " << dbSequences.size() << endl;
  for(int j = 0 ; j < dbSequences.size() ; j++)
    dumpDBSequence(ss, p,  dbSequences[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpDBSequence(stringstream &ss, string p, DBSequencePtr &DbSequence)
{
  p += "DBSequence.";
  dumpIdentifiableParamContainer(ss, p, *DbSequence);
  ss << p << "length :    " << DbSequence->length << endl;
  ss << p << "accession : " << DbSequence->accession << endl;
  ss << p << "seq :       " << DbSequence->seq << endl;
  dumpSearchDatabase(ss, p, DbSequence->searchDatabasePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpPeptides(stringstream &ss, string p, std::vector< pwiz::identdata::PeptidePtr> &peptides)
{
  ss << p << "peptides.size() : " << peptides.size() << endl;
  for(int j = 0 ; j < peptides.size() ; j++)
    dumpPeptide(ss, p, peptides[j]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpPeptide(stringstream &ss, string p,  pwiz::identdata::PeptidePtr &peptide)
{
  p += "Peptide.";
  dumpIdentifiableParamContainer(ss, p, *peptide);
  ss << p << "peptideSequence : " << peptide->peptideSequence << endl;
  dumpModificationsVector(ss, p, peptide->modification);
  dumpSubstitutionModificationsVector(ss, p, peptide->substitutionModification);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSubstitutionModificationsVector(stringstream &ss, string p,  vector<SubstitutionModificationPtr> &substitutionModification)
{
  ss << p << "substitutionModification.size() : " << substitutionModification.size() << endl;
  for(int i = 0 ; i < substitutionModification.size() ; i++)
    dumpSubstitutionModification(ss, p, substitutionModification[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSubstitutionModification(stringstream &ss, string p,  SubstitutionModificationPtr &substitutionModification)
{
  p += "SubstitutionModification.";
  ss << p << "originalResidue       " << substitutionModification->originalResidue << endl;
  ss << p << "replacementResidue    " << substitutionModification->replacementResidue << endl;
  ss << p << "location              " << substitutionModification->location << endl;
  ss << p << "avgMassDelta          " << substitutionModification->avgMassDelta << endl;
  ss << p << "monoisotopicMassDelta " << substitutionModification->monoisotopicMassDelta << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpModificationsVector(stringstream &ss, string p, std::vector<pwiz::identdata::ModificationPtr> &modifications)
{
  ss << p << "modifications.size() : " << modifications.size() << endl;
  for(int i = 0 ; i < modifications.size() ; i++)
    dumpModification(ss, p, modifications[i]);
}

void PWizInterface::dumpModification(stringstream &ss, string p, pwiz::identdata::ModificationPtr &modification)
{
  p += "Modification.";
  dumpParamContainer(ss, p, *modification);
  ss << p << "location :              " << modification->location << endl;
  ss << p << "avgMassDelta :          " << modification->avgMassDelta << endl;
  ss << p << "monoisotopicMassDelta : " << modification->monoisotopicMassDelta << endl;
  ss << p << "residues[j] : ";
  for(int j = 0 ; j < modification->residues.size() ; j++)
    ss << modification->residues[j] << " ; ";
  ss << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpProteinsVector(stringstream &ss, string p, std::vector<ProteinPtr> &proteinPtrs)
{
  ss << p << "proteins.size() : " << proteinPtrs.size() << endl;
  for(int i = 0 ; i < proteinPtrs.size() ; i++)
    dumpProtein(ss, p, proteinPtrs[i]);
}

void PWizInterface::dumpProtein(stringstream &ss, string p, ProteinPtr &proteinPtrs)
{
  p += "Protein.";
  dumpParamContainer(ss, p, *proteinPtrs);
  ss << p << "id :       " << proteinPtrs->id << endl;
  ss << p << "sequence : " << proteinPtrs->sequence << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpRetentionTimesVector(stringstream &ss, string p, std::vector<RetentionTime> &retentionTimes)
{
  ss << p << "retentionTimes.size() : " << retentionTimes.size() << endl;
  for(int i = 0 ; i < retentionTimes.size() ; i++)
    dumpRetentionTime(ss, p, retentionTimes[i]);
}

void PWizInterface::dumpRetentionTime(stringstream &ss, string p, RetentionTime &retentionTime)
{
  p += "RetentionTime.";
  dumpParamContainer(ss, p, retentionTime);
  dumpSoftware(ss, p, retentionTime.softwarePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSoftware(stringstream &ss, string p, pwiz::tradata::SoftwarePtr &software)
{
  p += "software.";
  dumpParamContainer(ss, p, *software);
  ss << p << "id :      " << software->id << endl;
  ss << p << "version : " << software->version << endl;
}
////////////////////////////////////////////////////////////////////////////////
void  PWizInterface::dumpEvidence(stringstream &ss, string p, Evidence &evidence)
{
  p += "Evidence.";
  dumpParamContainer(ss, p, evidence);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpPeptideEvidenceVector(stringstream &ss, string p, std::vector<PeptideEvidencePtr> &peptideEvidence)
{
  ss << p << "PeptideEvidence.size() : " << peptideEvidence.size() << endl;
  for(int i = 0 ; i < peptideEvidence.size() ; i++)
    dumpPeptideEvidence(ss, p, peptideEvidence[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpPeptideEvidence(stringstream &ss, string p, PeptideEvidencePtr &peptideEvidence)
{
  p += "PeptideEvidence.";
  dumpIdentifiableParamContainer(ss, p, *peptideEvidence);
  ss << p << "start :   " << peptideEvidence->start << endl;
  ss << p << "end :     " << peptideEvidence->end << endl;
  ss << p << "pre :     " << peptideEvidence->pre << endl;
  ss << p << "post :    " << peptideEvidence->post << endl;
  ss << p << "frame :   " << peptideEvidence->frame << endl;
  ss << p << "isDecoy : " << peptideEvidence->isDecoy << endl;
  dumpPeptide(ss, p, peptideEvidence->peptidePtr);
  dumpDBSequence(ss, p, peptideEvidence->dbSequencePtr);
  dumpTranslationTable(ss, p, peptideEvidence->translationTablePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpTranslationTableVector(stringstream &ss, string p, std::vector<TranslationTablePtr> &translationTables)
{
  ss << p << "TranslationTables.size() : " << translationTables.size() << endl;
  for(int i = 0 ; i < translationTables.size() ; i++)
    dumpTranslationTable(ss, p, translationTables[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpTranslationTable(stringstream &ss, string p, TranslationTablePtr &translationTablePtr)
{
  if(translationTablePtr == 0) return;
  p += "TranslationTable.";
  dumpIdentifiableParamContainer(ss, p, *translationTablePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationVector(stringstream &ss, string p, std::vector<SpectrumIdentificationPtr> &spectrumIdentification)
{
  ss << p << "SpectrumIdentification.size() : " << spectrumIdentification.size() << endl;
  for(int i = 0 ; i < spectrumIdentification.size() ; i++)
    dumpSpectrumIdentification(ss, p, spectrumIdentification[i]);
}

void PWizInterface::dumpSpectrumIdentification(stringstream &ss, string p, SpectrumIdentificationPtr &spectrumIdentification)
{
  p += "SpectrumIdentification.";
  dumpIdentifiable(ss, p, *spectrumIdentification);
  ss << p << "activityDate : " << spectrumIdentification->activityDate << endl;
  dumpSpectrumIdentificationProtocol(ss, p, spectrumIdentification->spectrumIdentificationProtocolPtr);
  dumpSpectrumIdentificationList(ss, p, spectrumIdentification->spectrumIdentificationListPtr);
  dumpSpectraDataVector(ss, p, spectrumIdentification->inputSpectra);
  dumpSearchDatabaseVector(ss, p, spectrumIdentification->searchDatabase);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpProteinDetectionProtocol(stringstream &ss, string p, std::vector<ProteinDetectionProtocolPtr> &proteinDetectionProtocol)
{
  p += "ProteinDetectionProtocol.";
  //ss << p << "proteinDetectionProtocolPtr : " << proteinDetectionProtocolPtr << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpProteinDetectionList(stringstream &ss, string p, ProteinDetectionList &proteinDetectionList)
{
  p += "ProteinDetectionList.";
  dumpIdentifiableParamContainer(ss, p, proteinDetectionList);
  //dumpProteinAmbiguityGroup(ss, p, proteinDetectionList.proteinAmbiguityGroup);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpProteinDetection(stringstream &ss, string p, ProteinDetection &proteinDetection)
{
  p += "ProteinDetection.";
  ss << p << "activityDate : " << proteinDetection.activityDate << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationProtocolVector(stringstream &ss, string p, std::vector<SpectrumIdentificationProtocolPtr> &spectrumIdentificationProtocol)
{
  ss << p << "SpectrumIdentificationProtocol.size() : " << spectrumIdentificationProtocol.size() << endl;
  for(int i = 0 ; i < spectrumIdentificationProtocol.size() ; i++)
    dumpSpectrumIdentificationProtocol(ss, p, spectrumIdentificationProtocol[i]);
}

void PWizInterface::dumpSpectrumIdentificationProtocol(stringstream &ss, string p, SpectrumIdentificationProtocolPtr &spectrumIdentificationProtocol)
{
  p += "SpectrumIdentificationProtocol.";
  dumpIdentifiable(ss, p, *spectrumIdentificationProtocol);
  dumpAnalysisSoftware(ss, p, spectrumIdentificationProtocol->analysisSoftwarePtr);
  dumpCVParam(ss, p, spectrumIdentificationProtocol->searchType);
  dumpParamContainer(ss, p, spectrumIdentificationProtocol->additionalSearchParams);
  dumpEnzymes(ss, p, spectrumIdentificationProtocol->enzymes);
  dumpParamContainer(ss, p, spectrumIdentificationProtocol->fragmentTolerance);
  dumpParamContainer(ss, p, spectrumIdentificationProtocol->parentTolerance);
  dumpParamContainer(ss, p, spectrumIdentificationProtocol->threshold);
  dumpDatabaseTranslation(ss, p, spectrumIdentificationProtocol->databaseTranslation);

  dumpSearchModificationVector(ss, p, spectrumIdentificationProtocol->modificationParams);
  dumpMassTableVector(ss, p, spectrumIdentificationProtocol->massTable);
  dumpFilterVector(ss, p, spectrumIdentificationProtocol->databaseFilters);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpFilterVector(stringstream &ss, string p, std::vector<FilterPtr> &filters)
{
  ss << p << "Filter.size() : " << filters.size() << endl;
  for(int i = 0 ; i < filters.size() ; i++)
    dumpFilter(ss, p, filters[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpFilter(stringstream &ss, string p, FilterPtr &filter)
{
  p += "Filter.";
  dumpParamContainer(ss, p, filter->filterType);
  dumpParamContainer(ss, p, filter->include);
  dumpParamContainer(ss, p, filter->exclude);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpMassTableVector(stringstream &ss, string p, std::vector<MassTablePtr> &massTables)
{
  ss << p << "MassTables.size() : " << massTables.size() << endl;
  for(int i = 0 ; i < massTables.size() ; i++)
    dumpMassTable(ss, p, massTables[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpMassTable(stringstream &ss, string p, MassTablePtr &massTable)
{
  if(massTable == 0) return;
  p += "MassTables.";
  ss << p << "id : " << massTable->id << endl;
  ss << p << "msLevel[j] : ";
  for(int j = 0 ; j < massTable->msLevel.size() ; j++)
    ss << massTable->msLevel[j] << " ; ";
  ss << endl;
  dumpResiduesVector(ss, p, massTable->residues);
  dumpAmbiguousResiduesVector(ss, p, massTable->ambiguousResidue);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAmbiguousResiduesVector(stringstream &ss, string p, std::vector<AmbiguousResiduePtr> &ambiguousResidues)
{
  ss << p << "residues.size() : " << ambiguousResidues.size() << endl;
  for(int i = 0 ; i < ambiguousResidues.size() ; i++)
    dumpAmbiguousResidue(ss, p, ambiguousResidues[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAmbiguousResidue(stringstream &ss, string p, AmbiguousResiduePtr &ambiguousResidue)
{
  p += "AmbiguousResidue.";
  ss << p << "code : " << ambiguousResidue->code << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpResiduesVector(stringstream &ss, string p, std::vector<ResiduePtr> &residues)
{
  ss << p << "residues.size() : " << residues.size() << endl;
  for(int i = 0 ; i < residues.size() ; i++)
    dumpResidue(ss, p, residues[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpResidue(stringstream &ss, string p, ResiduePtr &residue)
{
  p += "Residue.";
  ss << p << "code : " << residue->code << endl;
  ss << p << "mass : " << residue->mass << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSearchModificationVector(stringstream &ss, string p, std::vector<SearchModificationPtr> &searchModifications)
{
  ss << p << "SearchModifications.size() : " << searchModifications.size() << endl;
  for(int i = 0 ; i < searchModifications.size() ; i++)
    dumpSearchModification(ss, p, searchModifications[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSearchModification(stringstream &ss, string p, SearchModificationPtr &searchModifications)
{
  p += "SearchModification.";
  dumpParamContainer(ss, p, *searchModifications);
  ss << p << "fixedMod : " << searchModifications->fixedMod << endl;
  ss << p << "massDelta : " << searchModifications->massDelta << endl;
  ss << p << "searchModifications[j] : ";
  for(int j = 0 ; j < searchModifications->residues.size() ; j++)
    ss << searchModifications->residues[j] << " ; ";
  ss << endl;
  dumpCVParam(ss, p, searchModifications->specificityRules);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpDatabaseTranslation(stringstream &ss, string p, DatabaseTranslationPtr &databaseTranslation)
{
  if(databaseTranslation == 0) return;
  p += "databaseTranslation.";
  ss << p;
  for(int i = 0 ; i < databaseTranslation->frames.size() ; i++)
    ss << databaseTranslation->frames[i] << " ; ";
  ss << endl;

  dumpTranslationTableVector(ss, p, databaseTranslation->translationTable);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpEnzymes(stringstream &ss, string p, Enzymes &enzymes)
{
  ss << p << "enzymes.independent : " << enzymes.independent << endl;
  ss << p << "enzymes.enzymes.size() : " << enzymes.enzymes.size() << endl;
  for(int i = 0 ; i < enzymes.enzymes.size() ; i++)
    dumpEnzyme(ss, p, enzymes.enzymes[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpEnzyme(stringstream &ss, string p, EnzymePtr &enzyme)
{
  p += "Enzyme.";
  dumpIdentifiable(ss, p, *enzyme);
  ss << p << "enzyme.nTermGain : " << enzyme->nTermGain << endl;
  ss << p << "enzyme.cTermGain : " << enzyme->cTermGain << endl;
  ss << p << "enzyme.missedCleavages : " << enzyme->missedCleavages << endl;
  ss << p << "enzyme.minDistance : " << enzyme->minDistance << endl;
  ss << p << "enzyme.siteRegexp : " << enzyme->siteRegexp << endl;
  ss << p << "enzyme.terminalSpecificity : " << enzyme->terminalSpecificity << endl;
  dumpParamContainer(ss, p, enzyme->enzymeName);
//    proteome::Digestion::Specificity terminalSpecificity;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpInputs(stringstream &ss, string p, identdata::Inputs &inputs)
{
  p += "Inputs.";
  dumpSourceFile(ss, p, inputs.sourceFile);
  dumpSearchDatabaseVector(ss, p, inputs.searchDatabase);
  dumpSpectraDataVector(ss, p, inputs.spectraData);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSourceFile(stringstream &ss, string p, std::vector<identdata::SourceFilePtr> &sourceFile)
{
  p += "SourceFile.";
  ss << p << "size() : " << sourceFile.size() << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSearchDatabaseVector(stringstream &ss, string p, std::vector<SearchDatabasePtr> &searchDatabases)
{
  ss << p << "SearchDatabaseVector.size() : " << searchDatabases.size() << endl;
  for(int i = 0 ; i < searchDatabases.size() ; i++)
    dumpSearchDatabase(ss, p, searchDatabases[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSearchDatabase(stringstream &ss, string p, SearchDatabasePtr &searchDatabase)
{
  p += "SearchDatabase.";
  dumpIdentifiableParamContainer(ss, p, *searchDatabase);
  ss << p << "location :              " << searchDatabase->location << endl;
  ss << p << "version :               " << searchDatabase->version << endl;
  ss << p << "releaseDate :           " << searchDatabase->releaseDate << endl;
  ss << p << "numDatabaseSequences :  " << searchDatabase->numDatabaseSequences << endl;
  ss << p << "numResidues :           " << searchDatabase->numResidues << endl;
  dumpCVParam(ss, p, searchDatabase->fileFormat);
  dumpParamContainer(ss, p, searchDatabase->databaseName);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectraDataVector(stringstream &ss, string p, std::vector<SpectraDataPtr> &spectraData)
{
  ss << p << "SpectraData.size() : " << spectraData.size() << endl;
  for(int i = 0 ; i < spectraData.size() ; i++)
    dumpSpectraData(ss, p, spectraData[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectraData(stringstream &ss, string p, SpectraDataPtr &spectraData)
{
  p += "SpectraData.";
  dumpIdentifiable(ss, p, *spectraData);
  ss << p << "location :              " << spectraData->location << endl;
  dumpCVParam(ss, p, spectraData->fileFormat);
  dumpCVParam(ss, p, spectraData->spectrumIDFormat);
  for(int i = 0 ; i < spectraData->externalFormatDocumentation.size() ; i++)
    ss << p << "externalFormatDocumentation : " << spectraData->externalFormatDocumentation[i] << endl;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationList(stringstream &ss, string p, SpectrumIdentificationListPtr &spectrumIdentificationList)
{
  p += "spectrumIdentificationList.";
  dumpIdentifiable(ss, p, *spectrumIdentificationList);
  ss << p << "numSequencesSearched : " << spectrumIdentificationList->numSequencesSearched << endl;
  dumpMeasureVector(ss, p, spectrumIdentificationList->fragmentationTable);
  dumpSpectrumIdentificationResultVector(ss, p, spectrumIdentificationList->spectrumIdentificationResult);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationResultVector(stringstream &ss, string p, std::vector<SpectrumIdentificationResultPtr> &spectrumIdentificationResult)
{
  ss << p << "SpectrumIdentificationItem.size() : " << spectrumIdentificationResult.size() << endl;
  for(int i = 0 ; i < spectrumIdentificationResult.size() ; i++)
    dumpSpectrumIdentificationResult(ss, p, spectrumIdentificationResult[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationResult(stringstream &ss, string p, SpectrumIdentificationResultPtr &spectrumIdentificationResult)
{
  p += "spectrumIdentificationResult.";
  dumpIdentifiableParamContainer(ss, p, *spectrumIdentificationResult);
  ss << p << "spectrumID : " << spectrumIdentificationResult->spectrumID << endl;
  dumpSpectraData(ss, p, spectrumIdentificationResult->spectraDataPtr);
  dumpSpectrumIdentificationItemVector(ss, p, spectrumIdentificationResult->spectrumIdentificationItem);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationItemVector(stringstream &ss, string p, std::vector<SpectrumIdentificationItemPtr> &spectrumIdentificationItems)
{
  ss << p << "SpectrumIdentificationItem.size() : " << spectrumIdentificationItems.size() << endl;
  for(int i = 0 ; i < spectrumIdentificationItems.size() ; i++)
    dumpSpectrumIdentificationItem(ss, p, spectrumIdentificationItems[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpSpectrumIdentificationItem(stringstream &ss, string p, SpectrumIdentificationItemPtr &spectrumIdentificationItem)
{
  p += "SpectrumIdentificationItem.";
  dumpIdentifiableParamContainer(ss, p, *spectrumIdentificationItem);
  ss << p << "chargeState : " << spectrumIdentificationItem->chargeState << endl;
  ss << p << "experimentalMassToCharge : " << spectrumIdentificationItem->experimentalMassToCharge << endl;
  ss << p << "calculatedMassToCharge : " << spectrumIdentificationItem->calculatedMassToCharge << endl;
  ss << p << "rank : " << spectrumIdentificationItem->rank << endl;
  ss << p << "passThreshold : " << spectrumIdentificationItem->passThreshold << endl;
  dumpPeptide(ss, p, spectrumIdentificationItem->peptidePtr);
  dumpMassTable(ss, p, spectrumIdentificationItem->massTablePtr);
  dumpSample(ss, p, spectrumIdentificationItem->samplePtr);
  dumpPeptideEvidenceVector(ss, p, spectrumIdentificationItem->peptideEvidencePtr);
  dumpIonTypeVector(ss, p, spectrumIdentificationItem->fragmentation);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIonTypeVector(stringstream &ss, string p, std::vector<IonTypePtr> &ionType)
{
  ss << p << "IonType.size() : " << ionType.size() << endl;
  for(int i = 0 ; i < ionType.size() ; i++)
    dumpIonType(ss, p, ionType[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIonType(stringstream &ss, string p, IonTypePtr &ionType)
{
  p += "IonType.";
  dumpCVParam(ss, p, *ionType);
  ss << p << "charge : " << ionType->charge << endl;
  ss << p << "index : ";
  for(int i = 0 ; i < ionType->index.size() ; i++)
    ss << ionType->index[i] << " ; ";
  ss << endl;
  dumpFragmentArrayVector(ss, p, ionType->fragmentArray);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpFragmentArrayVector(stringstream &ss, string p, std::vector<FragmentArrayPtr> &fragmentArray)
{
  ss << p << "FragmentArray.size() : " << fragmentArray.size() << endl;
  for(int i = 0 ; i < fragmentArray.size() ; i++)
    dumpFragmentArray(ss, p, fragmentArray[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpFragmentArray(stringstream &ss, string p, FragmentArrayPtr &fragmentArray)
{
  p += "FragmentArray.";
  ss << p << "values : ";
  for(int i = 0 ; i < fragmentArray->values.size() ; i++)
    ss << fragmentArray->values[i] << " ; ";
  ss << endl;
  dumpMeasure(ss, p, fragmentArray->measurePtr);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpMeasureVector(stringstream &ss, string p, std::vector<MeasurePtr> &measure)
{
  ss << p << "measure.size() : " << measure.size() << endl;
  for(int i = 0 ; i < measure.size() ; i++)
    dumpMeasure(ss, p, measure[i]);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpMeasure(stringstream &ss, string p, MeasurePtr &measure)
{
  p += "Measure.";
  dumpIdentifiableParamContainer(ss, p, *measure);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpAnalysisData(stringstream &ss, string p, AnalysisData &analysisData)
{
  p += "AnalysisData.";
  ss << p << "spectrumIdentificationList.size() : " << analysisData.spectrumIdentificationList.size() << endl;
}
////////////////////////////////////////////////////////////////////////////////
//struct PWIZ_API_DECL IdentifiableParamContainer : public ParamContainer
//{
//    IdentifiableParamContainer(const std::string& id_ = "",
//                 const std::string& name_ = "");
//    virtual ~IdentifiableParamContainer() {}
//
//    std::string id;
//    std::string name;
//
//    virtual bool empty() const;
//};
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpIdentifiableParamContainer(stringstream &ss, string p, IdentifiableParamContainer &identifiableParamContainer)
{
  p += "IdentifiableParamContainer.";
  ss << p << "id : " << identifiableParamContainer.id << endl;
  ss << p << "name : " << identifiableParamContainer.name << endl;
  dumpParamContainer(ss, p, identifiableParamContainer);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpParamContainer(stringstream &ss, string p, ParamContainer &paramContainer)
{
  p += "ParamContainer.";
  dumpParamGroup(ss,p,paramContainer.paramGroupPtrs);
  dumpCVParamsVector(ss, p, paramContainer.cvParams);
  dumpUserParam(ss, p, paramContainer.userParams);
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpParamGroup(stringstream &ss, string p, std::vector<ParamGroupPtr> &paramGroupPtrs)
{
  p += "ParamGroup.";
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::dumpUserParam(stringstream &ss, string p, std::vector<UserParam> &userParams)
{
  p += "UserParam.";
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Find the index of a given ID in pwiz data structure for controlled vocabulary
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::findIndex(SpectrumPtr &spectrum, int id)
{
  for(int i = 0 ; i < spectrum->cvParams.size() ; i++) {
    if(spectrum->cvParams[i].cvid == id)
      return i;
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
// Get the scan from the ID sting
////////////////////////////////////////////////////////////////////////////////
bool PWizInterface::getScan(string &str, int &scan, char *scanStr)
{
  scan = -1;
  // test if id contains something
  if(!str.length())
    return false;

  // hold strings
  vector<string> aux;
  stringSplit(str, aux, " \t");

  for(int i = 0 ; i < aux.size() ; i++) {

    vector<string> aux2;
    stringSplit(aux[i], aux2, "=");
    if(aux2.size() == 2) {
      if(aux2[0].compare(scanStr) == 0) {
        scan = getInt(aux2[1].c_str());
        return true;
      }
    }
  }

  // if ID exists, get
  return false;
}
////////////////////////////////////////////////////////////////////////////////
// Read Peptide (PWiz)
////////////////////////////////////////////////////////////////////////////////
/*int PWizInterface::readPeptideData(PeptideID_pepXml &peptide, const string &fn)
{
  //stringstream ss;
  string empty;
  int ret;

  try {
    PeptideID_pepXml peptide(fn);
    //dumpPeptideID(ss, empty, peptide);
    ret = getPeptideIDSize(peptide);
  } catch(...) {
    return 0;
  }
  return ret;
} */
////////////////////////////////////////////////////////////////////////////////
// Convert IdentData (PWiz) to psmSet (SPS)
////////////////////////////////////////////////////////////////////////////////
/*
m_spectrumFile        -> File name where spectrum is
m_scanNum             -> Spectrum scan number
m_annotation         -> Peptide annotation. Note that pepXML stores modifications very differently than we do, this will require reading in a couple different fields in pepXML and building a valid SpecNets annotation.
m_protein             -> Not totally necessary, but nice if they've got that info
m_charge              -> VERY IMPORTANT!
m_score               -> If they've got some sort of scoring, this is useful.
m_pValue              -> Same as m_score, nice if they've got it.
m_isDecoy             -> I don't know whether they actually flag this, but it would be really nice to use that info if they do.

m_spectrumFile       = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.
                       spectrumIdentificationResult.SpectraData.Identifiable.name +

                       Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.
                       spectrumIdentificationResult.SpectraData.location

m_scanNum.    = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.spectrumID

m_annotation = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.PeptideEvidence.Peptide.peptideSequence  +
               Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.PeptideEvidence.Peptide.modifications

m_protein    = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.PeptideEvidence.DBSequence.IdentifiableParamContainer.ParamContainer.CVParam.value when

                Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.PeptideEvidence.DBSequence.IdentifiableParamContainer.ParamContainer.CVParam.cvid :   1001088 -> MS_protein_description

m_charge.     = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.chargeState

m_score.      = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
               SpectrumIdentificationItem.rank
m_pValue     =
m_isDecoy.    = Ident.AnalysisCollection.SpectrumIdentification.spectrumIdentificationList.spectrumIdentificationResult.
                SpectrumIdentificationItem.PeptideEvidence.isDecoy
*/
////////////////////////////////////////////////////////////////////////////////
// identification section
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
int PWizInterface::acquireSearchModification(SearchModsData &search, std::vector<SpectrumIdentificationProtocolPtr> &protocols)
{
  for(int i = 0 ; i < protocols.size() ; i++) {
    // get search mods list
    std::vector<SearchModificationPtr> &sm = protocols[i]->modificationParams;
    for(int j = 0 ; j < sm.size() ; j++) {

      SearchModData md;
      md.massShift = sm[j]->massDelta;

      for(int k = 0 ; k < sm[j]->residues.size() ; k++) {
        string aux;
        aux = sm[j]->residues[k];
        md.residues.push_back(aux);
      }

      if(sm[j]->fixedMod)
        search.fixedMods.push_back(md);
      else
        search.variable.push_back(md);
    }
  }
}

/*
void dumpSearchModification(stringstream &ss, string p, SearchModificationPtr &searchModifications)
{
  p += "SearchModification.";
  dumpParamContainer(ss, p, *searchModifications);
  ss << p << "fixedMod : " << searchModifications->fixedMod << endl;
  ss << p << "massDelta : " << searchModifications->massDelta << endl;
  ss << p << "searchModifications[j] : ";
  for(int j = 0 ; j < searchModifications->residues.size() ; j++)
    ss << searchModifications->residues[j] << " ; ";
  ss << endl;
  dumpCVParam(ss, p, searchModifications->specificityRules);
}

Ident.AnalysisCollection.SpectrumIdentification.SpectrumIdentificationProtocol.SearchModification.fixedMod : 1
Ident.AnalysisCollection.SpectrumIdentification.SpectrumIdentificationProtocol.SearchModification.massDelta : 57
Ident.AnalysisCollection.SpectrumIdentification.SpectrumIdentificationProtocol.SearchModification.searchModifications[j] : C ;
*/
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populateModification(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, pwiz::identdata::Modification & modification, vector<string> &vec, bool outputFixed)
{
  // check location
  int loc = modification.location;
  if(loc < 0 || loc >= vec.size())
    return 0;

  // get the mass value
  double mass = modification.monoisotopicMassDelta;

  // check if it is nTerm
  if(loc == 0) {
    stringstream res;
    if(mass != 0.0) {
      res << "[";
      res << mass;
      res << "]";

      vec[loc] = res.str();
    }
    return 1;
  }

  // get the residues for the location
  stringstream ss;
  if(modification.residues.size())
    for(int j = 0 ; j < modification.residues.size() ; j++)
      ss << modification.residues[j];
  else
    ss << vec[loc];

  // check if is a fixed mod
  if(searchModsData.fixed(mass, modification.residues) && !outputFixed)
    return 0;

  // calculate replacement string
  stringstream res;
  if(mass != 0.0) {
    res << "(";
    res << ss.str();
    res << ",";
    res << mass;
    res << ")";
  } else {
    res << ss.str();
  }

  // apply
  vec[loc] = res.str();

  //substitutionModification.originalResidue;
  //substitutionModification.replacementResidue;

  //substitutionModification.avgMassDelta;
  //substitutionModification.monoisotopicMassDelta;
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::splitSequence(vector<string> &vec, string &str)
{
  vec.resize(str.length());
  for(int i = 0 ; i < str.length() ; i++)
    vec[i] = str[i];
}
////////////////////////////////////////////////////////////////////////////////
void PWizInterface::joinSequence(vector<string> &vec, string &str)
{
  str = "";
  for(int i = 0 ; i < vec.size() ; i++)
    str += vec[i];
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populatePeptideEvidence(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, PeptideEvidence &peptideEvidence, IdData &data, bool outputFixed)
{
  // DECOY flag
  bool isDecoy = peptideEvidence.isDecoy;
  string peptideAnnotation;

  // ANNOTATION
  if(peptideEvidence.peptidePtr) {

    pwiz::identdata::Peptide &peptide = *(peptideEvidence.peptidePtr);
    // get the peptide
    peptideAnnotation = peptide.peptideSequence;

    // split the peptide in chunks of 1 aa
    vector<string> workingPeptide;
    splitSequence(workingPeptide, peptideAnnotation);
    // add string for nterm
    string aaa;
    workingPeptide.insert(workingPeptide.begin(), aaa);

    // the the subsitution mods
     std::vector<pwiz::identdata::ModificationPtr> &modifications = peptide.modification;
    // process the modifications
    for(int m = 0 ; m < modifications.size() ; m++) {
      populateModification(psmSet, searchModsData, *(modifications[m]), workingPeptide, outputFixed);
    }
    // put the peptide back together
    joinSequence(workingPeptide, peptideAnnotation);
    // PROTEIN
    data.m_protein = "";
    DBSequencePtr & dbSequence = peptideEvidence.dbSequencePtr;
    // get protein name, if it exists
    std::vector<CVParam> &cvParams = dbSequence->cvParams;
    for(int j = 0 ; j < cvParams.size() ; j++)
      if(cvParams[j].cvid == 1001088) {
        data.m_protein = cvParams[j].value;
        break;
      }

    // store the psm entry
    psmPtr psm(new PeptideSpectrumMatch);
    psm->m_isDecoy      = isDecoy;
    psm->m_annotation   = peptideAnnotation;
    psm->m_charge       = data.m_charge;
    psm->m_score        = data.m_score;
    psm->m_scanNum      = data.m_scan;
    psm->m_pValue       = data.m_score;
    psm->m_protein      = data.m_protein;
    psm->m_spectrumFile = data.m_spectrumFile;

    // add to psmset
    psmSet.push_back(psm);
  }
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populateSpectrumIdentificationItem(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationItem &spectrumIdentificationItem, IdData &data, bool outputFixed)
{
  // CHARGE
  data.m_charge = spectrumIdentificationItem.chargeState;
  // SCORE
  data.m_score = spectrumIdentificationItem.rank;

  std::vector<PeptideEvidencePtr> &peptideEvidences =  spectrumIdentificationItem.peptideEvidencePtr;

  for(int l = 0 ; l < peptideEvidences.size() ; l++)
    populatePeptideEvidence(psmSet, searchModsData, *(peptideEvidences[l]) , data, outputFixed);
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populateSpectrumIdentificationResult(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationResult & spectrumIdentificationResult, bool outputFixed)
{
  IdData data;

  // SCAN number
  data.m_scan = -1;
  int scan;
  if(getScan(spectrumIdentificationResult.spectrumID, scan, "scan"))
    data.m_scan = scan;
  else if(getScan(spectrumIdentificationResult.spectrumID, scan, "index"))
    data.m_scan = scan;

  SpectraDataPtr &spectraData = spectrumIdentificationResult.spectraDataPtr;
  data.m_spectrumFile = spectraData->location;

  std::vector<SpectrumIdentificationItemPtr> & spectrumIdentificationItems = spectrumIdentificationResult.spectrumIdentificationItem;

  for(int k = 0 ; k < spectrumIdentificationItems.size() ; k++)
     populateSpectrumIdentificationItem(psmSet, searchModsData, *(spectrumIdentificationItems[k]) , data, outputFixed);
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populateFromSpectrumIdentification(PeptideSpectrumMatchSet &psmSet, SearchModsData &searchModsData, SpectrumIdentificationListPtr &spectrumIdentificationListPtr, bool outputFixed)
{
  std::vector<SpectrumIdentificationResultPtr> spectrumIdentificationResults = spectrumIdentificationListPtr->spectrumIdentificationResult;

  for(int j = 0 ; j < spectrumIdentificationResults.size() ; j++)
    populateSpectrumIdentificationResult(psmSet, searchModsData, *(spectrumIdentificationResults[j]), outputFixed);
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::populateFromIdent(PeptideSpectrumMatchSet &psmSet, IdentData& idd, SearchModsData &searchModsData, bool outputFixed)
{
  AnalysisCollection &analysisCollection = idd.analysisCollection;
  std::vector<SpectrumIdentificationPtr> &spectrumIdentification = analysisCollection.spectrumIdentification;

  for(int i = 0 ; i < spectrumIdentification.size() ; i++) {
    // get mods
    populateFromSpectrumIdentification(psmSet, searchModsData, spectrumIdentification[i]->spectrumIdentificationListPtr, outputFixed);
  }

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::convertIdentData2psmSet(vector<IdentDataPtr> &iddList, PeptideSpectrumMatchSet &psmSet, vector<SearchModsData> &searchModsData, bool outputFixed)
{
  try  {
    for (size_t i=0; i < iddList.size(); ++i) {
      IdentData& idd = *iddList[i];
      // process the data
      if(searchModsData.size() > i)
        populateFromIdent( psmSet, idd, searchModsData[i], outputFixed);
    }

  } catch (exception& e) {
    cerr << e.what() << endl;
    return 0;
  }
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::acquaireMods(vector<SearchModsData> &searchModsData)
{
  searchModsData.resize(iddList.size());

  try  {
    for (size_t i=0; i < iddList.size(); ++i) {
      IdentData& idd = *iddList[i];
      AnalysisProtocolCollection &analysisProtocolCollection = idd.analysisProtocolCollection;
      std::vector<SpectrumIdentificationProtocolPtr> &protocols = analysisProtocolCollection.spectrumIdentificationProtocol;
      // get search mods list
      acquireSearchModification(searchModsData[i], protocols);
    }
  } catch (exception& e) {
    cerr << e.what() << endl;
    return 0;
  }
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
double PWizInterface::findpercursormz(SpectrumPtr &spectrum)
{
  for(int j = 0 ; j < spectrum->precursors.size() ; j++) {
    for(int k = 0 ; k < spectrum->precursors[j].selectedIons.size() ; k++) {
      for(int l = 0 ; l < spectrum->precursors[j].selectedIons[k].cvParams.size() ; l++) {
        if(spectrum->precursors[j].selectedIons[k].cvParams[l].cvid == MS_selected_ion_m_z) {
          double ret = atof(spectrum->precursors[j].selectedIons[k].cvParams[l].value.c_str());
          return ret;
        }
      }
    }
  }
  return -1.0;
}
////////////////////////////////////////////////////////////////////////////////
double PWizInterface::findPercursorIntensity(SpectrumPtr &spectrum)
{
  for(int j = 0 ; j < spectrum->precursors.size() ; j++) {
    for(int k = 0 ; k < spectrum->precursors[j].selectedIons.size() ; k++) {
      for(int l = 0 ; l < spectrum->precursors[j].selectedIons[k].cvParams.size() ; l++) {
        if(spectrum->precursors[j].selectedIons[k].cvParams[l].cvid == MS_peak_intensity) {
          double ret = atof(spectrum->precursors[j].selectedIons[k].cvParams[l].value.c_str());
          return ret;
        }
      }
    }
  }
  return -1.0;
}
////////////////////////////////////////////////////////////////////////////////
double PWizInterface::findRetentionTime(SpectrumPtr &spectrum)
{
  for(int i = 0 ; i < spectrum->scanList.scans.size() ; i++) {
    for(int j = 0 ; j < spectrum->scanList.scans[i].cvParams.size() ; j++) {
      //cout << "cvParam " <<  spectrum->scanList.scans[i].cvParams[j].cvid <<"   value " <<  spectrum->scanList.scans[i].cvParams[j].value << endl;
      if(spectrum->scanList.scans[i].cvParams[j].cvid == MS_scan_start_time) {
        double ret = atof(spectrum->scanList.scans[i].cvParams[j].value.c_str());
        return ret;
      }
    }
  }
  return -1.0;
}

//getRetentionTime(spectrum->scanList.scans[0])

string getRetentionTime(const Scan& scan)
{
    ostringstream oss;
    oss << "PT" << scan.cvParam(MS_retention_time).timeInSeconds() << "S";
    return oss.str();
}
////////////////////////////////////////////////////////////////////////////////
// Convert MSData (PWiz) to specset (SPS)
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::convertMSData2SpecSet(MSData &msd, SpecSet &specSet, int msl)
{
  SpectrumList& spectrumList = *msd.run.spectrumListPtr;
  SpectrumPtr spectrum;
  //Read m/z and intensity values from the spectra
  const bool getBinaryData = true;
  size_t numSpectra = spectrumList.size();

  specSet.resize(0);

  for (int i = 0 ; i < numSpectra ; ++i) {
    // define specnets spectrum object
    specnets::Spectrum spect;
    spect.resize(0);

    ///////// Get the peak list ///////
    // Get spectrum binary data
  	spectrum = spectrumList.spectrum(i, getBinaryData);
  	// Define pairs data holder
  	vector<MZIntensityPair> pairs;
  	// Get MZ intensity pairs
  	spectrum->getMZIntensityPairs(pairs);
    // Add the peaks
  	for (vector<MZIntensityPair>::const_iterator it = pairs.begin(), end = pairs.end(); it!=end; ++it) {
      // add a peak
      spect.push_back(TwoValues<float>( (float)it->mz, (float)it->intensity));

      //cout << it->mz << " ; " << it->intensity << endl;
    }

    //cout << "------------ spectrum " << i << "----------" << endl;
    //for(int j = 0 ; j < spect.size() ; j++) {
    //  cout << j << " (" << spect[j][0] << " ; " << spect[j][1] << ")" << endl;
    //}

    ///////// Data ///////
    int index;

    // MS Level
    if((index = findIndex(spectrum, MS_ms_level)) >= 0) {
      spect.msLevel             = getInt(spectrum->cvParams[index].value.c_str());
      if(spect.msLevel < msl) continue;
    }

    // parent M/Z
    //if((index = findIndex(spectrum, MS_base_peak_m_z)) >= 0)
    //  spect.parentMZ            = getFloat(spectrum->cvParams[index].value.c_str());
    //else
    //if((index = findIndex(spectrum, MS_base_peak)) >= 0)
    //  spect.parentMZ            = getFloat(spectrum->cvParams[index].value.c_str());

    spect.parentMZ = findpercursormz(spectrum);
    //cout <<   "parentMZ : " <<   spect.parentMZ << endl;


    // parent charge
    if(spect.msLevel == 1)
      spect.parentCharge = 1;
    else
    if((index = findIndex(spectrum, MS_charge_state)) >= 0)
      spect.parentCharge        = getInt(spectrum->cvParams[index].value.c_str());
    else
      try {
        for(int j = 0 ; j < spectrum->precursors.size() ; j++) {
          for(int k = 0 ; k < spectrum->precursors[j].selectedIons.size() ; k++) {
            for(int l = 0 ; l < spectrum->precursors[j].selectedIons[k].cvParams.size() ; l++) {
              if(spectrum->precursors[j].selectedIons[k].cvParams[l].cvid == MS_charge_state) {
                spect.parentCharge = getInt(spectrum->precursors[j].selectedIons[k].cvParams[l].value.c_str());
              }
            }
          }
        }
      } catch (...) {
        //can't find charge...just ignore it
      }
      
    // base peak intensity
    //if((index = findIndex(spectrum, MS_base_peak_intensity)) >= 0)
    //if((index = findIndex(spectrum, MS_intensity_of_precursor_ion)) >= 0)
    //  spect.precursor_intensity = getInt(spectrum->cvParams[index].value.c_str());
      
    spect.precursor_intensity =   findPercursorIntensity(spectrum);
    //cout << "precursor_intensity = " << spect.precursor_intensity << endl;


    // Retention time
    //spect.retention_time = getRetentionTime(spectrum->scanList.scans[0]); //findRetentionTime(spectrum);

    spect.retention_time = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

    //cout << "retention_time = " << spect.retention_time << "  " << findRetentionTime(spectrum) <<endl;
    //cout << "retention_time = " << spect.retention_time << "  " << getRetentionTime(spectrum->scanList.scans[0]) <<endl;


    // scan number
    int scan;
    spect.scan = -1;
    if(getScan(spectrum->id, scan, "scan"))
      spect.scan = scan;
    
    // File name
    size_t  found = inputFilename.find_last_of("/\\");
    spect.fileName  = inputFilename.substr(found + 1);

    // parent mass
    if(spect.parentCharge > 0)
      spect.parentMass          =  spect.parentMZ * spect.parentCharge - AAJumps::massHion * (spect.parentCharge - 1.0);
    else
      spect.parentMass          =  spect.parentMZ;

    // MS Fragmentation model
    spect.msFragType            = specnets::Spectrum::FragType_CID;
    if((index = findIndex(spectrum, MS_electron_transfer_dissociation)) >= 0)
      spect.msFragType          = specnets::Spectrum::FragType_ETD;


    if((index = findIndex(spectrum, MS_high_energy_collision_induced_dissociation)) >= 0) {
      spect.msFragType          = specnets::Spectrum::FragType_HCD;
    } else {
      try {
        for(int j = 0 ; j < spectrum->precursors.size() ; j++) {
              for(int k = 0 ; k < spectrum->precursors[j].activation.cvParams.size() ; k++) {
                if(spectrum->precursors[j].activation.cvParams[k].cvid == MS_high_energy_collision_induced_dissociation) {
                  spect.msFragType          = specnets::Spectrum::FragType_HCD;
                }
            }
        }
      } catch (...) {
        //can't find charge...just ignore it
      }
    }

    //cout << "------------ spectrum " << i << "----------" << endl;
    //for(int j = 0 ; j < spect.size() ; j++) {
    //  cout << j << " (" << spect[j][0] << " ; " << spect[j][1] << ")" << endl;
    //}

    // store the spectrum in the specset
    specSet.push_back(spect);

    //cout << "------------ spectrum " << i << "----------" << endl;
    //for(int j = 0 ; j < specSet[specSet.size()-1].size() ; j++) {
    //  cout << j << " (" << specSet[specSet.size()-1][j][0] << " ; " << specSet[specSet.size()-1][j][1] << ")" << endl;
    //}

  }

  //cout << "==================================================" << endl;
  //for(int i = 0 ; i < specSet.size() ; i++) {
  //  cout << "------------ spectrum " << i << "----------" << endl;
  //  for(int j = 0 ; j < specSet[i].size() ; j++) {
  //    cout << j << " (" << specSet[i][j][0] << " ; " << specSet[i][j][1] << ")" << endl;
  //  }
  //}

  return 1;
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::loadDataUsingPWiz(const string &inName, SpecSet &specs, int msl)
{
	FullReaderList readers;
	inputFilename = inName;

  try {
    MSDataFile msd(inName, &readers);
    convertMSData2SpecSet(msd, specs, msl);
    //dumpMSData(msd, msl);


	} catch (exception& e) {
    cout << "Error loadding file: " << e.what() << endl;
    return 0;
  }

  FilenameManager mngr(inName);
  for (unsigned int i = 0; i < specs.size(); i++)
  {
    if (specs[i].fileName.length() == 0)
    {
      specs[i].fileName = mngr.getFilenameWithExtension();
    }
  }
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
int PWizInterface::loadDataUsingPWiz(PeptideSpectrumMatchSet &psmSet, vector<SearchModsData> &searchModsData, bool outputFixed)
{
	//populate an MSData object from an MS data filepath
  try {
    convertIdentData2psmSet(iddList, psmSet, searchModsData, outputFixed);
	} catch (exception& e) {
    return 0;
  }
  return 1;
} // -lpwiz_data_identdata -lpwiz_data_identdata_version -lpwiz_utility_chemistry -lpwiz_data_proteome
////////////////////////////////////////////////////////////////////////////////
bool PWizInterface::openIdent(string &inName)
{
	pwiz::identdata::DefaultReaderList readers;
  pwiz::identdata::Reader::Config readerConfig;

  IterationListenerRegistry* pILR = 0;
  readerConfig.iterationListenerRegistry = pILR;

  inputFilename = inName;
	//populate an MSData object from an MS data filepath
  try {
    readers.read(inName, iddList, readerConfig);
    //convertIdentData2psmSet(iddList, psmSet, searchModsData, outputFixed);
	} catch (exception& e) {
    cout << "Error loadding file: " << e.what() << endl;
    return false;
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////
//#else
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
