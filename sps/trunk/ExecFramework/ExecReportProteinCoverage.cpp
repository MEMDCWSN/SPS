////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "ExecReportProteinCoverage.h"
#include "Logger.h"

#include "aminoacid.h"
#include "db_fasta.h"
#include "utils.h"

#include "copyright.h"

#define CSV_SEP ';'


using namespace specnets;
using namespace std;
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
const char * shtmlhead =
  "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"\n"
  "  \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"

  "<HTML xmlns=\"http://www.w3.org/1999/xhtml\">\n"

  "<head>\n"

  "  <meta http-equiv=\"Content-Type\" content=\"text/shtml; charset=ISO-8859-1\" />\n"
//   "  <meta http-equiv=\"PRAGMA\" content=\"NO-CACHE\" />\n"
//   "  <meta http-equiv=\"CACHE-CONTROL\" content=\"NO-STORE, NO-CACHE, MUST-REVALIDATE, POST-CHECK=0, PRE-CHECK=0\" />\n"
//   "  <meta http-equiv=\"EXPIRES\" content=\"01 Jan 1970 00:00:00 GMT\" />\n"
;


const char * shtmlhead2 =

  "  <script type=\"JavaScript\">\n"
  "  <!--hide\n"
  "  function modalWin(filename, name)\n"
  "  {\n"
  "    if (! name)\n"
  "      name = \"name\";\n"

  "    if (window.showModalDialog)\n"
  "    {\n"
  "      window.showModalDialog(filename, name, \"dialogWidth:640px; dialogHeight:480px\");\n"
  "    }\n"
  "    else\n"
  "    {\n"
  "      window.open(filename, name, 'width=640, height=480, toolbar=no, directories=no, status=no, menubar=no, scrollbars=no, resizable=no, modal=yes');\n"
  "    }\n"
  "  }\n"

  "  function modalForm(form, name)\n"
  "  {\n"
  "    if (! name)\n"
  "      name = \"name\";\n"

  "    if (window.focus)\n"
  "    {\n"
  "      modalWin(\"\", name);\n"
  "      form.target = name;\n"
  "    }\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script src=\"../../js/pager.js\" language=\"javascript\" type=\"text/javascript\"></script>\n"
  "  <noscript>\n"
  "    <meta http-equiv=\"refresh\" content=\"10\">\n"
  "  </noscript>\n"

  "  <script language=\"JavaScript\">\n"
  "  <!--\n"
  "  var sURL = unescape(window.location.pathname);\n"
  "  function doLoad()\n"
  "  {\n"
  "      setTimeout( \"refresh()\", 10*1000 );\n"
  "  }\n"

  "  function refresh()\n"
  "  {\n"
  "      window.location.href = sURL;\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script language=\"JavaScript1.1\">\n"
  "  <!--\n"
  "  function refresh()\n"
  "  {\n"
  "      window.location.replace( sURL );\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"

  "  <script language=\"JavaScript1.2\">\n"
  "  <!--\n"
  "  function refresh()\n"
  "  {\n"
  "      window.location.reload(true);\n"
  "  }\n"
  "  //-->\n"
  "  </script>\n"
  "</head>\n";
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
  ExecReportProteinCoverage::ExecReportProteinCoverage(void) : indexBase(1)
  {
    m_name = "ExecReportProteinCoverage";
    m_type = "ExecReportProteinCoverage";
  }

  // -------------------------------------------------------------------------
  ExecReportProteinCoverage::ExecReportProteinCoverage(const ParameterList & inputParams) :
    ExecBase(inputParams), indexBase(1)
  {
    m_name = "ExecReportProteinCoverage";
    m_type = "ExecReportProteinCoverage";
  }
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::validateParams(std::string & error)
{
  m_isValid = false;

  VALIDATE_PARAM_EXIST("OUTDIR");
  VALIDATE_PARAM_EXIST("INPUT_FASTA");
  VALIDATE_PARAM_EXIST("SPS_CONTIG_NAMES");
  VALIDATE_PARAM_EXIST("SPS_CONTIG_MATCHES_IDX");
  VALIDATE_PARAM_EXIST("SPS_CONTIG_MATCHES");
  VALIDATE_PARAM_EXIST("SPS_CONTIG_SPECTRA");
  VALIDATE_PARAM_EXIST("CSPS_CONTIG_MATCHES_IDX");
  VALIDATE_PARAM_EXIST("CSPS_CONTIG_MATCHES");
  VALIDATE_PARAM_EXIST("CSPS_CONTIG_SPECTRA");
  VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");

  m_isValid = true;
  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::merge(void)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
vector<ExecBase*> const & ExecReportProteinCoverage::split(int numSplit)
{
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
ExecBase * ExecReportProteinCoverage::clone(const ParameterList & input_params) const
{
  return new ExecReportProteinCoverage(input_params);
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::loadOutputData(void)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::saveInputData(std::vector<std::string> & filenames)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::saveOutputData(void)
{
  // Process each file
  for(int i = 0 ; i < m_protein.size() ; i++)     {
    // write output html file
    if(m_protein[i].page.length())
      outputFile(m_protein[i].oFileName, m_protein[i].page );
    // write output csv file
    if(m_protein[i].csvFile.length())
      outputFile(m_protein[i].oFileNameCsv, m_protein[i].csvFile );
  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::genheader(std::string & shtml)
{
  string pathPrefix = ".";
  if (m_params.exists("HTML_DEFS"))
    pathPrefix = m_params.getValue("HTML_DEFS");

  string pageTitle = "";
  if (m_params.exists("REPORT_TITLE"))
    pageTitle = m_params.getValue("REPORT_TITLE");


  shtml = shtmlhead;

  shtml += "  <title>UCSD Computational Mass Spectrometry Website</title>\n";
  shtml += "  <link href=\"";
  shtml += pathPrefix;
  shtml += "/styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n";
  shtml += "  <link rel=\"shortcut icon\" href=\"";
  shtml += pathPrefix;
  shtml += "/images/favicon.ico\" type=\"image/icon\" />\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/scripts/util.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/scripts/download.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/scripts/render.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/scripts/inspect.js\" language=\"javascript\" type=\"text/javascript\"></script>\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/js/mootools.js\" type=\"text/javascript\"></script>\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/js/slimbox.js\" type=\"text/javascript\"></script>\n";
  shtml += "  <link href=\"";
  shtml += pathPrefix;
  shtml += "/css/slimbox.css\" rel=\"stylesheet\" type=\"text/css\" media=\"screen\" />\n";
  shtml += "  <script src=\"";
  shtml += pathPrefix;
  shtml += "/js/sorttable.js\"></script>\n";

  shtml += shtmlhead2;
  shtml += '\n';

  shtml += "  <div id=\"rb_logos\">\n";
  shtml += "    <h3 align=center>";
  shtml += pageTitle;
//  shtml += "<img src=\"";
//  shtml += pathPrefix;
//  shtml += "/images/pagelogo_ccms_left.jpg\" border=\"0\"/>";
  shtml += "</h3>\n";
  shtml += "  </div>\n";
}

void ExecReportProteinCoverage::genFooter(std::string & page)
{
//  page += "<div class=\"footer\">";
//  page += "<p>";
//  page += REPORT_FOOTER_TEXT;
//  page += "</p>";
//  page += "</div>";
}

void ExecReportProteinCoverage::genTableFooter(std::string & page)
{
  page += "<table><tr><td>";
  page += "<p>";
  page += REPORT_FOOTER_TEXT;
  page += "</p>";
  page += "</td></tr></table>";
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
char *ExecReportProteinCoverage::readFile(std::string &iFileName, int &length)
{
  char * buffer;

  ifstream is;
  is.open( iFileName.c_str() , ios::binary );
  if(is.fail())
    return NULL;

  // get length of file:
  is.seekg (0, ios::end);
  length = is.tellg();
  is.seekg (0, ios::beg);

  // allocate memory:
  buffer = new char [length];

  // read data as a block:
  is.read (buffer,length);
  is.close();

  return buffer;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::outputFile(std::string &oFileName, std::string &page)
{
  // get the output directory
  string path = m_params.getValue("OUTDIR");
  rtrim(path);

  // Open output file for writing
  string aux = path;
  aux += "/";
  aux += oFileName;
  FILE *f = fopen(aux.c_str(), "wb");
  if(!f)
    return -1;

  fwrite(page.c_str(), page.length(), 1, f);
  fclose(f);
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs)
{
  // First file: homglue matches
	for(int i = 0 ; i < contigMatches.size() ; i++)
    // if equal to -1, then has no match
	  if(contigMatches[i][0] != -1) {
      // create strcuture for current contig
      ContigMatchData contigMatchData;
      // Store contig index
      contigMatchData.contigIndex = i;
	    // find curreint protein
	    std::map<int, std::vector<ContigMatchData> >::iterator it;
	    it = contigs.find(contigMatches[i][0]);
	    // if key exists, insert the new value
	    if(it != contigs.end()) {
	      it->second.push_back(contigMatchData);
	    } else {
	      // if it doesn't exist, create a new entry.
        // create the vector
        vector<ContigMatchData> aux;
        aux.push_back(contigMatchData);
       // Insert CSPS contig (index) info indexed by protein index
        contigs.insert(pair<int, vector<ContigMatchData> > (contigMatches[i][0], aux));
	    }
  	}

  // Second file: CSPS contig matches
	for(int i = 0 ; i < contigMatchesIndex.size() ; i++)
    // test if there is any data to store
	  if(contigMatchesIndex[i].size() > 0) {

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, i);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig)
    		for(int j = 0 ; j < contigMatchesIndex[i].size() ; j++) {
  		    PairContigProteinMassIdx pairItem;
  		    pairItem.contigMassIdx = (int)contigMatchesIndex[i][j][0];
  		    pairItem.proteinMassIdx = (int)contigMatchesIndex[i][j][1];
  		    contig->pair.push_back(pairItem);
  		  }
	  }

	// 3rd file: Contig spectra
	for(int i = 0 ; i < contigSpectra.size() ; i++)
	  if(contigSpectra[i].size() > 0) {

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, i);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig)
    		for(int j = 0 ; j < contigSpectra[i].size() ; j++) {
  		    float mass = contigSpectra[i][j][0];
  		    contig->contigMass.push_back(mass);
  		  }
	  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::dumpContig(std::map<int, std::vector<ContigMatchData> > &contigs)
{
	DEBUG_MSG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
	for(int i = 0 ; i < contigs.size() ; i++) {
	  for(int k = 0 ; k < contigs[i].size() ; k++) {
      DEBUG_MSG("Protein index:  " << i);
      DEBUG_MSG("Contig index: " << contigs[i][k].contigIndex);
      DEBUG_MSG("Masses: ");
      stringstream aux;
      aux << "# " << contigs[i][k].contigMass.size() << " - ";
      for(int j = 0 ; j < contigs[i][k].contigMass.size() ; j++)
        aux << contigs[i][k].contigMass[j] << " ; ";
      DEBUG_MSG(aux);
      stringstream aux2;
      aux2 << "(protIdx ; ContIdx) = ";
      aux2 << "# " << contigs[i][k].pair.size() << " - ";
      for(int j = 0 ; j < contigs[i][k].pair.size() ; j++)
        aux2 << "(" << contigs[i][k].pair[j].proteinMassIdx << " ; " << contigs[i][k].pair[j].contigMassIdx << ") ; ";
      DEBUG_MSG(aux2);
    }
	}
	DEBUG_MSG(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
ContigMatchData *ExecReportProteinCoverage::getContigData( std::map<int, std::vector<ContigMatchData> > &contig, int contigIndex)
{
  std::map<int, vector<ContigMatchData> >::iterator it;
  // search for contig linked to given contig index
  for( it = contig.begin() ; it != contig.end() ; it++)
    // Cycle thru contig vector (for each protein)
    for(int i = 0 ; i < (*it).second.size() ; i++)
      if( (*it).second[i].contigIndex == contigIndex)
        // return pointer to struct
        return &((*it).second[i]);

  return NULL;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::getContigMassAtProteinIndex(ContigMatchData &contig, int proteinIndex, float &contigMassAtIndex)
{
  for(int i = 0 ; i < contig.pair.size() ; i++) {
    if(contig.pair[i].proteinMassIdx == proteinIndex)
      if(contig.pair[i].contigMassIdx < contig.contigMass.size()) {
        contigMassAtIndex = contig.contigMass[contig.pair[i].contigMassIdx];
        return true;
      } else {
        return false;
      }
  }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::getContigNames(string &contigNamesFileName)
{
  int length;
  // read the file
  char *buffer = readFile(contigNamesFileName, length);
  if(!buffer) {
    ERROR_MSG("ERROR: couldn't open file " << contigNamesFileName);
    return -1;
  }

  // get references into vector
  stringSplit(string(buffer), m_contigNames, "\n\r");
  // delete the buffer
  delete [] buffer;
  // return ok
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
/*int ExecReportProteinCoverage::getProteinReferences(string &refIndex)
{
  int length;
  // read the file
  char *buffer = readFile(refIndex, length);
  if(!buffer) {
    ERROR_MSG("ERROR: couldn't open file " << refIndex);
    return -1;
  }

  // get references into vector
  stringSplit(string(buffer), m_references, "\n\r");
  // delete the buffer
  delete [] buffer;
  // return ok
  return 1;
} */
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::processProteinsFile(void)
{
  if(m_fasta.size() <= 0) return;

  for(int i = 0 ; i < m_fasta.size() ; i++) {

    // Protein sequence info
    PdProteinInfo protein;
    // set protein name
    protein.proteinName = m_fasta.getID(i);
    // get the sequence
    string sequence = m_fasta.getSequence(i);
    // cycle thru the sequence
    for(int j = 0 ; j < sequence.size() ; j++) {
      aaCell cell;
      // the aa sequence
      cell.aa = sequence[j];
      // just one
      cell.colspan = 0;
      // its position
      cell.startPosition = j;
      protein.proteinDetail.processedAA.push_back(cell);
    }
    // set protein size
    protein.proteinDetail.startPosition = 0;
    protein.proteinDetail.endPosition = protein.proteinDetail.processedAA.size();
    // set internal protein index value
    protein.ProteinIndex = i;
    // set protein detail file name
    protein.oFileName = "protein_details.";
    protein.oFileName += parseInt(protein.ProteinIndex + indexBase);
    // set CSV filename
    protein.oFileNameCsv = protein.oFileName;
    // html file extension
    protein.oFileName += ".html";
    // csv file extension
    protein.oFileNameCsv += ".txt";
    // add the protein to the protein list
    m_protein.insert(pair<int, PdProteinInfo>( protein.ProteinIndex, protein));
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::processContigs(int ProteinIndex, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source)
{
  // get peak tolerance
   float peakTolerance = getFloat(m_params.getValue("TOLERANCE_PEAK").c_str());

  // get protein data reference
  PdProteinInfo & protein = m_protein[ProteinIndex];

  // Get csps contig information. If not found for this protein, NULL is returned.
  std::map<int, vector<ContigMatchData> >::iterator it2;
  it2 = source.find(protein.ProteinIndex);

  // Add CSPS contig information (if exists)
  if(it2 == source.end())
    return;

  // There may be several contigs for this protein. Must cycle thru them
  for(int i = 0 ; i < (*it2).second.size() ; i++) {

    // get csps contig span and location
    int size = (*it2).second[i].pair.size();
    int minimunIndex = (*it2).second[i].pair[0].proteinMassIdx;
    int maximunIndex = (*it2).second[i].pair[size-1].proteinMassIdx;

    // Define the data structure to hold it and set some initial data
    ProteinDetaisContigInfo contig;
    contig.startPosition = minimunIndex;
    contig.endPosition = maximunIndex;


    // Loop thru contig size to generate info cells
    for(int j = 0 ; j < size-1 ; j++) {

      // get the corresponding contig indexes
      int contigIndex = (*it2).second[i].pair[j].contigMassIdx;
      int contigIndexAfter = (*it2).second[i].pair[j+1].contigMassIdx;

      // get the correspondif protein indexes
      int proteinIndex = (*it2).second[i].pair[j].proteinMassIdx;
      int proteinIndexAfter = (*it2).second[i].pair[j+1].proteinMassIdx;
      // Get the mass at this point

      if(contigIndex >= (*it2).second[i].contigMass.size() || contigIndex >= (*it2).second[i].contigMass.size()) {
        DEBUG_MSG("Contig #" << j << ": Invalid mass index: " << contigIndex << " > " << (*it2).second[i].contigMass.size());
        return;
      }

      float massCsps = (*it2).second[i].contigMass[contigIndex];
      // and at the point after
      float massAfter = (*it2).second[i].contigMass[contigIndexAfter];

      string sequence;

      // get AA sequence for interval
      if( (contigIndex >= protein.proteinDetail.startPosition ) && (contigIndexAfter < protein.proteinDetail.endPosition) )
        for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < protein.proteinDetail.processedAA.size()) ; k++)
          sequence += protein.proteinDetail.processedAA[k].aa;

      // get masses for sequence
      vector<float> masses;
      getMasses((char*)sequence.c_str(), masses);
      float mass = 0.0;
      for(int k = 0 ; k < masses.size() ; k++)
        mass += masses[k];

      string result = sequence;

      // massAux used to calculate mass difference
      float massAux = (massAfter - massCsps) - mass;

      if(fabs(massAux) > peakTolerance) {
        result = "(";
        result += sequence;
        result += ", ";
        result += parseFloat(massAux, 2);
        result += ")";
      }

      // cell to hold one data item
      aaCell cell;

      // save it
      cell.aa = result;
      // calculate the colspan value for this cell
      cell.colspan = (*it2).second[i].pair[j+1].proteinMassIdx - (*it2).second[i].pair[j].proteinMassIdx - 1;
      // and the start position
      cell.startPosition = (*it2).second[i].pair[j].proteinMassIdx;

      // add the cell to the list
      contig.processedAA.push_back(cell);
    }

    // set the contig name
    contig.name = parseInt((*it2).second[i].contigIndex);

    // set the base 1 index
    contig.base1Idx = (*it2).second[i].contigIndex + 1;

    // add the csps contig info to the csps contig list
    target.insert(pair<int, ProteinDetaisContigInfo > ((*it2).second[i].contigIndex, contig));
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ExecReportProteinCoverage::getIntFromSeqName(string seq)
{
  vector<string> aux;
  stringSplit(seq, aux, ":");
  if(aux.size() > 1)
    return aux[1];
  return "";
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ExecReportProteinCoverage::getContigName(int i)
{
	if( (i>= 0) && i < m_contigNames.size() )
		return m_contigNames[i];
  return string("--");
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::populateContigNames(PdProteinDetail &contig)
{
  PdProteinDetail::iterator it;
  for(it = contig.begin() ; it != contig.end() ; it++) {
    int contigID = it->first;
   	string name = getContigName(contigID);
    it->second.name = (name.length() > 0 ? name : parseInt(contigID + 1) );
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::generateOutputContig(int &i, int proteinIdx, vector<int> &vectorID, std::string &page, int cellPerLine, PdProteinDetail &contig, bool link)
{
  // get protein size
  int proteinSize =  m_protein[proteinIdx].proteinDetail.endPosition;

  // Add CSPS contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = contig[vectorID[j]];

    // Check if sequence in the range we are outputing now
    if( (contigSequence.startPosition < i + cellPerLine) &&
        (contigSequence.endPosition   > i              )    ) {

      // Write the contig id and link
      if(link) {
        page += "<tr><th align=\"right\"><a href=\"contig.";
        page += getIntFromSeqName(contigSequence.name);
        if(indexBase == 1)
          page += ".0";
        page += ".html\" style=\"color: white\">";
        page += contigSequence.name;
        page += "</a></th>";

        //page += "<tr><th align=\"right\">";
        //page += contigSequence.name;
        //page += "</th>";
      } else {
        page += "<tr>";
        page += "<td class=\"rc2\" width=\"100px\">";
        page += "CSPS ";
        //page += contigSequence.name;
        page += parseInt(contigSequence.base1Idx);
        page += "</TD>";
      }

      // find first cell to output
      int l = 0;
      while( (l < contigSequence.processedAA.size()) &&
             (i > contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan) )
        l++;

      // cycle thru
      for(int k = i ; (k < i+cellPerLine) && (k < proteinSize) ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence.startPosition)
          page += "<td class=\"rc2\" />\n";
        else if( (l >= contigSequence.processedAA.size()) )
          page += "<td class=\"rc2\" />\n";
        // otherwise, the content
        else {

          page += "<td ";
          // page += " class=\"rc1\" style=\"background-color: transparent; border: solid 0 #060; border-left-width:1px;border-right-width:1px;\" ";

          int border = 0;

          string outputString = contigSequence.processedAA[l].aa;

          // Calc colspan nr of occupied cells
          int colspan = contigSequence.processedAA[l].colspan + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan >= i + cellPerLine) {
            colspan -= contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence.processedAA[l].colspan + 1) / 2 )
              outputString = "";
          }
          // beggining
          if(contigSequence.processedAA[l].startPosition < i) {
            colspan -= i - contigSequence.processedAA[l].startPosition;
            border++;
            if(colspan <= (contigSequence.processedAA[l].colspan + 1) / 2 )
              outputString = "";
          }

          page += " class=\"rca";
          page += parseInt(border);
          page += "\" ";

          if(colspan > 1) {
            page += " colspan=\"";
            page += parseInt(colspan);
            page += '"';
            k += colspan-1;
          }
          page += '>';
          page +=  outputString;
          page += "</td>";
          l++;
        }
      }
      page += "</tr>\n";
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::generateOutputContigCsv(int &i, int proteinIdx, vector<int> &vectorID, std::string &page, int cellPerLine, PdProteinDetail &contig, bool link)
{
  // Add CSPS contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = contig[vectorID[j]];

    // Check if sequence in the range we are outputing now
    if( (contigSequence.startPosition < i + cellPerLine) &&
        (contigSequence.endPosition   > i              )    ) {

      // Write the contig id
      if(link)
        // sps contig: output name
        page += contigSequence.name;
      else
        // csps contig: output config index 1 based
        parseInt(contigSequence.base1Idx);
      // separator
      page += CSV_SEP;

      // find first cell to output
      int l = 0;
      while( (l < contigSequence.processedAA.size()) &&
             (i > contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan) )
        l++;

      // cycle thru
      for(int k = i ; k < i+cellPerLine ; k++) {

        // if start position is lower than current position, output an empty cell
        if(k < contigSequence.startPosition) {
          page += CSV_SEP;
          page += CSV_SEP;
        } else if( (l >= contigSequence.processedAA.size()) ) {
          page += CSV_SEP;
          page += CSV_SEP;
        // otherwise, the content
        } else {

          int border = 0;

          string outputString = contigSequence.processedAA[l].aa;

          // Calc colspan nr of occupied cells
          int colspan = contigSequence.processedAA[l].colspan + 1;
          // careful with split cells at the beggining or end -- end
          if(contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan >= i + cellPerLine) {
            colspan -= contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan - i - cellPerLine + 1;
            border += 2;
            if(colspan < (contigSequence.processedAA[l].colspan + 1) / 2 )
              outputString = "";
          }
          // begining
          if(contigSequence.processedAA[l].startPosition < i) {
            colspan -= i - contigSequence.processedAA[l].startPosition;
            border++;
            if(colspan <= (contigSequence.processedAA[l].colspan + 1) / 2 )
              outputString = "";
          }

          if(!(border & 0x01))
            page += '|';

          if(colspan == 1) {
            // cell content
            page += CSV_SEP;
            page += outputString;
          } else {

            if((outputString.length() > 0) && (outputString[0] == '(')) {
              // outputing (xx, yy), and empty cells following
              page += CSV_SEP;
              page += outputString;
              // until the end of cell
              while(--colspan && (k < i+cellPerLine)) {
                page += CSV_SEP;
                page += CSV_SEP;
                k++;
              }

            } else {
              // outputing A.B.C.D
              int cells = 0;
              while( (cells < colspan-1) && (k < i+cellPerLine)) {
                page += CSV_SEP;
                page += outputString[cells];
                page += CSV_SEP;
                page += '.';
                k++;
                cells++;
              }
              page += CSV_SEP;
              page += outputString[cells];

            }

          }
          // cell end
          page += CSV_SEP;
          // right border content (separator)
          if( (!(border & 0x02)) && (true) && (true) ) ;
//            page += '|';

          l++;
        }
      }
      page += '\n';
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ExecReportProteinCoverage::getOrder(PdProteinDetail &contig, int i, int size, vector<int> &order)
{
  // temporary vector to hold contig start position per contig
  vector<ContigOrdering> aux;
  // contig iterator
  PdProteinDetail::iterator it;
  // cycle thru all contigs
  for( it = contig.begin() ; it != contig.end() ; it++ ) {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = it->second;
    // Check for empty contigs
    if(!contigSequence.processedAA.size()) continue;
    // Check if sequence in the range we are outputing now
    if( (contigSequence.startPosition < i + size) &&
        (contigSequence.endPosition   > i              )    ) {
      // find first cell to output
      int l = 0;
      while( (l < contigSequence.processedAA.size()) &&
             (i > contigSequence.processedAA[l].startPosition + contigSequence.processedAA[l].colspan) )
        l++;

      // create ordering cell with contig info
      ContigOrdering orderingCell;
      orderingCell.contigIndex = it->first;
      orderingCell.startIndex = contigSequence.startPosition; //contigSequence.processedAA[l].startPosition;
      orderingCell.endIndex = contigSequence.endPosition;
      aux.push_back(orderingCell);
    }
  }

  // Order
  sort(aux.begin(), aux.end());
  // set order in order vector
  for(int k = 0 ; k < aux.size() ; k++)
    order.push_back(aux[k].contigIndex);
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::generateOutputFile2(int proteinIdx, std::string &page, int cellPerLine)
{
  // test for proteins with no contig mapped to them, and do not generate them
  if((m_protein[proteinIdx].cspsDetails.size() == 0) && (m_protein[proteinIdx].spsDetails.size() == 0))
    return 1;

  // details
  genheader(page);

  page += "<body>\n";
  page += "  <div id=\"bodyWrapper\">\n";
  page += "  <div id=\"textWrapper\">\n";
  page += "  <table class=\"result\" width=\"100%\" style=\"border-spacing: 0px\">\n";
  page += "    <tr>\n";
  page += "      <td colspan=\"0\"><h2><i>";
  page += m_protein[proteinIdx].proteinName;
  page += "</i></h2><hr></td>\n";
  page += "    </tr>\n";
  page += "  </table>\n";

  page += "<table><tr><td><a href=\"";
  page += m_protein[proteinIdx].oFileNameCsv;
  page += "\">Protein coverage as Excel-ready format (TXT file)</a></td></tr></table>";


  // Get the protein sequence. From this, we know the lenght
  ProteinDetaisContigInfo & proteinSequence = m_protein[proteinIdx].proteinDetail;

  // general position indexer
  int i = 0;


  // Keep it under protein size
  while(i < proteinSequence.processedAA.size()) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    vector<int> spsID, cspsID;
    // get the csps contig indexes
    getOrder(m_protein[proteinIdx].cspsDetails, i, cellPerLine, cspsID);
    // get the sps contig indexes
    getOrder(m_protein[proteinIdx].spsDetails, i, cellPerLine, spsID);


    // if we are starting a new table, add table beggining
    page += "<table class=\"result2\" width=\"100%\" style=\"background-color: #CCCFFF\">\n";
    page += "<tr>";
    page += "<td class=\"rc2\" width=\"100px\" align=\"right\">";
    page += parseInt(i+1);
    page += "&nbsp;</td>";


    // output protein
    for(int j = i ; ( j < i + cellPerLine ) && ( j < proteinSequence.processedAA.size() ) ; j++) {

      page += "<td align=\"center\" ";
      // if an empty cell, it's a separator column. It should be 1 pixel wide
      if( proteinSequence.processedAA[j].aa.length() == 0 )
        page += "class=\"rh2\" ";
      else {
        page += "class=\"rh1\" ";
      }
      page += ">";
      // The AA from the protein sequence
      page +=  proteinSequence.processedAA[j].aa;
      // header cell terminator
      page += "</td>";
    }
    // Header row terminator
    page += "</tr>\n";


    // Add CSPS contig information (if exists)
    generateOutputContig(i, proteinIdx, cspsID, page, cellPerLine, m_protein[proteinIdx].cspsDetails, false);

    // Add SPS contig information (if exists)
    generateOutputContig(i, proteinIdx, spsID, page, cellPerLine, m_protein[proteinIdx].spsDetails, true);

    // HTML table terminator
    page += "  </table><br>\n";

    i += cellPerLine;
  }

  page += "</div></div>";

  genFooter(page);

  page += "</body>";

	return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ExecReportProteinCoverage::generateOutputFileCsv(int proteinIdx, std::string &page, int cellPerLine)
{
  // test for proteins with no contig mapped to them, and do not generate them
  if((m_protein[proteinIdx].cspsDetails.size() == 0) && (m_protein[proteinIdx].spsDetails.size() == 0))
    return 1;

  // Get the protein sequence. From this, we know the lenght
  ProteinDetaisContigInfo & proteinSequence = m_protein[proteinIdx].proteinDetail;

  // general position indexer
  int i = 0;

  // Keep it under protein size
  while(i < proteinSequence.processedAA.size()) {

    // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
    vector<int> spsID, cspsID;
    // get the csps contig indexes
    getOrder(m_protein[proteinIdx].cspsDetails, i, cellPerLine, cspsID);
    // get the sps contig indexes
    getOrder(m_protein[proteinIdx].spsDetails, i, cellPerLine, spsID);


    // if we are starting a new table, add table beggining
    page += '\n';
    page += parseInt(i+1);
    page += CSV_SEP;

    // output protein
    for(int j = i ; ( j < i + cellPerLine ) && ( j < proteinSequence.processedAA.size() ) ; j++) {

      page += CSV_SEP;
      // The AA from the protein sequence
      page +=  proteinSequence.processedAA[j].aa;
      // header cell terminator
      page += CSV_SEP;
    }
    // Header row terminator
    page += '\n';


    // Add CSPS contig information (if exists)
    generateOutputContigCsv(i, proteinIdx, cspsID, page, cellPerLine, m_protein[proteinIdx].cspsDetails, false);

    // Add SPS contig information (if exists)
    generateOutputContigCsv(i, proteinIdx, spsID, page, cellPerLine, m_protein[proteinIdx].spsDetails, true);

    i += cellPerLine;
  }

	return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::loadInputData(void)
{
  // Protein reference file
  string refIndex = m_params.getValue("INPUT_FASTA");
	// get the references
  if(m_fasta.Load(refIndex.c_str()) == 0)
    return false;

  // Get the Contig Names
  string contigNamesFileName = m_params.getValue("SPS_CONTIG_NAMES");
  rtrim(contigNamesFileName);
  if(getContigNames(contigNamesFileName) < 0)
    return false;


  // auxiliary var to hold return values
  int ret;


 // **** SPS Contigs ****  //

  // Get SPS Contig matches index
  string spsContigMatchesIndexFileName = m_params.getValue("SPS_CONTIG_MATCHES_IDX");
  rtrim(spsContigMatchesIndexFileName);
	ret = m_spsContigMatchesIndex.loadPklBin(spsContigMatchesIndexFileName.c_str());
  if(ret < 0)
    return false;

  // Get SPS Contig matches
  string spsContigMatchesFileName = m_params.getValue("SPS_CONTIG_MATCHES");
  rtrim(spsContigMatchesFileName);
	ret = Load_binArray(spsContigMatchesFileName.c_str(), m_spsContigMatches);
  if(ret < 0)
    return false;

  // Get SPS Contig spectra
  string spsContigSpectraFileName = m_params.getValue("SPS_CONTIG_SPECTRA");
  rtrim(spsContigSpectraFileName);
	ret = m_spsContigSpectra.loadPklBin(spsContigSpectraFileName.c_str());
  if(ret < 0)
    return false;




 // **** CSPS Contig ****  //

  // Get CSPS Contig matches index
  string cspsContigMatchesIndexFileName = m_params.getValue("CSPS_CONTIG_MATCHES_IDX");
  rtrim(cspsContigMatchesIndexFileName);
	ret = m_cspsContigMatchesIndex.loadPklBin(cspsContigMatchesIndexFileName.c_str());
  if(ret < 0)
    return false;

  // Get CSPS Contig matches
  string cspsContigMatchesFileName = m_params.getValue("CSPS_CONTIG_MATCHES");
  rtrim(cspsContigMatchesFileName);
	ret = Load_binArray(cspsContigMatchesFileName.c_str(), m_cspsContigMatches);
  if(ret < 0)
    return false;

  // Get CSPS Contig spectra
  string cspsContigSpectraFileName = m_params.getValue("CSPS_CONTIG_SPECTRA");
  rtrim(cspsContigSpectraFileName);
	ret = m_cspsContigSpectra.loadPklBin(cspsContigSpectraFileName.c_str());
  if(ret < 0)
    return false;


  // csps DEBUG
  if(Logger::getLogger().getReportingLevel() >= 9) {

  	DEBUG_MSG("getContigMatchesIndex() has load " << m_cspsContigMatchesIndex.size());
  	for(int i = 0 ; i < m_cspsContigMatchesIndex.size() ; i++) {
    	DEBUG_MSG( i);
  		for(int j = 0 ; j < m_cspsContigMatchesIndex[i].size() ; j++) {
  	  	DEBUG_MSG("   " << i << "," << j << " --> " << m_cspsContigMatchesIndex[i][j][0]
  	  	<< ", " << m_cspsContigMatchesIndex[i][j][1]);
  	  }
  	}

   	DEBUG_MSG("getContigMatches() has load " << m_cspsContigMatches.size());
  	for(int i = 0 ; i < m_cspsContigMatches.size() ; i++) {
  		DEBUG_MSG(i << " --> " << m_cspsContigMatches[i][0] <<", " << m_cspsContigMatches[i][1]);
  	}

  	DEBUG_MSG("getContigSpectra() has load " << m_cspsContigSpectra.size());
  	for(int i = 0 ; i < m_cspsContigSpectra.size() ; i++) {
    	DEBUG_MSG(i);
  		for(int j = 0 ; j < m_cspsContigSpectra[i].size() ; j++){
  	  	DEBUG_MSG("   " << i << "," << j << " --> " << m_cspsContigSpectra[i][j][0] <<
  	  	", " << m_cspsContigSpectra[i][j][1]);
  	  }
  	}


    // sps
  	DEBUG_MSG("getContigMatchesIndex() has load " << m_spsContigMatchesIndex.size());
  	for(int i = 0 ; i < m_spsContigMatchesIndex.size() ; i++) {
    	DEBUG_MSG(i);
  		for(int j = 0 ; j < m_spsContigMatchesIndex[i].size() ; j++) {
  	  	DEBUG_MSG("   " << i << "," << j << " --> " << m_spsContigMatchesIndex[i][j][0] <<
  	  	", " << m_spsContigMatchesIndex[i][j][1]);
  	  }
  	}

    DEBUG_MSG("getContigMatches() has load " << m_spsContigMatches.size());
  	for(int i = 0 ; i < m_spsContigMatches.size() ; i++) {
  		DEBUG_MSG(i << " --> " << m_spsContigMatches[i][0] <<", " << m_spsContigMatches[i][1]);
  	}

  	DEBUG_MSG("getContigSpectra() has load " << m_spsContigSpectra.size());
  	for(int i = 0 ; i < m_spsContigSpectra.size() ; i++) {
    	DEBUG_MSG(i);
  		for(int j = 0 ; j < m_spsContigSpectra[i].size() ; j++){
  	  	DEBUG_MSG("   " << i << "," << j << " --> " << m_spsContigSpectra[i][j][0]
  	  	<< ", " << m_spsContigSpectra[i][j][1]);
  	  }
  	}
  }
// DEBUG ^^^

	return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportProteinCoverage::invoke(void)
{
  DEBUG_MSG("Entering  ExecReportProteinCoverage::invoke()");

  // Process proteins file
  processProteinsFile();

  // Process csps files
  processContigFiles(m_cspsContigMatchesIndex, m_cspsContigMatches, m_cspsContigSpectra, m_cspsContigInfo);
  if(Logger::getLogger().getReportingLevel() >= 9)
    dumpContig(m_cspsContigInfo);  // DEBUG

  // Process sps files
  processContigFiles(m_spsContigMatchesIndex, m_spsContigMatches, m_spsContigSpectra, m_spsContigInfo);
  if(Logger::getLogger().getReportingLevel() >= 9)
    dumpContig(m_spsContigInfo); // DEBUG

// get page construction parameter
  int aaPerLine = 20;
  if (m_params.exists("AA_PER_LINE"))
    aaPerLine = getInt(m_params.getValue("AA_PER_LINE").c_str());

  // Process each file
  for(int i = 0 ; i < m_protein.size() ; i++) {

    // second step csps processing
    //DEBUG_MSG("Entering  processContigs(" << i << ", ...)");
    processContigs(i, m_protein[i].cspsDetails, m_cspsContigInfo);

    // second step sps processing
    //DEBUG_MSG("Entering  processContigs(" << i << ", ...)");
    processContigs(i, m_protein[i].spsDetails, m_spsContigInfo);

    // Replace sps contig numbers by their names
    //DEBUG_MSG("Entering  populateContigNames(" << i << ", ...)");
		populateContigNames(m_protein[i].spsDetails);

    // Generate output HTML file
    //DEBUG_MSG("Entering  generateOutputFile2(" << i << ", ...)");
    generateOutputFile2(i, m_protein[i].page, aaPerLine);

    // Generate output CSV file
    //DEBUG_MSG("Entering  generateOutputFileCsv(" << i << ", ...)");
    generateOutputFileCsv(i, m_protein[i].csvFile, aaPerLine);

  }

  DEBUG_MSG("Exiting  ExecReportProteinCoverage::invoke()");

  // Exit program
  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////

