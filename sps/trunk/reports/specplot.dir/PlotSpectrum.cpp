///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <fstream>

#include "PlotSpectrum.h"
#include "util.h"
#include "spectrum_scoring.h"

#include "Defines.h"
#include "Timer.h"

///////////////////////////////////////////////////////////////////////////////
// Defines used in this file



#define PEAK_COLOR_NOT_ANNOTATED      "grey"
#define PEAK_COLOR_ANNOTATED_MULTIPLE "black"
#define PEAK_COLOR_ANNOTATED_SINGLE   "black"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {


///////////////////////////////////////////////////////////////////////////////
// Colors used:
//  [0]  label (y or b)
//  [1]  dashed upper line
//  [2]  mass dots
//  [3]  mass segment
//  [4]  peptides
//  [5]  vertical dashed line
//  [6]  annotations on peaks
//  [7]  annotation conneting lines
//  [8]  peaks when more then one peptide used

char *spectrumColorsB[][9] = {
  { "#0000C0", "#0000C0", "#0000C0", "#0000C0", "#0000C0", "#BFBFFF", "#0000C0", "#9090D0", "#0000C0" } ,
  { "#C00000", "#C00000", "#C00000", "#C00000", "#C00000", "#FFBFBF", "#C00000", "#D09090", "#C00000" } ,
  { "#00C000", "#00C000", "#00C000", "#00C000", "#00C000", "#BFFFBF", "#00C000", "#90D090", "#00C000" } ,
  { "#C0C000", "#C0C000", "#C0C000", "#C0C000", "#C0C000", "#D0D090", "#C0C000", "#D0D090", "#C0C000" } ,
  { "#00C0C0", "#00C0C0", "#00C0C0", "#00C0C0", "#00C0C0", "#00FFFF", "#00C0C0", "#90D0D0", "#00C0C0" }
};

char *spectrumColorsY[][9] = {
  { "#4900BF", "#4900BF", "#4900BF", "#4900BF", "#4900BF", "#70BFFF", "#4900BF", "#B89FDF", "#4900BF" } ,
  { "#BF4900", "#BF4900", "#BF4900", "#BF4900", "#BF4900", "#FDD8BF", "#BF4900", "#DF9FB8", "#BF4900" } ,
  { "#49BF00", "#49BF00", "#49BF00", "#49BF00", "#49BF00", "#D8FDBF", "#49BF00", "#9FDFB8", "#49BF00" } ,
  { "#BFBF49", "#BFBF49", "#BFBF49", "#BFBF49", "#BFBF49", "#F8FDBF", "#BFBF49", "#FFDF9F", "#BFBF49" } ,
  { "#49BFBF", "#49BFBF", "#49BFBF", "#49BFBF", "#49BFBF", "#D8FDBF", "#49BFBF", "#9FDFDF", "#49BFBF" }
};

///////////////////////////////////////////////////////////////////////////////
// Constructors and destructor
///////////////////////////////////////////////////////////////////////////////
PlotSpectrum::PlotSpectrum(void) :
  m_spectrum(NULL),
  m_range_min(-1.0),
  m_range_max(-1.0),
  m_trim(false),
  m_annotationByIntensity(false),
  m_offsetN(AAJumps::massHion),
  m_offset(AAJumps::massH2O + AAJumps::massHion)
{
}
///////////////////////////////////////////////////////////////////////////////
PlotSpectrum::~PlotSpectrum(void)
{
}
///////////////////////////////////////////////////////////////////////////////
// Annotation definition method
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::defineMainAnnotations(void)
{
  for(int i = 0 ; i < m_model.probs.size() ; i++)
    if(m_model.probs[i].isMainIon)
      if(m_model.probs[i].isNTerm) {
        m_mainIonN = m_model.probs[i].label;
        m_offsetN  = m_model.probs[i].massOffset;
      } else {
        m_mainIon = m_model.probs[i].label;
        m_offset  = m_model.probs[i].massOffset;
      }
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::defineAnnotations(void)
{
  for(int i = 0 ; i < m_model.probs.size() ; i++)
    if(m_model.probs[i].isBreakIon)
      if(m_model.probs[i].isNTerm) {
        ionsSearchDotsB.push_back(m_model.probs[i].name);
        ionsSearchSegmentsB.push_back(m_model.probs[i].name);
        if(m_model.probs[i].charge == 1)
          ionsSearchVerticalB.push_back(m_model.probs[i].name);

      } else {
        ionsSearchDotsY.push_back(m_model.probs[i].name);
        ionsSearchSegmentsY.push_back(m_model.probs[i].name);
        if(m_model.probs[i].charge == 1)
          ionsSearchVerticalY.push_back(m_model.probs[i].name);
      }

  //////////////////////////////////////////////////////////////////////////////
  // b ions anotation section

  // Draw b-ion upper dots
  //ionsSearchDotsB.push_back("b");
  //ionsSearchDotsB.push_back("b++");

  // Draw b-ion upper segments
  //ionsSearchSegmentsB.push_back("b");
  //ionsSearchSegmentsB.push_back("b++");

  // Draw b-ion vertical lines
  //ionsSearchVerticalB.push_back("b");


  //////////////////////////////////////////////////////////////////////////////
  // y ions anotation section

  // Draw y-ion upper dots
  //ionsSearchDotsY.push_back("y");
  //ionsSearchDotsY.push_back("y++");

  // Draw b-ion upper segments
  //ionsSearchSegmentsY.push_back("y");
  //ionsSearchSegmentsY.push_back("y++");

  // Draw b-ion vertical lines
  //ionsSearchVerticalY.push_back("y");
}
//////////////////////////////////////////////////////////////////////////////
// Define peak annotations section
void PlotSpectrum::definePeakAnnotations(vector<DrawingData> &ionLabels)
{
  // specify peaks to get -- (prefix - prefix/full search - show - label color - callout color - use B series color schema )
  //
  // [0] search string - search string used to search for annotation
  // [1] prefix/full search - if true, search string is used as a prefix. If false, the annotation must be exactly equal to the search string
  // [2] show - what to show in the label for this annotation
  // [3] color index (in color table) used to draw the label. -1 for black
  // [4] color index (in color table) used to draw the callout line. -1 for black
  // [5] color schema used. True to use B series, false to use Y series
  string label, name;
  int colorIdx1, colorIdx2;
  bool nTerm;

  for(int i = 0 ; i < m_model.probs.size() ; i++)
    if(m_model.probs[i].hasLabel) {
      name  = m_model.probs[i].name;
      label = (m_model.probs[i].label.length() ? m_model.probs[i].label : m_model.probs[i].name);
      colorIdx1 = m_model.probs[i].isIF ? 6 : -1;
      colorIdx2 = m_model.probs[i].isIF ? 7 : -1;
      nTerm = m_model.probs[i].isNTerm;

      ionLabels.push_back( DrawingData(name.c_str(), false, label.c_str(), colorIdx1, colorIdx2, nTerm) );
    }
/*
  // y series
  ionLabels.push_back( DrawingData("y",       false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y++",     false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z",     false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z1",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z2",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z3",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z4",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z5",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z6",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z7",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z8",    false, "y", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("y+z9",    false, "y", 6,   7, false  ) );

  // b series
  ionLabels.push_back( DrawingData("b",       false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b++",     false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z",     false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z1",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z2",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z3",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z4",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z5",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z6",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z7",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z8",    false, "b", 6,   7, true   ) );
  ionLabels.push_back( DrawingData("b+z9",    false, "b", 6,   7, true   ) );

  // a series
  ionLabels.push_back( DrawingData("a",       false, "a", 6,   7, false  ) );
  ionLabels.push_back( DrawingData("a++",     false, "a", 6,   7, false  ) );

  // P series
  ionLabels.push_back( DrawingData("P++",     false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z",     false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z1",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z2",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z3",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z4",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z5",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z6",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z7",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z8",    false, "P", -1, -1, false  ) );
  ionLabels.push_back( DrawingData("P+z9",    false, "P", -1, -1, false  ) );

  //ionLabels.push_back( DrawingData("P++-H2O", false, "P", -1, -1, false  ) );
  //ionLabels.push_back( DrawingData("P++-NH3", false, "P", -1, -1, false  ) );
*/
}
///////////////////////////////////////////////////////////////////////////////
// Auxiliary methods
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::getColors(int index, char **cb[], char **cy[])
{
  index %= 5;
  *cb = spectrumColorsB[index];
  *cy = spectrumColorsY[index];
}
///////////////////////////////////////////////////////////////////////////////
char *PlotSpectrum::getColor(int index, int position, bool b)
{
  index %= 5;
  if(b)
    return spectrumColorsB[index][position];
  return spectrumColorsY[index][position];
}
///////////////////////////////////////////////////////////////////////////////
// Tests if a given mass point is within parameters
bool PlotSpectrum::testMass(double &value)
{
  return ((value >= m_mzLowerLimit) && (value <= m_mzUpperLimit));
}
///////////////////////////////////////////////////////////////////////////////
// Tests if a given mass set is within parameters
bool PlotSpectrum::testMass(double &left, double &right)
{
  return ((right >= m_mzLowerLimit) && (left <= m_mzUpperLimit));
}
///////////////////////////////////////////////////////////////////////////////
// Tests if a given mass set is within parameters.
// Special case for unkown mzUpperLimit, but known user limit, if available
bool PlotSpectrum::testMass2(double &mass)
{
  return ((mass >= m_mzLowerLimit) && (m_range_max > 0.0 ? mass <= m_range_max : true));
}
///////////////////////////////////////////////////////////////////////////////
// Drawing procedure entry point with no parameters
///////////////////////////////////////////////////////////////////////////////
int PlotSpectrum::drawExec(void)
{
  // Check for spectrum
  if(m_spectrum == NULL)
    return ERROR;

  /////////////////////////////////////////////////////////////////////////////
  // set locations

  // Set the renderer location
  m_rendererObject->setRendererLocation(m_rendererLocation);
  // set the font files location
  m_rendererObject->setFontLocation(m_fontLocation);


  /////////////////////////////////////////////////////////////////////////////
  // Load annotation model

	// Define annotation file location
	string annotationFile = m_annotationModelDirectory;
	annotationFile += '/';
	// define correct annotation model
	if(m_readModelFromData) {
	  switch(m_spectrum->msFragType) {
	  case Spectrum::FragType_PRM:
	    //annotationFile += ANNOTATION_MODEL_PRM;
	    annotationFile += ANNOTATION_MODEL_CID;
	    m_massShift -= 1.0;
	    break;
	  case Spectrum::FragType_ETD:
	    annotationFile += ANNOTATION_MODEL_ETD;
	    break;
    case Spectrum::FragType_CID:
	    annotationFile += ANNOTATION_MODEL_CID;
	    break;
    case Spectrum::FragType_HCD:
	    annotationFile += ANNOTATION_MODEL_CID;
	    break;
    default:
	    annotationFile += DEFAULT_ANNOTATION_MODEL;
	    break;
    }
	} else
  	annotationFile += m_annotationModel;
	// Set annotatioln file
	m_model.LoadModel((char*)annotationFile.c_str());

  /////////////////////////////////////////////////////////////////////////////
  // output filename

  // set output file name
  if(m_fileOutDir.length() > 0) {
    // Set output directory, if defined
    m_rendererImage.outputFileName = m_fileOutDir;
    if(m_fileOutDir[m_fileOutDir.length()-1] != '/')
      m_rendererImage.outputFileName += '/';
  }

  if(m_fileNameSuffix.length()) {

    int found = m_fileName.find_last_of(".");
    // filename
    string fn = m_fileName.substr(0,found);
    fn += m_fileNameSuffix;
    fn += m_fileName.substr(found);

    // add filename suffix
    m_rendererImage.outputFileName += fn;

  } else
    // Add the filename
    m_rendererImage.outputFileName += m_fileName;



  /////////////////////////////////////////////////////////////////////////////
  // draw procedure

  // Init draw object
  setDrawProperties();

  // Draw procedure
  if(plotSpectrum() == ERROR)
    return ERROR;

  /////////////////////////////////////////////////////////////////////////////
  // define broken x axis
  /*
  // find first and second greates peaks
  double first, second;
  first = second = 0.0;

  for(int i = 0 ; i < m_spectrum->size() ; i++) {
    float aux = (*m_spectrum)[i][1];
    if(aux > first) {
      second = first;
      first = aux;
    } else if(aux > second)
      second = aux;
  }

  //cout << m_spectrumIndex << "    " << second / first << endl;

  // break axis if first greater then twice the second
  if(second / first < 0.50) {
    second *= 1.2; // gap of 20%
    first  *= 0.9; // gap of 10%
    // set broken axis vars
    m_rendererObject->breakAxisY(second, first);
  }
  */

  /////////////////////////////////////////////////////////////////////////////
  // Do the drawing

  Timer_c timer;

  m_rendererObject->execute();

  if(m_verbose)
    cout << "Time used by renderer object: " << timer.stop() << endl;
    

  // Dump masses (DEBUG)
  //dumpPeptideAndMasses();

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
// Draw the object - main draw procedure
///////////////////////////////////////////////////////////////////////////////
int PlotSpectrum::plotSpectrum(void)
{
  // defines main annotation ions
  defineMainAnnotations();
  // Sets the annotations vectors
  defineAnnotations();

  // initialize spectrumData data structure from Spectrum object
  initSpectrumData();

  // Amino acid masses object
  // note that the 'max allowed jumps' parameter has to be set, but we're not using it for anything.
  AAJumps jumps(1);
  // if a amino acids file was specified, use it
  if(m_aminoAcidsFile.length()) {
    string aux = m_aminoAcidsFileDirectory;
    aux += '/';
    aux += m_aminoAcidsFile;
//    cout << " loading jumps at " << aux << endl;
    jumps.loadJumps(aux.c_str(), true);
  }

  //psmPtr psm(new PeptideSpectrumMatch);
  psmPtr psm;
  if(m_spectrum->psmList.size()) {
    // get first PSM
    psm = m_spectrum->psmList.front();
    //Set spectrum
    psm->m_spectrum = m_spectrum;
  } else {
    psm = (psmPtr)new PeptideSpectrumMatch;
    //Set spectrum
    psm->m_spectrum = m_spectrum;
    //Set PSM on spectrum (probably not necessary, but might make things easier
    m_spectrum->psmList.push_back(psm);
  }


  // Annotate the spectrum for all peptides
  for(int i = 0 ; i < m_peptide.size() ; i++) {

    if(m_peptide[i].length() == 0 || m_peptide[i].compare("*..*") == 0) {
      m_peptide.erase(m_peptide.begin() + i);
      i--;
      continue;
    }

    // translate annotation style, if needed
    // check for inspect style flag
    if(m_annotationStyle == ANNOTATION_STYLE_INSPECT) {
      // auxiliary string
      string aux;
      // check leading *.
      if(m_peptide[i].size() == 1)
        aux = "*.";
      else if(m_peptide[i][1] != '.')
        aux = "*.";
      // add peptide
      aux += m_peptide[i];
      // check trailing .*
      if(m_peptide[i].size() == 1)
        aux += ".*";
      else if(m_peptide[i][m_peptide[i].size() - 2] != '.')
        aux += ".*";
      // where to store specnets peptide
      string peptideSpecnets;
      // translate
      psm->inspectToSpecNets(aux, peptideSpecnets);
      // put translated peptide back
      m_peptide[i] = peptideSpecnets;

    } else {

      if(m_peptide[i].size() > 1 && m_peptide[i][1] == '.') {
        m_peptide.erase(m_peptide.begin() + i);
        i--;
        continue;
      }
    }

    // Annotate one peptide
    annotate(m_peptide[i], psm);

    // Compute masses needed for drawing
    bool sequenceOK = jumps.getPRMandSRMMasses(m_peptide[i], m_prmMasses[i], m_srmMasses[i], m_peptideMass[i]);

    if(!sequenceOK) {
      cerr << "Invalid sequence: " << m_peptide[i] << endl;
      return ERROR;
    }

//    cout << endl << " ... prm masses" << endl;
//
//    for(int j = 0 ; j < m_prmMasses[i].size() ; j++)
//      cout << m_prmMasses[i][j] << endl;
//
//    cout << endl << " ... srm masses" << endl;
//
//    for(int j = 0 ; j < m_srmMasses[i].size() ; j++)
//      cout << m_srmMasses[i][j] << endl;
//
//    cout << endl << " ... peptide masses " << endl;
//
//    for(int j = 0 ; j < m_peptideMass.size() ; j++)
//      cout << m_peptideMass[j] << endl;

    // Add annotation data to internal representation
    addAnnotationData(i);
  }

  // Write annotation file, if it was specified
  if(m_annotationsOutputFile.length()) {
    buildAnnotationOutput();
    writeAnnotationFile();
  }


  // Calculate M/Z and Intensity limits
  calcultateLimits();

  // Initialize the graph
  initializeGraph();

  // Calculates where to draw the y and b labels for the various peptides
  calcDrawPosition();

  // Set the graph title
  setGraphTitle();

  // Draw upper peptide annotations
  drawAnnotationPeptides();

  // Draw spectrum labels
  drawPeakLabelsMain();

  // Draw spectrum
  drawSpectrumPeaks();
}
///////////////////////////////////////////////////////////////////////////////
// Draw Upper peptide annotations
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawAnnotationPeptides(void)
{
  // Build data structures for the upper mass dots and vertical lines
  for(int i = 0 ; i < m_peptide.size() ; i++) {
      // Build annotation vector (prm)
    buildAnnotationVector(i, m_peptide[i], m_peptideMass[i], m_prmMasses[i], peptideAnnotations[i], false);
    // Build annotation vector (srm)
    buildAnnotationVector(i, m_peptide[i], m_peptideMass[i], m_srmMasses[i], peptideAnnotationsReverse[i], true);
  }

  // draw annotations
  for(int index = 0 ; index < peptideAnnotations.size() ; index++)
    drawAnnotations(index);
}
///////////////////////////////////////////////////////////////////////////////
// Draw peak labels
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawSpectrumPeaks(void)
{
  // Prepare peaks for output
  preparePeaksToOutput();

  // Draw the spectrum lines
//  if(m_peptide.size() == 1)
    drawSpectrumLinesFast();
//  else
//    drawSpectrumLines();
}
///////////////////////////////////////////////////////////////////////////////
// Draw peak labels
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawPeakLabelsMain(void)
{
  // Aquire peak labels
  vector<DrawingData> ionLabels;
  // Get peak annotations
  definePeakAnnotations(ionLabels);
  // Get the peaks
  aquirePeakLabels(ionLabels);

  // find placement for the labels
  placePeakLabels();
  // draw the peak labels
  drawPeakLabels();
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::initSpectrumData(void)
{
  // cycle tru all SPectrum peaks
  for(int i = 0 ; i < m_spectrum->size() ; i++) {

    // get data
    m_spectrumData[i].peakMass =      (*m_spectrum)[i][0];
    m_spectrumData[i].peakIntensity = (*m_spectrum)[i][1];
  }
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::addAnnotationData(int index)
{
  // The annotation data item structure
  SpectrumAnnotationData anno;
  // Peptide index is the same to all elements
  anno.peptideIndex = index;
  // Cycle thru annotations

//  int size = m_spectrum->size();
//  for(register int j = 0 ; j <  size; j++)
//    if(m_spectrum->annotation[j].first != 0) {
//      anno.charge       = m_spectrum->annotation[j].first->charge;
//      anno.seriesIndex  = m_spectrum->annotation[j].second;
//      anno.annotation   = m_spectrum->annotation[j].first->name;
//      anno.prob         = m_spectrum->annotation[j].first->prob;
//      m_spectrumData[j].annotations.push_back(anno);
//    }


  int size = m_spectrum->size();
  psmPtr currAnnotation = m_spectrum->psmList.front();

  //cout << "Annotations: " << currAnnotation->m_peakAnnotations.size() << endl;

  for(register int j = 0 ; j <  size && j < currAnnotation->m_peakAnnotations.size() ; j++)
    if(currAnnotation->m_peakAnnotations[j].first != 0) {
      anno.charge       = currAnnotation->m_peakAnnotations[j].first->charge;
      anno.seriesIndex  = currAnnotation->m_peakAnnotations[j].second;
      anno.annotation   = currAnnotation->m_peakAnnotations[j].first->name;
      anno.prob         = currAnnotation->m_peakAnnotations[j].first->prob;
      m_spectrumData[j].annotations.push_back(anno);

      //cout << "       Anno: " << anno.annotation << endl;
    }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::calcDrawPosition(void)
{
  // space for both y a b peptide lines
  float peptidesSpace = 3.5;
  // pixels to map transformation factor
  float factor =  m_rendererImage.fontHeight * m_rendererImage.pixelHeight / m_rendererImage.zoom;
  // space between title and peptide
  float interleave = m_rendererImage.zoom >= 1.0 ? 2.0 : 1.0;
  // top margin position (adds space for title)
  float topFactor = (m_rendererImage.fileType.compare("eps") ? 1.5 : 0.5);
  m_topMarginPosition = 1.0 - m_titlePos * topFactor * factor;

  //cout << "-------------------------------------------" << endl;
  //cout << "fontHeight: " << m_rendererImage.fontHeight << endl;
  //cout << "pixelHeight: " << m_rendererImage.pixelHeight << endl;
  //cout << "Factor: " << factor << endl;
  //cout << "m_topMarginPosition: " << m_topMarginPosition << endl;

  // calculate positions
  for(int index = 0 ; index < m_peptide.size() ; index++) {
    float bPos = m_topMarginPosition - factor * (interleave  + (index * peptidesSpace));
    float yPos = m_topMarginPosition - factor * (interleave + 1.5  + (index * peptidesSpace));
    // store positions
    m_headerPosition[index] = make_pair(bPos, yPos);

    //cout << "bPos: " << bPos << endl;
    //cout << "yPos: " << yPos << endl;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::initializeGraph(void)
{
  // Calculate the margins
  float marginLeft, marginRight, marginTop, marginBottom;
  // size occupied by a peptide (b and y lines)
  float peptidesSpace = m_peptide.size() * 3.5;
  // space between title and peptide
  float interleave = m_rendererImage.zoom >= 1.0 ? 2.0 : 1.0;

  marginLeft    = 12.0;
  marginRight   = 3.0;
  // Different values for different image sizes
  marginTop     = m_titlePos * 2.0 + peptidesSpace + interleave;
  marginBottom  = (m_rendererImage.fileType.compare("eps") ? 3.0 : 6.0); //m_rendererImage.zoom >= 1.0 ? 6.0 : 4.0;

  // Calculate title offsets
  m_titleOffsetX = 0;
  m_titleOffsetY = peptidesSpace + interleave;

  // Set the margins
  m_rendererObject->setMarginLeft(marginLeft);
  m_rendererObject->setMarginRight(marginRight);
  m_rendererObject->setMarginTop(marginTop);
  m_rendererObject->setMarginBottom(marginBottom);

  // set graph size in pixels
  graphSizePixelsX = m_rendererImage.imageSizeX - (int)(marginLeft - marginRight);
  graphSizePixelsY = m_rendererImage.imageSizeY - (int)(marginTop - marginBottom);

  // save space on graph for one label above the higher peak
  // pixels to map transformation factor
  float factor =  m_rendererImage.fontHeight * m_rendererImage.pixelHeight / m_rendererImage.zoom;
  // space for label in graph units
  float space = 2.5 * factor;
  // graph space in graph units
  float graphSpace = (1 - (marginTop - marginBottom) * factor) / m_intensityLimit;
  // space for label percentage of data in graph units
  float labelSpace = space / graphSpace;


  //cout << "m_intensityLimit: " << m_intensityLimit << endl;

  m_intensityLimit += labelSpace;

  //cout << "factor:         " << factor << endl;
  //cout << "space:          " << space << endl;
  //cout << "graphSpace:     " << graphSpace << endl;
  //cout << "labelSpace:     " << labelSpace << endl;
  //cout << "m_intensityLimit: " << m_intensityLimit << endl;


  // Graph coverage
  m_rendererObject->setXrange(m_mzLowerLimit, m_mzUpperLimit);
  m_rendererObject->setYrange(0, m_intensityLimit);

  // Border location
  m_rendererObject->setBorder(3);

  // Labels if zoom >= 1
  if( m_rendererImage.zoom >= 1.0 ) {
    m_rendererObject->setXlabel("Mass / Charge (m/z)");
    m_rendererObject->setYlabel("Intensity");
  }

  // set number of labels on the scales
  float f1 = 7.0;
  float f2 = 8.0;

  //f1 = ( graphSpace * factor < f1 ? : f1 );
  float graphSize = graphSpace * m_intensityLimit / (space * 2.0);
  f2 = fmax(fmin(graphSize, f2), 3.0 );

   // Calculate intervals
  double msdX = m_rendererObject->calculateInterval(m_mzUpperLimit - m_mzLowerLimit, f1, m_rendererImage.zoom);
  double msdY = m_rendererObject->calculateInterval(m_intensityLimit, f2, m_rendererImage.zoom);
  // Set intervals
  m_rendererObject->setXinterval(msdX);
  m_rendererObject->setYinterval(msdY);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::setGraphTitle(void)
{
  // Set the title
  if(m_titlePos) {
    //stringstream aux;
    //aux << "Consensus spectrum " << m_spectrumIndexGlobal << " (Contig " << m_contig << ")";
    // Title object for plot object
    RendererTitle title;
    title.label     = m_title.c_str(); //aux.c_str();
    title.offsetX   = m_titleOffsetX;
    title.offsetY   = m_titleOffsetY;
    title.fontName  = m_rendererImage.fontName;
    title.fontSize  = int(m_rendererImage.fontSize + 2);
    // Do it
    m_rendererObject->setTitle(title);
    //m_rendererObject->setTitle(aux.str().c_str());
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawAnnotations(int index)
{
  // constants used
  //double bShift = AAJumps::massHion + m_massShift;
  //double yShift = AAJumps::massH2O + AAJumps::massHion + m_massShift;
  double bShift = m_offsetN + m_massShift;
  double yShift = m_offset + m_massShift;

  // Get the colors used for this index
  char **cb, **cy;
  getColors(index, &cb, &cy);

  // draw the "b" label
  drawUpperGuideLabel(m_mainIonN, m_headerPosition[index].first, cb[0]);
  // draw the "y" label
  drawUpperGuideLabel(m_mainIon, m_headerPosition[index].second, cy[0]);

  // draw the upper dashed line
  float aux1 = max(m_mzLowerLimit, bShift);
  float aux2 = min(m_mzUpperLimit, m_peptideMass[index] + bShift);
  drawUpperGuideLine(aux1, aux2, m_headerPosition[index].first, cb[1]);
  // draw the lower dashed line
  aux1 = max(m_mzLowerLimit, yShift);
  aux2 = min(m_mzUpperLimit, m_peptideMass[index] + yShift);
  drawUpperGuideLine(aux1, aux2, m_headerPosition[index].second, cy[1]);



  // Data used in the following operations
  ParamsData1 paramsData(&ionsSearchDotsB, &peptideAnnotations[index], &m_peptide[index]);
  paramsData.shift        = bShift;
  paramsData.peptideMass  = m_peptideMass[index];
  paramsData.yPosition    = m_headerPosition[index].first;

  // Draw b-ion upper labels
  paramsData.color  = cb[4];
  drawUpperLabelsAll(paramsData);

  // Draw b-ion upper dots
  paramsData.ions   = &ionsSearchDotsB;
  paramsData.color  = cb[2];
  drawUpperMassDots(paramsData);

  // Draw b-ion upper segments
  paramsData.ions   = &ionsSearchSegmentsB;
  paramsData.color  = cb[3];
  drawUpperSegments(paramsData);

  // Draw b-ion vertical lines
  paramsData.ions   = &ionsSearchVerticalB;
  paramsData.color  = cb[5];
  drawVerticalLines(paramsData);


  // Update parameter data for next operations
  paramsData.shift              = yShift;
  paramsData.yPosition          = m_headerPosition[index].second;
  paramsData.peptideAnnotation  = &peptideAnnotationsReverse[index];

  // Draw y-ion upper labels
  paramsData.color  = cy[4];
  drawUpperLabelsAll(paramsData);

  // Draw y-ion mass dots
  paramsData.ions   = &ionsSearchDotsY;
  paramsData.color  = cy[2];
  drawUpperMassDots(paramsData);

  // Draw y-ion upper segments
  paramsData.ions   = &ionsSearchSegmentsY;
  paramsData.color  = cy[3];
  drawUpperSegments(paramsData);

  // Draw y-ion vertical lines
  paramsData.ions   = &ionsSearchVerticalY;
  paramsData.color  = cy[5];
  drawVerticalLines(paramsData);
}
///////////////////////////////////////////////////////////////////////////////
// Draw spectrum lines using traditional method (slow)
///////////////////////////////////////////////////////////////////////////////
void  PlotSpectrum::drawSpectrumLines()
{
  // Draw spectrum lines using traditional method (slow)
  RendererLine line;
  line.start.x.type = RENDERER_COORD_NONE;
  line.start.y.type = RENDERER_COORD_NONE;
  line.end.x.type   = RENDERER_COORD_NONE;
  line.end.y.type   = RENDERER_COORD_NONE;

  line.width = 2;
  line.type  = 1;
  // Cycle thru spectrum items
  for(int j = 0 ; j < peaksOutput.size() ; j++) {
    // fill the line object
    line.start.x.value = peaksOutput[j].mass;
    line.start.y.value = 0;
    line.end.x.value = peaksOutput[j].mass;
    line.end.y.value = peaksOutput[j].intensity;
    line.color = peaksOutput[j].color;
    // Draw the line
    m_rendererObject->drawLine(line);
  }
}
///////////////////////////////////////////////////////////////////////////////
// Draw the graph using binary array method (fast)
///////////////////////////////////////////////////////////////////////////////
void  PlotSpectrum::drawSpectrumLinesFast(void)
{
  vector<RendererData> plotData;
  RendererData  rendererData;

  rendererData.lineWidth   = 2;
  rendererData.lineType    = 1;
  rendererData.dataType    = RENDERER_DATA_TYPE_DOUBLE;
  rendererData.dataRecords = 0;

  int size = peaksOutput.size();
  // Cycle thru spectrum items
  for(register int j = 0 ; j < size ; j++) {
    double peak = peaksOutput[j].mass;
    if(testMass(peak)) {
      if(rendererData.lineColor.compare(peaksOutput[j].color)) {
        if(j)
          plotData.push_back(rendererData);
        rendererData.dataDouble.clear();
        rendererData.lineColor = peaksOutput[j].color;
        rendererData.dataRecords = 0;
      }
      rendererData.dataDouble.push_back(make_pair<double,double>(peak, peaksOutput[j].intensity));
      rendererData.dataRecords++;
    }
  }

  // Add last peaks
  plotData.push_back(rendererData);

  // Draw the graph
  m_rendererObject->drawData(plotData);

/*
  vector<RendererData> allPlotData;
  RendererData  rendererDataBlack, rendererDataGrey;

  rendererDataBlack.lineWidth   = 2;
  rendererDataBlack.lineType    = 1;
  rendererDataBlack.dataType    = RENDERER_DATA_TYPE_DOUBLE;
  rendererDataBlack.lineColor   = PEAK_COLOR_ANNOTATED_SINGLE;
  rendererDataBlack.dataRecords = 0;

  rendererDataGrey.lineWidth    = 2;
  rendererDataGrey.lineType     = 1;
  rendererDataGrey.dataType     = RENDERER_DATA_TYPE_DOUBLE;
  rendererDataGrey.lineColor    = PEAK_COLOR_NOT_ANNOTATED;
  rendererDataGrey.dataRecords  = 0;

  int size = m_spectrumData.size();
  // Cycle thru spectrum items
  for(register int j = 0 ; j < size ; j++) {
    double peak = m_spectrumData[j].peakMass;
    if(testMass(peak)) {
      // Spectrum line color depends on the presence of annotation
      if(m_spectrumData[j].annotations.size() == 0) {
        rendererDataGrey.dataDouble.push_back(make_pair<double,double>(peak, m_spectrumData[j].peakIntensity));
        rendererDataGrey.dataRecords++;
      } else {
        rendererDataBlack.dataDouble.push_back(make_pair<double,double>(peak, m_spectrumData[j].peakIntensity));
        rendererDataBlack.dataRecords++;
      }
    }
  }

  // Add grey and black graph peaks
  allPlotData.push_back(rendererDataGrey);
  allPlotData.push_back(rendererDataBlack);

  // Draw the graph
  m_rendererObject->drawData(allPlotData); */
}
///////////////////////////////////////////////////////////////////////////////
// Spectra annotation procedure
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::annotate(string &peptide, psmPtr psm)
{
  // set the peptide for charge determination
  psm->m_annotation = peptide;
  //Calculate charge and set on psm:
  if(psm->m_charge <= 0)
     psm->setChargeByAnnotation();
  // ions to retrieve
	string allIons("all");
  // mass shift
  float aux = m_massShift;
  // Annotate the spectrum
  psm->annotate(peptide, allIons, m_model, aux, aux, m_peakMassTol);

}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::preparePeaksToOutput(void)
{
  // Elemt structure used
  spectrumOutputElem elem;

  for(int j = 0 ; j < m_spectrumData.size() ; j++) {
    // Get peak m/z value to speed up
    double mass = m_spectrumData[j].peakMass;
    // ignore peaks out of range
    if(testMass(mass)) {
      // store mass and intensity
      elem.mass       = mass;
      elem.intensity  = m_spectrumData[j].peakIntensity;

      // peak coloring. If the user didn't specify peptides, color all in black
      if(!m_peptide.size()) {
        elem.color      = PEAK_COLOR_ANNOTATED_MULTIPLE;
        elem.annotated  = false;

      } else {
        // If peptides were specified, spectrum line color depends on the presence of annotations
        switch(m_spectrumData[j].annotations.size()) {
        case 0:
          elem.color      = PEAK_COLOR_NOT_ANNOTATED;
          elem.annotated  = false;
          break;
        case 1: {
          // Find annotated peak with greater probability
          int maxProbIndex = -1;
          float maxProb = -1.0;
          for(int k = 0 ; k < m_spectrumData[j].annotations.size() ; k++)
            if(m_spectrumData[j].annotations[k].prob > maxProb) {
              maxProbIndex = k;
              maxProb = m_spectrumData[j].annotations[k].prob;
            }
          if(maxProbIndex >= 0) {
            // determine if we have a B series
            bool ccc = false;
            if(m_spectrumData[j].annotations[maxProbIndex].annotation.length())
              if(m_spectrumData[j].annotations[maxProbIndex].annotation[0] == 'b')
                ccc = true;
            // First is annotation index, second is color index 8 =
            elem.color    = getColor(m_spectrumData[j].annotations[maxProbIndex].peptideIndex, 8, ccc);
          }
          elem.annotated  = true;
          break;
        }
        default:
          elem.color      = PEAK_COLOR_ANNOTATED_MULTIPLE;
          elem.annotated  = true;
          break;
        }
      }
      // Store the peak
      peaksOutput.push_back(elem);
    }
  }
  sort(peaksOutput.begin(), peaksOutput.end());
}
///////////////////////////////////////////////////////////////////////////////
// Calcultate M/Z and intensity upper limits based on spectrum input data
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::calcultateLimits(void)
{
  // Initialize M/Z lower limit variable
  m_mzLowerLimit = 0.0;

  // Initialize intensity variable
  m_intensityLimit = 0.0;

  // use peak list size in an auxiliary variable for speed purposes
  int peakListSize = m_spectrumData.size();

  // Cycle thru intensity values to find max
  for(int i = 0 ; i < peakListSize ; i++) {
    double mz = m_spectrumData[i].peakMass;
    // test if mz is between the range window - min value
    if(m_range_min > 0.0)
      if(mz < m_range_min)
        continue;
    // test if mz is between the range window - max value
    if(m_range_max > 0.0)
      if(mz > m_range_max)
        continue;
    double aux = m_spectrumData[i].peakIntensity;
    if(aux >  m_intensityLimit)
      m_intensityLimit = aux;
  }


  // Size added by 5% to fit one label.
  //m_intensityLimit *= 1.05;

  // if the user has defined the m/z range upper limit, use it. If not, calculate it.
  if(m_range_max > 0.0) {

    m_mzUpperLimit = m_range_max;

  // use value -2 to specify parent_mz
  } else  {

    // Initialize m/z upper limit variable
    m_mzUpperLimit = 0.0;

    // Cycle thru m/z values to find max
    for(int i = 0 ; i < peakListSize ; i++) {
      double aux = m_spectrumData[i].peakMass;
      if(aux >  m_mzUpperLimit)
        m_mzUpperLimit = aux;
    }
    // The upper limit must be corrected according to pre-calculated peptide mass.
    for(int i = 0 ; i < m_peptideMass.size() ; i++)
      if(m_peptideMass[i] > m_mzUpperLimit)
        m_mzUpperLimit = m_peptideMass[i];

    // Added 5% for space purposes.
    m_mzUpperLimit *= 1.05;
  }

  if(m_peptide.size() < 1 && m_trim)
    m_mzUpperLimit = m_spectrum->parentMass;


  // if the user has defined the m/z range lower limit, use it.
  if(m_range_min > 0.0)
    m_mzLowerLimit = m_range_min;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideLabel(string theLabel, double coord, string color)
{
  drawLabel(theLabel, 0,      RENDERER_COORD_GRAPH,
                      coord,  RENDERER_COORD_SCREEN,
                      -2.0, 0,
                      RENDERER_OFFSET_RIGHT,
                      color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideLabelSegments(string theLabel, double x, double y, string color)
{
  drawLabel(theLabel, x,      RENDERER_COORD_NONE,
                      y,      RENDERER_COORD_SCREEN,
                      0.0,    0.5,
                      RENDERER_OFFSET_CENTER,
                      color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideLine(double x1, double x2, double y1, string color)
{
  drawLine(  1, 0,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideDots(double x, double y, int size, string color)
{
  // if type is EPS, use specific method
  if(m_fileFormat.compare("eps") == 0) {
    drawUpperGuideDotsEps(x, y, size, color);
    return;
  }


  drawLine(  size, 1,
             x, RENDERER_COORD_NONE,
             y, RENDERER_COORD_SCREEN,
             x, RENDERER_COORD_NONE,
             y, RENDERER_COORD_SCREEN,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideDotsEps(double x, double y, int size, string color)
{
  double newY1 = y - (size == 2 ? 0.0025 : 0.005);
  double newY2 = y + (size == 2 ? 0.0025 : 0.005);
  int newSize = (size == 2 ? 8 : 16);

  drawLine(  newSize, 1,
             x,     RENDERER_COORD_NONE,
             newY1, RENDERER_COORD_SCREEN,
             x,     RENDERER_COORD_NONE,
             newY2, RENDERER_COORD_SCREEN,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideSegments(double x1, double x2, double y1, string color)
{
  drawLine(  2, 1,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideVertical(double x1, double y1, string color)
{
  drawLine(  1, 0,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x1, RENDERER_COORD_NONE,
             0,  RENDERER_COORD_NONE,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperGuideVertical2(double x1, double y1, double x2, double y2, string color)
{
  drawLine(  1, 0,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_NONE,
             color);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
int PlotSpectrum::peakAnnotated(vector<DrawingData> &series, string &annotation)
{
  for(int j = 0 ; j < series.size() ; j++) {
    // If it is a prefix, just compare the first characters
    if(series[j].isPrefix) {
      if(annotation.compare(0, series[j].ion.length(), series[j].ion) == 0)
        return j;
    } else
      // Otherwise, compare the whole string
      if(annotation.compare(series[j].ion) == 0)
        return j;
  }
  return -1;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
int PlotSpectrum::findAnnotatedIndex(PeptideAnnotation &annotationData, vector<string> &series, bool chargeOne)
{
  // Annotations vector size
  register int annotationsSize = annotationData.annotations.size();
  // Search for the annotation
  for(register int j = 0 ; j < series.size() ; j++)
    for(register int i = 0 ; i < annotationsSize ; i++)
      if(annotationData.annotations[i].annotation.compare(series[j]) == 0)
        if( (!chargeOne) || (annotationData.annotations[i].charge == 1))
          return i;
  return -1;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
string PlotSpectrum::getMassLabel(string &peptide, int &position, bool reverseLabels)
{
  string aux1 = "()[]{}";
  // get the 1st residue
  string aa, b;
  aa = peptide[position];
  position += reverseLabels ? -1 : 1;
  // if its a delimiter, find the closing match
  if(aa.find_first_of(aux1) != string::npos)
    do {
      b = peptide[position];
      position += reverseLabels ? -1 : 1;
      // if we're moving backwards, we need to reverse the sequence.
      aa = reverseLabels ? b + aa : aa + b;
    } while(b.find_first_of(aux1) == string::npos);
  return aa;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::processUpperLabel(string &label)
{
  // check for empty labels
  if(!label.length())
    return;

  // check for mass item. if it is the case, return an empty label
  if(label[0] == '[') {
    label = "";
    return;
  }

  // check for AAs with a mass shift
  if(label[0] == '(') {
    string aux;
    for(int i = 0 ; i < label.length() ; i++) {
      if(label[i] == ',' || label[i] == ')') {
        aux += ')';
        break;
      }
      aux += label[i];
    }
    label = aux;
  }
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::processUpperLabel2(string &label)
{
  // check for empty labels
  if(!label.length())
    return;

  // store the position after each tag
  size_t lastPosition = 0;
  int index = 0;

  // return string accumulator
  stringstream aux;

  // check for mass item. if it is the case, return an empty label
  if(label[lastPosition] == '[') {
    // initialize walker
    size_t currentPosition = lastPosition + 1;
    // find end
    while(label[currentPosition] != ']' && currentPosition < label.size())
      currentPosition++;
    // get contents
    string aux2 = label.substr(lastPosition+1, currentPosition - lastPosition - 1);
    // translate to a float
    float value = getFloat(aux2.c_str());
    // build return value
    aux << '[';
    // set precision
    aux.unsetf(ios::floatfield);
    // set return value
    aux << std::setprecision(1) << fixed << value;
    // add closing brackets
    aux << ']';
    // set the return value
    label = aux.str();
  }

  // check for AAs with a mass shift
  if(label[lastPosition] == '(') {
    // initialize walker
    size_t currentPosition = lastPosition + 1;
    // find comma or end
    while(label[currentPosition] != ',' && label[currentPosition] != ')' && currentPosition < label.size())
      currentPosition++;
    // get sub sequence contents
    aux << '(' << label.substr(lastPosition+1, currentPosition - lastPosition - 1);

    lastPosition = currentPosition + 1;

    if(label[currentPosition] == ',') {
      // find end
      while(label[currentPosition] != ')' && currentPosition < label.size())
        currentPosition++;
      // get contents
      string aux2 = label.substr(lastPosition, currentPosition - lastPosition );
      // translate to a float
      float value = getFloat(aux2.c_str());
      // add comma
      aux << ',';
      // set precision to 1 decimal place
      aux.unsetf(ios::floatfield);
      if(value < 100.0)
        aux << std::setprecision(1);
      else
        aux << std::setprecision(0);
      // add value
      aux << fixed << value;
    }
    // closing brackets
    aux << ')';
    // set return value
    label = aux.str();
  }
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperLabelsAll(ParamsData1 &paramsData)
{
  // Speedup auxliliary var
  vector<PeptideAnnotation> &annot = *paramsData.peptideAnnotation;
  // size of annotations vector
  int masses = annot.size();

  // cycle on ion masses
  for(int i = 0 ; i < masses ; i++) {

    // get data
    double previousMass = annot[i].previousMass;
    double currentMass  = annot[i].mass;

    // only drawn if in range
    if(testMass(previousMass, currentMass)) {

      // get the peptide item at this mass position
      string aa = annot[i].peptideItem;

      // process annotation for zoom < 0.8
      if(m_rendererImage.zoom < 0.8)
        processUpperLabel(aa);
      else
        processUpperLabel2(aa);

      // The label position.
      float aux1 = max(m_mzLowerLimit, previousMass + paramsData.shift);
      float aux2 = min(m_mzUpperLimit, currentMass  + paramsData.shift);
      float labelPosition = (aux1 + aux2) / 2.0;

      // Draw the label
      drawUpperGuideLabelSegments(aa, labelPosition, paramsData.yPosition, paramsData.color);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperMassDots(ParamsData1 &paramsData)
{
  // Speedup auxliliary var
  vector<PeptideAnnotation> &annot = *paramsData.peptideAnnotation;
  // cycle speed up
  int masses = annot.size();

  // draw initial dot
  drawUpperGuideDots(paramsData.shift, paramsData.yPosition, 2, paramsData.color);

  // cycle on ion masses
  for(int i = 0 ; i < masses ; i++) {

    // get data
    double currentMass  = annot[i].mass;

    // only drawn if in range
    if(testMass(currentMass)) {

      // Find annotation
      int index = findAnnotatedIndex(annot[i], *(paramsData.ions) );

      // Annotation status
      bool annotated = (index >= 0);

      // Dot size dependes on the presence of annotation.
      int dotSize = annotated ? 5 : 2;

      // draw the dot
      drawUpperGuideDots(currentMass + paramsData.shift, paramsData.yPosition, dotSize, paramsData.color);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawUpperSegments(ParamsData1 &paramsData)
{
  // Is previous mass annotated?
  bool previous = true;
  // Speedup auxliliary var
  vector<PeptideAnnotation> &annot = *paramsData.peptideAnnotation;

  // start mass at 0
  int masses = annot.size();
  // Cycle on masses
  for(int i = 0 ; i < masses ; i++) {
    // Find annotation
    int index = findAnnotatedIndex(annot[i], *(paramsData.ions) );
    // Annotation status
    bool annotated = (index >= 0);

    // get data
    double previousMass = annot[i].previousMass;
    double currentMass  = annot[i].mass;
    // only drawn if in range
    if(testMass(previousMass, currentMass))
      // Upper segment if two consecutive peaks are annotated, or if previous and last
      if(previous && (annotated || i == masses-1) ) {
        double x1 = fmax(m_mzLowerLimit, previousMass + paramsData.shift);
        double x2 = fmin(m_mzUpperLimit, currentMass  + paramsData.shift);
        drawUpperGuideSegments(x1 ,x2 , paramsData.yPosition, paramsData.color);
      }

    // Save current item annotation status for next iteration
    previous = annotated;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawVerticalLines(ParamsData1 &paramsData)
{
   // Speedup auxliliary var
  vector<PeptideAnnotation> &annot = *paramsData.peptideAnnotation;
 // size of mass vector
  int masses = annot.size();

  for(int i = 0 ; i < masses ; i++) {
    // get data
    double currentMass  = annot[i].mass;
    // only drawn if in range
    if(testMass(currentMass)) {
      // Find annotation
      //int index = findAnnotatedIndex(paramsData.ions, i+1);
      int index = findAnnotatedIndex(annot[i], *(paramsData.ions), true );
      // Annotation status : -1 -> doesn't exist
      if(index >= 0)
        // Vertical line on peak, if annotated and if charge equals 1
        if(annot[i].annotations[index].charge == 1) {
    //      drawUpperGuideVertical(massesVector[i] + shift, topPosition, color);
          double x1 = paramsData.shift + currentMass;
          double y1 = paramsData.yPosition;
          double x2 = annot[i].annotations[index].peakMass;
          double y2 = annot[i].annotations[index].peakIntensity;
          drawUpperGuideVertical2(x1, y1, x2, y2, paramsData.color);
        }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::buildAnnotationVector(int peptideIndex, string &peptide, float peptideMass, vector<float> &masses, vector<PeptideAnnotation> &peptideAnnotation, bool reverse)
{
  // size of mass vector
  int massesCount = masses.size();
  // start mass at 0
  double previousMass = 0.0;
  // position in peptide string - used for parsing
  //int position = 0;
  int position = reverse ? peptide.size()-1 : 0;

  // cycle thru all peptide masses
  for(int i = 0 ; i < massesCount ; i++) {
    // build annotation item
    PeptideAnnotation annotation;
    annotation.mass             = masses[i];
    annotation.previousMass     = previousMass;
    annotation.peptideItem      = getMassLabel(peptide, position, reverse);

    // cycle thru all the peaks
    for(int j = 0 ; j < m_spectrumData.size() ; j++)
      // cycle thru all the annotations for each peak
      for(int k = 0 ; k < m_spectrumData[j].annotations.size() ; k++)
        // if the series index is the one we are looking for...
        if( (m_spectrumData[j].annotations[k].seriesIndex == i+1) &&
            (m_spectrumData[j].annotations[k].peptideIndex == peptideIndex) ) {
          // Store the data
          AnnotationData  annotationData;
          annotationData.peakMass       = m_spectrumData[j].peakMass;
          annotationData.peakIntensity  = m_spectrumData[j].peakIntensity;
          annotationData.charge         = m_spectrumData[j].annotations[k].charge;
          annotationData.seriesIndex    = m_spectrumData[j].annotations[k].seriesIndex;
          annotationData.annotation     = m_spectrumData[j].annotations[k].annotation;
          annotationData.prob           = m_spectrumData[j].annotations[k].prob;
          // Store annotation for this mass item
          annotation.annotations.push_back(annotationData);
        }
    // Store mass item info in list
    peptideAnnotation.push_back(annotation);
   // Save current mass for next iteration
    previousMass = masses[i];
  }

  // final item (reversed only)
  if(reverse) {
    // build annotation item
    PeptideAnnotation annotation;
    annotation.mass             = peptideMass;
    annotation.previousMass     = previousMass;
    annotation.peptideItem      = getMassLabel(peptide, position, reverse);
    peptideAnnotation.push_back(annotation);
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::aquirePeakLabels(vector<DrawingData> &ions)
{
   // cycle tru all the peaks
  for(int i = 0 ; i < m_spectrumData.size() ; i++) {

    // get data
    double currentMass  = m_spectrumData[i].peakMass;

    // only considered if in range
    if(testMass2(currentMass)) {

      // Check for annotations, given the annotation list
      for(int j = 0 ; j < m_spectrumData[i].annotations.size() ; j++) {

        // Check for the annotation, given the annotation list
        int ionIndex = peakAnnotated(ions,  m_spectrumData[i].annotations[j].annotation);

        if(ionIndex > -1) {
          // peptide index
          int peptideIdx = m_spectrumData[i].annotations[j].peptideIndex;
          // color schema used
          bool ccc          = ions[ionIndex].useBColors;

          int pi = m_spectrumData[i].annotations[j].peptideIndex;

          string series = (m_spectrumData[i].annotations[j].seriesIndex < m_prmMasses[pi].size() ?
          parseInt(m_spectrumData[i].annotations[j].seriesIndex) : "");

          //cout <<   m_spectrumData[i].annotations[j].seriesIndex << " < " << m_prmMasses[pi].size() << endl;

          // Build peak annotation item
          PeakLabelItem  peak;
          peak.mass         = currentMass;
          peak.intesity     = m_spectrumData[i].peakIntensity;
          peak.prob         = m_spectrumData[i].annotations[j].prob;
          peak.annotation   = ions[ionIndex].show;
          peak.subscript    = (m_spectrumData[i].annotations[j].charge == 1 ? "" : parseInt(m_spectrumData[i].annotations[j].charge));
          peak.superscript  = series;
          peak.color        = (ions[ionIndex].colorIdx == -1 ? "black" : getColor(peptideIdx, ions[ionIndex].colorIdx, ccc));
          peak.lineColor    = (ions[ionIndex].lineColorIdx == -1 ? "black" : getColor(peptideIdx, ions[ionIndex].lineColorIdx, ccc));
          peak.placed       = false;

          // Store the peak
          peakLabelList.push_back(peak);
        }
      }
    }
  }
  // sort peak list - rank from higher to lower intensity
  sort(peakLabelList.begin(), peakLabelList.end());
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
bool PlotSpectrum::checkPlacementPosition(ParamsData2 &params)
{
  register int xx = params.currentX;
  register int yy = params.currentY;
  register int lx = params.labelSizeX;
  register int ly = params.labelSizeY;

  // Check map boundaries. If outside, can't do it
  if( (xx < 0) || (xx >= params.matX - lx) ||
      (yy < 0) || (yy >= params.matY - ly) )
    return false;

  register const vector< vector<int> > & mm = params.matrix;

  // check label placement position against matrix
  for(register int j = 0 ; j < lx ; j++)
    for(register int k = 0 ; k < ly ; k++)
      // if occupied by something else, can't use this position
      if(mm[yy + j][xx + k] != 0)
        return false;
  return true;
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::checkPlacementPositions(ParamsData2 &params)
{
  // constants used
  // placement contains the relative coordinates for label possible positions to be tested.
  float placementX[]    = { 0,0,0,0,0.5,-0.5, 1,-1,0.5,-0.5,1,-1,1.5,-1.5,2,-2,0.5,-0.5,1,-1,1.5,-1.5,2,-2,1,-1,2,-2 };
  float placementY[]    = { 0,1,2,3,  3,   3, 3, 3,   2,  2,2, 2,   2,  2,2, 2,  1,   1,1, 1,  1,   1,1, 1,0, 0,0, 0 };
  int placementElems  = 28;
  int labelSize       = 26;

  // Constants used by the function to speed up: transformation to graph coordinates.
  float factorX = (m_mzUpperLimit - m_mzLowerLimit) / params.matX;
  float factorY = m_intensityLimit / params.matY;

  // Adjust label size according to graph size and zoom
  params.labelSizeX = (int)(labelSize / m_rendererImage.zoom);
  params.labelSizeY = (int)(labelSize / m_rendererImage.zoom);

  // cycle thru possible placement positions
  for(int position = 0 ; position < placementElems ; position++) {

    // calculate map position based on displacement
    int currentX = (int)(params.x1 + placementX[position] * params.labelSizeX);
    int currentY = (int)(params.y1 + placementY[position] * params.labelSizeY);

    // Center position on map, used for map marking and checking purposes
    params.currentX = currentX - params.labelSizeX / 2;
    params.currentY = currentY;

    // check label placement position against matrix. Returns 'true' if it can placed at given position.
    if(checkPlacementPosition(params)) {

      // speed up
      register int i   = params.i;
      register int cx  = params.currentX;
      register int cy  = params.currentY;
      register int lsx = params.labelSizeX;
      register int lsy = params.labelSizeY;
      register int mx  = params.matX;
      register int my  = params.matY;
      register vector< vector<int> > & mm = params.matrix;

      // calculate positions in graph
      peakLabelList[i].x = currentX * factorX + m_mzLowerLimit;
      peakLabelList[i].y = currentY * factorY;

      // mark in matrix
      for(register int j = 0 ; (j <  lsy) && (j < my - cy) ; j++)
        for(register int k = 0 ; (k <  lsx) && (k < mx - cx) ; k++)
          mm[cy + j][cx + k] = 2;

      // mark label as placed.
      peakLabelList[i].placed = true;

      // If position is right above the peak, no need to draw the conneting line.
      peakLabelList[i].drawConnectionLine = position;
      return;
    }
  }
  // if no position found, mark as 'not placed' -> default by aquirePeakLabels()
  //peakLabelList[params.i].placed = false;
}
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::placePeakLabels(void)
{
  // define matrix dimensions
  int matX = graphSizePixelsX;
  int matY = graphSizePixelsY;

  // Constants used by the function to speed up
  float factorX = matX / (m_mzUpperLimit - m_mzLowerLimit);
  float factorY = matY / m_intensityLimit;

  // define screen matrix, initialized with 0s
  vector< vector<int> > matrix(matY, vector<int>(matX,0));

  // mark peaks in pixel map
  for(int i = 0 ; i < m_spectrumData.size() ; i++) {

    // get data
    double currentMass = m_spectrumData[i].peakMass;

    // Ignore out of range peaks
    if(testMass(currentMass)) {

       // calc positions in pixel map. Adjust to mz Lower Limit.
      int x1 = (int)((currentMass - m_mzLowerLimit) * factorX);
      int y1 = (int)(m_spectrumData[i].peakIntensity * factorY);
      for(int j = 0 ; j < y1 ; j++)
        // mark in map as "graph"
        matrix[j][x1] = 1;
    }
  }

  // Data structure used for parameter passing
  ParamsData2 params(matrix);
  params.matX = matX;
  params.matY = matY;

  // cycle thru labels
  for(params.i = 0 ; params.i < peakLabelList.size() ; params.i++) {

    // get label coordinates (pixel), based on stored data. Adjust to mz Lower Limit.
    params.x1 = (int)((peakLabelList[params.i].mass - m_mzLowerLimit) * factorX);
    params.y1 = (int)(peakLabelList[params.i].intesity * factorY);

    // find a placement position for label
    checkPlacementPositions(params);
  }

  // debug :)
//  for(int y = 0 ; y < matY ; y++) {
//    for(int x = 0 ; x < matX ; x++)
//      cout << matrix[y][x];
//    cout << endl;
//  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::drawPeakLabels(void)
{
  // Build annotations, starting from higher intesity to lower (previously ordered)
  for(int i = 0 ; i < peakLabelList.size() ; i++) {

  // Only draw placed labels -- the ones a valid position was found
    if(peakLabelList[i].placed) {

      // Draw the connecting line, if needed
      if(peakLabelList[i].drawConnectionLine) {
        drawLine(  1,  0,
                   peakLabelList[i].mass,      RENDERER_COORD_NONE,
                   peakLabelList[i].intesity,  RENDERER_COORD_NONE,
                   peakLabelList[i].x,         RENDERER_COORD_NONE,
                   peakLabelList[i].y,         RENDERER_COORD_NONE,
                   peakLabelList[i].lineColor);
      }

      // Set the label text
      string labelToPrint = peakLabelList[i].annotation;
      // Add superscript and subscript prefix
      if(peakLabelList[i].superscript.length() && peakLabelList[i].subscript.length())
        labelToPrint += "@";
      // Add superscript.
      if(peakLabelList[i].superscript.length()) {
        labelToPrint += "^{";
        labelToPrint += peakLabelList[i].superscript;
        labelToPrint += "}";
      }
      // Add subscript
      if(peakLabelList[i].subscript.length()) {
        labelToPrint += "_{";
        labelToPrint += peakLabelList[i].subscript;
        labelToPrint += "}";
      }
      // Draw the label
      drawLabel(labelToPrint,
                peakLabelList[i].x, RENDERER_COORD_NONE,
                peakLabelList[i].y, RENDERER_COORD_NONE,
                0.0, 1.0,
                RENDERER_OFFSET_CENTER,
                peakLabelList[i].color);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::buildAnnotationOutput(void)
{
  int ident[] = {10, 10, 15, 8, 9, 10, 10};

  // Headers
  oAnnotFileContent << "#################   Header   ############################################# " << endl;
  // Peptide indexes and sequences   Peptide <number>: <string>
  for(int i = 0 ; i < m_peptide.size() ; i++)
    oAnnotFileContent << "Peptide " << i << ':' << m_peptide[i] << endl;
  // Title
  oAnnotFileContent << "Title " << m_title << endl;
  // Zoom
  oAnnotFileContent << "Zoom " << getZoom() << endl;
  // Range
  oAnnotFileContent << "Range " << m_range_min <<':' << m_range_max << endl;
  // Format
  oAnnotFileContent << "Format " << getFileFormat() << endl;
  // Separator
  oAnnotFileContent << "########################################################################## " << endl;
  oAnnotFileContent << "# mass         height           type    order   charge   peptide   prob    " << endl;
  oAnnotFileContent << "#------------------------------------------------------------------------- " << endl;

  // build ordering array
  vector<orderData> od;
	for(int i = 0 ; i < m_spectrumData.size() ; i++) {
	  orderData o;
	  o.index = i;
	  o.data  = m_spectrumData[i].peakIntensity;
	  od.push_back(o);
  }

  if(m_annotationByIntensity)
    sort(od.begin(), od.end());

  // Peak list
	for(unsigned int i = 0 ; i < m_spectrumData.size() ; i++) {
	  // get current index
	  int k = od[i].index;
    // Annotation auxiliary var
    stringstream oAnnoAux;
    // clear auxliary streingstream object
    oAnnoAux.clear();
    // Add peak mass and intensity data
    oAnnoAux.fill(' ');
    oAnnoAux.width(ident[0]);
    oAnnoAux << left << m_spectrumData[k].peakMass;
    oAnnoAux.fill(' ');
    oAnnoAux.width(ident[1]);
    oAnnoAux << right << m_spectrumData[k].peakIntensity;
    // Add the annotations
    if(m_spectrumData[k].annotations.size()) {
    	for(unsigned int j = 0 ; j < m_spectrumData[k].annotations.size() ; j++) {
    	  // Add the peak data
    	  oAnnotFileContent << oAnnoAux.str();
    	  // Add the specific annotation data
        oAnnotFileContent.fill(' ');
        oAnnotFileContent.width(ident[2]);
        oAnnotFileContent << right << m_spectrumData[k].annotations[j].annotation;
        oAnnotFileContent.fill(' ');
        oAnnotFileContent.width(ident[3]);
        oAnnotFileContent << right << parseInt(m_spectrumData[k].annotations[j].seriesIndex);
        oAnnotFileContent.fill(' ');
        oAnnotFileContent.width(ident[4]);
        oAnnotFileContent << right << parseInt(m_spectrumData[k].annotations[j].charge);
        oAnnotFileContent.fill(' ');
        oAnnotFileContent.width(ident[5]);
        oAnnotFileContent << right << parseInt(m_spectrumData[k].annotations[j].peptideIndex);
        oAnnotFileContent.fill(' ');
        oAnnotFileContent.width(ident[6]);
        oAnnotFileContent << right << parseFloat(m_spectrumData[k].annotations[j].prob, 4);
        oAnnotFileContent << endl;
      }
    } else {
  	  // Add the peak data
  	  oAnnotFileContent << oAnnoAux.str();
      oAnnotFileContent << endl;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::writeAnnotationFile(void)
{
  string fullFileName = m_fileOutDir;
  if(fullFileName[fullFileName.size() - 1] != '/' && fullFileName[fullFileName.size() - 1] != '\\')
    fullFileName += '/';
  fullFileName += m_annotationsOutputFile;

  fstream annotFile(fullFileName.c_str(), fstream::out | fstream::trunc | fstream::binary);
  annotFile << oAnnotFileContent.str() << endl;
  annotFile.close();
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::readAnnotationFile(void)
{

}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotSpectrum::dumpPeptideAndMasses(void)
{
  for(int i = 0 ; i < m_peptide.size() ; i++) {

    const vector<PeptideAnnotation> &pa = peptideAnnotations[i];
    const vector<PeptideAnnotation> &rpa = peptideAnnotationsReverse[i];

    cout << " -----------------------------" << endl;
    cout << " - PRM masses" << endl;

    // dump direct
    for(int j = 0 ; j < pa.size() ; j++){
      cout << pa[j].peptideItem << '\t' << '\t' << pa[j].mass << '\t' << pa[j].mass -  pa[j].previousMass << endl;
    }

    cout << endl;
    cout << " -----------------------------" << endl;
    cout << " - SRM masses" << endl;

    // dump reverse
    for(int j = 0 ; j < rpa.size() ; j++){
      cout << rpa[j].peptideItem << '\t' << '\t' << rpa[j].mass << '\t' << rpa[j].mass -  rpa[j].previousMass << endl;
    }

    cout << endl;

  }

}
///////////////////////////////////////////////////////////////////////////////
}; //namespace specnets
///////////////////////////////////////////////////////////////////////////////
