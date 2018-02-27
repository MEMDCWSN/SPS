///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_SPECTRUM_H__
#define __PLOT_SPECTRUM_H__
///////////////////////////////////////////////////////////////////////////////

#include "PlotBase.h"

#include "spectrum.h"
#include "aminoacid.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  AUxiliary data structures used
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Internal spectrum representation with several annotations
  /*! \brief Internal spectrum representation with several annotations
   */
struct SpectrumAnnotationData {
  int     peptideIndex;
  int     charge;
  int     seriesIndex;
  float   prob;
  string  annotation;
};
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Internal spectrum item representation
   */
struct SpectrumItem {
  double                          peakMass;
  double                          peakIntensity;
  vector<SpectrumAnnotationData>  annotations;
};

// unsigned is the peak index, considering mass incresing ordering
typedef map<unsigned, SpectrumItem> SpectrumData;

////////////////////////////////////////////////////////////////////////////////
// Label data container for label drawing and placement over the graph
  /*! \brief abel data container for label drawing and placement over the graph
   */
struct PeakLabelItem {
  // peak mass value
  /*! \brief peak mass value
   */
  double mass;
  // peak intensity value
  /*! \brief peak intensity value
   */
  double intesity;
  // Annotation probability
  /*! \brief Annotation probability
   */
  float  prob;
  // label coordinates after processing
  /*! \brief label coordinates after processing
   */
  double x,y;
  // Is label placed?
  /*! \brief Is label placed?
   */
  bool placed;
  // Does the label need a connection line?
  /*! \brief Does the label need a connection line?
   */
  bool drawConnectionLine;
  // label to print
  /*! \brief label to print
   */
  string annotation;
  // superscript label part
  /*! \brief superscript label part
   */
  string superscript;
  // subscript label part
  /*! \brief subscript label part
   */
  string subscript;
  // color
  /*! \brief color
   */
  string color;
  // line color
  /*! \brief line color
   */
  string lineColor;

  // used to order items - higher to lower intensity, in this case
  /*! \brief used to order items - higher to lower intensity, in this case
   */
  bool operator<(const PeakLabelItem &o) const {return intesity < o.intesity;};

};
////////////////////////////////////////////////////////////////////////////////
// Auxiliary data structure used to simplify data drawing and searching
  /*! \brief Auxiliary data structure used to simplify data drawing and searching
   */
struct DrawingData {
  // ion to search for
  /*! \brief ion to search for
   */
  string  ion;
  // specifies if ion is used as a prefix or an exact string to search for
  /*! \brief specifies if ion is used as a prefix or an exact string to search for
   */
  bool    isPrefix;
  // string to draw
  /*! \brief string to draw
   */
  string  show;
  // color to draw
  /*! \brief color to draw
   */
  int     colorIdx;
  // Color to draw conneting line
  /*! \brief Color to draw conneting line
   */
  int     lineColorIdx;
  // use B series color schema if true (use y if false)
  /*! \brief use B series color schema if true (use y if false)
   */
  bool    useBColors;

  /*! \brief Constructor
   */
  DrawingData(const char *i, bool p, const char *s, int c, int lc, bool bc) :
    ion(i), isPrefix(p), show(s), colorIdx(c), lineColorIdx(lc), useBColors(bc)
    {};

};
////////////////////////////////////////////////////////////////////////////////
// Annotation data container per mass item
// Used to speed up annotations searchs
  /*! \brief Annotation data container per mass item
   */
struct AnnotationData {
  double  peakMass;
  double  peakIntensity;
  int     charge;
  int     seriesIndex;
  float   prob;
  string  annotation;
};
// Mass item data container. Used as an alternative data structure to speed up searches
  /*! \brief Mass item data container. Used as an alternative data structure to speed up searches
   */
struct PeptideAnnotation {
  string                  peptideItem;
  double                  mass;
  double                  previousMass;
  vector<AnnotationData>  annotations;
};
////////////////////////////////////////////////////////////////////////////////
// Data structure used to simplify and speed up parameter passing
  /*! \brief Data structure used to simplify and speed up parameter passing
   */
struct ParamsData1 {
  vector<string>            *ions;
  vector<PeptideAnnotation> *peptideAnnotation;
  string                    *peptide;

  double                    shift;
  double                    yPosition;
  double                    peptideMass;
  string                    color;

  ParamsData1(vector<string> *i, vector<PeptideAnnotation> *a, string *p) : ions(i), peptideAnnotation(a), peptide(p) {};
};
////////////////////////////////////////////////////////////////////////////////
// Data structure used to simplify and speed up parameter passing
  /*! \brief Data structure used to simplify and speed up parameter passing
   */
struct ParamsData2 {
  vector< vector<int> >&matrix;
  int matX;
  int matY;
  int x1;
  int y1;
  int i;
  int currentX;
  int currentY;
  int labelSizeX;
  int labelSizeY;

  ParamsData2(vector< vector<int> >&m) : matrix(m) {};
};
////////////////////////////////////////////////////////////////////////////////
// spectrumOutputElem contains information about a spectrum item to be printed
  /*! \brief spectrumOutputElem contains information about a spectrum item to be printed
   */
struct spectrumOutputElem {
  string color;
  float  mass;
  float  intensity;
  bool   annotated;

  bool operator<(const spectrumOutputElem &o) const {return intensity > o.intensity;};
};
////////////////////////////////////////////////////////////////////////////////
// ordering data
  /*! \brief ordering data
   */
struct orderData {
  int index;
  double data;

  bool operator<(const orderData &o) const {return data > o.data;};
};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
  /*! /class PlotSpectrum PlotSpectrum.h
    \brief Execution class for spectrum
   */
class PlotSpectrum : public PlotBase {


 protected:


  //////////////////////////////////////////////////////////////////////////////
  // Data properties

  // Intensity upper limit
  /*! \brief Intensity upper limit
   */
  double m_intensityLimit;

  // M/Z limits
  /*! \brief M/Z limits
   */
  double m_mzUpperLimit;
  /*! \brief  M/Z limits
   */
  double m_mzLowerLimit;


  bool m_trim;

  //////////////////////////////////////////////////////////////////////////////
  // Auxiliary data structures used for processing

  // Header positions (2 lines per record)
  // first is b peptide Y axis position
  // second is y peptide Y axis position
  /*! \brief Header positions (2 lines per record)
   */
  map<int, pair<double,double> >  m_headerPosition;

  //prefix masses, i.e. b ions
  /*! \brief prefix masses, i.e. b ions
   */
  map<int, vector<float> >        m_prmMasses;
  //suffix masses, i.e. y ions
  /*! \brief suffix masses, i.e. y ions
   */
  map<int, vector<float> >        m_srmMasses;
  // peptide global masses
  /*! \brief peptide global masses
   */
  map<int, float>                 m_peptideMass;

  // Holds the spectrum data with all the annotations per peak
  /*! \brief Holds the spectrum data with all the annotations per peak
   */
  SpectrumData                    m_spectrumData;

  // processed annotations for module usage, to speed up searches
  /*! \brief processed annotations for module usage, to speed up searches
   */
  map<int, vector<PeptideAnnotation> >  peptideAnnotations;
  /*! \brief processed annotations for module usage, to speed up searches
   */
  map<int, vector<PeptideAnnotation> >  peptideAnnotationsReverse;

  // vector to hold peak information (for peak label generation)
  /*! \brief vector to hold peak information (for peak label generation)
   */
  vector<PeakLabelItem>           peakLabelList;

  // Vector to hold peak information to be output
  /*! \brief Vector to hold peak information to be output
   */
  vector<spectrumOutputElem>      peaksOutput;


  //////////////////////////////////////////////////////////////////////////////
  // Annotation model

  /*! \brief Annotation model
   */
  MS2ScoringModel m_model;


  //////////////////////////////////////////////////////////////////////////////
  // Annotation definition containers

  // Mass dots
  /*! \brief Mass dots
   */
  vector<string> ionsSearchDotsB, ionsSearchDotsY;
  // Mass segments
  /*! \brief Mass segments
   */
  vector<string> ionsSearchSegmentsB, ionsSearchSegmentsY;
  // Vertical lines
  /*! \brief Vertical lines
   */
  vector<string> ionsSearchVerticalB, ionsSearchVerticalY;



  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data

  // Data spectrum for output
  /*! \brief Data spectrum for output
   */
  Spectrum *m_spectrum;

  // Spectrum to draw
  /*! \brief Spectrum to draw
   */
  int     m_spectrumIndex;

  // Spectrum scan number
  /*! \brief Spectrum scan number
   */
  string  m_spectrumScan;

  //
  /*! \brief
   */
  string  m_spectrumInfo;


  //////////////////////////////////////////////////////////////////////////////
  // Peptides data

  // Peptide sequences
  /*! \brief Peptide sequences
   */
  vector<string>  m_peptide;

  //////////////////////////////////////////////////////////////////////////////
  // input file location and name

  // Annotation output data
  /*! \brief Annotation output data
   */
  stringstream oAnnotFileContent;

  // Annotations output file name
  /*! \brief Annotations output file name
   */
  string m_annotationsOutputFile;
  // Annotations input file name
  /*! \brief Annotations input file name
   */
  string m_annotationsInputFile;

  // annotation order
  /*! \brief annotation order
   */
  bool m_annotationByIntensity;

  //////////////////////////////////////////////////////////////////////////////
  // Graph & axis variables

  // Specified m/z axis range
  /*! \brief pecified m/z axis range
   */
  double  m_range_min;
  /*! \brief pecified m/z axis range
   */
  double  m_range_max;

  // Graph size, in pixels. Used in label placement
  /*! \brief Graph size, in pixels. Used in label placement
   */
  int graphSizePixelsX, graphSizePixelsY;

  // main ion representations for N-Term
  /*! \brief main ion representations for N-Term
   */
  string m_mainIon, m_mainIonN;

  // mass offsets
  /*! \brief mass offsets
   */
  double m_offset, m_offsetN;


  //! \name MODIFIERS

  //////////////////////////////////////////////////////////////////////////////
  // The following methods contain low level calls to the renderer object

  // Used to draw the "y" and "b" labels
  /*! \brief Draws the "y" and "b" labels

    Draws the leading labels for the y and b series on top of the image

   @sa drawLabel()
   */
  void drawUpperGuideLabel(string label, double coord, string color);

  // Draws a upper line (dashed)
  /*! \brief Draws the upper dashed lines.

   Draws the upper base dashed lines under the peptide sequence at the top of the image

   @sa drawLine()
   */
  void drawUpperGuideLine(double x1, double x2, double y1, string color);

  // Draw an upper amino-acids (or peptide items)
  /*! \brief Draws one of the upper aminoacis (breaks)
   @sa drawLabel()
   */
  void drawUpperGuideLabelSegments(string label, double x, double y, string color);

  // Draw a upper mass dots
  /*! \brief Draws one of the upper mass dots (breaks)

    Draws the upper line breaks under the annotated amino acids

   @sa drawUpperGuideDotsEps
   */
  void drawUpperGuideDots(double x, double y, int size, string color);

  // Draw a upper mass dots for EPS format
  /*! \brief Draws one of the upper mass dots (breaks)

    Draws the upper

   @sa drawLine()
   */
  void drawUpperGuideDotsEps(double x, double y, int size, string color);

  // Draw an upper mass line segment
  /*! \brief Draw on of the upper mass line segment (breaks) for EPS file format
   @sa drawLine()
   */
  void drawUpperGuideSegments(double x1, double x2, double y1, string color);

  // Draws the vertical line from the dot to the bottom of the graph (depracated)
  /*! \brief Draw one of the vertical lines (from break to peak)
   @sa drawLine()
   */
  void drawUpperGuideVertical(double x1, double y1, string color);

  // Draws the vertical line from the dot to the top of the peak
  /*! \brief Draw one of the vertical lines (from break to peak)
   @sa drawLine()
   */
  void drawUpperGuideVertical2(double x1, double y1, double x2, double y2, string color);


  //////////////////////////////////////////////////////////////////////////////
  // Methos used to build auxiliary (annotation) data structures

  // Main annotations draw procedure
  /*! \brief Draws an annotation peptide
   @sa buildAnnotationVector()
   @sa drawAnnotations()
   */
  void drawAnnotationPeptides(void);

  // initialize spectrumData data structure from Spectrum object
  /*! \brief Initializes spectrum data
   */
  void initSpectrumData(void);

  // Add annotation data to spectrumData structure
  /*! \brief Populate the spectrumData structure with annotation data
    @param spectrum index in specset
   */
  void addAnnotationData(int index);

  // Annotate spectra using annotation module
  /*! \brief Annotates a spectrum
    @param Peptide
    @param psmPtr object
   */
  void annotate(string &peptide, psmPtr);

  // Builds annotation vector used for drawing
  /*! \brief Build annotation vector used for drawing
    @param int Peptide index in specset
    @param string & Peptide
    @param float peptide mass
    @param vector<float> vector of masses (breaks)
    @param vector<PeptideAnnotation> vector of peptide annotations
    @param bool reverse flag. If true, annotation is drawn in reverse
    @param psmPtr object
    @sa psmPtr->setChargeByAnnotation();
    @sa psmPtr->annotate();
   */
  void buildAnnotationVector(int peptideIndex, string &peptide, float peptideMass, vector<float> &masses, vector<PeptideAnnotation> &peptideAnnotation, bool reverse);

  // Sets annotations vectors
  /*! \brief Populate the annotation vectors
   */
  void defineAnnotations(void);

  // Sets annotations vector for peak labels
  /*! \brief Defines the peaks annotation (read from annotation model)
    @param vector<DrawingData> &ionLabels
   */
  void definePeakAnnotations(vector<DrawingData> &ionLabels);

  // Parses the peptide string, and returns the next item
  /*! \brief Gets a single mass element from a peptide.
    @param string &peptide - the peptide string
    @param int &position - the current position (will be changed)
    @param bool reverse - direction
    @return string the retrieved mass label
   */
  string getMassLabel(string &peptide, int &position, bool reverse);

  // Checks if a peak is annotated
  /*! \brief Checks if a peak is annotated.
    @param vector<DrawingData> &series -
    @param string &annotation -
    @return The annotation index in the series if found, -1 otherwise.
   */
  int  peakAnnotated(vector<DrawingData> &series, string &annotation);

  // Checks if a peak is annotated (2)
  /*! \brief Checks if a peak is annotated.
    @param vector<DrawingData> &series -
    @param int peakIndex -
    @return The annotation index in the series if found, -1 otherwise.
   */
  int  peakAnnotated(vector<DrawingData> &series, int peakIndex);

  // Checks if a peak is annotated (3)
  /*! \brief Checks if a peak is annotated.
    @param PeptideAnnotation &annotation -
    @param vector<string> &series -
    @param bool chargeOne = false -
    @return The annotation index in the series if found, -1 otherwise.
   */
  int  findAnnotatedIndex(PeptideAnnotation &annotation, vector<string> &series, bool chargeOne = false);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to draw graph elements

  // Draws upper peptide
  /*! \brief Draws upper peptide.
    @sa drawUpperGuideLabelSegments()
    @sa processUpperLabel()
    @sa processUpperLabel2()
   */
  void drawUpperLabelsAll(ParamsData1 &);

  // Draws mass dots
  /*! \brief Draws mass dots.
    @sa drawUpperGuideDots()
    @sa findAnnotatedIndex()
   */
  void drawUpperMassDots(ParamsData1 &);

  // Draws upper mass segments
  /*! \brief Draws upper mass segments.
   @safindAnnotatedIndex()
   */
  void drawUpperSegments(ParamsData1 &);

  // Draws vertical lines on annotated peaks
  /*! \brief Draws vertical lines on annotated peaks.
   @sa findAnnotatedIndex()
   @sa drawUpperGuideVertical2()
   */
  void drawVerticalLines(ParamsData1 &);

  // Draws annotations
  /*! \brief Draws annotations.
   @param int index
   @sa drawUpperGuideLabel()
   @sa drawUpperGuideLine()
   @sa drawUpperLabelsAll()
   @sa drawUpperMassDots()
   @sa drawUpperSegments()
   @sa drawVerticalLines()
   */
  void drawAnnotations(int index);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to process and draw peak labels

  // Draw peak labels entry point
  /*! \brief Top level for drawing peak labels.
   @sa definePeakAnnotations()
   @sa aquirePeakLabels()
   @sa placePeakLabels()
   @sa drawPeakLabels()
   */
  void drawPeakLabelsMain(void);

  // process label for lower zoom levels -- removing of mass values and keeping AAs
  /*! \brief processes label for lower zoom levels by removing of mass values and keeping AAs
   */
  void processUpperLabel(string &label);
  // process label for upper zoom levels -- rouding of mass values to 1 decimal places
  /*! \brief process label for upper zoom levels; rouding of mass values to 1 decimal places
   */
  void processUpperLabel2(string &label);


  // aquire main peak labels
  /*! \brief aquire main peak labels
   */
  void defineMainAnnotations(void);

  // Build peak list to print labels
  /*! \brief Build peak list to print labels
   @sa testMass2()
   @sa parseInt()
   */
  void aquirePeakLabels(vector<DrawingData> &ions);

  // Place peak labels so that they don't overlap
  /*! \brief Place peak labels so that they don't overlap
   @sa testMass()
   @sa checkPlacementPositions()
   */
  void placePeakLabels(void);

  // Draw the peaks labels
  /*! \brief Draw the peaks labels
   @sa drawLine()
   @sa drawLabel()
   */
  void drawPeakLabels(void);

  // Checks if a specific position is already taken (peak label placement)
  /*! \brief Checks if a specific position is already taken (peak label placement)
   */
  bool checkPlacementPosition(ParamsData2 &params);

  // Checks all the possible positions for a peak label
  /*! \brief Checks all the possible positions for a peak label
   @sa checkPlacementPosition()
   */
  void checkPlacementPositions(ParamsData2 &params);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to process and draw peaks

  // Draw peaks entry point
  /*! \brief Entry point for peaks drawing
   @sa preparePeaksToOutput()
   @sa drawSpectrumLinesFast()
   */
  void drawSpectrumPeaks(void);

  // Builds a structure based on peakColors. This structure is used to produce output
  /*! \brief Builds a structure based on peakColors. This structure is used to produce output
   */
  void preparePeaksToOutput(void);

  // Draws the spectrum lines sending "drawLine" commands to the draw engine. Slow.
  /*! \brief
   @sa m_rendererObject->drawLine()
   */
  void drawSpectrumLines(void);

  // Draws the spectum lines using a binary array data method. Fast.
  /*! \brief Draws the spectum lines using a binary array data method.
   @sa drawData()
   */
  void drawSpectrumLinesFast(void);


  //////////////////////////////////////////////////////////////////////////////
  // Auxiliary methods (misc)

  // Gets a color set used for rendering, given a peptide index
  /*! \brief
   */
  void getColors(int index, char **cb[], char **cy[]);

  // Get a single color at a specified protein index for a specified position. Flag means b or y reference
  /*! \brief
   @param int color index
   @param int color position
   @param bool b or y series
   @return pointer to color as string
   */
  char *getColor(int index, int position, bool b = true);

  // Tests if a given mass point is within parameters
  /*! \brief Tests if a given mass point is within m/z parameters
   @param double mass to test
   @return test result
   */
  bool testMass(double &value);

  // Tests if a given mass set is within parameters
  /*! \brief Tests if a given mass set is within m/z parameters
   @param double lower mass value
   @param double upper mass value
   @return test result
   */
  bool testMass(double &left, double &right);

  // Tests if a given mass point is within parameters, when peptide mass is not known.
  /*! \brief Tests if a given mass point is within parameters, when peptide mass is not known.
   @param double mass value to be tested
   @return test result
   */
  bool testMass2(double &value);


  //////////////////////////////////////////////////////////////////////////////
  // Annotation file output methods

  // Add annotation data to output stream
  /*! \brief Add annotation data to output stream
   */
  void buildAnnotationOutput(void);

  // Write the annotation data file
  /*! \brief Write the annotation data file
   */
  void writeAnnotationFile(void);

  // Read the annotation data file
  /*! \brief Read the annotation data file
   */
  void readAnnotationFile(void);


  //////////////////////////////////////////////////////////////////////////////
  // Initialization methods

  // Calculates image m/z upper and lower bounds
  /*! \brief Calculate image m/z upper and lower bounds
   */
  void calcultateLimits(void);

  // Initialize several graph parameters
  /*! \brief Initialize several image parameters
   @sa m_rendererObject->calculateInterval()
   @sa m_rendererObject->setXinterval()
   @sa m_rendererObject->setYinterval()
   @sa m_rendererObject->setXlabel()
   @sa m_rendererObject->setYlabel()
   @sa m_rendererObject->setXrange()
   @sa m_rendererObject->setYrange()
   @sa m_rendererObject->setMarginLeft()
   @sa m_rendererObject->setMarginRight()
   @sa m_rendererObject->setMarginTop()
   @sa m_rendererObject->setMarginBottom()
   */
  void initializeGraph(void);

  // Sets the image title
  /*! \brief Sets the image title
   @sa m_rendererObject->setTitle()
   */
  void setGraphTitle(void);

  // Calculate b and y Y coordinate positions, and image height
  /*! \brief Calculate b and y Y coordinate positions, and image height
   */
  void calcDrawPosition(void);


  //////////////////////////////////////////////////////////////////////////////
  // Main draw routine

  /*! \brief Main image drawing routine
   @sa defineMainAnnotations()
   @sa defineAnnotations()
   @sa initSpectrumData()
   @sa loadJumps()
   @sa inspectToSpecNets()
   @sa getPRMandSRMMasses()
   @sa addAnnotationData()
   @sa buildAnnotationOutput()
   @sa writeAnnotationFile()
   @sa calcultateLimits()
   @sa initializeGraph()
   @sa calcDrawPosition()
   @sa setGraphTitle()
   @sa drawAnnotationPeptides()
   @sa drawPeakLabelsMain()
   @sa drawSpectrumPeaks()
   */
  int plotSpectrum(void);

  //@}


  //////////////////////////////////////////////////////////////////////////////
  // Debug

  //! \name DEBUG
  //@{
  /*! \brief Outputs peptide and masses to the screen
   */
  void dumpPeptideAndMasses(void);
  //@}


 public:


  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{
  /*! \brief The constructor
   */
  PlotSpectrum(void);
  //@}

  //! \name DESTRUCTOR
  //@{
  ~PlotSpectrum(void);
  //@}

  //////////////////////////////////////////////////////////////////////////////
  // Variable initializers

  //! \name ACCESSORS
  //@{

  // clear object
  /*! \brief Clears all data.
  */
  virtual void clear(void)
  {
    PlotBase::clear();

    m_headerPosition.clear();
    m_prmMasses.clear();
    m_srmMasses.clear();
    m_peptideMass.clear();
    peptideAnnotations.clear();
    peptideAnnotationsReverse.clear();
    peakLabelList.clear();
    peaksOutput.clear();
    ionsSearchDotsB.clear();
    ionsSearchDotsY.clear();
    ionsSearchSegmentsB.clear();
    ionsSearchSegmentsY.clear();
    ionsSearchVerticalB.clear();
    ionsSearchVerticalY.clear();
    m_peptide.clear();
  };

  //////////////////////////
  // Spectrum
  /*! \brief Sets a spectrum object by pointer.
   */
  virtual void setSpectrum(Spectrum *s)               {m_spectrum                 = s;    };
  /*! \brief Specified spectrum index
   */
  virtual void setSpectrumIndex(int idx)              {m_spectrumIndex            = idx;  };
  /*! \brief Specifies spectrum scan
   */
  virtual void setSpectrumScan(string &ss)            {m_spectrumScan             = ss;   };
  virtual void setSpectrumInfo(string ss)             {m_spectrumInfo             = ss;   };

  virtual void setTrim(bool t)                        {m_trim                     = t;    };

  // Peptide
  /*! \brief Specifies the peptide
   */
  virtual void addPeptide(string &p)                  {m_peptide.push_back(p);            };

  // File & directory
  /*! \brief Specifies annotation output file name.
   */
  virtual void setAnnotationOutputFile(string a)      {m_annotationsOutputFile    = a;    };
  virtual void setAnnotationInputFile(string a)       {m_annotationsInputFile     = a;    };

  // misc
  /*! \brief Sort annotationns by intensity instead of m/z.
   */
  virtual void setAnnotationByIntensity(void)         {m_annotationByIntensity = true;};

  // Graph
  /*! \brief Specifies maximun m/z range to display
   */
  virtual void setRangeMax(float r)                   {m_range_max                = r;    };
  /*! \brief Specifies minimun m/z range to display
   */
  virtual void setRangeMin(float r)                   {m_range_min                = r;    };
  //@}


  //////////////////////////////////////////////////////////////////////////////
  // Default draw object
  //! \name MODIFIERS
  //@{
  /*! \brief Generates the image.
   */
  virtual int drawExec(void);
  //@}

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
