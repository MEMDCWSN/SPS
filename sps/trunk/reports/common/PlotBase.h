///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_BASE_H__
#define __PLOT_BASE_H__
///////////////////////////////////////////////////////////////////////////////

#include "RendererBase.h"
#include "Defines.h"
#include "ReportDefines.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
#define ANNOTATION_STYLE_SPECNETS   1
#define ANNOTATION_STYLE_INSPECT    2

#define DEFAULT_ANNOTATION_FILE_LOCATION "."




///////////////////////////////////////////////////////////////////////////////
// PlotBase class
//
// Base class for all plot object classes. Defines base methods for drawing
//
 /*! \brief PlotBase class

   Base class for all plot object classes. Defines base methods for drawing

   */
class PlotBase {

 protected:

  // debug state
  /*! \brief debug state
  */
   int m_debug;

  // EPS CORRECTION
  /*! \brief EPS CORRECTION
  */
  bool m_epsCorrection;

  // use minimun intensity
  /*! \brief use minimun intensity
  */
  bool m_useMinimunIntensity;

  /*! \brief Renderer object

    Base class pointer for renderer object. Should be initialized as a parameter on object contruction.
  */
  RendererBase   *m_rendererObject;


  /*! \brief Renderer owner flag

  is this object owner of the renderer object?
  */
  bool            m_ownRenderer;

  // Image properties
  /*! \brief Image properties
  */
  RendererImage   m_rendererImage;

  // Image contents

  // axis sizes
  /*! \brief Image X axis size
  */
  double m_axisSizeX;
  /*! \brief Image Y axis size
  */
  double m_axisSizeY;


  bool m_verbose;

  /*! \brief Title

  Specifies if the title is printed
  */
  int     m_titlePos;


  /*! \brief Title contents

  Specifies the title
  */
  string  m_title;

  // Top margin position
  /*! \brief Top margin position
  */
  double m_topMarginPosition;

  // Title offsets
  /*! \brief Title offset X
  */
  double m_titleOffsetX;
  /*! \brief Title offset Y
  */
  double m_titleOffsetY;

  // Mass shift value
  /*! \brief Mass shift value
  */
  double  m_massShift;

  // Peak mass tolerance
  /*! \brief Peak mass tolerance
  */
  double  m_peakMassTol;

  // use peak mass tolerance.
  /*! \brief use peak mass tolerance flag
  */
  bool    m_PeakMassTolSet;

  // use individual peak mass tolerance
  /*! \brief use individual peak mass tolerance flag
  */
  bool    m_useIndividualPeakMassTol;

  // peptide annotation style
  /*! \brief peptide annotation style

  Specifies the peptide annotation style (inspect of specnets)
  */
  int m_annotationStyle;


  //////////////////////////////////////////////////////////////////////////////
  // input file location and name

  // Amino-acids file
  /*! \brief Amino-acids file name
  */
  string  m_aminoAcidsFile;

  //  Amino-acids file location
  /*! \brief Amino-acids file location
  */
  string  m_aminoAcidsFileDirectory;

  // Annotation model used
  /*! \brief Annotation model file name
  */
  string  m_annotationModel;

  // Annotation model file location
  /*! \brief Annotation model file location
  */
  string  m_annotationModelDirectory;

  // if true, the annotation model is acquired from the file
  /*! \brief Annotation model use flag

     if true, the annotation model is acquired from the file
  */
  bool m_readModelFromData;

  // Renderer directory
  /*! \brief External renderer location
  */
  string m_rendererLocation;

  // Font file directory
  /*! \brief Font file location
  */
  string m_fontLocation;

  /*! \brief Input file name (spectra file)
  */
  string m_inputSpectraFilename;



  // Image file attributes

  // File name
  /*! \brief Output image file name
  */
  string m_fileName;

  // output directory
  /*! \brief Output image file location
  */
  string m_fileOutDir;

  // File prefix for file name formatting
  /*! \brief File prefix
  */
  string m_filePrefix;

  // File suffix for file name formatting
  /*! \brief File name suffix
  */
  string m_fileNameSuffix;

  // File extension (defines file format)
  /*! \brief Output file extension
  */
  string m_fileFormat;

  // File encoding
  /*! \brief Output file encoding
  */
  string m_encoding;

  // output target
  /*! \brief Output target
  */
  string m_target;

  // generated image stored internally
  /*! \brief Image container
  */
  string m_image;



  /*! \brief Defines pixel size considering the image size is 1
  */
  virtual void initializePixelSize(void);


  // Draws a generic label
  /*! \brief Draws a generic label
  */
  virtual void drawLabel(string theLabel,
                          double x1, RendererCoordType coordX1,
                          double y1, RendererCoordType coordY1,
                          double offSetX,
                          double offSetY,
                          RendererLabelOffsetType  location,
                          string color,
                          RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

  // Draws a generic line
  /*! \brief Draws a generic line
  */
  virtual void drawLine(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

  // draw a generic arrow
  /*! \brief Draws a generic arrow
  */
  virtual void drawArrow(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

 public:


  // Default constructor

  //! \name CONSTRUCTORS
  //@{
  /*! \brief Default constructor

  */
  PlotBase();

    //@}

  // Class Destructor.

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
  */
  ~PlotBase();

  // debug method
  /*! \brief Debug mode
  */
  virtual void setDebug(int d)                  {m_debug = d;};

  // use eps correction
  /*! \brief Uses EPS coordinate correction
  */
  virtual void setEpsCorrection(bool c)         {m_epsCorrection = c;};

  // use minimun intensity value
  /*! \brief Defines minimun intensity value
  */
  virtual void useMinimunIntentity(bool u)      {m_useMinimunIntensity = u;};

  // Parameter setting methods
  /*! \brief Sets output filename
  */
  virtual void setFileName(const string &fn)    {m_fileName   = fn;};

  /*! \brief Sets output filename suffix
  */
  virtual void setFnSuffix(const string &fns)   {m_fileNameSuffix = fns;};

  /*! \brief Sets output file directory
  */
  virtual void setFileOutDir(const string &od)  {m_fileOutDir = od;};

  /*! \brief Sets file preffix
  */
  virtual void setFilePrefix(const string &fp)  {m_filePrefix = fp;};

  /*! \brief Sets the file format
  */
  virtual void setFileFormat(const string &ff)  {m_fileFormat = ff;};

  /*! \brief Sets the file encoding
  */
  virtual void setEncoding(const string &e)     {m_encoding   = e; };

  /*! \brief Sets the result image target
  */
  virtual void setTarget(const string &t)       {m_target     = t; };

  /*! \brief Sets the zoom value
  */
  virtual void setZoom(const double &z)         {m_rendererImage.zoom = z;};

  /*! \brief Sets the output image dimensions
  */
  virtual void setImageDimensions(const int x, const int y)     {m_rendererImage.imageSizeX = x;  m_rendererImage.imageSizeY  = y; initializePixelSize();};

  /*! \brief Sets the output image height, in pixels
  */
  virtual void setImageHeight(const int y)                      {m_rendererImage.imageSizeY = y; initializePixelSize();};

  /*! \brief Sets the output image width, in pixels
  */
  virtual void setImageWidth(const int x)                       {m_rendererImage.imageSizeX = x; initializePixelSize();};

  /*! \brief Sets the pixel size
  */
  virtual void setPixelSize(const double &pw, const double &ph) {m_rendererImage.pixelWidth = pw; m_rendererImage.pixelHeight = ph;};

  /*! \brief Sets the maximun image height
  */
  virtual void setImageStretchHeight(const int y)               {m_rendererImage.imageStretchLimitY = y;};

  /*! \brief Sets the maximun image width
  */
  virtual void setImageStretchWidth(const int x)                {m_rendererImage.imageStretchLimitX = x;};

  /*! \brief Sets the line thickness
  */
  virtual void setLineThickness(const int x)                    {m_rendererImage.lineThickness = x;};

  /*! \brief Sets the name of the font to use in text
  */
  virtual void setFontName(const string &fn)    {m_rendererImage.fontName    = fn;};

  /*! \brief Sets the font height
  */
  virtual void setFontHeight(const int fh)      {m_rendererImage.fontHeight  = fh;};

  /*! \brief Sets the font width
  */
  virtual void setFontWidth(const int fw)       {m_rendererImage.fontWidth   = fw;};

  /*! \brief Sets the font size
  */
  virtual void setFontSize(const double &fs)    {m_rendererImage.fontSize    = fs;};

  /*! \brief Sets the font color
  */
  virtual void setFontColor(const int fc)       {m_rendererImage.fontColor   = fc;};

  /*! \brief Specifies if a title will be rendered
  */
  virtual void setTitlePresence(bool tp)                {m_titlePos = tp;}; // True for no title; false to use the title

  /*! \brief Sets the title
  */
  virtual void setTitle(const string &t)                {m_title    = t;};

  /*! \brief Specifies the top margin position
  */
  virtual void setTopMarginPosition(const double &tmp)  {m_topMarginPosition = tmp;};

  // Mass
  /*! \brief Sets the Peak Mass Tolerance
  */
  virtual void setPeakMassTol(float m)                {m_peakMassTol              = m;
                                                       m_PeakMassTolSet           = true; };
  /*! \brief Specifies if individual Peak Mass Tolerance is to be used
  */
  virtual void useIndividualPeakMassTol(bool u)       {m_useIndividualPeakMassTol = u;    };

  /*! \brief Sets the mass shift value
  */
  virtual void setMassShift(float m)                  {m_massShift                = m;    };


  // Peptide
  /*! \brief Sets the peptide annotation style
  */
  virtual void setAnnotatinStyle(int s)               {m_annotationStyle          = s;    };


  // File & directory

   /*! \brief Sets the aminoacids filename
   */
  virtual void setAminoAcidsFile(string f)            {m_aminoAcidsFile           = f;    };

   /*! \brief Sets the aminoacids file location
   */
  virtual void setAminoAcidsFileDirectory(string l)   {m_aminoAcidsFileDirectory  = l;    };

   /*! \brief Sets the annotation model filename
   */
  virtual void setAnnotationModel(string f)           {m_annotationModel          = f;    };

   /*! \brief Sets the annotation model file location
   */
  virtual void setAnnotationModelDirectory(string l)  {m_annotationModelDirectory = l;    };

   /*! \brief Specifies if the model used is load from file
   */
  virtual void setModelFromFile(void)                 {m_readModelFromData        = false;};

   /*! \brief Sets the renderer location in the file system
   */
  virtual void setRendererLocation(string l)          {m_rendererLocation         = l;    };

   /*! \brief Sets the font file location
   */
  virtual void setFontLocation(string l)              {m_fontLocation             = l;    };


   /*! \brief Sets the input spectra file name
   */
  virtual void setInputSpectraFilename(string fn)     {m_inputSpectraFilename     = fn;   };

  // renderer object delegation


  // Parameter access methods

   /*! \brief Gets the file format
   */
  virtual string &getFileFormat(void)           {return m_fileFormat;};

   /*! \brief Gets the zoom value
   */
  virtual double  getZoom(void)                 {return m_rendererImage.zoom;};

   /*! \brief Gets a reference to the internal image
   */
  virtual const string &getImage(void) const {return m_image;};

  // Object initializer
   /*! \brief initializes the renderer object
   */
  virtual void initialize(void);

    /*! \brief setss verbosity flag
   */
  virtual void setVerbose(bool v) {m_verbose = v;};

  // clear the object
    /*! \brief Clears the renderer object
   */
  virtual void clear(void) {if(m_rendererObject) m_rendererObject->clear();};

  // Initial image properties
   /*! \brief Calculates the draw properties based on the current options
   */
  virtual void setDrawProperties(void);

  // Define renderer
   /*! \brief Defines the renderer
   */
  virtual void setDrawMethod(RendererType);

   /*! \brief Default draw method
   */
  virtual void draw(void);

   /*! \brief Draw method, using a file as target.
   */
  virtual void draw(char *filename);

   /*! \brief Draw execution
   */
  virtual int drawExec(void) = 0;

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
