///////////////////////////////////////////////////////////////////////////////
#ifndef __RENDERER_BASE_H__
#define __RENDERER_BASE_H__
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;

typedef enum {RENDERER_TYPE_NONE, RENDERER_TYPE_GNU} RendererType;

typedef enum { RENDERER_DATA_ORDER_NONE, RENDERER_DATA_ORDER_FRONT, RENDERER_DATA_ORDER_BACK } RendererDataOrder;


typedef enum { RENDERER_DATA_TYPE_INT, RENDERER_DATA_TYPE_DOUBLE } RendererDataType;

typedef enum { RENDERER_COORD_NONE, RENDERER_COORD_FIRST, RENDERER_COORD_SECOND, RENDERER_COORD_SCREEN, RENDERER_COORD_GRAPH } RendererCoordType;

typedef enum { RENDERER_OFFSET_RIGHT, RENDERER_OFFSET_CENTER, RENDERER_OFFSET_LEFT } RendererLabelOffsetType;


///////////////////////////////////////////////////////////////////////////////
// Renderer Coordinate
//
// Stores information about a coordinate: value and coordinate system used
//
 /*! \brief RendererCoordinate class

   Stores information about a coordinate: value and coordinate system used
   
   */
struct RendererCoordinate {
  /*! \brief Coordinate value
   */
  double        value;
  /*! \brief Coordinate type
   */
  RendererCoordType type;
};
///////////////////////////////////////////////////////////////////////////////
// Point class
//
// Stores a 2D point coordinates
//
 /*! \brief RendererPoint class

   Stores a 2D point coordinates
   
   */
struct RendererPoint {
  /*! \brief X and Y Coordinates
   */
  RendererCoordinate   x, y;
};
///////////////////////////////////////////////////////////////////////////////
// Line class
//
// Stores line info
//

 /*! \brief Line class
 
   Stores line info
   
   */
struct RendererLine {

 /*! \brief Start point coordinates
   */
  RendererPoint     start;

 /*! \brief End point coordinates
   */
  RendererPoint     end;

  /*! \brief Line type
   */
  int               type;

  /*! \brief Line width
   */
  int               width;

  /*! \brief Line color
   */
  string            color;

  /*! \brief Drawing order
   */
  RendererDataOrder order;

  /*! \brief Default constructor
   */
  RendererLine() : order(RENDERER_DATA_ORDER_NONE) {};

};
///////////////////////////////////////////////////////////////////////////////
// RendererLabel class
//
// Stores label info
//
 /*! \brief RendererLabel class
 
   Stores label info
   
   */
struct RendererLabel {
  /*! \brief Label contents
   */
  string                    label;
  /*! \brief 
   */
  RendererPoint             position;
  /*! \brief Label location
   */
  double                    offsetX;
  /*! \brief Label X offset
   */
  double                    offsetY;
  /*! \brief Label Y offset
   */
  RendererLabelOffsetType   location;
  /*! \brief Font color
   */
  string                    color;
  // drawing order
  /*! \brief drawing order
   */
  RendererDataOrder         order;

  /*! \brief Default constructor
   */
  RendererLabel() : order(RENDERER_DATA_ORDER_NONE) {};

};
///////////////////////////////////////////////////////////////////////////////
// RendererTitle class
//
// Stores title info
//
 /*! \brief RendererTitle class
 
   Stores title info
   
   */
struct RendererTitle {
  /*! \brief Title string
   */
  string    label;
  /*! \brief X offset
   */
  double    offsetX;
  /*! \brief Y offset
   */
  double    offsetY;
  /*! \brief Font color
   */
  string    color;
  /*! \brief Font name
   */
  string    fontName;
  /*! \brief Font size
   */
  double    fontSize;
};

///////////////////////////////////////////////////////////////////////////////
// RendererData class
//
// Plot data container. Contains data to be ploted.
//
 /*! \brief RendererData class
 
   Plot data container. Contains data to be ploted.
   
   */
struct RendererData {
  /*! \brief Data title.
   */
  string                              title;        // Data title.
  /*! \brief Data type: may be int or double
   */
  RendererDataType                    dataType;     // Data type: may be int or double
  /*! \brief # of records
   */
  int                                 dataRecords;  // # of records
  /*! \brief Data records for double data type
   */
  std::vector<pair<double,double> >   dataDouble;   // Data records for double data type
  /*! \brief Data records for int data type
   */
  std::vector<pair<int,int> >         dataInt;      // Data records for int data type
  /*! \brief width of data lines
   */
  int                                 lineWidth;    // width of data lines
  /*! \brief type of data lines
   */
  int                                 lineType;     // type of data lines
  /*! \brief color of data lines
   */
  string                              lineColor;    // color of data lines
};
///////////////////////////////////////////////////////////////////////////////
// Renderer Area - specifies a drawing area
///////////////////////////////////////////////////////////////////////////////
 /*! \brief RendererArea class
 
   Specifies a drawing area
   
   */
struct RendererArea {

  /////////////////////////////////////////////////////////////////////////////
  // grapth specific margins and location

  // borders
  /*! \brief borders
   */
  int    m_border;

  // image margins
  /*! \brief image margins
   */
  double m_marginLeft, m_marginRight, m_marginBottom, m_marginTop;

  // area size
  /*! \brief Use area size
   */
  bool m_useSize;
  /*! \brief area size
   */
  double m_sizeX, m_sizeY;

  // area location
  /*! \brief Use area location
   */
  bool m_useOrigin;
  /*! \brief area location
   */
  double m_originX, m_originY;


  //////////////////////////////////////////////////////////////////////////////
  // axis data

  // axis labels
  /*! \brief axis labels
   */
  string m_xLabel, m_yLabel;

  // x and y start and end values
  /*! \brief x and y start and end values
   */
  bool   m_useXstart, m_useXend, m_useYstart, m_useYend;
  /*! \brief x and y start and end values
   */
  double m_startX, m_endX;
  /*! \brief x and y start and end values
   */
  double m_startY, m_endY;

  // data for x tics
  /*! \brief Use data for x tics
   */
  bool   m_useTicsX;
  /*! \brief Data for y tics
   */
  double m_ticsX;
  /*! \brief Extended x tics data
   */
  vector<pair<string,string> > m_ticsExtendedX;
  // data for y tics
  /*! \brief Use data Y tics
   */
  bool   m_useTicsY;
  /*! \brief Data Y tics
   */
  double m_ticsY;
  /*! \brief Extended y tics data
   */
  vector<pair<string,string> > m_ticsExtendedY;


  //////////////////////////////////////////////////////////////////////////////
  // data on area

  // general labels
  /*! \brief general labels
   */
  vector<RendererLabel> m_labels;

  // arbitrary lines
  /*! \brief arbitrary lines
   */
  vector<RendererLine>  m_lines;

  // arbitrary arrow
  /*! \brief arbitrary arrows
   */
  vector<RendererLine> m_arrows;

  // data to renderer
  /*! \brief datasets to render
   */
  vector<vector<RendererData> > m_data;


  //////////////////////////////////////////////////////////////////////////////
  // Methods

  // constructor
  /*! \brief constructor
   */
  RendererArea() {reset();};

  /*! \brief Reset object
   */
  void reset(void)
  {
    m_useXstart   = m_useXend     = false;
    m_useYstart   = m_useYend     = false;
    m_border      = 0;
    m_marginLeft  = m_marginRight = m_marginBottom = m_marginTop = 0.0;
    m_useSize     = false;
    m_useOrigin   = false;
    m_xLabel      = ""; m_yLabel  = "";
    m_useTicsX    = m_useTicsY    = false;
    m_ticsExtendedX.clear(); m_ticsExtendedY.clear();
    m_labels.clear();
    m_lines.clear();
    m_arrows.clear();
    m_data.clear();
  };

};

///////////////////////////////////////////////////////////////////////////////
// RendererImage class
//
// Image data container
//
 /*! \brief RendererImage class
 
   Image data container
   
   */
struct RendererImage {

  // Image dimensions
  /*! \brief Image dimensions X
   */
  int     imageSizeX;
  /*! \brief Image dimensions Y
   */
  int     imageSizeY;

  //Image stretch limits. Used to limit automatic image stretching
  /*! \brief Image stretch limits. Used to limit automatic image stretching
   */
  int     imageStretchLimitX;
  /*! \brief Image stretch limits. Used to limit automatic image stretching
   */
  int     imageStretchLimitY;

  // Zoom value for image
  /*! \brief Zoom value for image
   */
  double  zoom;

  // Pixel size. Typically 1.0
  /*! \brief Pixel size. Typically 1.0
   */
  double  pixelWidth;
  /*! \brief Pixel size. Typically 1.0
   */
  double  pixelHeight;

  // line thickness
  /*! \brief line thickness
   */
  int     lineThickness;

  // Font name, size and default color
  /*! \brief Font name
   */
  string  fontName;
  /*! \brief Font height
   */
  int     fontHeight;
  /*! \brief Font width
   */
  int     fontWidth;
  /*! \brief Font size
   */
  double  fontSize;
  /*! \brief Font color
   */
  int     fontColor;

  // File name used to write the .png, .jpg, ...
  /*! \brief File name used to write the image to
   */
  string  outputFileName;

  // File type: png, eps, jpg, ...
  /*! \brief File type: png, eps
   */
  string  fileType;

  /*! \brief Default constructor
   */
  RendererImage() : zoom(1.0), lineThickness(1), fileType("png") {};

};
///////////////////////////////////////////////////////////////////////////////
// PlotBase class
//
// Low level drawing wrapper base class.
//
 /*! \brief RendererBase class

   Provides base funtionality for image rendering classes.
   
   */
class RendererBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // renderer data

  /*! \brief Defines the renderer location. Defaults to '.'
   */
  string    m_rendererLocation;

  /*! \brief Defines the font files location. Defaults to '.'
   */
  string    m_fontLocation;

  /*! \brief defines if LD_LIBRARY_PATH is set to the same directory as to the same directory the executable is located in
   */
  bool      m_setLibraryPath;

  //////////////////////////////////////////////////////////////////////////////
  // extra image properties

  // title data
  /*! \brief Defines if the title is to be rendered
   */
  bool m_useTitle;

  /*! \brief Defines the title string
   */
  RendererTitle m_title;


  //////////////////////////////////////////////////////////////////////////////
  // Areas data

  /*! \brief areas in store to draw
   */
  vector<RendererArea>  m_areas;

  /*! \brief current area being rendered
   */
  RendererArea          m_areaCurrent;

  // stores data for x broken axis
  //bool    m_breakYAxis;
  //double  m_yAxisStop, m_yAxisRestart;


  //////////////////////////////////////////////////////////////////////////////
  // auxiliary methods

  /*! \brief Debug flag
   */
  int m_debug;

  /*! \brief EPS correction flag
   */
  bool m_epsCorrection;


  /*! \brief Method used by derived classes. May need to be overloaded
   */
  virtual void outputVector(vector<string> &) {};


 public:

  /*! \brief Image properties
   */
  RendererImage   m_rendererImage;


  //////////////////////////////////////////////////////////////////////////////
  // constructor and destructor

  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Default constructor
   */
  RendererBase();

    //@}

  // Class Destructor.

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
   */
  ~RendererBase();

  // Initialization method
  /*! \brief Initialization method
   */
  virtual void initialize(void);

  /*! \brief Clear renderer object
   */
  virtual void clear(void)                      {m_areas.clear();m_areaCurrent.reset();};

  /*! \brief Set debug mode
      @param d Debug level
   */
  virtual void setDebug(int d)                  {m_debug = d;};

  // use eps correction
  /*! \brief use eps correction
      Used for correcting positional errors when gnuplot renders an EPS image
    @param c True to use EPS correction
   */
  virtual void setEpsCorrection(bool c)         {m_epsCorrection = c;};

  //
  /*! \brief Execute (render) image
   */
  virtual int execute(void)                                 {};

  // Internal method to calculate auxiliary values used in label position calculation
  /*! \brief Internal method to calculate auxiliary values used in label position calculation
   */
  virtual double calculateInterval(double distance, double factor, double zoom);

  // add a new drawing area
  /*! \brief add a new drawing area
   */
  virtual void addArea(void)                                {m_areas.push_back(m_areaCurrent); m_areaCurrent.reset();};

  // set the renderer executable location
  /*! \brief set the renderer executable location
    @param l Location in the filesystem (excluding filename)
   */
  virtual void setRendererLocation(string &l)               {m_rendererLocation = l;};
  // set the font files location
  /*! \brief set the font files location
    @param f Location in the filesystem (excluding filename)
   */
  virtual void setFontLocation(string &f)                   {m_fontLocation = f;};
  // use the LD_LIBRARY_PATH
  /*! \brief use the LD_LIBRARY_PATH variable
   */
  virtual void setLibraryPath(void)                         {m_setLibraryPath = true;};
  // do not use the LD_LIBRARY_PATH
  /*! \brief do not use the LD_LIBRARY_PATH
   */
  virtual void clrLibraryPath(void)                         {m_setLibraryPath = false;};


  //sets the title
  /*! \brief sets the title
    @param t Title
   */
  virtual void setTitle(const RendererTitle &t)             {m_title = t;m_useTitle = true;};

  // Draw a single line
  /*! \brief Draw a single line
   */
  virtual void drawLine(const RendererLine &line)           {m_areaCurrent.m_lines.push_back(line);};
  // Draw a single arrow
  /*! \brief Draw a single arrow
   */
  virtual void drawArrow(const RendererLine &arrows)        {m_areaCurrent.m_arrows.push_back(arrows);};
  // Draw a curve
  /*! \brief Draw a curve
   */
  virtual void drawCurve(void)                              {};
  // Draw graph axis
  /*! \brief Draw graph axis
   */
  virtual void drawAxis(void)                               {};
  // Draw a label
  /*! \brief Draw a label
   */
  virtual void drawLabel(void)                              {};

  // Draw data
  /*! \brief Draw data
   */
  virtual void drawData(vector<RendererData> &data)         {m_areaCurrent.m_data.push_back(data);};

  // Define margins
  /*! \brief Define margins
   */
  virtual void setMarginLeft(const double &margin)          {m_areaCurrent.m_marginLeft    = margin;};
  /*! \brief Define margins
   */
  virtual void setMarginRight(const double &margin)         {m_areaCurrent.m_marginRight   = margin;};
  /*! \brief Define margins
   */
  virtual void setMarginBottom(const double &margin)        {m_areaCurrent.m_marginBottom  = margin;};
  /*! \brief Define margins
   */
  virtual void setMarginTop(const double &margin)           {m_areaCurrent.m_marginTop     = margin;};
  /*! \brief Define margins
   */
  virtual void setMarginAll(const double &margin)           {m_areaCurrent.m_marginLeft    = m_areaCurrent.m_marginRight =
                                                             m_areaCurrent.m_marginBottom  = m_areaCurrent.m_marginTop   = margin;};
  // define origin and size
  /*! \brief define image size
   */
  virtual void setSize(double ox, double oy)                {m_areaCurrent.m_sizeX   = ox; m_areaCurrent.m_sizeY   = oy; m_areaCurrent.m_useSize   = true;};
  /*! \brief define image origin
   */
  virtual void setOrigin(double sx, double sy)              {m_areaCurrent.m_originX = sx; m_areaCurrent.m_originY = sy; m_areaCurrent.m_useOrigin = true;};

  // Set Border
  /*! \brief Set Border
   */
  virtual void setBorder(const int &border)                   {m_areaCurrent.m_border = border;};

  // define axis range - values covered by axis
  /*! \brief define axis range - values covered by axis
   */
  virtual void setXrange(const double &sx, const double &ex)  {m_areaCurrent.m_startX = sx; m_areaCurrent.m_endX = ex;};
  /*! \brief define axis range - values covered by axis
   */
  virtual void setXrange(const double &sx)                    {m_areaCurrent.m_startX = sx;};
  /*! \brief define axis range - values covered by axis
   */
  virtual void setYrange(const double &sy, const double &ey)  {m_areaCurrent.m_startY = sy; m_areaCurrent.m_endY = ey;};
  /*! \brief define axis range - values covered by axis
   */
  virtual void setYrange(const double &sy)                    {m_areaCurrent.m_startY = sy;};

  // define a broken axis
  /*! \brief define a broken axis
   */
  virtual void breakAxisY(const double until, const double restart);

  // Set axis labels
  /*! \brief Set axis labels
   */
  virtual void setXlabel(const string &l)                   {m_areaCurrent.m_xLabel = l;};
  /*! \brief Set axis labels
   */
  virtual void setXlabel(const char *l)                     {m_areaCurrent.m_xLabel = l;};
  /*! \brief Set axis labels
   */
  virtual void setYlabel(const string &l)                   {m_areaCurrent.m_yLabel = l;};
  /*! \brief Set axis labels
   */
  virtual void setYlabel(const char *l)                     {m_areaCurrent.m_yLabel = l;};

  // Set an arbitrary label
  /*! \brief Set an arbitrary label
   */
  virtual void addLabel(const RendererLabel &l)             {m_areaCurrent.m_labels.push_back(l);};

  // set axis tics
  /*! \brief set axis tics
   */
  virtual void setXinterval(const double &v)                  {m_areaCurrent.m_ticsX = v; m_areaCurrent.m_useTicsX = true;};
  /*! \brief set axis tics
   */
  virtual void setYinterval(const double &v)                  {m_areaCurrent.m_ticsY = v; m_areaCurrent.m_useTicsY = true;};

  /*! \brief set axis tics
   */
  virtual void setXinterval(vector<pair<string,string> > &v)  {m_areaCurrent.m_ticsExtendedX = v; m_areaCurrent.m_useTicsX = true;};
  /*! \brief set axis tics
   */
  virtual void setYinterval(vector<pair<string,string> > &v)  {m_areaCurrent.m_ticsExtendedY = v; m_areaCurrent.m_useTicsY = true;};

  /*! \brief Disable axis tics
   */
  virtual void disableXinterval(void)                         {m_areaCurrent.m_useTicsX = false;};
  /*! \brief Disable axis tics
   */
  virtual void disableYinterval(void)                         {m_areaCurrent.m_useTicsY = false;};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
