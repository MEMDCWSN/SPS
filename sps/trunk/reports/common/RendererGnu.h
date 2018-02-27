///////////////////////////////////////////////////////////////////////////////
#ifndef __RENDERER_GNU_H__
#define __RENDERER_GNU_H__
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>


#include "RendererBase.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////


using namespace std;

///////////////////////////////////////////////////////////////////////////////
 /*! \brief RendererGnu class

   Provides image rendering funtionality based on GnuPlot.
   
   */
class RendererGnu : public RendererBase {

  // Where to hold GnuRenderer command stream (command stack)
  /*! \brief Where to hold GnuRenderer command stream (command stack)
   */
  stringstream   m_commandStack, m_gnuplotQueue;

  // Stores gnuRenderer command file name
  /*! \brief Stores gnuRenderer command file name
   */
  string          m_gnuplotCommandFileName;

  // Keeps track of terminator
  /*! \brief Keeps track of terminator
   */
  bool  m_terminatorIssued;

  // Output helpers section
  /*! \brief Output helpers section
   */
  virtual void outputVector(stringstream &out, vector<string> &);
  /*! \brief Outputs a vector of strings
   */
  virtual void outputDoubleVector(stringstream &out, vector<pair<string,string> > &);

  // Sets GnuRenderer Command File name
  /*! \brief  Sets GnuRenderer Command File name
   */
  virtual void setGnuplotCommandFileName(void);
  // Generate gnuplot command queue
  /*! \brief Generate gnuplot command queue
   */
  virtual void createGnuplotCommandList(void);
  // Writes command file for GnuRenderer
  /*! \brief Writes command file for GnuRenderer
   */
  virtual void writeGnuplotFile(void);
  // invokes GnuRenderer
  /*! \brief invokes GnuRenderer
   */
  virtual void callGnuplot(void);
  virtual void callGnuplot2(void);
  // Removes gnuRenderer command file
  /*! \brief gnuRenderer command file
   */
  virtual void deleteGnuplotCommandFile(void);

  // Axis range tail
  /*! \brief  Axis range tail
   */
  virtual void setAxisRange(stringstream &out, const double &RangeStart, const double &RangeEnd);
  virtual void setAxisRange(stringstream &out, const double &RangeStart);

  // Method for partial arrow and line drawing methods
  /*! \brief Method for partial arrow and line drawing methods
   */
  virtual void drawArrowInternal(stringstream &out, const RendererLine &line);

  // output a single coordinate
  /*! \brief output a single coordinate
   */
  virtual void drawCoordinate(stringstream &out, const RendererCoordinate &coordinate);

  // output item order string
  /*! \brief output item order string
   */
  virtual void outputOrder(stringstream &out, const RendererDataOrder &order);


  // output a single point
  /*! \brief output a single point
   */
  virtual void drawPoint(stringstream &out, const RendererPoint &RendererPoint);

  // get file format command for gnuplot
  /*! \brief get file format command for gnuplot
   */
  virtual string getFileFormat(void);



  // sets the title
  /*! \brief sets the title
   */
  virtual void setTitle(stringstream &out, const char *title);
  /*! \brief sets the title
   */
  virtual void setTitle(stringstream &out, const RendererTitle &title);

  // Clear command stack
  /*! \brief Clears command stack
   */
  virtual void clear(stringstream &out);

  // Draw a single line
  /*! \brief Draws a single line
   */
  virtual void drawLine(stringstream &out, const RendererLine &line);
  // Draw an arrow
  /*! \brief Draws an arrow
   */
  virtual void drawArrow(stringstream &out, const RendererLine &line);
  // Draw a label
  /*! \brief Draws a label
   */
  virtual void drawLabel(stringstream &out);
  // Draw a curve
  /*! \brief Draws a curve
   */
  virtual void drawCurve(stringstream &out);
  // Draw graph axis
  /*! \brief Draws graph axis
   */
  virtual void drawAxis(stringstream &out);

  // Draw data
  /*! \brief Draw data
   */
  virtual void drawData(stringstream &out, vector<RendererData> &data);

  // Set margins
  /*! \brief  Set margin left
   */
  virtual void setMarginLeft(stringstream &out, const double &margin, bool src = false);
  /*! \brief  Set margin right
   */
  virtual void setMarginRight(stringstream &out, const double &margin, bool src = false);
  /*! \brief  Set margin bottom
   */
  virtual void setMarginBottom(stringstream &out, const double &margin, bool src = false);
  /*! \brief  Set margin top
   */
  virtual void setMarginTop(stringstream &out, const double &margin, bool src = false);
  /*! \brief  Set margins
   */
  virtual void setMarginAll(stringstream &out, const double &margin, bool src = false);

  // Set origin and size
  /*! \brief Set area origin
   */
  virtual void setOrigin(stringstream &out, double &x, double &y);
  /*! \brief Set area size
   */
  virtual void setSize(stringstream &out, double &x, double &y);

  // Set draw area
  /*! \brief Set draw area
   */
  virtual void writeDataSection(stringstream &out, RendererArea &area);

  // Set Border
  /*! \brief Set Border
   */
  virtual void setBorder(stringstream &out, const int &border);

  // Set axis range
  /*! \brief Set axis X range
   */
  virtual void setXrange(stringstream &out, const double &, const double &);
  /*! \brief Set axis X range
   */
  virtual void setXrange(stringstream &out, const double &);
  /*! \brief Set axis Y range
   */
  virtual void setYrange(stringstream &out, const double &, const double &);
  /*! \brief Set axis Y range
   */
  virtual void setYrange(stringstream &out, const double &);

  // Set axis labels
  /*! \brief Set axis labels
   */
  virtual void setXlabel(stringstream &out, const string &);
  /*! \brief Set axis labels
   */
  virtual void setXlabel(stringstream &out, const char *);
  /*! \brief Set axis labels
   */
  virtual void setYlabel(stringstream &out, const string &);
  /*! \brief Set axis labels
   */
  virtual void setYlabel(stringstream &out, const char *);

  /*! \brief Set labels
   */
  virtual void setLabel(stringstream &out, const RendererLabel &label);


  // set axis intervals
  /*! \brief set axis intervals
   */
  virtual void setXinterval(stringstream &out);
  /*! \brief set axis intervals
   */
  virtual void setYinterval(stringstream &out);

  /*! \brief set axis intervals
   */
  virtual void setXinterval(stringstream &out, const double &);
  /*! \brief set axis intervals
   */
  virtual void setYinterval(stringstream &out, const double &);

  /*! \brief set axis intervals
   */
  virtual void setXinterval(stringstream &out, vector<pair<string,string> > &);
  /*! \brief set axis intervals
   */
  virtual void setYinterval(stringstream &out, vector<pair<string,string> > &);

  /*! \brief set axis intervals
   */
  virtual void setXintervalOut(stringstream &out);
  /*! \brief set axis intervals
   */
  virtual void setYintervalOut(stringstream &out);



 public:


  // Constructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Default constructor
   */
  RendererGnu();

    //@}

  // Class Destructor.

    //! \name DESTRUCTOR
    //@{

  /*! \brief Default destructor
   */
  ~RendererGnu();

  // initialize object
  /*! \brief initialize rendering object
   */
  virtual void initialize(void);

  // Execute commands in the stack
  /*! \brief  Execute commands in the stack
   */
  virtual int execute(void);

  // Clear command stack
  /*! \brief Clear command stack
   */
  virtual void clear(void)       {m_commandStack.str("");m_gnuplotQueue.str("");RendererBase::clear();};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
