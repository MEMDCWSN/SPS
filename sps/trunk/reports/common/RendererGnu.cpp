////////////////////////////////////////////////////////////////////////////////
//#if defined(__MINGW32__) || defined(__CYGWIN__)
//#include <windows.h>
//#endif

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>


#include "RendererGnu.h"
#include "cstream.h"
#include "utils.h"
#include "Specific.h"

////////////////////////////////////////////////////////////////////////////////
// Defines used in this file

#define UNIX_DIRECTORY_SEPARATOR    '/'
#define GNU_PLOT_DATA_EXTENSION     ".gp"
#define DIRECTORY_SEPARATOR         '/'

#if  defined(__MINGW32__) || defined(__CYGWIN__)
#define GNU_PLOT_COMMAND            "gnuplot.exe "
#else
#define GNU_PLOT_COMMAND            "gnuplot "
#endif




//#if defined(__MINGW32__) || defined(__CYGWIN__)

//#else
  int (* gnuplot_main)(int, char **);
  int pdll;
  //extern "C" int gnuplot_main(int argc, char **argv);
//#endif




////////////////////////////////////////////////////////////////////////////////
namespace specnets {
////////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////////
//Constructor
RendererGnu::RendererGnu()
{
  initialize();

  gnuplot_main = NULL;
  pdll = (int)NULL;

//#if defined(__MINGW32__) || defined(__CYGWIN__)

  /*pdll = LoadLibrary("gnuplot.dll");

  if (! pdll)  {
    cerr << "error: cannot load 'gnuplot.dll': " << GetLastError() << endl;
    return;
  }

  gnuplot_main = (int (*)(int, char **)) (GetProcAddress(pdll , "gnuplot_main"));

  if (! gnuplot_main) {
    cerr << "error: cannot find 'gnuplot_main'" << endl;
  }
  */
//#endif
}

// Destructor
RendererGnu::~RendererGnu()
{
//#if defined(__MINGW32__) || defined(__CYGWIN__)
  /*if ( pdll != NULL )
    FreeLibrary( pdll );
    */
//#endif
}
////////////////////////////////////////////////////////////////////////////////
// Initialize object
void RendererGnu::initialize(void)
{
  m_terminatorIssued = false;
  //m_commandStack.str(std::string());
  m_commandStack.clear();
  //m_gnuplotQueue.str(std::string());
  m_gnuplotQueue.clear();
}
////////////////////////////////////////////////////////////////////////////////
// File management section
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setGnuplotCommandFileName(void)
{
  // GluRenderer command file name is output file name with ".gp" appended
  m_gnuplotCommandFileName = m_rendererImage.outputFileName;
  m_gnuplotCommandFileName += GNU_PLOT_DATA_EXTENSION;
}
////////////////////////////////////////////////////////////////////////////////
string RendererGnu::getFileFormat(void)
{
  // default return string
  string ret = "png";

  // format EPS
  if(m_rendererImage.fileType.compare("eps") == 0)
    ret = "postscript eps";

  // Return the value
  return ret;
}
////////////////////////////////////////////////////////////////////////////////
// Generates command stack for gnuplot
void RendererGnu::createGnuplotCommandList(void)
{
  bool eps = (m_rendererImage.fileType.compare("eps") == 0);

  // Define nan
  m_gnuplotQueue << "nan=0\n";

  // if format is eps, we need to specify the prolog file location
  if(eps) {
    m_gnuplotQueue << "set loadpath ";
    m_gnuplotQueue << "\"" << m_rendererLocation << "\"\n";
  }

  // Define output type
  m_gnuplotQueue << "set term ";
  m_gnuplotQueue << getFileFormat();
  m_gnuplotQueue << " enhanced ";
  if(!eps)
    m_gnuplotQueue <<"crop ";

  // Font
  if(!eps) {
    m_gnuplotQueue << "font \"";
    m_gnuplotQueue << m_fontLocation;
    m_gnuplotQueue << UNIX_DIRECTORY_SEPARATOR;
    m_gnuplotQueue << m_rendererImage.fontName;
    m_gnuplotQueue << "\" ";
    m_gnuplotQueue << m_rendererImage.fontSize * (eps ? 2.0 : 1.0);
  }

  // Image size
  if(eps) {
     m_gnuplotQueue << " color";
  }// else {
    m_gnuplotQueue << " size ";
    m_gnuplotQueue << m_rendererImage.imageSizeX * m_rendererImage.zoom / (eps ? 100.0 : 1.0);
    m_gnuplotQueue << ",";
    m_gnuplotQueue << m_rendererImage.imageSizeY * m_rendererImage.zoom / (eps ? 100.0 : 1.0);
  //}
  m_gnuplotQueue << "\n";

  // File name
  m_gnuplotQueue << "set output \"";
  m_gnuplotQueue << m_rendererImage.outputFileName.c_str();
  m_gnuplotQueue << "\"\n";


  if(m_areas.size()) {
    // Multiplot allows for more than one plot array
    m_gnuplotQueue << "set multiplot layout " << m_areas.size() << ",1 ";
    if(m_title.label.length())
       m_gnuplotQueue << "title \"" << m_title.label << "\"";
     m_gnuplotQueue << "\n";
  } else {
    // draw image title
    if(m_useTitle)
      setTitle(m_gnuplotQueue, m_title);
  }

   m_gnuplotQueue << "set size 1.0, 1.0\n";
   m_gnuplotQueue << "set origin 0.0, 0.0\n";

  // autoscale fix
  m_gnuplotQueue << "set autoscale fix\n";

  // write data stack
  for(int i = 0 ; i < m_areas.size() ; i++) {
    writeDataSection(m_commandStack, m_areas[i]);
    setXintervalOut(m_commandStack);
    setYintervalOut(m_commandStack);
    m_commandStack << "unset key\n";
    m_commandStack << "unset title\n";
    m_commandStack << "unset label\n";
    m_commandStack << "unset arrow\n";
  }

  // write current data item
  writeDataSection(m_commandStack, m_areaCurrent);

  // Add the command list to the gnuPlot command list
  m_gnuplotQueue << m_commandStack.str();

  // Output terminator, if needed
  //if(!m_terminatorIssued)
  //  m_gnuplotQueue << "plot \"-\" binary record=0 format=\"%double%double\" with impulses lt 1 lw 2 lc rgb \"grey\" notitle\n";

//  m_gnuplotQueue << "replot\n";
  m_gnuplotQueue << "unset multiplot\n";

  m_gnuplotQueue << "unset arrow\n";
  m_gnuplotQueue << "unset label\n";
  //gnuplotCmdFile << "unset output";
}
////////////////////////////////////////////////////////////////////////////////
// Generates command stack about a data section
void RendererGnu::writeDataSection(stringstream &out, RendererArea &area)
{
  // global margins
  setMarginLeft(out, area.m_marginLeft);
  setMarginRight(out, area.m_marginRight);
  setMarginBottom(out, area.m_marginBottom);
  setMarginTop(out, area.m_marginTop);
  setBorder(out, area.m_border);

  //cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
  //cout << "m_marginLeft   " << area.m_marginLeft << endl;
  //cout << "m_marginRight  " << area.m_marginRight << endl;
  //cout << "m_marginBottom " << area.m_marginBottom << endl;
  //cout << "m_marginTop    " << area.m_marginTop << endl;
  //cout << "m_border       " << area.m_border<< endl;

  //cout << "m_originX      " << area.m_originX << endl;
  //cout << "m_originY      " << area.m_originY << endl;
  //cout << "m_sizeX        " << area.m_sizeX << endl;
  //cout << "m_sizeX        " << area.m_sizeX << endl;

  double shift = 0.0, increase = 0.0;
  if(m_rendererImage.fileType.compare("eps") == 0 && m_epsCorrection) {
    shift = 0.005 * (1 - area.m_originX * 1.7);
    increase = area.m_sizeX * 0.005;
  }

  double originX = area.m_originX - shift;
  double sizeX = area.m_sizeX + increase;

  // graph origin in image
  if(area.m_useOrigin)
    setOrigin(out, originX, area.m_originY);
    //setOrigin(out, area.m_originX, area.m_originY);

  // graph area in image
  if(area.m_useSize)
    setSize(out, sizeX, area.m_sizeY);
    //setSize(out, area.m_sizeX, area.m_sizeY);

  // axis data

  // x tics
  if(area.m_useTicsX)
    if(area.m_ticsExtendedX.size())
      setXinterval(out, area.m_ticsExtendedX);
    else
      setXinterval(out, area.m_ticsX);
  else
    setXintervalOut(out);

  // y tics
  if(area.m_useTicsY)
    if(area.m_ticsExtendedY.size())
      setYinterval(out, area.m_ticsExtendedY);
    else
      setYinterval(out, area.m_ticsY);
  else
    setYintervalOut(out);

  // labels
  if(area.m_xLabel.length())
    setXlabel(out, area.m_xLabel);

  if(area.m_yLabel.length())
    setYlabel(out, area.m_yLabel);

  // X axis range
  if(area.m_startX >= 0.0 && area.m_endX >= 0.0)
    setXrange(out, area.m_startX, area.m_endX);
  else if(area.m_startX >= 0.0)
    setXrange(out, area.m_startX);

  // Y axis range
  if(area.m_startY >= 0.0 && area.m_endY >= 0.0)
    setYrange(out, area.m_startY, area.m_endY);
  else if(area.m_startY >= 0.0)
    setYrange(out, area.m_startY);


  // write lines
  for(int i = 0 ; i < area.m_lines.size() ; i++)
    drawLine(out, area.m_lines[i]);

  // write arrows
  for(int i = 0 ; i < area.m_arrows.size() ; i++)
    drawArrow(out, area.m_arrows[i]);

  // write labels
  for(int i = 0 ; i < area.m_labels.size() ; i++)
    setLabel(out, area.m_labels[i]);

  // write data
  for(int i = 0 ; i < area.m_data.size() ; i++)
    drawData(out, area.m_data[i]);

  if(!m_terminatorIssued) {
    out << "plot [";
    out << m_areaCurrent.m_startX;
    out << ":";
    out << m_areaCurrent.m_endX;
    out << "] [";
    out << m_areaCurrent.m_startY;
    out << ":";
    out << m_areaCurrent.m_endY;
    out << "] NaN notitle\n";
  }
    //out << "plot \"-\" binary record=0 format=\"%double%double\" with impulses lt 1 lw 2 lc rgb \"grey\" notitle\n";
}
////////////////////////////////////////////////////////////////////////////////
// Writes command file for GnuRenderer
void RendererGnu::writeGnuplotFile(void)
{
  // Declare command file object
  ofstream gnuplotCmdFile;

  // Open the command file to write to
  gnuplotCmdFile.open(m_gnuplotCommandFileName.c_str(), fstream::in | fstream::out | fstream::trunc | fstream::binary);

  // check if file open succeeded
  if(!gnuplotCmdFile.is_open()) {
    cerr << "Failed to create gnuplot command file" << endl;
    return;
  }

  // send command list to file
  gnuplotCmdFile << m_gnuplotQueue.str();

  // Close the file
  gnuplotCmdFile.close();
}
////////////////////////////////////////////////////////////////////////////////
// invokes GnuRenderer
void RendererGnu::callGnuplot(void)
{
  // defines the command line string
  string aux, pars;

  // if we shoud set the library path
  //if(m_setLibraryPath) {
  //  aux += "LD_LIBRARY_PATH=";
  //  aux += m_rendererLocation;
  //  aux += ':';
  //  aux += m_rendererLocation;
  //  aux += DIRECTORY_SEPARATOR;
  //  aux += "lib";
  //  aux += ':';
  //  aux += m_rendererLocation;
  //  aux += DIRECTORY_SEPARATOR;
  //  aux += "..";
  //  aux += DIRECTORY_SEPARATOR;
  //  aux += "lib";
  //  aux += ' ';
  //}

  aux += m_rendererLocation;
  aux += DIRECTORY_SEPARATOR;
  aux += GNU_PLOT_COMMAND;
  
  string aux2 = m_gnuplotCommandFileName;
  //char *beg = "./";
  //if(aux2.compare(0, 2, beg) == 0)
  //  aux2 = aux2.substr(2);
  
  pars  = '"';
  pars += aux2;
  pars += '"';

  // if we shoud set the library path
  if(m_setLibraryPath) {
    string variableValue;
    variableValue  = m_rendererLocation;
    variableValue += ':';
    variableValue += m_rendererLocation;
    variableValue += DIRECTORY_SEPARATOR;
    variableValue += "lib";
    variableValue += ':';
    variableValue += m_rendererLocation;
    variableValue += DIRECTORY_SEPARATOR;
    variableValue += "..";
    variableValue += DIRECTORY_SEPARATOR;
    variableValue += "lib";

    addEnvironmentVariable(aux, "LD_LIBRARY_PATH", variableValue, false);
  }

  // Call GnuRenderer
  int ret = spsSystem( aux.c_str() , pars.c_str() );
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::callGnuplot2(void)
{
  // declare cstream object
  //cstream splot;
  // write gnu command stream to cstream output
  //splot << m_gnuplotQueue.str();
  // write stream to stdin
  //cout <<  splot.dup2(stdin) << endl;

  // prepara parameter list for gnuplot
  char *caux[3];
  caux[0] = NULL;
  caux[1] = (char *)m_gnuplotCommandFileName.c_str();
  caux[2] = NULL;

  // call gnuplot
  if(gnuplot_main)
    gnuplot_main(2, caux);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::deleteGnuplotCommandFile(void)
{
  remove( m_gnuplotCommandFileName.c_str() );
}
////////////////////////////////////////////////////////////////////////////////
// Execute commands
////////////////////////////////////////////////////////////////////////////////
int RendererGnu::execute(void)
{
  // Sets Gnuplot command file name
  setGnuplotCommandFileName();
  // Generate gnuplot command queue
  createGnuplotCommandList();
  // Create Gnuplot command file
  writeGnuplotFile();
  // Execute it
//#if defined(__MINGW32__) || defined(__CYGWIN__)
//  callGnuplot2();
//#else
  // Execute gnuplot method
  callGnuplot();
//#endif
  // clean up
  if(!m_debug)
    deleteGnuplotCommandFile();

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
// Command building section
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::clear(stringstream &out)
{
  out.clear();
  m_terminatorIssued = false;
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::outputVector(stringstream &out, vector<string> &data)
{
  for(int i = 0 ; i < data.size() ; i++) {
    // i == 0 means it's the first item, so no comma needed.
    if(i != 0) out << ',';
    out << '\"' << data[i] << '\"';
  }
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::outputDoubleVector(stringstream &out, vector<pair<string,string> > &data)
{
  for(int i = 0 ; i < data.size() ; i++) {
    // i == 0 means it's the first item, so no comma needed.
    if(i != 0) out << ',';
    out << '\"' << data[i].first << "\" " << data[i].second;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Component methods
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setTitle(stringstream &out, const char *title)
{
  out << "set title \"" << title << "\"\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setTitle(stringstream &out, const RendererTitle &title)
{
  out << "set title \"" << title.label << "\"";
  //out << "set label \"" << title.label << "\"";
  out << " " << title.offsetX << "," << title.offsetY;
  //out << " at screen " << 0.3 + title.offsetX << "," << 0.95 + title.offsetY;
  if(m_rendererImage.fileType.compare("eps") != 0)
    if(title.fontName.length() > 0) {
      out << " font \"";
      out << m_fontLocation << UNIX_DIRECTORY_SEPARATOR << title.fontName;
      out << "," << title.fontSize << "\"";
    }
  out << "\n";
}
////////////////////////////////////////////////////////////////////////////////
// Draw a single line
void RendererGnu::drawLine(stringstream &out, const RendererLine &line)
{
  drawArrowInternal(out, line);
  out << " nohead \n";
}
////////////////////////////////////////////////////////////////////////////////
// Draw an arrow
void RendererGnu::drawArrow(stringstream &out, const RendererLine &line)
{
  drawArrowInternal(out, line);
  float size = m_rendererImage.pixelWidth * 7;
  out << " head filled size screen " << size << ",20,45 ";
  out << " \n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::drawArrowInternal(stringstream &out, const RendererLine &line)
{
  out << "set arrow from ";
  drawPoint(out, line.start);
  out << " to ";
  drawPoint(out, line.end);
  out << " lt " << line.type;
  out << " lw " << line.width;
  if(line.color.length())
    out <<" lc rgb \"" << line.color << "\"";
  outputOrder(out, line.order);
}
////////////////////////////////////////////////////////////////////////////////
// Output a point - pair of coordinates
void RendererGnu::drawPoint(stringstream &out, const RendererPoint &RendererPoint)
{
  drawCoordinate(out, RendererPoint.x);
  out << ", ";
  drawCoordinate(out, RendererPoint.y);
}
////////////////////////////////////////////////////////////////////////////////
// Output coordinate
void RendererGnu::drawCoordinate(stringstream &out, const RendererCoordinate &coordinate)
{
  switch(coordinate.type) {
    case RENDERER_COORD_FIRST:
      out << "first";
      break;
    case RENDERER_COORD_SECOND:
      out << "second";
      break;
    case RENDERER_COORD_SCREEN:
      out << "screen";
      break;
    case RENDERER_COORD_GRAPH:
      out << "graph";
      break;
    case RENDERER_COORD_NONE:
    default:
      break;
  }
  out << " " << coordinate.value;
}
////////////////////////////////////////////////////////////////////////////////
// output order
void RendererGnu::outputOrder(stringstream &out, const RendererDataOrder &order)
{
  switch(order) {
  case RENDERER_DATA_ORDER_NONE:
    break;
  case RENDERER_DATA_ORDER_FRONT:
    out << " front";
    break;
  case RENDERER_DATA_ORDER_BACK:
    out << " back";
    break;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Draw a label
void RendererGnu::drawLabel(stringstream &out)
{
}
////////////////////////////////////////////////////////////////////////////////
// Draw a curve
void RendererGnu::drawCurve(stringstream &out)
{
}
////////////////////////////////////////////////////////////////////////////////
// Draw graph axis
void RendererGnu::drawAxis(stringstream &out)
{
}
////////////////////////////////////////////////////////////////////////////////
// Draw data
void RendererGnu::drawData(stringstream &out, vector<RendererData> &data)
{
  out << "\n";

  // header section
  if(data.size() == 0) {
    out << "plot \"-\" binary record=0";
  } else {
    // Output the several data point group headers
    for(int i = 0 ; i < data.size() ; i++) {

      // lead command
      if(i == 0)
        out << "plot ";
      else
        out << ", ";

      // number of records
      int dataSize = data[i].dataRecords;

      // line specs
      out << "\"-\" binary record=";
      out << dataSize;
      out << " ";
      out << "format=";

      m_terminatorIssued = true;

      switch(data[i].dataType) {
      case RENDERER_DATA_TYPE_INT:
        out << "\"%int%int\" ";
        break;
      case RENDERER_DATA_TYPE_DOUBLE:
        out << "\"%double%double\" ";
        break;
      }
      out << "with impulses ";

      // line type
      out << "lt " << data[i].lineType << " ";

      // line width
      out << "lw " << data[i].lineWidth << " ";

      // line color
      if(data[i].lineColor.length())
        out << "lc rgb \"" << data[i].lineColor << "\" ";

      // title
      if(data[i].title.length() == 0)
        out << "notitle";
      else
        out << "title \"" << data[i].title << "\"";
    }


    // data section
    out << "\n";
    for(int i = 0 ; i < data.size() ; i++) {

      // number of records
      int dataSize = data[i].dataRecords;

      switch(data[i].dataType) {

      case RENDERER_DATA_TYPE_INT:

        for(int j = 0 ; (j < data[i].dataInt.size()) && (j < dataSize) ; j++) {
          out.write((char*)&data[i].dataInt[j].first, sizeof(int));
          out.write((char*)&data[i].dataInt[j].second, sizeof(int));
        }

        break;

      case RENDERER_DATA_TYPE_DOUBLE:

        for(int j = 0 ; (j < data[i].dataDouble.size() ) && (j < dataSize) ; j++) {
          out.write((char*)&data[i].dataDouble[j].first, sizeof(double));
          out.write((char*)&data[i].dataDouble[j].second, sizeof(double));
        }

        break;
      }

    }
  }
}
////////////////////////////////////////////////////////////////////////////////
// Set the margins
void RendererGnu::setMarginTop(stringstream &out, const double &margin, bool scr)
{
  out << "set tmargin ";
  if(scr)
    out << "at screen ";
  out << margin << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setMarginBottom(stringstream &out, const double &margin, bool scr)
{
  out << "set bmargin ";
  if(scr)
    out << "at screen ";
  out << margin << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setMarginLeft(stringstream &out, const double &margin, bool scr)
{
  out << "set lmargin ";
  if(scr)
    out << "at screen ";
  out << margin << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setMarginRight(stringstream &out, const double &margin, bool scr)
{
  out << "set rmargin ";
  if(scr)
    out << "at screen ";
  out << margin << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setMarginAll(stringstream &out, const double &margin, bool scr)
{
  setMarginTop(out, margin, scr);
  setMarginBottom(out, margin, scr);
  setMarginLeft(out, margin, scr);
  setMarginRight(out, margin, scr);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setBorder(stringstream &out, const int &border)
{
  out << "set border " << border << "\n";
}
////////////////////////////////////////////////////////////////////////////////
// Set axis range
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setXrange(stringstream &out, const double &rangeStart, const double &rangeEnd)
{
  out << "set xrange";
  setAxisRange(out, rangeStart, rangeEnd);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setXrange(stringstream &out, const double &rangeStart)
{
  out << "set xrange";
  setAxisRange(out, rangeStart);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setYrange(stringstream &out, const double &rangeStart, const double &rangeEnd)
{
  out << "set yrange";
  setAxisRange(out, rangeStart, rangeEnd);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setYrange(stringstream &out, const double &rangeStart)
{
  out << "set yrange";
  setAxisRange(out, rangeStart);
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setAxisRange(stringstream &out, const double &rangeStart, const double &rangeEnd)
{
  out << " [" << rangeStart << ":" << rangeEnd <<"]\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setAxisRange(stringstream &out, const double &rangeStart)
{
  out << " [" << rangeStart << ":]\n";
}
////////////////////////////////////////////////////////////////////////////////
// Set axis labels
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setXlabel(stringstream &out, const string &str)
{
  setXlabel(out, str.c_str());
}

void RendererGnu::setXlabel(stringstream &out, const char *str)
{
  out << "set xlabel \"" << str << "\"\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setYlabel(stringstream &out, const string &str)
{
  setYlabel(out, str.c_str());
}

void RendererGnu::setYlabel(stringstream &out, const char *str)
{
  out << "set ylabel \"" << str << "\"\n";
}
////////////////////////////////////////////////////////////////////////////////
// set axis intervals
void RendererGnu::setXinterval(stringstream &out)
{
  out << "set xtics axis nomirror tc rgb \"dark-gray\"\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setYinterval(stringstream &out)
{
  out << "set ytics axis nomirror tc rgb \"dark-gray\"\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setXinterval(stringstream &out, const double &val)
{
  out << "set xtics axis nomirror tc rgb \"dark-gray\" " << val << "\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setYinterval(stringstream &out, const double &val)
{
  out << "set ytics axis nomirror tc rgb \"dark-gray\" " << val << "\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setXinterval(stringstream &out, vector<pair<string,string> > &data)
{
  out << "set xtics (";
  outputDoubleVector(out, data);
  out << ") out nomirror\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setYinterval(stringstream &out, vector<pair<string,string> > &data)
{
  out << "set ytics (";
  outputDoubleVector(out, data);
  out << ") out nomirror\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setXintervalOut(stringstream &out)
{
  //out << "set xtics out\n";
  out << "unset xtics\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setYintervalOut(stringstream &out)
{
  //out << "set ytics out\n";
  out << "unset ytics\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setOrigin(stringstream &out, double &x, double &y)
{
  out << "set origin " << x << ", " << y << "\n";
}
/*----------------------------------------------------------------------------*/
void RendererGnu::setSize(stringstream &out, double &x, double &y)
{
  out << "set size " << x << ", " << y << "\n";
}
////////////////////////////////////////////////////////////////////////////////
void RendererGnu::setLabel(stringstream &out, const RendererLabel &label)
{
  // label command and content
  out << "set label \"";
  out << label.label;
  out << "\"";

  // label position
  out << " at ";
  drawPoint(out, label.position);

  // label offset
  out << " offset ";
  out << label.offsetX;
  out << ',';
  out << label.offsetY;
  out << " ";

  // offset position
  switch(label.location) {
  case RENDERER_OFFSET_RIGHT:
    out << "right";
    break;
  case RENDERER_OFFSET_CENTER:
    out << "center";
    break;
  case RENDERER_OFFSET_LEFT:
    out << "left";
    break;
  default:
    break;
  }

  // order
  outputOrder(out, label.order);

  // Color
  if(label.color.length())
    out << " tc rgb " << label.color;
  out << "\n";
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
////////////////////////////////////////////////////////////////////////////////
