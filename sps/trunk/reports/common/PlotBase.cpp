///////////////////////////////////////////////////////////////////////////////
#include "PlotBase.h"
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "PlotBase.h"
#include "RendererGnu.h"
#include "utils.h"
#include "ReportBase64.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
// Constructors and destructor
PlotBase::PlotBase() :
  m_rendererObject(NULL),
  m_ownRenderer(false),
  m_titlePos(1),
  m_target("file"),
  m_fileFormat("png"),
  m_fileOutDir("."),
  m_fileName("file.png"),
  m_peakMassTol(0.45),
  m_PeakMassTolSet(false),
  m_useIndividualPeakMassTol(false),
  m_massShift(0),
  m_annotationModel(DEFAULT_ANNOTATION_MODEL),
  m_annotationModelDirectory(DEFAULT_ANNOTATION_FILE_LOCATION),
  m_readModelFromData(true),
  m_annotationStyle(ANNOTATION_STYLE_SPECNETS),
  m_debug(0),
  m_epsCorrection(false),
  m_useMinimunIntensity(false),
  m_verbose(false)
{
  // Call base class initizalizer
  initialize();
}

PlotBase::~PlotBase()
{
  if(m_rendererObject != NULL && m_ownRenderer)
    delete m_rendererObject;
}
///////////////////////////////////////////////////////////////////////////////
// initialize object
void PlotBase::initialize(void)
{
  // Set image size
  m_rendererImage.imageSizeX  = 800;
  m_rendererImage.imageSizeY  = 600;

  // set default stretch limit to NO_STRETCH
  m_rendererImage.imageStretchLimitX = -2;
  m_rendererImage.imageStretchLimitY = -2;

  // set pixel size
  initializePixelSize();

  // set the zoom
  m_rendererImage.zoom = 1.0;

  // Font name, size and default color
  m_rendererImage.fontName    = "VeraBd.ttf";
  m_rendererImage.fontHeight  = 12;
  m_rendererImage.fontWidth   = m_rendererImage.fontHeight * 2 / 3;
  m_rendererImage.fontColor   = 0x00;

  m_rendererImage.fontSize    = roundDouble(m_rendererImage.fontHeight / 1.55 , 1);

  if(m_rendererObject)
    m_rendererObject->initialize();
}
///////////////////////////////////////////////////////////////////////////////
void PlotBase::initializePixelSize(void)
{
  // set pixel size
  m_rendererImage.pixelWidth  = 1.0 / m_rendererImage.imageSizeX;
  m_rendererImage.pixelHeight = 1.0 / m_rendererImage.imageSizeY;
}
///////////////////////////////////////////////////////////////////////////////
// Define renderer
///////////////////////////////////////////////////////////////////////////////
void PlotBase::setDrawMethod(RendererType rendererType)
{
  if(m_rendererObject != NULL && m_ownRenderer)
    delete m_rendererObject;

  switch(rendererType) {
  case RENDERER_TYPE_NONE:
    m_rendererObject = NULL;
    m_ownRenderer = false;
    break;
  case RENDERER_TYPE_GNU:
    m_rendererObject = new RendererGnu();
    m_ownRenderer = true;
    break;
  default:
    m_rendererObject = NULL;
    m_ownRenderer = false;
  }
}
///////////////////////////////////////////////////////////////////////////////
// Initialize draw object properties
///////////////////////////////////////////////////////////////////////////////
void PlotBase::setDrawProperties(void)
{
  // set the file format
  m_rendererImage.fileType = m_fileFormat;

  m_rendererObject->m_rendererImage = m_rendererImage;

  m_rendererObject->setDebug(m_debug);

  m_rendererObject->setEpsCorrection(m_epsCorrection);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotBase::drawLabel( string theLabel,
                          double x1, RendererCoordType coordX1,
                          double y1, RendererCoordType coordY1,
                          double offSetX,
                          double offSetY,
                          RendererLabelOffsetType location,
                          string color,
                          RendererDataOrder order)

{
  RendererLabel label;

  label.label               = theLabel;
  label.position.x.value    = x1;
  label.position.x.type     = coordX1;
  label.position.y.value    = y1;
  label.position.y.type     = coordY1;
  label.color               = "\"" + color + "\"";
  label.offsetX             = offSetX;
  label.offsetY             = offSetY;
  label.location            = location;
  label.order               = order;
  // Draw the label
  m_rendererObject->addLabel(label);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotBase::drawArrow(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order)
{
  // Line object
  RendererLine arrow;

  arrow.width          = lineWidth;
  arrow.type           = lineType;
  arrow.start.x.value  = x1;
  arrow.start.x.type   = coordX1;
  arrow.start.y.value  = y1;
  arrow.start.y.type   = coordY1;
  arrow.end.x.value    = x2;
  arrow.end.x.type     = coordX2;
  arrow.end.y.value    = y2;
  arrow.end.y.type     = coordY2;
  arrow.color          = color;
  arrow.order          = order;
  // Draw the line
  m_rendererObject->drawArrow(arrow);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotBase::drawLine(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order)
{
  // Line object
  RendererLine line;

  line.width          = lineWidth;
  line.type           = lineType;
  line.start.x.value  = x1;
  line.start.x.type   = coordX1;
  line.start.y.value  = y1;
  line.start.y.type   = coordY1;
  line.end.x.value    = x2;
  line.end.x.type     = coordX2;
  line.end.y.value    = y2;
  line.end.y.type     = coordY2;
  line.color          = color;
  line.order          = order;
  // Draw the line
  m_rendererObject->drawLine(line);
}
///////////////////////////////////////////////////////////////////////////////
// Draw
///////////////////////////////////////////////////////////////////////////////
void PlotBase::draw(char *fileName)
{
  // set the file name and execute the draw method
  m_fileName = fileName;
  draw();
}
///////////////////////////////////////////////////////////////////////////////
void PlotBase::draw(void)
{
  // execute the drawexec method
  if(drawExec() == ERROR)
    return;

  // if a file is requested, do nothing
  if(m_target.compare("file") == 0)
    return;

  // open the file for reading (binary)
  int length;
  char * buffer;

  string filename;
  // compose filename
  if(m_fileOutDir.length() > 0) {
    // Set output directory, if defined
    filename = m_fileOutDir;
    if(m_fileOutDir[m_fileOutDir.length()-1] != '/')
      filename += '/';
  }
  // Add the filename
  filename += m_fileName;

  ifstream is;
  is.open (filename.c_str(), ios::binary );

  // check for file open
  if(!is.is_open())
    return;

  // get length of file:
  is.seekg (0, ios::end);
  length = is.tellg();
  is.seekg (0, ios::beg);

  // allocate memory:
  buffer = new char [length+1];

  // read data as a block:
  is.read(buffer, length);
  is.close();
  // terminator
  buffer[length] = 0;


  // if uuencoded, encode
  if(m_encoding.compare("uu64") == 0) {
    m_image = base64_encode(reinterpret_cast<const unsigned char*>(buffer), length);
  } else {
    m_image = (char *)buffer;
  }


  // the output to cout...
  if( m_target.compare("cout") == 0) {

    if(m_encoding.compare("uu64") == 0) {
      cout << m_image;
    } else {
      cout.write((char *)buffer, length);
    }
  }

  // check is the image is just to be stored in the object
  if( m_target.compare("internal") == 0) {
  }


  // remove the file
  remove( filename.c_str() );
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
