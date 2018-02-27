///////////////////////////////////////////////////////////////////////////////
#include "RendererBase.h"
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <math.h>
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
///////////////////////////////////////////////////////////////////////////////
RendererBase::RendererBase(void)
{
  initialize();
}
// Destructor
RendererBase::~RendererBase(void)
{
}
///////////////////////////////////////////////////////////////////////////////
void RendererBase::initialize(void)
{
  m_rendererLocation  = '.';
  m_fontLocation      = '.';
  m_useTitle          = false;
  m_setLibraryPath    = true;
}
///////////////////////////////////////////////////////////////////////////////
double RendererBase::calculateInterval(double distance, double factor, double zoom)
{
  // distance according to zoom and factor
  double temp = distance / zoom / factor;
  // logaritmic scale
  double mag  = floor(log10(temp));
  // back to linear
  double powX  = pow(10.0, mag);
  //
  double msd  = unsigned(temp / powX + 0.5);

  // stepping
  if(msd > 5.0)
    msd = 10.0;
  else if(msd > 2.0)
    msd = 5.0;
  else if (msd > 1.0)
    msd = 2.0;

  return msd * powX;
}
///////////////////////////////////////////////////////////////////////////////
void RendererBase::breakAxisY(const double until, const double restart)
{
  // recalculate y axis tics for lower section
  double msdY = calculateInterval(until - m_areaCurrent.m_startY, 8.0, m_rendererImage.zoom);
  setYinterval(msdY);

  // store data & keep current
  m_areas.push_back(m_areaCurrent);
  // location of last object
  int loc = m_areas.size() - 1;

  // recalculate y axis tics for upper section
  msdY = calculateInterval(m_areaCurrent.m_endY - restart, 8.0, m_rendererImage.zoom);
  setYinterval(msdY);


  // space beetwen the 2 parts
  double gap = 0.005;

  // image size relative to total image
  double size = (m_rendererImage.imageSizeY - m_areaCurrent.m_marginTop - m_areaCurrent.m_marginBottom) / m_rendererImage.imageSizeY;

  double disp = m_areaCurrent.m_marginTop / m_rendererImage.imageSizeY;

  double breakPoint = 0.9;

  // adjust stored area (top)
  m_areas[loc].m_useSize  = true;
  m_areas[loc].m_sizeX    = 1.0;
  m_areas[loc].m_sizeY    = 1 - breakPoint + gap - disp;

  m_areas[loc].m_useOrigin = true;
  m_areas[loc].m_originX  = 0.0;
  m_areas[loc].m_originY  = breakPoint - gap - disp - m_areas[loc].m_sizeY;

  m_areas[loc].m_startY   = restart;

  m_areas[loc].m_useTicsX = false;
  m_areas[loc].m_border  &= 0x0E;
  m_areas[loc].m_border   = 0x0F;
  m_areas[loc].m_xLabel   = "";


  // adjust current area (bottom)
  //setSize(1.0, ratio * size  - gap);
  //setOrigin(0.0, 1.0 - ratio + gap - disp1);
  setSize(1.0, breakPoint);
  setOrigin(0.0, 0.0);
  m_areaCurrent.m_endY     = until;
  m_areaCurrent.m_useTicsX = true;
  m_areaCurrent.m_border  &= 0x0B;
  m_areaCurrent.m_border   = 0x0F;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
