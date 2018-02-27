///////////////////////////////////////////////////////////////////////////////
#ifndef __LABEL_PLACER_H__
#define __LABEL_PLACER_H__
///////////////////////////////////////////////////////////////////////////////

#include "PlotBase.h"

#include "spectrum.h"
#include "aminoacid.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////

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
//
////////////////////////////////////////////////////////////////////////////////
  /*! \brief Spectrum label class. Used to draw the top labels in specplot
   */
class SpectrumLabel : public PlotBase  {


  // vector to hold peak information (for peak label generation)
  /*! \brief vector to hold peak information (for peak label generation)
   */
  vector<PeakLabelItem>           peakLabelList;



  /*! \brief Pleaces the peak labels
   */
  void placePeakLabels(void);

  // Checks if a specific position is already taken (peak label placement)
  /*! \brief Checks if a specific position is already taken (peak label placement)
   */
  bool checkPlacementPosition(ParamsData2 &params);
  // Checks all the possible positions for a peak label
  /*! \brief Checks all the possible positions for a peak label
   */
  void checkPlacementPositions(ParamsData2 &params);


  // mz lower and upper limits
  /*! \brief mz lower limit
   */
  double m_mzLowerLimit;
  /*! \brief mz upper limit
   */
  double m_mzUpperLimit;

  
 protected:

  // Constructors and Destructor
  /*! \brief Constructors
   */
  SpectrumLabel(RendererBase *renderer, double ll, double ul) :
    m_mzLowerLimit(ll), m_mzUpperLimit(ul)
    {m_rendererObject = renderer; m_ownRenderer = false;};

  /*! \brief Destructor
   */
  ~SpectrumLabel() {};


  // Default draw object entry point
  /*! \brief Default draw object entry point
   */
  virtual void drawExec(void);


};
///////////////////////////////////////////////////////////////////////////////
}; //namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
