///////////////////////////////////////////////////////////////////////////////
#include "SpectrumLabel.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
bool SpectrumLabel::checkPlacementPosition(ParamsData2 &params)
{
  // Check map boundaries. If outside, can't do it
  if( (params.currentX < 0) || (params.currentX >= params.matX - params.labelSizeX) ||
      (params.currentY < 0) || (params.currentY >= params.matY - params.labelSizeY) )
    return false;

  // check label placement position against matrix
  for(int j = 0 ; j < params.labelSizeX ; j++)
    for(int k = 0 ; k < params.labelSizeY ; k++)
      // if occupied by something else, can't use this position
      if(params.matrix[params.currentY + j][params.currentX + k] != 0)
        return false;
  return true;
}
///////////////////////////////////////////////////////////////////////////////
void SpectrumLabel::checkPlacementPositions(ParamsData2 &params)
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

      // calculate positions in graph
      peakLabelList[params.i].x = currentX * factorX + m_mzLowerLimit;
      peakLabelList[params.i].y = currentY * factorY;

      // mark in matrix
      for(int j = 0 ; (j <  params.labelSizeY) && (j < params.matY - params.currentY) ; j++)
        for(int k = 0 ; (k <  params.labelSizeX) && (k < params.matX - params.currentX) ; k++)
          params.matrix[params.currentY + j][params.currentX + k] = 2;

      // mark label as placed.
      peakLabelList[params.i].placed = true;

      // If position is right above the peak, no need to draw the conneting line.
      peakLabelList[params.i].drawConnectionLine = position;
      return;
    }
  }
  // if no position found, mark as 'not placed' -> default by aquirePeakLabels()
  //peakLabelList[params.i].placed = false;
}
///////////////////////////////////////////////////////////////////////////////
void SpectrumLabel::placePeakLabels(void)
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
};
///////////////////////////////////////////////////////////////////////////////
