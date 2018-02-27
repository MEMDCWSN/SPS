///////////////////////////////////////////////////////////////////////////////
#include "ContigSpectrum.h"
#include "aminoacid.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
void ContigSpectrum::addGreenArrow(double x1, double y1, double x2, double y2)
{
  drawArrow( 2, 1,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_GRAPH,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_GRAPH,
             "#00A000",
             RENDERER_DATA_ORDER_FRONT);
}
///////////////////////////////////////////////////////////////////////////////
void ContigSpectrum::addRedArrow(double x1, double y1, double x2, double y2)
{
  drawArrow( 2, 1,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_GRAPH,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_GRAPH,
             "#A00000",
             RENDERER_DATA_ORDER_FRONT);
}
///////////////////////////////////////////////////////////////////////////////
void ContigSpectrum::addRedLabel(double x1, double y1, string &label)
{
      // Draw the label
      drawLabel(label,
                x1, RENDERER_COORD_NONE,
                y1, RENDERER_COORD_GRAPH,
                0.0, 1.0,
                RENDERER_OFFSET_CENTER,
                "#A00000");
}
///////////////////////////////////////////////////////////////////////////////
void ContigSpectrum::buildVector(void)
{
  if(!abData)
    return;
  //const vector<pair< vector<int>, vector<double> > >  &abData = (*m_abinfo)[m_contigIndex].second;

  for(int i = 0 ; i < abData->size() ; i++)
    for(int j = 0 ; j < (*abData)[i].first.size() ; j++)
      if( (*abData)[i].first[j] == m_spectrumIndex) {
        abMassItem item;
        item.mass = (*abData)[i].second[j];
        item.idxInContig = i;
        m_abMasses.push_back(item);
        //m_abMasses.push_back( (*abData)[i].second[j]);
      }
}
///////////////////////////////////////////////////////////////////////////////
int ContigSpectrum::find(double &value)
{
  int a = 0;
  int b = m_abMasses.size() - 1;
  int m = (a + b) / 2;

  while(a < b) {

    if( fabs(m_abMasses[m].mass - value) < m_peakMassTol )
      return m_abMasses[m].idxInContig;

    else if(value < m_abMasses[m].mass)
      b = m - 1;

    else // if(m_abMasses[m] < value)
      a = m + 1;

    m = (a + b) / 2;
  }

  // final item
  if(m >= 0 && m < m_abMasses.size())
    if( fabs(m_abMasses[m].mass - value) < m_peakMassTol )
      return m_abMasses[m].idxInContig;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
void ContigSpectrum::draw(double limit, double factor)
{
  if(!m_spectrum || !m_deNovoIntervals)
    return;

  double peak, intensity;

  vector<RendererData> plotData;
  RendererData  rendererData;
  string color;
  double first, second;
  bool firstFound = true;
  int  lastIndex = 0;

  rendererData.lineWidth   = m_rendererImage.lineThickness; //4;
  rendererData.lineType    = 1;
  rendererData.dataType    = RENDERER_DATA_TYPE_DOUBLE;
  rendererData.dataRecords = 0;

  int size = m_spectrum->size();

  // the two initial and final bar
  rendererData.lineColor = "#0000C0";
  intensity = limit * factor;

  peak = 0.0;
  rendererData.dataDouble.push_back(make_pair<double,double>(peak * factor, intensity));
  rendererData.dataRecords++;

  peak = fmax(m_spectrum->parentMass, (*m_spectrum)[size-1][0]) - AAJumps::massMH;
  //peak = m_spectrum->parentMass;
  rendererData.dataDouble.push_back(make_pair<double,double>(peak * factor, intensity));
  rendererData.dataRecords++;


  //plotData.push_back(rendererData);
  //rendererData.dataRecords = 0;

  buildVector();

  for(register int j = 0 ; j < size ; j++) {

    peak = (*m_spectrum)[j][0];

    //cout << peak << endl;
    int index = find(peak);

    if((index != -1) && (index < m_deNovoIntervals->size())) {
      if(firstFound) {
        first = peak;
        firstFound = false;
      } else {
        second = peak;

        double starPeakDiff = second - first;
        double contPeakDiff = (*m_deNovoIntervals)[index] - (*m_deNovoIntervals)[lastIndex];
        double diff = starPeakDiff - contPeakDiff;

        //cout << "starPeakDiff: " << second << " - " << first << " = " << starPeakDiff << endl;
        //cout << "contPeakDiff: " << (*m_deNovoIntervals)[index] << " - " << (*m_deNovoIntervals)[lastIndex] << " = " << contPeakDiff << endl;
        //cout << "diff: " << starPeakDiff << " - " << contPeakDiff << " = " << diff << endl;

        double auxX = first * factor;
        double auxY = second * factor;

        // define peak mass tolerace for this peak
        double peakMassTol = m_peakMassTol;
        if(m_useIndividualPeakMassTol)
          m_spectrum->getTolerance(index);

        if( fabs( diff ) < peakMassTol )
          addGreenArrow(auxX, 0.0, auxY, 0.0);
        else {
          addRedArrow(auxX, 0.0, auxY, 0.0);

          double lx = (second + first) / 2.0;
          double ly = 0.0;
          string label = "[";
          if(contPeakDiff < starPeakDiff) label += '+';
          label += parseFloat(diff, 2);
          label += ']';
          addRedLabel(lx * factor, ly, label);
        }

        first = second;
      }
      lastIndex = index;
      color = "#C00000";
    } else
      color = "#000000";


    if(rendererData.lineColor.compare(color)) {
      //if(j)
        plotData.push_back(rendererData);
      rendererData.dataDouble.clear();
      rendererData.lineColor = color;
      rendererData.dataRecords = 0;
    }

    intensity = (*m_spectrum)[j][1];

    if(m_useMinimunIntensity && (intensity < limit * 0.2))
      intensity = limit * 0.2;

    rendererData.dataDouble.push_back(make_pair<double,double>(peak * factor, intensity));
    rendererData.dataRecords++;
  }

  // Add last peaks
  plotData.push_back(rendererData);

  // Draw the graph
  m_rendererObject->drawData(plotData);
}
///////////////////////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////////////////////
