///////////////////////////////////////////////////////////////////////////////
#include "ContigSequence.h"
#include "Tokenizer.h"

#include <limits>

///////////////////////////////////////////////////////////////////////////////
using namespace spsReports;
using namespace std;

namespace specnets {

///////////////////////////////////////////////////////////////////////////////
// Drawing procedure entry point with no parameters
///////////////////////////////////////////////////////////////////////////////
int ContigSequence::drawExec(void)
{
  // Check for spectrum
  if( (m_sequence == NULL) || (m_breaks == NULL) )
    return ERROR;

  drawDashedLines();

  drawTopArrows();

  drawPeptideLabel();

  drawPeptide();

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::addDashedSegment(double x1, double y1, double x2, double y2)
{
  drawLine(  1, 0,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_SCREEN,
             m_colour);
}
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::addDashedSegmentBent(double x1, double y1, double x2, double y2)
{
  double disp2 = (double)m_rendererImage.fontHeight / (double)m_rendererImage.imageSizeY / m_rendererImage.zoom * 1.0;
  double disp3 = (m_mzUpperLimit - m_mzLowerLimit) * 0.01;

  double x11 = x1;
  double y11 = y1;
  double x21 = x1 + disp3;
  double y21 = y1 - disp2;

  double x12 = x1 + disp3;
  double y12 = y1 - disp2;
  double x22 = x2 + disp3;
  double y22 = y2 + disp2;

  double x13 = x2 + disp3;
  double y13 = y2 + disp2;
  double x23 = x2;
  double y23 = y2;


  drawLine(  1, 0,
             x11, RENDERER_COORD_NONE,
             y11, RENDERER_COORD_SCREEN,
             x21, RENDERER_COORD_NONE,
             y21, RENDERER_COORD_SCREEN,
             m_colour);


  drawLine(  1, 0,
             x12, RENDERER_COORD_NONE,
             y12, RENDERER_COORD_SCREEN,
             x22, RENDERER_COORD_NONE,
             y22, RENDERER_COORD_SCREEN,
             m_colour);

  drawLine(  1, 0,
             x13, RENDERER_COORD_NONE,
             y13, RENDERER_COORD_SCREEN,
             x23, RENDERER_COORD_NONE,
             y23, RENDERER_COORD_SCREEN,
             m_colour);
}
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::drawDashedLines(void)
{
  // nothing to draw if previous breaks is null
  if(!m_previousBreaks) return;
  // go thru all breaks
  for(int i = 0 ; i < m_breaks->size() ; i++)

    for(int j = 0 ; j < m_previousBreaks->size() ; j++) {

      //cout << "sequence: " << m_offset + (*m_breaks)[i] << " <-> " << m_offsetPrevious + (*m_previousBreaks)[j];

      if( fabs( m_offset + (*m_breaks)[i] - m_offsetPrevious - (*m_previousBreaks)[j]) < m_peakMassTol)

        //cout << "...ok";

        // draw the vertical dashed line
        if(m_bendLine)
          addDashedSegmentBent( m_offset + (*m_breaks)[i], m_yPosition, m_offsetPrevious + (*m_previousBreaks)[j], m_yPositionPrevious );
        else
          addDashedSegment( m_offset + (*m_breaks)[i], m_yPosition, m_offsetPrevious + (*m_previousBreaks)[j], m_yPositionPrevious );

      //cout << endl;
    }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::addRedArrow(double x1, double y1, double x2, double y2)
{
  drawArrow( 2, 1,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_SCREEN,
             m_colour);
}
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::drawTopArrows(void)
{
  if(!m_breaks->size())
    return;

  double first, second;

  first = (*m_breaks)[0];
  for(int i = 1 ; i < m_breaks->size() ; i++) {
    second = (*m_breaks)[i];
    addRedArrow(m_offset + first, m_yPosition, m_offset + second, m_yPosition);
    first = second;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ContigSequence::cleanSequence(string &iSequence)
{
  string aux;
  string valid("ARNDCEQGHILKMFPSTWYV0123456789,.-[]()");
  for(int i = 0 ; i < iSequence.length() ; i++)
    if(valid.find_first_of(iSequence[i]) != string::npos)
      aux += iSequence[i];
  iSequence = aux;
}
///////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ContigSequence::breakSequence(vector<string> &oSequenceMapping, string &iSequence, bool omitValues, bool useParantesis)
{
  string aa("ARNDCEQGHILKMFPSTWYV");
  string aux;
  string cell;

  oSequenceMapping.clear();

  // store the position after each tag
  size_t lastPosition = 0;

  while(lastPosition < iSequence.size()) {

    cell = "";

    // in case of an AA
    if(aa.find_first_of(iSequence[lastPosition]) != string::npos) {

      cell = iSequence[lastPosition];

      lastPosition++;

    // in case of mass
    } else if(iSequence[lastPosition] == '[') {

      // sequence initiator
      string strValue = "[";

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find end
      while(currentPosition < iSequence.size() && iSequence[currentPosition] != ']')
        currentPosition++;

      // get contents
      double value;
      string mass = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
      strValue += roundMasses(mass, 1, value);

      // update last position
      lastPosition = currentPosition + 1;

      // sequence terminator
      strValue += ']';

      if(!omitValues && value > 30.0)
       cell = strValue;
      else
        cell = '.';

    } else if(iSequence[lastPosition] == '(') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find comma
      while(currentPosition < iSequence.size() && iSequence[currentPosition] != ',' && iSequence[currentPosition] != ')')
        currentPosition++;

      // get sub sequence contents
      string cellAux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);

      lastPosition = currentPosition + 1;

      bool massDeviation = false;

      if(iSequence[currentPosition] == ',') {
        // find end
        while(currentPosition < iSequence.size() && iSequence[currentPosition] != ')')
          currentPosition++;

        // add the comma
        if(!omitValues)
          cellAux += ',';

        // get contents
        if(!omitValues) {
          string mass = iSequence.substr(lastPosition, currentPosition - lastPosition);
          double value;
          cellAux += roundMasses(mass, 1, value);
        }

        // update last position
        lastPosition = currentPosition + 1;

        // we had a comma and a mass deviation
        bool massDeviation = true;
      }

      // initial
      if(useParantesis || massDeviation)
        cell = '(';

      // contents
      cell += cellAux;

      // terminator
      if(useParantesis || massDeviation)
        cell += ')';

    } else { // ERROR
      lastPosition++;
    }

    // if we're just drawing dots, letit be a dot.
    if(m_allDots)
      cell = '.';

    // add the cell to the list
    oSequenceMapping.push_back(cell);
  }
}
///////////////////////////////////////////////////////////////////////////////
string ContigSequence::roundMasses(string &strValue, int decimals, double &value)
{
  std::ostringstream ss;

  value = getFloat(strValue.c_str());
  ss << std::fixed << std::setprecision(decimals) << value;
  std::string s = ss.str();
  if(decimals > 0 && s[s.find_last_not_of('0')] == '.') {
      s.erase(s.size() - decimals + 1);
  }
  return s;
}
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::drawPeptide(void)
{
  vector<string> items;

  cleanSequence(*m_sequence);

  breakSequence(items, *m_sequence, m_rendererImage.zoom < 0.9);
  //StringExplode(*m_sequence, items, true, "|");

  if(items.size() < 1 || m_breaks->size() < 2)
    return;

  for(int i = 1 ; i < m_breaks->size() && i <= items.size() ; i++) {

    drawLabel(items[i-1],
              m_offset + ((*m_breaks)[i] + (*m_breaks)[i-1]) / 2, RENDERER_COORD_NONE,
              m_yPosition, RENDERER_COORD_SCREEN,
              0,
              1,
              RENDERER_OFFSET_CENTER,
              "#000000");
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void ContigSequence::drawPeptideLabel(void)
{
  if(!m_drawLabel)
    return;

  drawLabel(m_label,
            0,   RENDERER_COORD_GRAPH,
            m_yPosition,  RENDERER_COORD_SCREEN,
            0, 1.0,
            RENDERER_OFFSET_RIGHT,
            "#000000");
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
