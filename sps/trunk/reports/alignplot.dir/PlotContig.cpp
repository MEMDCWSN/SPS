///////////////////////////////////////////////////////////////////////////////
#include <limits>
#include <iostream>
#include <algorithm>

#include "PlotContig.h"
#include "Tokenizer.h"


///////////////////////////////////////////////////////////////////////////////
#define PEAK_COLOR_ANNOTATED_SINGLE "#cccccc"
#define PEAK_COLOR_NOT_ANNOTATED    "#eeeeee"
///////////////////////////////////////////////////////////////////////////////
using namespace spsReports;
using namespace std;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
// Class constructors and destructor
///////////////////////////////////////////////////////////////////////////////
PlotContig::PlotContig() :
  m_star(NULL),
  m_abinfo(NULL),
  m_contigIndex(-1),
  m_referenceOffset(0.0),
  m_homologOffset(0.0),
  m_deNovoOffset(0.0),
  m_userOffset(0.0),
  m_reverse(false)
{
}
///////////////////////////////////////////////////////////////////////////////
PlotContig::~PlotContig()
{
}
///////////////////////////////////////////////////////////////////////////////
// Test mass boundaries
///////////////////////////////////////////////////////////////////////////////
bool PlotContig::testMass(float mass)
{
  return true;
}
///////////////////////////////////////////////////////////////////////////////
// initializer methods
///////////////////////////////////////////////////////////////////////////////
void PlotContig::setMassIntervals(string &strValues, vector<double> &values)
{
  vector<string> items;

  StringExplode(strValues, items, true, "|");

  for(int i = 0 ; i < items.size() ; i++) {
    double aux = getFloat(items[i].c_str());
    values.push_back(aux);
  }
}
///////////////////////////////////////////////////////////////////////////////
// Draw Exec
///////////////////////////////////////////////////////////////////////////////
void PlotContig::drawExec(void)
{
  // Check for spectrum
  if( (m_star == NULL) || (m_abinfo == NULL) )
    return;


  //dump_abruijn(cout, m_abinfo);

  /////////////////////////////////////////////////////////////////////////////
  // set locations

  // Set the renderer location
  m_rendererObject->setRendererLocation(m_rendererLocation);
  // set the font files location
  m_rendererObject->setFontLocation(m_fontLocation);


  // set output file name
  if(m_fileOutDir.length() > 0) {
    // Set output directory, if defined
    m_rendererImage.outputFileName = m_fileOutDir;
    if(m_fileOutDir[m_fileOutDir.length()-1] != '/')
      m_rendererImage.outputFileName += '/';
  }
  // Add the filename
  m_rendererImage.outputFileName += m_fileName;

  // reverse contig, if needed
  if(m_reverse)
    getReversedAbruijn();

  // Init draw object
  setDrawProperties();

  // Draw procedure
  plotContig();

  // Do the drawing
  m_rendererObject->execute();
}
///////////////////////////////////////////////////////////////////////////////
// Draw the object - main draw procedure
///////////////////////////////////////////////////////////////////////////////
void PlotContig::plotContig(void)
{
  // calc contig start
  calcStarInContig();

  // calculate offsets
  initData();

  // Calculates where to draw the labels for the peptides
  initDeNovo();

  // Calculates where to draw the labels for the peptides
  initUserIntervals();

  // Align spectra
  alignSpectra();

  // Calculate M/Z and Intensity limits
  calcultateLimits();

  // Initialize the graph
  initializeGraph();

  // Set the graph title
  setGraphTitle();

  // draw the sequences
  drawSequences();

  // draw vertical lines
  drawVerticalLines();

  // multiplot - store primary data
  m_rendererObject->addArea();

  // Draw spectra
  drawSpectra();
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::getReversedAbruijn(void)
{
  // test for needed files
  if(!m_abinfo) return;

  // contig in abruijn
  abContig_t &ab = (*m_abinfo)[m_contigIndex];
  // contig data in abruijn
  abContigData_t &cd = ab.second;
  abContigData_t nodes;

  // reverse the spectra - travel thru the contig nodes
  for(int i =  cd.size() - 1 ; i >= 0 ; i--) {
    // star IDs for node
    vector<int> IDs;
    // data for node
    vector<double> data;
    // cycle thru the data in the nodes
    for(int j = 0 ; j < cd[i].first.size() ; j++) {
      int id = cd[i].first[j];
      IDs.push_back(id);
      double parentMass = (*m_star)[id].parentMass;
      data.push_back(parentMass - AAJumps::massHion - cd[i].second[j]);
    }
    nodes.push_back(make_pair<vector<int>, vector<double> >(IDs, data));
  }

  ab.second = nodes;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::alignSpectra(void)
{
  // get the contig data section
  const vector<pair< vector<int>, vector<double> > >  &abData = (*m_abinfo)[m_contigIndex].second;

  // go thru all spectra
  for(int i = 0 ; i < m_offsets.size() ; i++) {

    // where to hold data
    vector<AlignmentData>  data;

    int id = m_offsets[i].id;

    // cycle thru all contig nodes
    for(int j = 0 ; j < abData.size() ; j++) {

      const vector<int> & nodeIds = abData[j].first;

      const vector<double> & nodeData = abData[j].second;

      // for each node, cycle thru all items
      for(int k = 0 ; k < nodeIds.size() ; k++) {

        if(nodeIds[k] == id) {

          AlignmentData item;
          item.nodeIndex    = j;
          item.abruijnData  = nodeData[k];
          item.seqsData     = (j < m_deNovoIntervals.size() ? m_deNovoIntervals[j] : -1.0);
          data.push_back(item);

        }
      }
    }

    // abruijn node index
    int   nodeIndex;
    // value in abruijn node
    float abruijnData;
    // value in seqs
    float seqsData;
    // number of diagonals if this node is aligned
    int   disp;

    // align spectrum

    // step 1: calculate the total number of unaligned nodes when aligning by each node
    for(int j = 0 ; j < data.size() ; j++) {

      // get displacement data to align on this node
      float disp = data[j].seqsData - data[j].abruijnData;

      // calculate number of oblic nodes if align by this nodes
      data[j].disp = ( data.size() ? 0 : std::numeric_limits<int>::max() );

      // check all nodes, and count how many are not alligned in this context
      for(int k = 0 ; k < data.size() ; k++)
        // if node is not aligned, increment data.disp in [j] context
        if( fabs(data[k].seqsData - data[k].abruijnData - disp) > m_peakMassTol)
          data[j].disp++;
    }

    // step 2: choose the node with least unaligned nodes
    int bestCount = data.size()+1;
    int bestIndex = -1;

    for(int j = 0 ; j < data.size() ; j++)

      if(data[j].disp < bestCount) {
        bestCount = data[j].disp;
        bestIndex = j;
      }

    // step 3: set displacement
    if(data.size())
      m_offsets[i].offset =  data[bestIndex].seqsData - data[bestIndex].abruijnData;
    else
      m_offsets[i].offset = 0.0;

    // debug info
    // cerr << "-------------- " << id << " -------------- " << endl;
    // for(int j = 0 ; j < data.size() ; j++)
    //   cerr << "(" << data[j].nodeIndex << " ; " << data[j].abruijnData << " ; " << data[j].seqsData << " ; " << data[j].disp << ") " << endl;
    // cerr << m_offsets[i].offset << endl;
  }

  // step 4: sort spectra from lower offset to greater
  sort(m_offsets.begin(), m_offsets.end());

  // step 5: get minimun offset
  double minOffset = m_offsets[0].offset;

  // step 6: set the deNovo sequence offset to be equal to the minimun offset, if negative
  if(minOffset < 0.0) {
    m_deNovoOffset = -minOffset;

    for(int i = 0 ; i < m_offsets.size() ; i++)
      m_offsets[i].offset -= minOffset;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::calcultateLimits(void)
{
  // Initialize M/Z lower limit variable
  m_mzLowerLimit = 0.0; //numeric_limits<double>::max();
  // Initialize m/z upper limit variable
  m_mzUpperLimit = 0.0;
  // Initialize intensity variable
  m_intensityLimit = 0.0;


  for(int i = 0 ; i < m_offsets.size() ; i++) {


    double offset = m_offsets[i].offset;

    int currentStar = m_offsets[i].id; //m_stars[id];

    Spectrum &sp = (*m_star)[currentStar];

    // use peak list size in an auxiliary variable for speed purposes
    int peakListSize = sp.size();

    // Cycle thru intensity values to find max m/z and intensity values
    for(int i = 0 ; i < peakListSize ; i++) {

      float intensity = sp[i][1];

      if(intensity >  m_intensityLimit)
        m_intensityLimit = intensity;
    }

    if(sp.parentMass + offset > m_mzUpperLimit)
      m_mzUpperLimit = sp.parentMass + offset;

    if(offset < m_mzLowerLimit)
      m_mzLowerLimit = offset;
    //m_mzLowerLimit = 0.0;
  }

  //cout << "Lower: " << m_mzLowerLimit << endl;
  //cout << "Upper: " << m_mzUpperLimit << endl;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::initializeGraph(void)
{
  // local lookups
  const unsigned font_height = 12;
  const unsigned font_width = font_height * 2 / 3;
  // number of contigs
  const unsigned nps = m_stars.size();
  // calculate image size
  double nzoom = m_rendererImage.zoom;

  double nx = m_rendererImage.imageSizeX;
  double ny = m_rendererImage.imageSizeY;
  // pixel width
  double npw = 1.0 / nx;
  // pixel height
  double nph = 1.0 / ny;
  // x axis size
  nxs = 1.0 - 8.0 * font_width * npw / nzoom - 3.0 * font_width * npw / nzoom;
  // y axis size
  nys = 1.0 - 12.0 * font_height * nph / nzoom - 2.0 * font_height * nph / nzoom;
  // x axis position
  nxp = 8.0 * font_width * npw / nzoom;
  // y axis position
  nyp = 2.0 * font_height * nph / nzoom;
  // top margin position
  ntp = 1.0 - 12.0 * font_height * nph / nzoom;

  // stretch graph
  if (nps * font_height / ny > nys) {

    ny = nps * font_height + nyp / nph + (1. - ntp) / nph;
    // pixel height
    nph = 1.0 / ny;
    // y axis size
    nys = 1.0 - 12.0 * font_height * nph / nzoom - 2.0 * font_height * nph / nzoom;
    // y axis position
    nyp = 2.0 * font_height * nph / nzoom;
    // top margin position
    ntp = 1.0 - 12.0 * font_height * nph / nzoom;
  }

  m_rendererImage.pixelHeight = nph;
  m_rendererImage.imageSizeY = (int)ny;

  // Calculate the margins
  float marginLeft, marginRight, marginTop, marginBottom;

  marginLeft    = 9.0;
  marginRight   = 3.0;
  marginTop     = 12.0;
  marginBottom  = 2.0;

  // Set the margins
  m_rendererObject->setMarginLeft(marginLeft);
  m_rendererObject->setMarginRight(marginRight);
  m_rendererObject->setMarginTop(marginTop);
  m_rendererObject->setMarginBottom(marginBottom);

  // set graph size in pixels
  //graphSizePixelsX = m_rendererImage.imageSizeX - (int)(marginLeft - marginRight);
  //graphSizePixelsY = m_rendererImage.imageSizeY - (int)(marginTop - marginBottom);

  // Calculate intervals
  double msdX = m_rendererObject->calculateInterval(m_mzUpperLimit - m_mzLowerLimit, 7.0, m_rendererImage.zoom);
  // Set intervals
  m_rendererObject->setXinterval(msdX);


  // Border
  m_rendererObject->setBorder(3);

  // y axis
  double one = m_intensityLimit / (double)m_stars.size();
  double half = one / 2.0;


  // m_rendererObject->setXintervalOut();
  // Y tics are spectra IDs
  vector<pair<string,string> > ytics;
  // Cycle thru all the spectra
  for(int i = 0 ; i < m_offsets.size() ; i++) {
    int loc = m_offsets.size() - 1 - i;
    ytics.push_back(make_pair<string,string>(parseInt(m_offsets[i].id + 1), parseFloat(one * loc + half, 2)));
  }
  m_rendererObject->setYinterval(ytics);

  // Graph coverage
  m_rendererObject->setXrange(m_mzLowerLimit, m_mzUpperLimit);
  m_rendererObject->setYrange(0, m_intensityLimit);
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::calcStarInContig(void)
{
  // contig in abruijn
  abContig_t &ab = (*m_abinfo)[m_contigIndex];
  m_stars = ab.first.first;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::initDeNovo(void)
{
  //const vector<pair< vector<int>, vector<double> > >  &data = (*m_abinfo)[m_contigIndex].second;

  //for(int i = 0 ; i < data.size() ; i++)
  //  for(int j = 0 ; j < data[i].first.size() ; j++)
  //    if(data[i].first[j] == m_stars[0])
  //      m_breaks.push_back(data[i].second[j]);


  Spectrum *sp = &(*m_seqs)[m_contigIndex];

  m_deNovoIntervals.clear();
  if(m_reverse) {

    for(int i = sp->size()-1 ; i >=0  ; i--) {
      double value =  sp->parentMass -  AAJumps::massHion  - (*sp)[i][0];
      m_deNovoIntervals.push_back(value);
    }

  } else {

    for(int i = 0 ; i < sp->size() ; i++)
      m_deNovoIntervals.push_back((*sp)[i][0]);

  }

  m_deNovoIntervals2 = m_deNovoIntervals;


  //if(m_deNovoIntervals2.size())
  //  return;

  //AAJumps jumps(1);
  //vector<float> aux;
  //jumps.getPRMMasses(m_denovo, aux);

  //m_deNovoIntervals2.clear();
  //m_deNovoIntervals2.push_back(0.0);
  //for(int i = 0 ; i < aux.size() ; i++)
  //  m_deNovoIntervals2.push_back( (double)aux[i] );
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::initUserIntervals(void)
{
  if(m_userIntervals.size())
    return;

  AAJumps jumps(1);
  vector<float> aux;
  jumps.getPRMMasses(m_user, aux);

  m_userIntervals.clear();
  m_userIntervals.push_back(0.0);
  for(int i = 0 ; i < aux.size() ; i++)
    m_userIntervals.push_back( (double)aux[i] );
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::initData(void)
{
  // cycly thru stars for this contig
  for(int i = 0 ; i < m_stars.size() ; i++) {
    // star index
    int star = m_stars[i];
    // get the spectrum object
    Spectrum *sp = &(*m_star)[star];
    // ---
    sp->addZPMpeaks(m_peakMassTol, 0, true, true);
    // store info
    ContigSpectrumOffsets offset;
    offset.id       = star;
    offset.offset   = 0.0;
    offset.reversed = getReverseState(star);

    //if(offset.reversed ^ m_reverse) {
    if((!( offset.reversed && m_reverse )) && ( offset.reversed || m_reverse )) {
      //cout << "reversing spectra: " << star << endl;
      reverseSpectra(sp);
    }

    m_offsets.push_back(offset);
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::reverseSpectra(Spectrum *sp)
{
  sp->reverse(0);
  /*
  sps::vector<TwoValues<float> > peakListReversed;

  for(int i = sp->size()-1 ; i >= 0 ; i--) {
    TwoValues<float> peak;
    peak[0] = sp->parentMass -  AAJumps::massHion  - (*sp)[i][0];
    peak[1] = (*sp)[i][1];
    peakListReversed.push_back(peak);
  }

  sp->peakList = peakListReversed;
  */
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::setGraphTitle(void)
{
  if(m_rendererImage.zoom < 1.0)
    return;
  // splot << "set title \"Protein Contig " << ssns[0] << "\" 0," << 3 + (cacid.size() - 1) * 2 << " font \"" << font_name << "," << font_size + 2 << "\"\n";

  // Set the title
  if(m_titlePos) {
    //stringstream aux;
    //aux << "Consensus spectrum " << m_spectrumIndexGlobal << " (Contig " << m_contig << ")";
    // Title object for plot object
    RendererTitle title;
    title.label     = m_title.c_str(); //aux.c_str();
    title.offsetX   = 0; //m_titleOffsetX;
    title.offsetY   = 3; //m_titleOffsetY;
    title.fontName  = m_rendererImage.fontName;
    title.fontSize  = int(m_rendererImage.fontSize + 2);
    // Do it
    m_rendererObject->setTitle(title);
    //m_rendererObject->setTitle(aux.str().c_str());
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::addRedArrow(double x1, double y1, double x2, double y2)
{
  drawArrow( 2, 1,
             x1, RENDERER_COORD_NONE,
             y1, RENDERER_COORD_SCREEN,
             x2, RENDERER_COORD_NONE,
             y2, RENDERER_COORD_SCREEN,
             "#A00000");
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::drawSequences(void)
{
  // distance between sequendes in the Y axis
  double disp = (double)m_rendererImage.fontHeight / (double)m_rendererImage.imageSizeY / m_rendererImage.zoom * 2;
  // start Y position (lower)
  double coord = ntp + disp;
  // start Y position for previous sequence
  double coordPrev = ntp;
  // labels are only drawn fully if zoom >= 0.8
  bool drawLabels = (m_rendererImage.zoom >= 0.8);

  // reversed spectra offset
  //double revOffset = (m_reverse ? AAJumps::massMH : 0.0);
  double revOffset = (m_reverse ? AAJumps::massH2O : 0.0);

  // offset auxiliary variable
  double offset;
  double offsetPrevious;

  // pointer to the previous sequence
  vector<double> *previous = NULL;

  // is the user sequence bent?
  bool bend = false;

  // deNovo sequence
  if(m_denovo.length()) {
    offset = m_deNovoOffset;
    offsetPrevious = 0.0;
    ContigSequence seq(m_rendererObject, &m_deNovoIntervals2, previous, &m_denovo, offset, offsetPrevious, "de Novo:", coord, coordPrev, drawLabels, m_mzLowerLimit, m_mzUpperLimit, "#A00000");
    seq.setPeakMassTol(m_peakMassTol);
    seq.setZoom(m_rendererImage.zoom);
    seq.setImageDimensions(m_rendererImage.imageSizeX, m_rendererImage.imageSizeY);
    seq.setPixelSize(m_rendererImage.pixelWidth, m_rendererImage.pixelHeight);
    seq.setFontHeight(m_rendererImage.fontHeight);
    seq.drawExec();
    previous = &m_deNovoIntervals2;
    coord += disp;
    coordPrev += disp;
    offsetPrevious = offset;
  }

  // homolog sequence
  if(m_homolog.length()) {
    offset = m_deNovoOffset + m_homologOffset + revOffset;
    ContigSequence seq(m_rendererObject, &m_homologIntervals, previous, &m_homolog, offset, offsetPrevious, "Homolog:", coord, coordPrev, drawLabels, m_mzLowerLimit, m_mzUpperLimit, "#A00000");
    seq.setPeakMassTol(m_peakMassTol);
    seq.setZoom(m_rendererImage.zoom);
    seq.setImageDimensions(m_rendererImage.imageSizeX, m_rendererImage.imageSizeY);
    seq.setPixelSize(m_rendererImage.pixelWidth, m_rendererImage.pixelHeight);
    seq.setFontHeight(m_rendererImage.fontHeight);
    seq.drawExec();
    previous = &m_homologIntervals;
    coord += disp;
    coordPrev += disp;
    bend = true;
    offsetPrevious = offset;
  }

  //cout << m_reference << endl;
  //for(int i = 0 ; i < m_referenceIntervals.size() ; i++)
  //  cout << m_referenceIntervals[i] << " ; " << endl;
  //cout << m_referenceOffset << endl;

  // reference sequence
  if(m_reference.length()) {
    offset = m_deNovoOffset + m_referenceOffset + revOffset;
    ContigSequence seq(m_rendererObject, &m_referenceIntervals, previous, &m_reference, offset, offsetPrevious, "Reference:", coord, coordPrev, drawLabels, m_mzLowerLimit, m_mzUpperLimit, "#A00000");
    seq.setPeakMassTol(m_peakMassTol);
    seq.setZoom(m_rendererImage.zoom);
    seq.setImageDimensions(m_rendererImage.imageSizeX, m_rendererImage.imageSizeY);
    seq.setPixelSize(m_rendererImage.pixelWidth, m_rendererImage.pixelHeight);
    seq.setFontHeight(m_rendererImage.fontHeight);
    seq.drawExec();
    previous = &m_referenceIntervals;
    coord += disp;
    coordPrev += disp;
    bend = true;
    offsetPrevious = offset;
  }

  // user sequence (on top)
  if(m_user.length()) {
    offset = m_deNovoOffset + m_userOffset;
    previous = &m_deNovoIntervals2;
    coordPrev = ntp + disp;
    ContigSequence seq(m_rendererObject, &m_userIntervals, previous, &m_user, offset, offsetPrevious, "User:", coord, coordPrev, drawLabels, m_mzLowerLimit, m_mzUpperLimit, "#4040F0", bend);
    seq.setPeakMassTol(m_peakMassTol);
    seq.setZoom(m_rendererImage.zoom);
    seq.setImageDimensions(m_rendererImage.imageSizeX, m_rendererImage.imageSizeY);
    seq.setPixelSize(m_rendererImage.pixelWidth, m_rendererImage.pixelHeight);
    seq.setFontHeight(m_rendererImage.fontHeight);
    seq.drawExec();
    coord += disp;
    coordPrev += disp;
  }

}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
bool PlotContig::checkContigItem(int spectrum, double mass)
{
  const vector<pair< vector<int>, vector<double> > >  &data = (*m_abinfo)[m_contigIndex].second;

  for(int i = 0 ; i < data.size() ; i++)
    for(int j = 0 ; j < data[i].first.size() ; j++)
      if(data[i].first[j] == spectrum)
        if(fabs(data[i].second[j] - mass) <= m_peakMassTol)
          return true;

  return false;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
bool PlotContig::checkContigItemAtPosition(int spectrum, int position, float &value)
{
  // get the contig data section
  const vector<pair< vector<int>, vector<double> > >  &data = (*m_abinfo)[m_contigIndex].second;

  if(data.size() <= position)
    return false;

  // get the node IDs list
  const vector<int> &nodeIds = data[position].first;

  // go thru the
  for(int j = 0 ; j < nodeIds.size() ; j++)
    if(nodeIds[j] == spectrum) {
      value = data[position].second[j];
      return true;
    }

  return false;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
bool PlotContig::getReverseState(int spectrum)
{
  //cout << " ----------------- " << spectrum << " ----------------- " << endl;

  // get the contig data section
  const pair< vector<int>, vector<int> >  &data = (*m_abinfo)[m_contigIndex].first;
  // get the node IDs list
  const vector<int> &nodeIds = data.first;
  const vector<int> &nodeRev = data.second;

  for(int i = 0 ; i < nodeIds.size() ; i++) {
    //cout << nodeIds[i] << " : " << nodeRev[i] << endl;
    if(nodeIds[i] == spectrum)
      return nodeRev[i];
  }
  return false;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::drawVerticalLines(void)
{
  // this vector stores the topmost positions. Lines should be drawn from this points down, and updated afterwards
  // first pair is x position, second is y positiob
  vector<pair<float, float> > topPoints;

  double disp =  (double)m_rendererImage.fontHeight / (double)m_rendererImage.imageSizeY / m_rendererImage.zoom * 2.0;
  double disp2 = (double)m_rendererImage.fontHeight / (double)m_rendererImage.imageSizeY / m_rendererImage.zoom * 1.0;
  double disp3 = (m_mzUpperLimit - m_mzLowerLimit) * 0.01;


  // fill the vector and draw the top dialgonal lines
  for(int i = 0 ; i < m_deNovoIntervals.size() ; i++) {

    double x1 = m_deNovoOffset + m_deNovoIntervals[i];
    double y1 = ntp + disp;
    double x2 = x1  + disp3;
    double y2 = ntp + disp - disp2;

    // draw the small top line
    drawLine(  1, 1,
               x1, RENDERER_COORD_NONE,
               y1, RENDERER_COORD_SCREEN,
               x2, RENDERER_COORD_NONE,
               y2, RENDERER_COORD_SCREEN,
               "#A00000");

    // store the topmost position as the next starting point
    topPoints.push_back(make_pair<float,float>(x2, y2));
  }

  // go thru all the star spectra and draw the vertical lines accordingly
  for(int i = 0 ; i < m_offsets.size() ; i++ ) {
  //for(int i = m_offsets.size()-1 ; i >=0  ; i-- ) {

    int id = m_offsets[i].id;

    int position = m_offsets.size() - i - 1;

    //go thru all the points of the star spectra
    for(int j = 0 ; j < m_deNovoIntervals.size() ; j++ ) {

      float value;

      if(checkContigItemAtPosition(id, j, value)) {

        double x1 = topPoints[j].first;
        double y1 = topPoints[j].second;

        double x2 = m_offsets[i].offset + value + disp3;
        double y2 = 1.0 / m_stars.size() * position * nys + nyp  + 0.01;

        // draw the vertical dashed line
        drawLine(  1, 0,
                   x1, RENDERER_COORD_NONE,
                   y1, RENDERER_COORD_SCREEN,
                   x2, RENDERER_COORD_NONE,
                   y2, RENDERER_COORD_SCREEN,
                   "#A00000");

        topPoints[j].first  = x2;
        topPoints[j].second = y2;

        x1 = x2 - disp3;
        y1 = y2 = 1.0 / m_stars.size() * position * nys + nyp  + 0.01;

        // draw small horizontal line near the
        drawLine(  1, 1,
                   x1, RENDERER_COORD_NONE,
                   y1, RENDERER_COORD_SCREEN,
                   x2, RENDERER_COORD_NONE,
                   y2, RENDERER_COORD_SCREEN,
                   "#A00000");

        //cerr << x1 << " ; " << x2 << " ; ";
      }
    }
    //cerr << endl;
  }
}
///////////////////////////////////////////////////////////////////////////////
void PlotContig::drawSpectra(void)
{
  for(int i = 0 ; i < m_offsets.size(); i++) {
    drawSpectrum(m_offsets[i].id, i);
    m_rendererObject->addArea();
  }
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
void PlotContig::drawSpectrum(int spectrumIndex, int spectrumCount)
{
  //cerr << " ------------ " << spectrumIndex << " ---------" << endl;

  // Set the margins
  m_rendererObject->setMarginLeft(0);
  m_rendererObject->setMarginRight(0);
  m_rendererObject->setMarginTop(0);
  m_rendererObject->setMarginBottom(0);
  // border just on the lower side
  m_rendererObject->setBorder(1);

  // get the spectrum object
  Spectrum *sp = &(*m_star)[spectrumIndex];

  //cerr << "size draw " << spectrumIndex << ":  " << sp->size() << endl;

  // number of peaks
  int size = sp->size();
  // number of items (used to calculate Y position)
  double nps = m_stars.size();
  // offset for this spectrum
  double offset = m_offsets[spectrumCount].offset;
  // position formula denominator
  double denom = (m_mzUpperLimit - m_mzLowerLimit) / nxs;

  double mzUpperLimit = sp->parentMass;

  //cerr << mzUpperLimit << " >> ";
  //for(int i = 0 ; i < sp->size() ; i++)
  //  cerr << (*sp)[i][0] << " ; ";
  //cerr << endl;

  int position = m_offsets.size() - spectrumCount - 1;

  // calculate positions and sizes
  double ox = ((*sp)[0][0] - m_mzLowerLimit + offset) / denom + nxp;
  //double ox = ((*sp)[0][0] - m_mzLowerLimit) / denom + nxp;
  double oy = 1.0 / nps * position * nys + nyp;
  //double sx = abs(m_mzUpperLimit - (*sp)[0][0] - offset) / denom;
  //double sx = abs(sp->parentMass * 1.05  - 0.0) / denom;
  //double sx = abs(m_mzUpperLimit - (*sp)[0][0]) / denom;
  double sx = fabs(mzUpperLimit - (*sp)[0][0]) / denom;
  double sy = (1.0 / nps) * nys;


  double disp3 = (m_mzUpperLimit - m_mzLowerLimit) * 0.005;
  double factor = (m_mzUpperLimit + disp3) / mzUpperLimit;
  //cout << "denom: " << denom << endl;
  //cout << "offset: " << offset << endl;
  //cout << "ox: " << ox << endl;
  //cout << "sx: " << sx << endl;
  //cout << "" << endl;

  // set origin
  m_rendererObject->setOrigin(ox, oy);
  // set size
  m_rendererObject->setSize(sx, sy);

  // draw the spectrum
  ContigSpectrum  cs(m_rendererObject, m_contigIndex, m_abinfo, sp, &m_deNovoIntervals, spectrumIndex);
  cs.setPeakMassTol(m_peakMassTol);
  cs.draw(m_intensityLimit, factor);
}
////////////////////////////////////////////////////////////////////////////////
// Dump the abruijn graph to the screen
void PlotContig::dump_abruijn(ostream &sout, abinfo_t *abruijn)
{
  if(!abruijn) return;

  sout << "----------------- abruijn graph----------------" << endl;

  // DUMP abruijn
  abinfo_t::iterator i0 = abruijn->begin();
  for(; i0 != abruijn->end() ; i0++) {
    sout << "Contig " << i0->first << endl;

    // output first vector in pair
    for(int i = 0 ; i < i0->second.first.first.size() ; i++) {
      sout << "     spectrum: " << i0->second.first.first[i] << "  --  flipped: " << i0->second.first.second[i] << endl;
    }

    sout << "     ------ " << endl;

    // output second pair -> vector of pairs of vectors
    // output first vector in pair
    std::vector<std::pair< vector<int>, vector<double> > >::iterator i1 = i0->second.second.begin();
    for(; i1 != i0->second.second.end() ; i1++) {
      sout << "     ------ " << endl;

      for(int i = 0 ; i < i1->first.size() ; i++) {
        sout << "     spectrum: " << i1->first[i] << "  --  peak mass: " << i1->second[i] << endl;
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
