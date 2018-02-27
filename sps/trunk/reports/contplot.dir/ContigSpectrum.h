///////////////////////////////////////////////////////////////////////////////
#ifndef __CONTIG_SPECTRUM_H__
#define __CONTIG_SPECTRUM_H__
///////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "RendererBase.h"
#include "ReportDefines.h"
#include "PlotBase.h"
#include "abruijn.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
//
struct abMassItem {
  double  mass;
  int     idxInContig;
};
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
	/*! \brief Class for displaying spectra in contig images

   Draws a spectrum in a contig and the arrows between peaks.

   */
class ContigSpectrum : public PlotBase {

   // spectrum to draw
	/*! \brief spectrum to draw
   */
  specnets::Spectrum *m_spectrum;

  // ABinfo data set
  //specnets::abinfo_t *m_abinfo;
	/*! \brief ABinfo data set
   */
  abContigData_t  *abData;

  // vector if masses from the abruijn for the star
	/*! \brief vector if masses from the abruijn for the star
   */
  vector<abMassItem> m_abMasses;

  // sps_seqs data, properly reversed, if the case
	/*! \brief sps_seqs data, properly reversed, if the case
   */
  vector<double> *m_deNovoIntervals;

  // Contig to draw
	/*! \brief Contig to draw
   */
  int m_contigIndex;

  // spectrum index
	/*! \brief spectrum index
   */
  int m_spectrumIndex;


  // build mass value list from abruijn graph masses
	/*! \brief build mass value list from abruijn graph masses
   */
  void buildVector(void);
  // find mass value in mass value list
	/*! \brief find mass value in mass value list
   */
  int find(double &value);

  // red horizontal line
	/*! \brief red horizontal line
   */
  void addRedArrow(double x1, double y1, double x2, double y2);
  // green horizontal line
	/*! \brief green horizontal line
   */
  void addGreenArrow(double x1, double y1, double x2, double y2);
  // add the red label with the mass shift value
	/*! \brief add the red label with the mass shift value
   */
  void addRedLabel(double x1, double y1, string &label);



 public:

	/*! \brief Constructor
   */
  ContigSpectrum(RendererBase *renderer, int ci, abContigData_t *a, specnets::Spectrum *s, vector<double> *seqs, int idx) :
    m_contigIndex(ci), m_spectrum(s), abData(a), m_spectrumIndex(idx), m_deNovoIntervals(seqs)
    {
      m_rendererObject              = renderer;
      m_ownRenderer                 = false;
      m_rendererImage.lineThickness = 4;
    };

	/*! \brief Destructor
   */
  ~ContigSpectrum() {};

	/*! \brief Draws the image
   */
  virtual void draw(double limit, double factor = 1.0);

	/*! \brief Main draw method entry point
   */
  virtual int  drawExec(void) {};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
