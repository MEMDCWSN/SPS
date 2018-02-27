///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_CONTIG_H__
#define __PLOT_CONTIG_H__
///////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "abruijn.h"
#include "ContigSequence.h"
#include "ContigSpectrum.h"
#include "PlotBase.h"
#include "aminoacid.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
// Auxiliary classes and methods

/*---------------------------------------------------------------------------*/
// spectrum info for spectra representation in a contig
  /*! \brief  spectrum info for spectra representation in a contig
   */
struct AlignmentData {
  // abruijn node index
  /*! \brief abruijn node index
   */
  int   nodeIndex;
  // value in abruijn node
  /*! \brief value in abruijn node
   */
  float abruijnData;
  // value in seqs
  /*! \brief value in seqs
   */
  float seqsData;
  // number of diagonals if this node is aligned
  /*! \brief number of diagonals if this node is aligned
   */
  int   disp;
};
/*---------------------------------------------------------------------------*/
  /*! \brief Contig spectrum offsets
   */
struct ContigSpectrumOffsets {

  // star id
  /*! \brief star id
   */
  int     id;
  // mass offset for the spectra
  /*! \brief mass offset for the spectra
   */
  double  offset;
  // parent mass + offset. The m/z mass limit for this spectra
  /*! \brief parent mass + offset. The m/z mass limit for this spectra
   */
  double  endingMass;
  // true if spectra needs to be reversed
  /*! \brief  true if spectra needs to be reversed
   */
  bool    reversed;
  // peak tolerance
  /*! \brief peak tolerance
   */
  double  tolerance;


  // used for sorting
  /*! \brief used for sorting
   */
  bool operator<(const ContigSpectrumOffsets &o) const {return (fabs(offset-o.offset)>tolerance ? offset < o.offset : endingMass < o.endingMass);};

};
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
  /*! \brief Contig plot data structure
   */
class PlotContig : public PlotBase {

  // ABinfo data set
  /*! \brief ABinfo data set
   */
  abinfo_t    *m_abinfo;

  // Data spectra loaded from file
  /*! \brief Data spectra loaded from file
   */
  SpecSet     *m_star;

  // Data spectra loaded from file
  /*! \brief Data spectra loaded from file
   */
  SpecSet     *m_seqs;

  // abruijn node data
  /*! \brief abruijn node data
   */
  abContigData_t  abData;

  // Contig to draw
  /*! \brief Contig to draw
   */
  int         m_contigIndex;

  double      m_mzLowerLimit;
  double      m_mzUpperLimit;
  double      m_intensityLimit;
  double      m_range_max;
  double      m_range_min;

  // contig reverse flag
  /*! \brief contig reverse flag
   */
  bool        m_reverse;

  // offset structure
  /*! \brief offset structure
   */
  vector<ContigSpectrumOffsets> m_offsets;


  // peptide
  /*! \brief peptide
   */
  string m_reference;
  /*! \brief peptide
   */
  string m_homolog;
  /*! \brief peptide
   */
  string m_denovo;
  /*! \brief peptide
   */
  string m_user;


  // star ids for the contig
  /*! \brief star ids for the contig
   */
  vector<int>  m_stars;

  // top draw positions
  /*! \brief top draw positions
   */
  vector<double> m_referenceIntervals;

  /*! \brief top draw positions
   */
  vector<double> m_homologIntervals;

  /*! \brief top draw positions
   */
  vector<double> m_deNovoIntervals, m_deNovoIntervals2;

  /*! \brief top draw positions
   */
  vector<double> m_userIntervals;


  // sequence offsets
  /*! \brief sequence offsets
   */
  double m_referenceOffset, m_homologOffset, m_deNovoOffset, m_userOffset;

  // sequence draw style
  /*! \brief sequence draw style
   */
  bool m_allDots;


  // graph data
  double nxs;
  double nys;
  double nxp;
  double nyp;
  double ntp;
  double np;
  double nm[4];

  // AAs
  /*! \brief AAs
   */
  AAJumps *jumps;

  // test if a mass value is within boundaries
  /*! \brief test if a mass value is within boundaries
   */
  bool testMass(float mass);

  //////////////////////////////////////////////////////////////////////////////
  // Upper level methods

  // contig drawing method
  /*! \brief contig drawing method
   */
  void plotContig(void);

  //////////////////////////////////////////////////////////////////////////////
  // Data initializers

  // calc offsets
  /*! \brief calc offsets
   */
  void initData(void);

  // set mass intervals
  /*! \brief set mass intervals
   */
  void setMassIntervals(string &strValues, vector<double> &values);

  // Set deNovo masses in a vector for easyer access
  /*! \brief Set deNovo masses in a vector for easyer access
   */
  void initDeNovo(void);

  // calculate user mass intervals based on sequence
  /*! \brief calculate user mass intervals based on sequence
   */
  void initMassIntervals(vector<double> &massIntervals, string &masses);

  // gets star ids based on contig id
  /*! \brief gets star ids based on contig id
   */
  void calcStarInContig(void);

  // Calculate graph limits
  /*! \brief Calculate graph limits
   */
  void calcultateLimits(void);

  // calculate areas sizes and origins
  /*! \brief calculate areas sizes and origins
   */
  void reverseSpectra(Spectrum *sp);

  // Align spectra, by calculating a start position
  /*! \brief Align spectra, by calculating a start position
   */
  void alignSpectra(void);


  //////////////////////////////////////////////////////////////////////////////
  // Draw object initializers

  // Initialize graph variables and y axis labels
  /*! \brief Initialize graph variables and y axis labels
   */
  void initializeGraph(void);

  // Sets the graph title
  /*! \brief Sets the graph title
   */
  void setGraphTitle(void);


  //////////////////////////////////////////////////////////////////////////////
  // Sequence processing & drawing

  // draws sequences
  /*! \brief draws sequences
   */
  void drawSequences(void);

  // draw a red arrow
  /*! \brief draw a red arrow
   */
  void addRedArrow(double x1, double y1, double x2, double y2);


  //////////////////////////////////////////////////////////////////////////////
  // Vertical lines drawing

  // draw vertical dashed lines
  /*! \brief draw vertical dashed lines
   */
  void drawVerticalLines(void);


  //////////////////////////////////////////////////////////////////////////////
  // Spectrum drawing

  // Draws the spectra graphs
  /*! \brief Draws the spectra graphs
   */
  void drawSpectra(void);

  // Draws a single spectrum
  /*! \brief Draws a single spectrum
   */
  void drawSpectrum(int spectrumIndex, int spectrumCount);


  //////////////////////////////////////////////////////////////////////////////
  // Contig access

  // find a mass item in a contig star
  /*! \brief find a mass item in a contig star
   */
  bool checkContigItem(int spectrum, double mass);

  // check spectrum reverse flag
  /*! \brief check spectrum reverse flag
   */
  bool getReverseState(int spectrum);

  // check for a star in a contig node given it's index
  /*! \brief check for a star in a contig node given it's index
   */
  bool checkContigItemAtPosition(int spectrum, int position, float &value);

  // get reverse abruijn data for node
  /*! \brief get reverse abruijn data for node
   */
  void getReversedAbruijn(void);

  // get abruijn data for node
  /*! \brief get abruijn data for node
   */
  void getDirectAbruijn(void);



 public:


  // Constructors and Destructor
  /*! \brief Constructors
   */
  PlotContig();
  /*! \brief Destructor
   */
  ~PlotContig();


  // Initialization methods
  /*! \brief Initialization methods
   */
  virtual void setStar(SpecSet *s)                { m_star              = s;};
  /*! \brief Initialization methods
   */
  virtual void setSeqs(SpecSet *s)                { m_seqs              = s;};
  /*! \brief Initialization methods
   */
  virtual void setAbinfo(abinfo_t *a)             { m_abinfo            = a;};
  /*! \brief Initialization methods
   */
  virtual void setContigIndex(int i)              { m_contigIndex       = i;};

  // peptide
  /*! \brief peptide
   */
  virtual void setReference(string &s)            {m_reference = s;};
  /*! \brief peptide
   */
  virtual void setHomolog(string &s)              {m_homolog = s;};
  /*! \brief peptide
   */
  virtual void setDenovo(string &s)               {m_denovo = s;};
  /*! \brief peptide
   */
  virtual void setUser(string &s)                 {m_user = s;};

  /*! \brief Peptide masses
   */
  virtual void setReferenceMass(string &s)        {setMassIntervals(s, m_referenceIntervals);};
  /*! \brief Peptide masses
   */
  virtual void setHomologMass(string &s)          {setMassIntervals(s, m_homologIntervals);};
  /*! \brief Peptide masses
   */
  virtual void setUserMass(string &s)             {setMassIntervals(s, m_userIntervals);};

  /*! \brief Peptide mass offset
   */
  virtual void setReferenceOffset(double o)       {m_referenceOffset =  o;};
  /*! \brief Peptide mass offset
   */
  virtual void setHomologOffset(double o)         {m_homologOffset =    o;};
  /*! \brief Peptide mass offset
   */
  virtual void setDenovoOffset(double o)          {m_deNovoOffset =     o;};
  /*! \brief Peptide mass offset
   */
  virtual void setUserOffset(double o)            {m_userOffset =       o;};

  /*! \brief Sets the reverse flag
   */
  virtual void setReverseFlag(void)               {m_reverse = true;};


  // Default draw object entry point
  /*! \brief Default draw object entry point
   */
  virtual int  drawExec(void);

  /*! \brief Dump the abruijn graph (debug)
   */
  void dump_abruijn(ostream &sout, abinfo_t *abruijn);


};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
