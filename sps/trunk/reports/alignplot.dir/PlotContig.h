///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_CONTIG_H__
#define __PLOT_CONTIG_H__
///////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "abruijn.h"
#include "ContigSpectrum.h"
#include "PlotBase.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
// Auxiliary classes and methods

/*---------------------------------------------------------------------------*/
// spectrum info for spectra representation in a contig
struct AlignmentData {
  // abruijn node index
  int   nodeIndex;
  // value in abruijn node
  float abruijnData;
  // value in seqs
  float seqsData;
  // number of diagonals if this node is aligned
  int   disp;
};
/*---------------------------------------------------------------------------*/
struct ContigSpectrumOffsets {

  // star id
  int     id;
  // mass offset for the spectra
  double  offset;
  // true if spectra needs to be reversed
  bool    reversed;


  // used for sorting
  bool operator<(const ContigSpectrumOffsets &o) const {return offset < o.offset;};

};
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
class PlotContig : public PlotBase {

  // ABinfo data set
  abinfo_t    *m_abinfo;

  // Data spectra loaded from file
  SpecSet     *m_star;

  // Data spectra loaded from file
  SpecSet     *m_seqs;

  // Contig to draw
  int         m_contigIndex;

  double      m_mzLowerLimit;
  double      m_mzUpperLimit;
  double      m_intensityLimit;
  double      m_range_max;
  double      m_range_min;

  // contig reverse flag
  bool        m_reverse;

  // offset structure
  vector<ContigSpectrumOffsets> m_offsets;


  // peptide
  string m_reference;
  string m_homolog;
  string m_denovo;
  string m_user;


  // star ids for the contig
  vector<int>  m_stars;

  // top draw positions
  vector<double> m_referenceIntervals;

  vector<double> m_homologIntervals;

  vector<double> m_deNovoIntervals, m_deNovoIntervals2;

  vector<double> m_userIntervals;


  // sequence offsets
  double m_referenceOffset, m_homologOffset, m_deNovoOffset, m_userOffset;


  // graph data
  double nxs;
  double nys;
  double nxp;
  double nyp;
  double ntp;
  double np;
  double nm[4];


  // test if a mass value is within boundaries
  bool testMass(float mass);

  //////////////////////////////////////////////////////////////////////////////
  // Upper level methods

  // contig drawing method
  void plotContig(void);

  //////////////////////////////////////////////////////////////////////////////
  // Data initializers

  // calc offsets
  void initData(void);

  // set mass intervals
  void setMassIntervals(string &strValues, vector<double> &values);

  // Set deNovo masses in a vector for easyer access
  void initDeNovo(void);

  // calculate user mass intervals based on sequence
  void initUserIntervals(void);

  // gets star ids based on contig id
  void calcStarInContig(void);

  // Calculate graph limits
  void calcultateLimits(void);

  // calculate areas sizes and origins
  void reverseSpectra(Spectrum *sp);

  // Align spectra, by calculating a start position
  void alignSpectra(void);


  //////////////////////////////////////////////////////////////////////////////
  // Draw object initializers

  // Initialize graph variables and y axis labels
  void initializeGraph(void);

  // Sets the graph title
  void setGraphTitle(void);


  //////////////////////////////////////////////////////////////////////////////
  // Sequence processing & drawing

  // draws sequences
  void drawSequences(void);

  // draw a red arrow
  void addRedArrow(double x1, double y1, double x2, double y2);


  //////////////////////////////////////////////////////////////////////////////
  // Vertical lines drawing

  // draw vertical dashed lines
  void drawVerticalLines(void);


  //////////////////////////////////////////////////////////////////////////////
  // Spectrum drawing

  // Draws the spectra graphs
  void drawSpectra(void);

  // Draws a single spectrum
  void drawSpectrum(int spectrumIndex, int spectrumCount);


  //////////////////////////////////////////////////////////////////////////////
  // Contig access

  // find a mass item in a contig star
  bool checkContigItem(int spectrum, double mass);

  // check spectrum reverse flag
  bool getReverseState(int spectrum);

  // check for a star in a contig node given it's index
  bool checkContigItemAtPosition(int spectrum, int position, float &value);

  // reverse the abruijn
  void getReversedAbruijn(void);



 public:


  // Constructors and Destructor
  PlotContig();
  ~PlotContig();


  // Initialization methods
  virtual void setStar(SpecSet *s)                { m_star              = s;};
  virtual void setSeqs(SpecSet *s)                { m_seqs              = s;};
  virtual void setAbinfo(abinfo_t *a)             { m_abinfo            = a;};
  virtual void setContigIndex(int i)              { m_contigIndex       = i;};

  // peptide
  virtual void setReference(string &s)            {m_reference = s;};
  virtual void setHomolog(string &s)              {m_homolog = s;};
  virtual void setDenovo(string &s)               {m_denovo = s;};
  virtual void setUser(string &s)                 {m_user = s;};

  virtual void setReferenceMass(string &s)        {setMassIntervals(s, m_referenceIntervals);};
  virtual void setHomologMass(string &s)          {setMassIntervals(s, m_homologIntervals);};
  virtual void setUserMass(string &s)             {setMassIntervals(s, m_userIntervals);};

  virtual void setReferenceOffset(double o)       {m_referenceOffset =  o;};
  virtual void setHomologOffset(double o)         {m_homologOffset =    o;};
  virtual void setDenovoOffset(double o)          {m_deNovoOffset =     o;};
  virtual void setUserOffset(double o)            {m_userOffset =       o;};

  virtual void setReverseFlag(void)               {m_reverse = true;};


  // Default draw object entry point
  virtual void drawExec(void);

  void dump_abruijn(ostream &sout, abinfo_t *abruijn);


};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
