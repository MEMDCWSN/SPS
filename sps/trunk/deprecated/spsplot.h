#ifndef SPSPLOT_H
#define SPSPLOT_H

#include <time.h>

#include <map>
#include <list>
#include <vector>
#include <string>
#include <utility>
#include <sstream>

#include "ion.h"
#include "batch.h"
#include "tuple.h"
#include "vector.h"
#include "abruijn.h"
#include "spectrum.h"
#include "db_fasta.h"
#include "aminoacid.h"
#include "uncertain.h"
#include "SpecSet.h"

#if defined(__linux__)
#include "shared.h"
#endif

#define STR(s) #s
#define XSTR(s) STR(s)


using namespace std;
using namespace __gnu_cxx;

namespace specnets {

typedef uncertain<double> uf;


/**
    Shotgun Protein Sequencing Plot.
*/

struct SpsPlot
{
  friend class SpecPlot;
  friend class ContPlot;

  // in
  enum bopt_t { eannot, everbose, eqtok, ektoq, eshift, estats, enotitle, ebopts };
  enum sopt_t
  {
    eproject, estars, ecomp, eseqs, epklbin, eindex, emp, emidx, efasta, erefindex, erefmp, erefmidx,
    efont, eoutdir, econtig, enodes, ep, espectrum, espectrumscan, eformat, ezoom, emgf, espectruminfo,
    eprefix, eaa, epeakmasstol, eparentmasstol, epkl, emzxml, eoutfile, elabel,
    einspect, eprm, estarsindex, epairs, ematchpairs, estarspairs, epartialoverlaps,
    erange, htmlDefs, ePageTitle, esopts
  };
  enum copt_t { ecluster, eclusterms, eclusterscan, epeptide, ecopts };

  bool bopt[ebopts];
  string sopt[esopts];
  sps::vector<string> copt[ecopts];

  static std::map<bopt_t, const char *> sbopt;
  static std::map<sopt_t, const char *> ssopt;
  static std::map<copt_t, const char *> scopt;
    
  void genHtmlHeader(ostream & shtml);
  void genHtmlHeader2(ostream & shtml);


  // dynamic call
  typedef pair< int (SpsPlot::*)(const sps::vector<string> &), sps::vector<string> > queue_element;
  typedef list<queue_element> queue_t;

  // info
  string spath;

#if defined(__linux__)
  SpsPlot() : cdata(sps::shared<Data>())
#else
  SpsPlot()
#endif
  {
    fill(bopt, bopt + sizeof(bopt) / sizeof(* bopt), false);

    sopt[eformat] = "png";
    sopt[ezoom] = "1";
#if defined(__linux__)
    sopt[efont] = "/usr/share/fonts/bitstream-vera";
#else
    sopt[efont] = "/usr/share/fonts/TTF/";
#endif
    sopt[eoutdir] = "html";
    sopt[eprefix] = "spectrum";
    sopt[epeakmasstol] = ".500001";
    sopt[eparentmasstol] = "1.000001";

    sopt[epartialoverlaps] = "1";

    cdata.ccolor.push_back("gray");
    cdata.ccolor.push_back("black");
  }

  int operator () ();

  int genheader(ostream &);
  void genFooter(ostream &);
  void genTableFooter(ostream &);

#if ! defined(SPECPLOT)
  int genhtml();
  int genindex(const sps::vector<string> &);
  int gencontig(unsigned nc);
  int stats();
#endif // SPECPLOT

  static int help(ostream & sout = cout);
  static int version(ostream & sout = cout);
  static int error(const string & a);
  
  
  void dumpAbruijn(char *);
  
  

private:
  struct Data;
  typedef map< const type_info *, pair<SpsPlot *, queue_element> > ccron_t;

  time_t ntime[2];
  static ccron_t ccron;

  static void sigext(int nsig = 0);
  static void sigimp(int nsig = 0);
  static void sigexp(int nsig = 0);

  struct Data
  {
    typedef sps::vector< sps::vector< sps::vector<unsigned> > > cclusterscan_t;
    typedef map< unsigned, list< pair<unsigned, unsigned> > > ccluster_t;
    typedef sps::vector< sps::vector<int> > cproteinidx_t;
    typedef sps::vector<SpecSet> cspecmap_t;
    typedef map<unsigned, string> cclusteridx_t;
    typedef map<unsigned, unsigned> crefindex_t;
    typedef map<unsigned, map<unsigned, pair<string, string> > > cpeptide_t;
    typedef map<unsigned, map<unsigned, pair<unsigned, string> > > csequence_t;
    typedef map< unsigned, pair<string, double[29]> > cstat_t;
    typedef map< unsigned, sps::vector<bool> > cmark_t;
    typedef map< unsigned, sps::tuple<string, string, unsigned> > cinspect_t;
    typedef hash_map< uf, hash_map< uf, sps::tuple<string, string, unsigned, unsigned, string> > > clabel_t;
    typedef list<string> ccolor_t;

    mutable SpecSet cspecset, cconsset;
    SpecSet cprotset, crefset, cmsset;
    cclusterscan_t cclusterscan;
    ccluster_t ccluster;
    cclusteridx_t cclusteridx;
    cspecmap_t cspecmap;
    mutable cproteinidx_t cproteinidx;
    cproteinidx_t creferenceidx;
    crefindex_t crefindex;
    mutable abinfo_t cabinfo;
    DB_fasta cfasta;
    hash_map<uf, char> caa;
    map<char, double> ciaa;
    cinspect_t cinspect;
    clabel_t clabel;
    ccolor_t ccolor;

    SpecSet cprm, cmatchpairs;
    sps::vector< sps::vector<unsigned> > cstarsindex;
    sps::vector<Results_PA> cpairs, cstarspairs;

    bool reference(unsigned nc) const
    {
      return crefindex.find(nc) != crefindex.end() && cproteinidx[nc][0] != creferenceidx[crefindex.find(nc)->second][0];
    }
  };

#if defined(__linux__)
  Data & cdata;
#else
  Data cdata;
#endif
};


/**
    Spectrum Plot.
*/

struct SpecPlot
{
  // in
  enum { erecurse, esuperpose, eannot, everbose, eshift, enoheader, enotitle, ebounds, ebopts };
  enum { epeakmasstol, econtig, espectrum, ezoom, eformat, eoutdir, eprefix, ehtml, efollow, eindex, efasta, efont, escan, eoutfile, erange, etitle, esopts };
  enum { ecluster, eclusterms, eclusterscan, ecopts };

  bool bopt[ebopts];
  string sopt[esopts];
  sps::vector<string> copt[ecopts];

  const SpsPlot::Data & cdata;
  SpecSet & cspecset;
  const SpsPlot::Data::cpeptide_t & cpeptide;

  // out
  ostream & splot, & shtml;

  SpsPlot::Data::cstat_t & cstat;
  SpsPlot::Data::cmark_t & cmark;

  SpsPlot::Data::cstat_t cstat_null;
  SpsPlot::Data::cmark_t cmark_null;

  void init()
  {
    fill(bopt, bopt + sizeof(bopt) / sizeof(* bopt), false);

    sopt[epeakmasstol] = ".500001";
    sopt[eoutdir] = ".";
    sopt[eprefix] = "spectrum";
  }

  SpecPlot(SpsPlot::Data & cdata, SpecSet & cspecset, ostream & splot, ostream & shtml, SpsPlot::Data::cpeptide_t & cpeptide)
  :
    cdata(cdata), cspecset(cspecset), splot(splot), shtml(shtml), cpeptide(cpeptide), cstat(cstat_null), cmark(cmark_null)
  {
    init();
  }

  SpecPlot(SpsPlot::Data & cdata, SpecSet & cspecset, ostream & splot, ostream & shtml, SpsPlot::Data::cpeptide_t & cpeptide, SpsPlot::Data::cstat_t & cstat)
  :
    cdata(cdata), cspecset(cspecset), splot(splot), shtml(shtml), cpeptide(cpeptide), cstat(cstat), cmark(cmark_null)
  {
    init();
  }

  SpecPlot(SpsPlot::Data & cdata, SpecSet & cspecset, ostream & splot, ostream & shtml, SpsPlot::Data::cpeptide_t & cpeptide, SpsPlot::Data::cmark_t & cmark)
  :
    cdata(cdata), cspecset(cspecset), splot(splot), shtml(shtml), cpeptide(cpeptide), cstat(cstat_null), cmark(cmark)
  {
    init();
  }

  int operator () ();

  friend ostream & operator << (ostream & sout, const SpecPlot & cplot)
  {
    return sout << "contig " << dec << cplot.sopt[econtig] << ", spectrum " << dec << cplot.sopt[espectrum] << hex;
  }
};


#if ! defined(SPECPLOT)
/**
    Contig Plot.
*/

struct ContPlot
{
  // in
  enum { erecurse, ebopts };
  enum { epeakmasstol, econtig, eoutdir, eformat, ezoom, eprefix, ehtml, efollow, emp, emidx, efasta, efont, erefmp, erefmidx, esopts };

  bool bopt[ebopts];
  string sopt[esopts];

  const SpsPlot::Data & cdata;

  // out
  ostream & splot, & shtml;

  SpsPlot::Data::cpeptide_t & cpeptide;
  SpsPlot::Data::csequence_t & csequence;

  SpsPlot::Data::cpeptide_t cpeptide_null;
  SpsPlot::Data::csequence_t csequence_null;

  void init()
  {
    fill(bopt, bopt + sizeof(bopt) / sizeof(* bopt), false);

    sopt[epeakmasstol] = ".500001";
    sopt[eoutdir] = ".";
  }

  ContPlot(SpsPlot::Data & cdata, ostream & splot, ostream & shtml)
  :
    cdata(cdata), splot(splot), shtml(shtml), cpeptide(cpeptide_null), csequence(csequence_null)
  {
    init();
  }

  ContPlot(SpsPlot::Data & cdata, ostream & splot, ostream & shtml, SpsPlot::Data::csequence_t & csequence)
  :
    cdata(cdata), splot(splot), shtml(shtml), cpeptide(cpeptide_null), csequence(csequence)
  {
    init();
  }

  ContPlot(SpsPlot::Data & cdata, ostream & splot, ostream & shtml, SpsPlot::Data::cpeptide_t & cpeptide)
  :
    cdata(cdata), splot(splot), shtml(shtml), cpeptide(cpeptide), csequence(csequence_null)
  {
    init();
  }

  int operator () ();

  friend ostream & operator << (ostream & sout, const ContPlot & cplot)
  {
    return sout << "contig " << dec << cplot.sopt[econtig] << hex;
  }
};
#endif // SPECPLOT


/**
    External utils.
*/

extern void dbg_handler(int nsig) __attribute__ ((weak));

extern void getMasses(const string & sequence, sps::vector<double> &cmass, sps::vector<string> &cname, const map<char, double> & caa, bool simple = false);


};  // namespece specnets

#endif
