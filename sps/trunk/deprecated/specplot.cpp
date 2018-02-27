#include <fcntl.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/file.h>

#include <cmath>
#include <list>
#include <limits>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

/*
#ifdef __GNUC__
#include <ext/hash_set>
#else
#include <hash_set>
#endif
*/

#include "spsplot.h"
#include "ion.h"
#include "label.h"
#include "range.h"
#include "aminoacid.h"
#include "uncertain.h"
#include "dbg_print.h"
#include "tuple.h"
#include "vector.h"
#include "iomanip.h"
#include "cstream.h"





using namespace std;
namespace specnets {


typedef list< range2d<double> > label_t;


int SpecPlot::operator () ()
{
  // process options
  unsigned nspectrum = 0;
  istringstream(sopt[espectrum]) >> nspectrum;

  if (! sopt[espectrum].empty())
    -- nspectrum;

  double npeakmasstol;
  istringstream(sopt[epeakmasstol]) >> npeakmasstol;

  double nzoom = 1.;
  istringstream(sopt[ezoom]) >> nzoom;

  char c;
  double nrange[3] = {0., 0., 0.};
  istringstream(sopt[erange]) >> nrange[0] >> c >> nrange[1];

  const unsigned npeptide = cpeptide.size();
  const unsigned ntitle = ! bopt[esuperpose] ? 1 : bopt[enoheader] ? 0 : npeptide;

  /**
      @note:  caa[0][peptide][spectra]
              caa[1][clusterfile][clusterscan]
  */
  map<unsigned, map<unsigned, string> > caa[2];

  // fetch
  string sc;

  // iterate peptides
  for (SpsPlot::Data::cpeptide_t::const_iterator ip = cpeptide.begin(); ip != cpeptide.end(); ++ ip)
  {
    int j = 0;

    // iterate spectra
    for (map<unsigned, pair<string, string> >::const_iterator is = ip->second.begin(); is != ip->second.end(); ++ is)
    {
      sc = is->second.first;

      if (sopt[espectrum].empty() || nspectrum == is->first)
        caa[0][ip->first][is->first] = is->second.second;
    }

    if (! caa[0][ip->first].size())
    {
      cerr << "error: invalid spectrum '" << sopt[espectrum] << "' (" << sopt[econtig] << ")" << endl;
      abort();

      return -1;
    }
  }

  // skip
  if (cpeptide.find(0) == cpeptide.end())
    return 0;

  /// peptide override
  // iterate spectra
  for (map<unsigned, pair<string, string> >::const_iterator is = cpeptide.find(0)->second.begin(); is != cpeptide.find(0)->second.end(); ++ is)
  {
    sc = is->second.first;

    stringstream sns;
    sns << "." << is->first + 1;

    ifstream speptide((sopt[eoutdir] + "/spectrum" + sns.str() + ".txt").c_str());

    if (speptide)
    {
      string s;

      getline(speptide, s);

      caa[0][3][is->first] = s;
    }

    if (! copt[ecluster].empty())
      // iterate cluster
      for (list< pair<unsigned, unsigned> >::const_iterator ie = cdata.ccluster.find(is->first)->second.begin(); ie != cdata.ccluster.find(is->first)->second.end(); ++ ie)
      {
        stringstream sns;

        if (cdata.cclusterscan.empty())
          sns << "." << is->first << "." << ie->first << "-" << ie->second;
        else
          sns << "." << is->first << "." << cdata.cclusterscan[ie->first][ie->second][0];

        ifstream speptide((sopt[eoutdir] + "/spectrum" + sns.str() + ".txt").c_str());

        if (speptide)
        {
          string s;

          getline(speptide, s);

          caa[1][ie->first][ie->second] = s;
        }
      }
  }

  map<unsigned, map<unsigned, sps::vector<double> > > cmass[2];
  map<unsigned, map<unsigned, sps::vector<string> > > cname[2];

  // iterate level
  for (unsigned nl = 0; nl < sizeof(caa) / sizeof(* caa); ++ nl)
    // iterate peptides
    for (map<unsigned, map<unsigned, string> >::iterator ip = caa[nl].begin(); ip != caa[nl].end(); ++ ip)
      // iterate spectra
      for (map<unsigned, string>::iterator is = ip->second.begin(); is != ip->second.end(); ++ is)
        getMasses(is->second, cmass[nl][ip->first][is->first], cname[nl][ip->first][is->first], cdata.ciaa, nzoom < 1.);

  // explicit colors
  map<string, unsigned> ccolor;

  for (SpsPlot::Data::ccolor_t::const_iterator i = cdata.ccolor.begin(); i != cdata.ccolor.end(); ++ i)
    ccolor[* i] = ccolor.size() - 1;

  stringstream sns1;

  if (! bopt[esuperpose])
    sns1 << "." << sc;

  const string ssns1 = sns1.str();
  const string sfilename = sopt[eprefix] + ssns1;
  const string slongfile = sopt[eoutdir] + "/" + sfilename;

  /// build indexes
  map< Ion::Type, list<const Ion *> > cion_t;
  map< Ion::Type, map<string, list<const Ion *> > > cion_ta;

  for (const Ion * p = Ion::ion; p != Ion::ion + sizeof(Ion::ion) / sizeof(* Ion::ion); ++ p)
    cion_t[p->type].push_back(p);

  for (const Ion * p = Ion::ion; p != Ion::ion + sizeof(Ion::ion) / sizeof(* Ion::ion); ++ p)
    cion_ta[p->type][p->acid].push_back(p);

  /// write script
  const unsigned font_height = 12;
  const unsigned font_width = font_height * 2 / 3;
  const string font_name = sopt[efont].empty() ? string("VeraBd.ttf") : sopt[efont] + "/VeraBd.ttf";

  double font_size;
  if (sopt[eformat] == "png")
    // png driver corresponding font
    font_size = font_height / 1.55;
  else
    // eps driver
    font_size = font_height;

  /**
      @note
      The screen size is: nx x ny
  */

  const double nx = 640.;
  const double ny = 480.;
  const double npw = 1. / nx; // pixel width
  const double nph = 1. / ny; // pixel height
  const double nxs = 1. - 10. * font_width * npw / nzoom - 3. * font_width * npw / nzoom; // x axis size
  double nys; // y axis size
  double ntp; // top margin position

  if (nzoom >= 1.)
  {
    nys = 1. - (ntitle * 2. + 6.) * font_height * nph / nzoom - 4. * font_height * nph / nzoom;
    ntp = 1. - (ntitle * 2. + 6.) * font_height * nph / nzoom;
  }
  else
  {
    nys = 1. - (ntitle * 2. + 4.) * font_height * nph / nzoom - 1. * font_height * nph / nzoom;
    ntp = 1. - (ntitle * 2. + 4.) * font_height * nph / nzoom;
  }

  // header position
  double np[] =
  {
    ntp + font_height / ny / nzoom * 3.5,
    ntp + font_height / ny / nzoom * 2.
  };

  splot << "nan=0\n";

  const unsigned nuo = unsigned(fabs(floor(log10(npeakmasstol))));

  shtml.precision(nuo);
  shtml.setf(ios::fixed, ios::floatfield);

  cout.precision(nuo);
  cout.setf(ios::fixed, ios::floatfield);

  // annotation
  ofstream sannot;

  // iterate peptides
  for (map<unsigned, string>::iterator ip = caa[0][0].begin(); ip != caa[0][0].end(); ++ ip)
  {
    list< pair<unsigned, unsigned> >::const_iterator ie;

    // associated cluster
    if (cdata.ccluster.find(nspectrum) != cdata.ccluster.end())
      ie = cdata.ccluster.find(nspectrum)->second.begin();

    // iterate cluster
    for (; copt[ecluster].empty() || cdata.ccluster.find(nspectrum) != cdata.ccluster.end() && ie != cdata.ccluster.find(nspectrum)->second.end(); ++ ie)
    {
      // explicit spectrum scan
      if (! sopt[escan].empty())
      {
        ostringstream sscan;
        sscan << ie->first << "-" << ie->second;

        if (sscan.str() != sopt[escan])
          continue;
      }

      // global spec clusteridx
      const unsigned ns = copt[ecluster].empty() ? ip->first : ie->second;

      // spectrum
      const Spectrum & spec = copt[ecluster].empty() ? cspecset.specs[ns] : cdata.cspecmap[ie->first].specs[ie->second];

      // pending ions
      map<Ion::Type, map<unsigned, map<unsigned, map<unsigned, pair<double, double> > > > > cionmh;

      // impulses
      unsigned nplotb[ccolor.size()];
      stringstream splotb[ccolor.size()];

      for (unsigned k = 0; k < sizeof(nplotb) / sizeof(* nplotb); ++ k)
        nplotb[k] = 0;

      // cluster peptide override
      if (! copt[ecluster].empty() && caa[1].find(ie->first) != caa[1].end() && caa[1][ie->first].find(ie->second) != caa[1][ie->first].end())
      {
        caa[0][3][ip->first].swap(caa[1][ie->first][ie->second]);
        cmass[0][3][ip->first].swap(cmass[1][ie->first][ie->second]);
        cname[0][3][ip->first].swap(cname[1][ie->first][ie->second]);
      }

      // iterate peptide inputs
      for (map<unsigned, map<unsigned, string> >::iterator ii = caa[0].begin(); ii != caa[0].end(); ++ ii)
      {
        // labels ordered by height
        multimap<double, sps::tuple<Ion::Type, unsigned, unsigned, const char * (*), double, string> > cion;

        // observed ions
        map<Ion::Type, map<unsigned, map<bool, set<double> > > > coion;

        if (ii->second.find(ip->first) != ii->second.end())
          splot.clear();
        else
          splot.setstate(ios::badbit);

        // clear stats
        cstat.erase(ns);
        cmark.erase(ns);

        // filename suffix
        stringstream sns2;

        if (copt[ecluster].empty())
          sns2 << "." << ns + 1;
        else
          sns2 << "." << ie->first << "." << ie->second;

        sns2 << "." << ii->first;

        const string ssns2 = sns2.str();

        // not superposed
        if (! bopt[esuperpose])
        {
          cionmh.clear();

          for (unsigned k = 0; k < sizeof(nplotb) / sizeof(* nplotb); ++ k)
            nplotb[k] = 0;

          for (unsigned k = 0; k < sizeof(splotb) / sizeof(* splotb); ++ k)
            splotb[k].str("");
        }

        // stats
        label_t label;
        double hmax = 0;
        double mmax = 0;

        for (unsigned k = 0; k < spec.peakList.size(); ++ k)
        {
          double nh = spec.peakList[k][1];

          if (hmax < nh)
            hmax = nh;
        }

        const double nshift = bopt[eshift] ? 0. : AAJumps::massHion;

        if (! bopt[esuperpose] || bopt[enoheader])
        {
          double nm = AAJumps::massH2O + nshift;

          for (unsigned k = 0; k < cmass[0][ii->first][ip->first].size(); ++ k)
            nm += cmass[0][ii->first][ip->first][k];

          mmax = spec.peakList.empty() ? nm : max<double>(nm, (* spec.peakList.rbegin())[0]);
        }
        else for (unsigned i = 0; i < ntitle; ++ i)
        {
          double nm = AAJumps::massH2O + nshift;

          for (unsigned k = 0; k < cmass[0][i][ip->first].size(); ++ k)
            nm += cmass[0][i][ip->first][k];

          mmax = spec.peakList.empty() ? nm : max<double>(mmax, max<double>(nm, (* spec.peakList.rbegin())[0]));
        }

        if (sopt[erange].empty())
          nrange[1] = mmax;

        nrange[2] = nrange[1] - nrange[0];

        label.push_back(range2d<double>(- numeric_limits<double>::max(), nrange[0] * nx * nxs * nzoom / nrange[2] + font_width * 2., - numeric_limits<double>::max(), numeric_limits<double>::max()));
        label.push_back(range2d<double>(nrange[1] * nx * nxs * nzoom / nrange[2] + font_width, numeric_limits<double>::max(), - numeric_limits<double>::max(), numeric_limits<double>::max()));

        if (! bopt[esuperpose] || ii->first == 0)
        {
          if (sopt[eformat] == "png")
            splot << "set term png enhanced crop font \"" << font_name << "\" " << font_size << " size " << nx * nzoom << ", " << ny * nzoom << '\n';
          else
            splot << "set term " << (sopt[eformat] == "eps" ? "postscript eps" : sopt[eformat] + " crop") << " enhanced font \"" << font_name << "\" " << font_size << " size " << 5. * nzoom << ", " << 3.5 * nzoom << '\n';

          if (sopt[eoutfile].empty())
            splot << "set output \"" << slongfile << ssns2 << "." << sopt[eformat] << "\"\n";
          else
            splot << "set output \"" << sopt[eoutdir] << "/" << sopt[eoutfile] << "\"\n";

          if (bopt[esuperpose])
            if (! bopt[everbose])
              if (sopt[eoutfile].empty())
                cout << slongfile << ssns2 << "." << sopt[eformat] << endl;
              else
                cout << sopt[eoutdir] << "/" << sopt[eoutfile] << endl;

          // annotation
          if (bopt[eannot])
          {
            if (sannot.open((slongfile + ssns2 + ".txt").c_str()), ! sannot)
            {
              cerr << "error: cannot open '" << slongfile << ssns2 << ".txt" << "': " << strerror(errno) << endl;
              return -2;
            }

            if (! bopt[everbose])
              cout << slongfile << ssns2 << ".txt" << endl;
          }

          if (nzoom >= 1.)
            if (! bopt[enotitle])
            {
              splot << "set title \"";

              if (! bopt[esuperpose])
                splot << "Consensus ";

              // implicit
              if (sopt[etitle].empty())
                if (copt[ecluster].empty())
                  splot << "Spectrum " << ns + 1;
                else if (cdata.cclusterscan.empty())
                  splot << "Spectrum Scan " << ie->first << "-" << ie->second;
                else
                  splot << "Spectrum Scan " << cdata.cclusterscan[ie->first][ie->second][0];
              // explicit
              else
                splot << sopt[etitle];

              if (! bopt[esuperpose])
                splot << " (Contig " << ssns1.substr(1) << ")";

              splot << "\" 0," << (ntitle * 2. + 3.) << " font \"" << font_name << "," << font_size + 2 << "\"\n";
            }

          splot << "set lmargin " << 10. << "\n";
          splot << "set rmargin " << 3. << "\n";

          if (nzoom >= 1.)
          {
            splot << "set tmargin " << (ntitle * 2. + 6.) << "\n";
            splot << "set bmargin " << 4. << "\n";
          }
          else
          {
            splot << "set tmargin " << (ntitle * 2. + 4.) << "\n";
            splot << "set bmargin " << 2. << "\n";
          }
          splot << "set autoscale fix\n";
          splot << "set xrange [" << nrange[0] << ":" << nrange[1] << "]\n";
          splot << "set yrange [" << 0 << ":]\n";
          splot << "set border 3\n";

          //splot << "unset xlabel\n";
          //splot << "unset ylabel\n";
          if (nzoom >= 1.)
          {
            splot << "set xlabel \"Mass / Charge (m/z)\"\n";
            splot << "set ylabel \"Intensity\"\n";
          }

          // intervals
          double xtemp = nrange[2] / nzoom / 7.;
          double xmag = floor(log10(xtemp));
          double xpow = pow(10., xmag);
          double xmsd = unsigned(xtemp / xpow + 0.5);

          if (xmsd > 5.)
            xmsd = 10.;
          else if (xmsd > 2.)
            xmsd = 5.;
          else if (xmsd > 1.)
            xmsd = 2.;

          double ytemp = hmax / nzoom / 8.;
          double ymag = floor(log10(ytemp));
          double ypow = pow(10., ymag);
          double ymsd = unsigned(ytemp / ypow + 0.5);

          if (ymsd > 5.)
            ymsd = 10.;
          else if (ymsd > 2.)
            ymsd = 5.;
          else if (ymsd > 1.)
            ymsd = 2.;

          splot << "set xtics axis nomirror tc rgb \"dark-gray\" " << xmsd * xpow << '\n';
          splot << "set ytics axis nomirror tc rgb \"dark-gray\" " << ymsd * ypow << '\n';
        }

        const char * scolor[][4] =
        {
          {
            "\"#0000C0\"",
            "\"#4900BF\"",
            "\"#BFBFFF\"",
            "\"#D8BFFF\""
          },
          {
            "\"#C00000\"",
            "\"#BF4900\"",
            "\"#FFBFBF\"",
            "\"#FDD8BF\""
          }
        };

        const char * (* pcolor);

        if (! bopt[esuperpose] || ii->first % 2 == 0)
          pcolor = scolor[0];
        else
          pcolor = scolor[1];

        double b = nshift;
        double y = AAJumps::massH2O + nshift;

        if (! bopt[enoheader])
        {
          splot << "set label \"b\" at graph 0, screen " << np[0] << " offset -2,0 right tc rgb " << pcolor[0] << "\n";
          splot << "set label \"y\" at graph 0, screen " << np[1] << " offset -2,0 right tc rgb " << pcolor[1] << "\n";

          for (int j = 0; j < cmass[0][ii->first][ip->first].size(); b += cmass[0][ii->first][ip->first][j], ++ j)
          {
            splot << "set arrow from " << b << ", screen " << np[0] << " to " << b << ", screen " << np[0] << " nohead lw 2 lc rgb " << pcolor[0] << "\n";

            if (j < cmass[0][ii->first][ip->first].size() - 1)
              splot << "set label \"" << cname[0][ii->first][ip->first][j] << "\" at " << b + cmass[0][ii->first][ip->first][j] / 2 << ", screen " << np[0] << " offset 0,0.5 center tc rgb " << pcolor[0] << " front\n";
          }

          for (int j = cmass[0][ii->first][ip->first].size(); j > 0; y += cmass[0][ii->first][ip->first][j - 1], j --)
          {
            splot << "set arrow from " << y << ", screen " << np[1] << " to " << y << ", screen " << np[1] << " nohead lw 2 lc rgb " << pcolor[1] << "\n";

            if (j > 0)
              splot << "set label \"" << cname[0][ii->first][ip->first][j - 1] << "\" at " << y + cmass[0][ii->first][ip->first][j - 1] / 2 << ", screen " << np[1] << " offset 0,0.5 center tc rgb " << pcolor[1] << " front\n";
          }

          if (sopt[eformat] == "png")
          {
            splot << "set arrow from " << nshift << ", screen " << np[0] << " to " << b << ", screen " << np[0] << " nohead lt 0 lc rgb " << pcolor[0] << "\n";
            splot << "set arrow from " << AAJumps::massH2O + nshift << ", screen " << np[1] << " to " << y << ", screen " << np[1] << " nohead lt 0 lc rgb " << pcolor[1] << "\n";
          }
          else
          {
            splot << "set arrow from " << nshift << ", screen " << np[0] << " to " << b << ", screen " << np[0] << " nohead lt 2 lc rgb " << pcolor[0] << "\n";
            splot << "set arrow from " << AAJumps::massH2O + nshift << ", screen " << np[1] << " to " << y << ", screen " << np[1] << " nohead lt 2 lc rgb " << pcolor[1] << "\n";
          }
        }

        if (! bopt[esuperpose] || ii->first == 0)
        {
          /// pre stats
          // m measured
          cstat[ns].second[1] = spec.parentMass;

          // m theoretical
          if (cdata.cinspect.find(ns) != cdata.cinspect.end())
            cstat[ns].second[2] = (y - nshift) +  AAJumps::massHion * cdata.cinspect.find(ns)->second.m2;
          else
            cstat[ns].second[2] = numeric_limits<double>::quiet_NaN();

          for (unsigned j = 0; j < cname[0][ii->first][ip->first].size(); ++ j)
            cstat[ns].first += cname[0][ii->first][ip->first][j];

          // m/z measured
          if (spec.parentCharge > 0.)
            cstat[ns].second[0] = (spec.parentMass + (spec.parentCharge - 1) * AAJumps::massHion) / spec.parentCharge;
          else
            cstat[ns].second[0] = spec.parentMass;

          // m/z theoretical
          if (cdata.cinspect.find(ns) != cdata.cinspect.end())
            cstat[ns].second[3] = cstat[ns].second[2] / cdata.cinspect.find(ns)->second.m2;
          else
            cstat[ns].second[3] = numeric_limits<double>::quiet_NaN();
        }

        // variables
        sps::vector<double> cline[2][2] =
        {
          {sps::vector<double>(cmass[0][ii->first][ip->first].size() + 2, 0.), sps::vector<double>(cmass[0][ii->first][ip->first].size() + 2, 0.)},
          {sps::vector<double>(cmass[0][ii->first][ip->first].size() + 2, 0.), sps::vector<double>(cmass[0][ii->first][ip->first].size() + 2, 0.)},
        };

        * cline[0][0].begin() = nshift;
        * cline[0][1].begin() = AAJumps::massH2O + nshift;
        * cline[0][0].rbegin() = b;
        * cline[0][1].rbegin() = y;

        * cline[1][0].begin() = nshift;
        * cline[1][1].begin() = AAJumps::massH2O + nshift;
        * cline[1][0].rbegin() = b;
        * cline[1][1].rbegin() = y;

        cmark[ns].resize(spec.peakList.size());

        // iterate peaks
        for (unsigned k = 0; k < spec.peakList.size(); ++ k)
        {
          bool bfound = false;
          const double nm = spec.peakList[k][0], nh = spec.peakList[k][1];

          // explicit
          const SpsPlot::Data::clabel_t::data_type::data_type * pl = 0;

          if (cdata.clabel.find(uf(nm, npeakmasstol)) != cdata.clabel.end())
            if (cdata.clabel.find(uf(nm, npeakmasstol))->second.find(uf(nh, npeakmasstol)) != cdata.clabel.find(uf(nh, npeakmasstol))->second.end())
              pl = & cdata.clabel.find(uf(nm, npeakmasstol))->second.find(uf(nh, npeakmasstol))->second;

          // iterate thru each possible ion type
          for (map< Ion::Type, list<const Ion *> >::iterator it = cion_t.begin(); it != cion_t.end(); ++ it)
          {
            const unsigned nheader = it->first / 3;

            // iterate thru each possible ion charge
            for (unsigned nz = 1; nz <= 2; ++ nz)
            {
              // get base mass
              const list<const Ion *> & cion_base = cion_ta[it->first][""];

              // no base mass for current ion type
              if (cion_base.empty())
                continue;

              double nm_base = cion_base.back()->mass + AAJumps::massHion * (nz - 1);
              double nm_peptide = nm_base - (AAJumps::massHion - nshift);
              unsigned norder = 0;
              string cname[0][ii->first];

              // stack ion masses
              switch (it->first)
              {
              case Ion::a:
              case Ion::b:
              case Ion::c:
                for (unsigned j = 0; j < cmass[0][ii->first][ip->first].size(); ++ j, ++ norder)
                {
                  if (nm_peptide + max(nm_base, 0.) + npeakmasstol > (nm * nz - AAJumps::massHion * (nz - 1)))
                    break;

                  nm_peptide += cmass[0][ii->first][ip->first][j];
                }
                break;

              case Ion::x:
              case Ion::y:
              case Ion::z:
                for (unsigned j = cmass[0][ii->first][ip->first].size(); j > 0; -- j, ++ norder)
                {
                  if (nm_peptide + max(nm_base, 0.) + npeakmasstol > (nm * nz - AAJumps::massHion * (nz - 1)))
                    break;

                  nm_peptide += cmass[0][ii->first][ip->first][j - 1];
                }
                norder = norder > 0 ? norder - 1 : 0;
                break;
              }

              // find matching b/y ion or isotopic
              list<const Ion *> cm = cion_t[it->first];

              // deterministic mass lookup
              for (list<const Ion *>::iterator i = cm.begin(), j; j = i, i != cm.end(); i = j)
              {
                ++ j;

                if (uf((* i)->mass, npeakmasstol) != uf((nm * nz - AAJumps::massHion * (nz - 1)) - nm_peptide + nm_base, npeakmasstol))
                  cm.erase(i);
              }

              bool bmatch[] = {false, false, false, false};

              if (! cm.empty())
                if (norder)
                  if (bopt[ebounds] || nm > 2. && nm < b - 2.)
                  {
                    bmatch[0] = true;

                    list<const Ion *>::iterator im;

                    for (im = cm.begin(); im != cm.end(); ++ im)
                      if ((* im)->acid == "")
                        break;

                    if (im != cm.end())
                      bmatch[1] = true;
                  }

              // explicit
              if (pl)
                bmatch[2] = true;

              // phos
              if (uf((cstat[ns].second[3] - 80.) / nz, npeakmasstol) == uf((nm * nz - AAJumps::massHion * (nz - 1)) - nm_peptide + nm_base, npeakmasstol))
                  bmatch[3] = true;
              else if (uf((cstat[ns].second[3] - 80. - AAJumps::massH2O) / nz, npeakmasstol) == uf((nm * nz - AAJumps::massHion * (nz - 1)) - nm_peptide + nm_base, npeakmasstol))
                  bmatch[3] = true;

              // implicit, explicit match or phos
              if (cdata.clabel.empty() && bmatch[0] || bmatch[2] || bmatch[3])
              {
                // pure implicit or explicit match
                if (bmatch[1] || bmatch[2])
                {
                  // pure implicit
                  if (bmatch[1])
                  {
                    // charge 1
                    if (nz == 1)
                    {
                      switch (it->first)
                      {
                      case Ion::b:
                      case Ion::y:
                        if (sopt[eformat] == "png")
                          splot << "set arrow from " << nm << ", graph 0 to " << nm << ", screen " << np[nheader] << " lc rgb " << pcolor[nheader + 2] << "lt 0 nohead back\n";
                        else
                          splot << "set arrow from " << nm << ", graph 0 to " << nm << ", screen " << np[nheader] << " lc rgb " << pcolor[nheader + 2] << "lt 2 nohead back\n";

                        // plain header segments
                        cline[0][nheader][norder] = nm_peptide;

                        if (norder >= 1)
                          if (double nm_prior = cline[0][nheader][norder - 1])
                            splot << "set arrow from " << nm_prior << ", screen " << np[nheader] << " to " << nm_peptide << ", screen " << np[nheader] << " nohead lc rgb " << pcolor[nheader] << '\n';
                        if (norder < cline[0][nheader].size() - 1)
                          if (double nm_prior = cline[0][nheader][norder + 1])
                            splot << "set arrow from " << nm_prior << ", screen " << np[nheader] << " to " << nm_peptide << ", screen " << np[nheader] << " nohead lc rgb " << pcolor[nheader] << '\n';
                        break;
                      }

                      switch (it->first)
                      {
                      case Ion::b:
                        // b ions
                        cstat[ns].second[12] += nh; // intensity

                        if (norder < cmass[0][ii->first][ip->first].size() - 1)
                          if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                          {
                            coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                            cstat[ns].second[16] += 1.; // observed
                          }
                        break;

                      case Ion::y:
                        // y ions
                        cstat[ns].second[13] += nh; // intensity

                        if (norder < cmass[0][ii->first][ip->first].size() - 1)
                          if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                          {
                            coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                            cstat[ns].second[17] += 1.; // observed
                          }
                        break;

                      case Ion::a:
                        // a ions
                        cstat[ns].second[14] += nh; // intensity

                        if (norder < cmass[0][ii->first][ip->first].size() - 1)
                          if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                          {
                            coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                            cstat[ns].second[19] += 1.; // observed
                          }
                        break;
                      }

                      // prm & stars
                      if (bopt[eshift])
                        // explained
                        if (! bfound)
                        {
                          cstat[ns].second[10] += nh; // intensity
                          cstat[ns].second[11] += 1.; // peaks
                        }
                    }
                    // charge 2
                    else if (nz == 2)
                      switch (it->first)
                      {
                      case Ion::b:
                      case Ion::y:
                        cstat[ns].second[15] += nh; // intensity

                        if (norder < cmass[0][ii->first][ip->first].size() - 1)
                          if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                          {
                            coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                            cstat[ns].second[18] += 1.; // observed
                          }
                        break;
                      }

                    switch (it->first)
                    {
                    case Ion::b:
                      cmark[ns][k] = true;

                      if (norder < cmass[0][ii->first][ip->first].size() - 1)
                        if (coion[it->first][3][bmatch[1]].find(norder) == coion[it->first][3][bmatch[1]].end())
                        {
                          coion[it->first][3][bmatch[1]].insert(norder); // observed ion
                          cstat[ns].second[23] += 1.; // observed
                        }
                      break;

                    case Ion::y:
                      cmark[ns][k] = true;

                      if (norder < cmass[0][ii->first][ip->first].size() - 1)
                        if (coion[it->first][3][bmatch[1]].find(norder) == coion[it->first][3][bmatch[1]].end())
                        {
                          coion[it->first][3][bmatch[1]].insert(norder); // observed ion
                          cstat[ns].second[24] += 1.; // observed
                        }
                      break;
                    }

                    // consecutive breaks
                    switch (it->first)
                    {
                    case Ion::b:
                    case Ion::y:
                      cline[1][nheader][norder] = nm_peptide;
                      break;
                    }
                  }

                  // implicit
                  if (cdata.clabel.empty())
                    if (cionmh[it->first][nz][norder][! bopt[esuperpose] ? 0 : ii->first].second < nh)
                      if (! bmatch[3])
                        cionmh[it->first][nz][norder][! bopt[esuperpose] ? 0 : ii->first] = make_pair(nm, nh);
                }
                // neutral loss ions, charge 1
                else if (nz == 1)
                {
                  list<const Ion *>::iterator im;

                  for (im = cm.begin(); im != cm.end(); ++ im)
                    if ((* im)->acid == "-H_20" || (* im)->acid == "-NH_3" || (* im)->acid == "(iso)")
                      break;

                  if (im != cm.end())
                  {
                    cstat[ns].second[14] += nh; // intensity

                    if (norder < cmass[0][ii->first][ip->first].size() - 1)
                      if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                      {
                        coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                        cstat[ns].second[19] += 1.; // observed
                      }
                  }
                }
                // neutral loss ions, charge 2
                else if (nz == 2)
                {
                  list<const Ion *>::iterator im;

                  for (im = cm.begin(); im != cm.end(); ++ im)
                    if ((* im)->acid == "-H_20" || (* im)->acid == "-NH_3" || (* im)->acid == "(iso)")
                      break;

                  if (im != cm.end())
                  {
                    cstat[ns].second[15] += nh; // intensity

                    if (norder < cmass[0][ii->first][ip->first].size() - 1)
                      if (coion[it->first][nz][bmatch[1]].find(norder) == coion[it->first][nz][bmatch[1]].end())
                      {
                        coion[it->first][nz][bmatch[1]].insert(norder); // observed ion
                        cstat[ns].second[18] += 1.; // observed
                      }
                  }
                }

                // phos
                if (bmatch[3])
                  sannot << nm << "\t" << nh << "\t" << Ion::cname[Ion::phos] << "\t\t\t\n";
                else
                  for (list<const Ion *>::iterator im = cm.begin(); im != cm.end(); ++ im)
                    sannot << nm << "\t" << nh << "\t" << Ion::cname[it->first] << "\t" << (*im)->acid << "\t" << norder << "\t" << nz << '\n';

                // ms/ms
                if (! bopt[eshift])
                  // explained
                  if (! bfound)
                  {
                    cstat[ns].second[10] += nh; // intensity
                    cstat[ns].second[11] += 1.; // peaks
                  }

                // multiple counted total
                if (bopt[ebounds] || nm > 2. && nm < b - 2.)
                {
                  cstat[ns].second[20] += nh; // intensity
                  cstat[ns].second[21] += 1.; // peaks
                }

                bfound = true;
              }
            }
          }

          if (! bfound)
            sannot << nm << "\t" << nh << '\n';

          // total
          if (bopt[ebounds] || nm > 2. && nm < b - 2.)
          {
            cstat[ns].second[8] += nh; // intensity
            cstat[ns].second[9] += 1.; // peaks
            cstat[ns].second[7] = cmass[0][ii->first][ip->first].size() - 2; // ions
          }

          // found
          if (bfound)
            // implicit or explicit with no color
            if (! pl)
            {
              splotb[ccolor["black"]].write((const char *)(const void *)(& nm), sizeof(nm));
              splotb[ccolor["black"]].write((const char *)(const void *)(& nh), sizeof(nh));
              ++ nplotb[ccolor["black"]];
            }
            // explicit
            else
            {
              const unsigned ncolor = ccolor[pl->m4];

              splotb[ncolor].write((const char *)(const void *)(& nm), sizeof(nm));
              splotb[ncolor].write((const char *)(const void *)(& nh), sizeof(nh));
              ++ nplotb[ncolor];
            }
          // not found
          else if (! bopt[esuperpose] || ii->first == 0)
          {
            splotb[ccolor["gray"]].write((const char *)(const void *)(& nm), sizeof(nm));
            splotb[ccolor["gray"]].write((const char *)(const void *)(& nh), sizeof(nh));
            ++ nplotb[ccolor["gray"]];
          }

          // explicit
          if (pl)
            cion.insert(make_pair(nh, sps::make_tuple(Ion::ctype[pl->m0], pl->m2, pl->m3, pcolor, nm, pl->m0 + pl->m1)));
        }

        if (bopt[esuperpose])
        {
          np[0] += font_height / ny / nzoom * 3.;
          np[1] += font_height / ny / nzoom * 3.;
        }

        unsigned nc[] = {0, 0, 0};

        for (unsigned j = 0; j < 2; ++ j)
          for (unsigned i = 1; i < cline[1][j].size(); nc[2] = max(nc[j], nc[2]), ++ i)
            if (cline[1][j][i - 1] && cline[1][j][i])
              ++ nc[j];
            else
              nc[j] = 0;

        cstat[ns].second[28] = nc[2];

        if (! bopt[esuperpose] || ii == -- caa[0].end())
        {
          unsigned nfound[] = {0, 0};

          // implicit
          if (cdata.clabel.empty())
            // iterate type
            for (map<Ion::Type, map<unsigned, map<unsigned, map<unsigned, pair<double, double> > > > >::iterator it = cionmh.begin(); it != cionmh.end(); ++ it)
              // iterate charge
              for (map<unsigned, map<unsigned, map<unsigned, pair<double, double> > > >::iterator ic = it->second.begin(); ic != it->second.end(); ++ ic)
                // iterate order
                for (map<unsigned, map<unsigned, pair<double, double> > >::iterator io = ic->second.begin(); io != ic->second.end(); ++ io)
                  // iterate layer
                  for (map<unsigned, pair<double, double> >::iterator il = io->second.begin(); il != io->second.end(); ++ il)
                  {
                    const unsigned nz = ic->first;
                    const Ion::Type type = it->first;
                    const unsigned no = io->first;
                    const unsigned nl = il->first;
                    const double nm = il->second.first;
                    const double nh = il->second.second;
                    const char * (* pcolor);

                    if (nl % 2 == 0)
                      pcolor = scolor[0];
                    else
                      pcolor = scolor[1];

                    cion.insert(make_pair(nh, sps::make_tuple(type, nz, no, pcolor, nm, Ion::cname[type])));
                  }

          // pixel / graph coordinate factor
          const double nf[] = {nx * nxs * nzoom / nrange[2], ny * nys * nzoom / hmax};

          // iterate height
          for (multimap<double, sps::tuple<Ion::Type, unsigned, unsigned, const char * (*), double, string> >::reverse_iterator ih = cion.rbegin(); ih != cion.rend(); ++ ih)
          {
            const unsigned nz = ih->second.m1;
            const Ion::Type type = ih->second.m0;
            const unsigned no = ih->second.m2;
            const double nm = ih->second.m4;
            const double nh = ih->first;
            const char * (* pcolor) = ih->second.m3;

            // non-empty label
            if (! ih->second.m5.empty())
            {
              const unsigned nlabel = ih->second.m5.substr(0, ih->second.m5.find('@')).size() + 1;

              // unclash labels
              for (unsigned rho = 0; rho < (nys * ny * nzoom / font_height); rho += 1)
              {
                label_t::iterator il;

                for (double phi = 0., npx = 0., npy = 0.; phi <= 3.1416; phi += 3.1416/4.)
                {
                  const range2d<double> sr((nm * nf[0] - npx) - font_width * (nlabel / 2), (nm * nf[0] - npx) + font_width * (nlabel / 2), (nh * nf[1] + npy), (nh * nf[1] + npy) + font_height * 2.);

                  for (il = label.begin(); il != label.end(); ++ il)
                    if (* il == sr)
                      break;

                  if (il == label.end())
                  {
                    label.push_back(sr);

                    // pixel to graph
                    const double ngx = npx / nf[0];
                    const double ngy = npy / nf[1];

                    // implicit
                    if (cdata.clabel.empty())
                      if (nz < 2)
                        splot << "set label \"" << ih->second.m5 << "@" << "^{" << no << "}\" at " << nm - ngx << ", " << nh + ngy << " center offset 0,1 front tc rgb " << pcolor[type / 3] << "\n";
                      else
                        splot << "set label \"" << ih->second.m5 << "@_{" << nz << "}^{" << no << "}\" at " << nm - ngx << ", " << nh + ngy << " center offset 0,1 front tc rgb " << pcolor[type / 3] << "\n";
                    // explicit
                    else
                      splot << "set label \"" << ih->second.m5 << "\" at " << nm - ngx << ", " << nh + ngy << " center offset 0,1 front tc rgb " << pcolor[type / 3] << "\n";

                    if (phi != 0.)
                      if (sopt[eformat] == "png")
                        splot << "set arrow from " << nm << ", " << nh << " to " << nm - ngx << ", " << nh + ngy << " lt 0 lw 2 lc rgb \"gray\"\n";
                      else
                        splot << "set arrow from " << nm << ", " << nh << " to " << nm - ngx << ", " << nh + ngy << " lt 2 lw 2 lc rgb \"gray\"\n";

                    break;
                  }

                  npx = cos(phi) * (font_width * rho);
                  npy = sin(phi) * (font_height * rho);
                }

                if (il == label.end())
                  break;
              }
            }

            switch (type)
            {
            case Ion::a:
            case Ion::b:
            case Ion::c:
              nfound[0] ++;
              break;

            case Ion::x:
            case Ion::y:
            case Ion::z:
              nfound[1] ++;
              break;
            }
          }

          cion.clear();

          // multiplot
          bool bfirst = true;
          for (SpsPlot::Data::ccolor_t::const_iterator i = cdata.ccolor.begin(); i != cdata.ccolor.end(); ++ i)
            if (splotb[ccolor[* i]].rdbuf()->in_avail())
            {
              if (bfirst)
                splot << "plot";
              else
                splot << ",";

              splot << " \"-\" binary record=" << nplotb[ccolor[* i]] << " format=\"%double%double\" with impulses lw 2 lt 1 lc rgb \"" << * i << "\" notitle";
              bfirst = false;
            }
          splot << "\n";

          for (SpsPlot::Data::ccolor_t::const_iterator i = cdata.ccolor.begin(); i != cdata.ccolor.end(); ++ i)
            if (splotb[ccolor[* i]].rdbuf()->in_avail())
              splot << splotb[ccolor[* i]].rdbuf() << sps::clear;

          splot << "unset arrow\n";
          splot << "unset label\n";
          splot << "unset output\n";

          // annotation
          if (bopt[eannot])
            sannot.close();

          if (ii == -- caa[0].end())
          {
            // contig number
            unsigned nc = 0;
            istringstream(sc) >> nc;
            -- nc;

            stringstream sns2;

            if (copt[ecluster].empty())
              sns2 << "." << ns + 1;
            else
              sns2 << "." << ie->first << "." << ie->second;

            const string ssns2 = sns2.str();

            shtml << "<tr align=\"center\" valign=\"middle\">\n";

            // index / scan
            if (! sopt[eindex].empty() || ! sopt[escan].empty())
              shtml << "<td rowspan=\"2\">\n";
            else
              shtml << "<td>";
            if (copt[ecluster].empty())
              shtml << ns + 1;
            else if (cdata.cclusterscan.empty())
              shtml << ie->first << "-" << ie->second;
            else
              shtml << cdata.cclusterscan[ie->first][ie->second][0];

            // contig
            if (! sopt[escan].empty())
              shtml << "<br>(" << sc << ")";
            shtml << "</td>\n";

            shtml << "<td><table><tr align=\"center\">\n";

            // spectrum image
            if (caa[0].find(3) != caa[0].end() && caa[0][3].find(ip->first) != caa[0][3].end() && ! caa[0][3][ip->first].empty())
              shtml << "<td><b>User</b></td>\n";
            else if (caa[0].find(2) != caa[0].end() && caa[0][2].find(ip->first) != caa[0][2].end() && ! caa[0][2][ip->first].empty())
              shtml << "<td><b>Reference</b></td>\n";
            else if (caa[0].find(1) != caa[0].end() && caa[0][1].find(ip->first) != caa[0][1].end() && ! caa[0][1][ip->first].empty())
              shtml << "<td><b>Homolog</b></td>\n";
            else
              shtml << "<td><b>De Novo</b></td>\n";
            shtml << "</tr>\n";

            shtml << "<tr>\n";

            // last peptide input
            ostringstream sns3;
            if (caa[0].find(3) != caa[0].end() && caa[0][3].find(ip->first) != caa[0][3].end() && ! caa[0][3][ip->first].empty())
              sns3 << ".3";
            else if (caa[0].find(2) != caa[0].end() && caa[0][2].find(ip->first) != caa[0][2].end() && ! caa[0][2][ip->first].empty())
              sns3 << ".2";
            else if (caa[0].find(1) != caa[0].end() && caa[0][1].find(ip->first) != caa[0][1].end() && ! caa[0][1][ip->first].empty())
              sns3 << ".1";
            else
              sns3 << ".0";
            const string ssns3 = sns3.str();

            shtml << "<td>";

            if (! sopt[escan].empty())
              shtml << "<a href=\"cluster." << sc << "." << sopt[espectrum] << ".html\"><img src=\"" << sfilename << ssns2 << ssns3 << "." << sopt[eformat] << "\" rel=\"lightbox\" style=\"border-style: none\"/></a>\n";
            else
              shtml << "<a href=\"l_" << sfilename << ssns2 << ssns3 << "." << sopt[eformat] << "\" rel=\"lightbox\"><img src=\"" << sfilename << ssns2 << ssns3 << "." << sopt[eformat] << "\" rel=\"lightbox\" style=\"border-style: none\"/></a>\n";

            shtml << "</td>";

          // Cluster button -- removed
            if (! sopt[efollow].empty())
              shtml << "<tr align=\"center\" valign=\"middle\"><td colspan=\"" << npeptide << "\"><input type=\"button\" value=\"Cluster\" onClick=\"window.location='" << sopt[efollow] << ssns1 << ssns2 << ssns3 << ".html'\"/></a></td></tr>\n";

            shtml << "</tr></table></td>\n";

            // peptide
            shtml << "<td>\n";
            shtml << "  <table>\n";
            if (caa[0].find(3) != caa[0].end() && caa[0][3].find(ip->first) != caa[0][3].end() && ! caa[0][3][ip->first].empty())
              shtml << "    <tr align=\"center\"><td><b>User</b></td></tr><tr align=\"center\"><td><a href=\"l_" << sfilename << ssns2 << ".3." << sopt[eformat] << "\" rel=\"lightbox\">" << caa[0][3][ip->first] << "</a></td></tr>\n";
            if (caa[0].find(2) != caa[0].end() && caa[0][2].find(ip->first) != caa[0][2].end() && ! caa[0][2][ip->first].empty())
              shtml << "    <tr align=\"center\"><td><b>Reference</b></td></tr><tr align=\"center\"><td><a href=\"l_" << sfilename << ssns2 << ".2." << sopt[eformat] << "\" rel=\"lightbox\">" << caa[0][2][ip->first] << "</a></td></tr>\n";
            if (caa[0].find(1) != caa[0].end() && caa[0][1].find(ip->first) != caa[0][1].end() && ! caa[0][1][ip->first].empty())
              shtml << "    <tr align=\"center\"><td><b>Homolog</b></td></tr><tr align=\"center\"><td><a href=\"l_" << sfilename << ssns2 << ".1." << sopt[eformat] << "\" rel=\"lightbox\">" << caa[0][1][ip->first] << "</a></td></tr>\n";
            shtml << "    <tr align=\"center\"><td><b>De Novo</b></td></tr><tr align=\"center\"><td><a href=\"l_" << sfilename << ssns2 << ".0." << sopt[eformat] << "\" rel=\"lightbox\">" << caa[0][0][ip->first] << "</a></td></tr>\n";
            shtml << "    <tr align=\"center\">\n";
            shtml << "      <td><br>\n";
//            shtml << "        <form method=\"POST\" action=\"/cgi-bin/spsplot.fcgi\">\n";
//            shtml << "          <script type=\"text/javascript\">common();</script>\n";
//            shtml << "          <input type=\"hidden\" name=\"--contig\" value=\"" << sopt[econtig] << "\" />\n";
//            shtml << "          <input type=\"hidden\" name=\"--spectrum\" value=\"";
//            if (copt[ecluster].empty())
//              shtml << ns + 1;
//            else if (cdata.cclusterscan.empty())
//              shtml << sopt[espectrum] << "." << ie->first << "-" << ie->second;
//            else
//              shtml << sopt[espectrum] << "." << cdata.cclusterscan[ie->first][ie->second][0];
//            shtml << "\" />\n";
//            shtml << "          <input type=\"text\" name=\"--peptide\" style=\"text-transform: uppercase; width:100%\" /><br>\n";
//            shtml << "          <input type=submit value=\"Update\"/>\n";
//            shtml << "        </form>\n";
            shtml << "      </td>\n";
            shtml << "    </tr>\n";
            shtml << "  </table>\n";
            shtml << "</td>\n";

            shtml << "<td>\n";
            shtml << spec.parentMass << "<br>";
            shtml << "</td>\n";
            shtml << "<td>\n";
            shtml << spec.parentCharge << "<br>";
            shtml << "</td>\n";
            shtml << "<td>\n";
            shtml << cstat[ns].second[23] * 100. / cstat[ns].second[7];
            shtml << "<br></td>\n";
            shtml << "<td>\n";
            shtml << cstat[ns].second[24] * 100. / cstat[ns].second[7];
            shtml << "<br></td>\n";
            shtml << "<td>\n";
            shtml << cstat[ns].second[10] * 100. / cstat[ns].second[8];
            shtml << "<br></td>\n";
            shtml << "</tr>\n";

            if (! sopt[eindex].empty() || ! sopt[escan].empty())
            {
              shtml << "<tr>\n";
              shtml << "<td colspan=\"0\" align=\"center\"><i>";

              if (! sopt[escan].empty())
              {
                if (! sopt[efasta].empty())
                  if (! cdata.reference(nc))
                    shtml << cdata.cfasta.getID(cdata.cproteinidx[nc][0]);
                  else
                    shtml << cdata.cfasta.getID(cdata.creferenceidx[cdata.crefindex.find(nc)->second][0]);
              }
              else if (! copt[ecluster].empty())
                shtml << cdata.cclusteridx.find(ie->first)->second;

              shtml << "</i></td></tr>\n";
            }

/// post stats
            cstat[ns].second[4] = cstat[ns].second[23] * 100. / cstat[ns].second[7];
            cstat[ns].second[5] = cstat[ns].second[24] * 100. / cstat[ns].second[7];
            cstat[ns].second[6] = cstat[ns].second[10] * 100. / cstat[ns].second[8];
          }
        }

        splot.clear();
      }

      if (copt[ecluster].empty())
        break;

      // cluster peptide override
      if (! copt[ecluster].empty() && caa[1].find(ie->first) != caa[1].end() && caa[1][ie->first].find(ie->second) != caa[1][ie->first].end())
      {
        caa[0][3][ip->first].swap(caa[1][ie->first][ie->second]);
        cmass[0][3][ip->first].swap(cmass[1][ie->first][ie->second]);
        cname[0][3][ip->first].swap(cname[1][ie->first][ie->second]);
      }
    }

    if (cpeptide.empty())
      break;
  }

  return 0;
}


void getMasses(const string & speptide, sps::vector<double> &cmass, sps::vector<string> &cname, const map<char, double> & caa, bool simple)
{
  istringstream cpeptide(speptide);

  enum {epeptide, eion} estate(epeptide);
  string::size_type ng;

  cmass.push_back(0);
  cname.push_back("");

  for (char nc; cpeptide >> nc; )
  {
    double nf = 0;

    switch (nc)
    {
    case '(':
      {
        estate = eion;

        string s;
        streampos p = cpeptide.tellg();
        getline(cpeptide, s, ')');
        cpeptide.seekg(p);

        ng = s.find(',');

        if (ng != string::npos)
          cname.back() += "(";
      }
      break;

    case ',':
      {
        ostringstream out;

        cpeptide >> nf;
        out << nf;

        cmass.back() += nf;

        if (! simple)
          cname.back() += string(1, nc) + out.str();
      }
      // reverse bracket
      continue;

    case '[':
      {
        ostringstream out;

        cpeptide >> nf;
        out << nf;

        cmass.back() += nf;

        if (! simple)
          cname.back() += string(1, nc) + out.str();
      }
      // reverse bracket
      continue;

    case ')':
      {
        estate = epeptide;

        if (ng != string::npos)
          cname.back() += ")";
      }
      break;

    case ']':
      {
        if (! simple)
          cname.back() += "]";
      }
      break;

    case '+':
      {
        ostringstream out;

        // robust double parsing
        const char * p[] = {&* (speptide.begin() + cpeptide.tellg()), 0};
        nf = strtod(p[0], & const_cast<char * &>(p[1]));
        cpeptide.ignore(p[1] - p[0]);

        cmass.pop_back();
        cname.pop_back();

        cmass.back() += nf;
        out << showpos << nf;

        if (! simple)
          cname.back() += out.str();
      }
      break;

    case '-':
      {
        ostringstream out;

        // robust double parsing
        const char * p[] = {&* (speptide.begin() + cpeptide.tellg()), 0};
        nf = strtod(p[0], & const_cast<char * &>(p[1]));
        cpeptide.ignore(p[1] - p[0]);

        cmass.pop_back();
        cname.pop_back();

        cmass.back() -= nf;
        out << showpos << nf;

        if (! simple)
          cname.back() += out.str();
      }
      break;

    default:
      {
        map<char, double>::const_iterator ia = caa.find(nc);
        if (ia != caa.end())
          cmass.back() += ia->second;

        cname.back() += string(1, nc);
      }
      break;
    }

    if (estate == epeptide)
    {
      cmass.push_back(0);
      cname.push_back("");
    }
  }
}


}; // namespece specnets
