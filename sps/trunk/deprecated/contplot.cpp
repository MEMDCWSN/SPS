#include <cmath>
#include <cctype>
#include <list>
#include <limits>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "spsplot.h"
#include "ion.h"
#include "tuple.h"
#include "utils.h"
#include "label.h"
#include "range.h"
#include "vector.h"
#include "iomanip.h"
#include "spectrum.h"
#include "aminoacid.h"
#include "uncertain.h"
#include "dbg_print.h"


using namespace std;
namespace specnets {


int ContPlot::operator () ()
{
  // process options
  double npeakmasstol;
  istringstream(sopt[epeakmasstol]) >> npeakmasstol;

  double nzoom = 1.;
  istringstream(sopt[ezoom]) >> nzoom;

  // explicit contig and explicit subcontig
  const bool bn[] = {! sopt[econtig].empty() && sopt[econtig].find(':') == string::npos, ! sopt[econtig].empty() && sopt[econtig].find('-') != string::npos};

  // start, end & temporary contig indexes
  unsigned nc[3] = {0, 0};

  // subcontig index
  unsigned nba = 0;

  if (! sopt[econtig].empty())
  {
    istringstream s(sopt[econtig]);

    s >> nc[0];
    nc[0] --;

    switch (s.get())
    {
    case '-':
      s >> nba;
      break;
    case ':':
      s >> nc[1];
      nc[1] --;
      break;
    }
  }

  stringstream sns;

  // explicit contig
  if (bn[0])
    sns << '.' << sopt[econtig];

  const string sfilename = sopt[eprefix] + sns.str();
  const string slongfile = sopt[eoutdir] + "/" + sfilename;

  // local lookups
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

  splot << "nan=0\n";

  // render spectra
  const unsigned nuo = unsigned(fabs(floor(log10(npeakmasstol))));

  shtml.precision(nuo);
  shtml.setf(ios::fixed, ios::floatfield);

  cout.precision(nuo);
  cout.setf(ios::fixed, ios::floatfield);

  // iterate contigs
  for (abinfo_t::iterator i = cdata.cabinfo.lower_bound(nc[0]); bn[0] ? i != cdata.cabinfo.end() : i != cdata.cabinfo.end() && i->first <= nc[1]; ++ i)
  {
    nc[0] = i->first;

    // skip
    if (cdata.cconsset.specs[nc[0]].peakList.empty())
      // explicit contig
      if (bn[0])
        break;
      else
        continue;

    stringstream sns;
    sns << nc[0] + 1;

    // reverse flag
    bool breverse;

    // cancel flag
    if (breverse = ! sopt[emp].empty() && cdata.cproteinidx[nc[0]][2])
      cdata.cproteinidx[nc[0]][2] = 0;

    // global spec index to component spec index map
    map<unsigned, unsigned> ccomponent;

    // iterate spectra
    for (unsigned i = 0; i < cdata.cabinfo[nc[0]].first.first.size(); ++ i)
    {
      // global spec index
      const unsigned ns = cdata.cabinfo[nc[0]].first.first[i];

      // map global spec index with component spec index
      ccomponent[ns] = i;

      // expand bounds
      cdata.cspecset.specs[ns].addZPMpeaks(npeakmasstol, 0., true);

      // reverse if requested
      if (cdata.cabinfo[nc[0]].first.second[i] != breverse)
        cdata.cspecset.specs[ns].reverse(0);

      // cancel flag
      if (cdata.cabinfo[nc[0]].first.second[i])
        cdata.cabinfo[nc[0]].first.second[i] = 0;
    }

    // reverse if requested
    if (breverse)
    {
      cdata.cconsset.specs[nc[0]].reverse(0);

      // iterate vertices
      for (unsigned nv = 0; nv < cdata.cabinfo[nc[0]].second.size(); ++ nv)
        // iterate a vertex
        for (unsigned i = 0, j = cdata.cabinfo[nc[0]].second[nv].second.size(); i < j; ++ i)
        {
          // global spec index
          const unsigned ns = cdata.cabinfo[nc[0]].second[nv].first[i];

          // consistency check
          if (cdata.cspecset.specs[nc[0]].peakList.empty())
          {
            ostringstream smessage;
            throw runtime_error(note(smessage, __PRETTY_FUNCTION__, __LINE__).str().c_str());
          }

          cdata.cabinfo[nc[0]].second[nv].second[i] = (* cdata.cspecset.specs[ns].peakList.rbegin())[0] - cdata.cabinfo[nc[0]].second[nv].second[i];
        }

      reverse(cdata.cabinfo[nc[0]].second.begin(), cdata.cabinfo[nc[0]].second.end());
    }

    typedef list< sps::tuple<unsigned, string, sps::vector<string> > > cacid_t;

    cacid_t cacid;

    // do novo
    cacid.push_back(sps::make_tuple(0U, string("De Novo"), sps::vector<string>(cdata.cconsset.specs[nc[0]].peakList.size())));

    // homolog
    if (! sopt[emp].empty())
    {
      if (cdata.cproteinidx[nc[0]][0] < cdata.cfasta.size() && ! cdata.cprotset.specs[nc[0]].peakList.empty())
      {
        string sprotein(cdata.cfasta[cdata.cproteinidx[nc[0]][0]]);

        cacid.push_back(sps::make_tuple(1U, string("Homolog"), sps::vector<string>(unsigned((* cdata.cprotset.specs[nc[0]].peakList.rbegin())[0]))));

        for (unsigned j = 0; j < cdata.cprotset.specs[nc[0]].size() - 1; ++ j)
          for (unsigned k = unsigned(cdata.cprotset.specs[nc[0]][j][1]), l = 0; k < unsigned(cdata.cprotset.specs[nc[0]][j + 1][1]); ++ k, ++ l)
            cacid.back().m2[unsigned(cdata.cprotset.specs[nc[0]][j][0])].push_back(sprotein[unsigned(cdata.cprotset.specs[nc[0]][j][1]) + l]);
      }
    }

    // reference
    if (! sopt[erefmp].empty() && cdata.reference(nc[0]))
    {
      nc[2] = cdata.crefindex.find(nc[0])->second;

      if (cdata.creferenceidx[nc[2]][0] < cdata.cfasta.size() && ! cdata.crefset.specs[nc[2]].peakList.empty())
      {
        string sprotein(cdata.cfasta[cdata.creferenceidx[nc[2]][0]]);

        cacid.push_back(sps::make_tuple(2U, string("Reference"), sps::vector<string>(unsigned((* cdata.crefset.specs[nc[2]].peakList.rbegin())[0]))));

        for (unsigned j = 0; j < cdata.crefset.specs[nc[2]].size() - 1; ++ j)
          for (unsigned k = unsigned(cdata.crefset.specs[nc[2]][j][1]), l = 0; k < unsigned(cdata.crefset.specs[nc[2]][j + 1][1]); ++ k, ++ l)
            cacid.back().m2[unsigned(cdata.crefset.specs[nc[2]][j][0])].push_back(sprotein[unsigned(cdata.crefset.specs[nc[2]][j][1]) + l]);
      }
    }

    // override
    sps::vector<double> cmass;
    sps::vector<string> cname;
    ifstream speptide((sopt[eoutdir] + "/contig." + sns.str() + ".txt").c_str());

    if (speptide)
    {
      string s;

      getline(speptide, s);
      getMasses(s, cmass, cname, cdata.ciaa);

      cacid.push_back(sps::make_tuple(3U, string("User"), sps::vector<string>()));
      cacid.back().m2 = cname;
    }

    // vertex magnitude & vertex component index map
    multimap< pair<unsigned, double> , unsigned> cvmvi;

    // iterate vertices
    for (unsigned nv = 0; nv < cdata.cabinfo[nc[0]].second.size() && nv < cdata.cconsset.specs[nc[0]].peakList.size(); ++ nv)
      // map vertex order
      cvmvi.insert(make_pair(make_pair(cdata.cabinfo[nc[0]].second[nv].second.size(), cdata.cconsset.specs[nc[0]].peakList[nv][1]), nv));


    // vertex/spectrum knots (unique keys)
    map<unsigned, map<double, unsigned> > cvsknot_map;
    map<unsigned, hash_map<uf, unsigned> > cvsknot_hash;

    // iterate vertices
    for (unsigned nv = 0; nv < cdata.cabinfo[nc[0]].second.size(); ++ nv)
      // iterate a vertex
      for (unsigned i = 0; i < cdata.cabinfo[nc[0]].second[nv].second.size(); ++ i)
      {
        // global spec index
        const unsigned ns = cdata.cabinfo[nc[0]].second[nv].first[i];

        // vertex mass
        const double nm = cdata.cabinfo[nc[0]].second[nv].second[i];

        cvsknot_map[ns].insert(make_pair(nm, nv));
        cvsknot_hash[ns].insert(make_pair(uf(nm, npeakmasstol), nv));
      }

    // spectrum correction
    map<unsigned, double> ccorr;

    struct recurse_t
    {
      enum {evertex, edones};

      set<unsigned> done[edones];

      unsigned & nc;
      abinfo_t & cabinfo;
      map<unsigned, double> & ccorr;
      map<unsigned, map<double, unsigned> > & cvsknot_map;

      recurse_t(unsigned & nc, abinfo_t & cabinfo, map<unsigned, double> & ccorr, map<unsigned, map<double, unsigned> > & cvsknot_map) : nc(nc), cabinfo(cabinfo), ccorr(ccorr), cvsknot_map(cvsknot_map) {}

      void vertex(unsigned nv, unsigned ns)
      {
        if (done[evertex].find(nv) != done[evertex].end())
          return;

        done[evertex].insert(nv);

        unsigned i = 0;

        // iterate a vertex
        for (; i < cabinfo[nc].second[nv].first.size(); ++ i)
          if (cabinfo[nc].second[nv].first[i] == ns)
            break;

        const double nm = cabinfo[nc].second[nv].second[i];

        list<unsigned> spec;

        // iterate a vertex
        for (unsigned j = 0; j < cabinfo[nc].second[nv].first.size(); ++ j)
        {
          const unsigned ns2 = cabinfo[nc].second[nv].first[j];

          if (ns2 == ns)
            continue;

          if (ccorr.find(ns2) != ccorr.end())
            continue;

          const double nm2 = cabinfo[nc].second[nv].second[j];

          ccorr[ns2] = ccorr[ns] + nm2 - nm;

          spec.push_back(ns2);
        }

        // iterate spectra
        for (list<unsigned>::iterator j = spec.begin(); j != spec.end(); ++ j)
          // iterate a spectrum
          for (map<double, unsigned>::iterator i = cvsknot_map[* j].begin(); i != cvsknot_map[* j].end(); ++ i)
          {
            if (i->second == nv)
              continue;

            vertex(i->second, * j);
          }
      }
    } recurse(nc[0], cdata.cabinfo, ccorr, cvsknot_map);

    if (! cvmvi.empty())
      recurse.vertex(cvmvi.rbegin()->second, cdata.cabinfo[nc[0]].second[cvmvi.rbegin()->second].first[0]);


    // graph boundaries (xmin, xmax, ymax)
    double nbound[3] = {numeric_limits<double>::max(), 0, 0};

    // iterate spectra
    for (unsigned i = 0; i < cdata.cabinfo[nc[0]].first.first.size(); ++ i)
    {
      // global spec index
      const unsigned ns = cdata.cabinfo[nc[0]].first.first[i];

      if (nbound[0] > (* cdata.cspecset.specs[ns].peakList.begin())[0] - ccorr[ns])
        nbound[0] = (* cdata.cspecset.specs[ns].peakList.begin())[0] - ccorr[ns];
      if (nbound[1] < (* cdata.cspecset.specs[ns].peakList.rbegin())[0] - ccorr[ns])
        nbound[1] = (* cdata.cspecset.specs[ns].peakList.rbegin())[0] - ccorr[ns];

      for (unsigned i = 0; i < cdata.cspecset.specs[ns].peakList.size(); ++ i)
        if (nbound[2] < cdata.cspecset.specs[ns].peakList[i][1])
          nbound[2] = cdata.cspecset.specs[ns].peakList[i][1];
    }

    // iterate spectra
    for (unsigned i = 0; i < cdata.cabinfo[nc[0]].first.first.size(); ++ i)
    {
      // global spec index
      const unsigned ns = cdata.cabinfo[nc[0]].first.first[i];

      // leftmost spectrum shift
      ccorr[ns] += nbound[0];
    }

    nbound[1] -= nbound[0];
    nbound[0] -= nbound[0];

    // component spec index ordering
    sps::vector<unsigned> corder;
    map<unsigned, unsigned> cdisorder;

    {
      // spectrum correction inverse
      multimap<double, unsigned> ccorr_i;

      for (map<unsigned, double>::iterator i = ccorr.begin(); i != ccorr.end(); ++ i)
        ccorr_i.insert(make_pair(i->second, i->first));

      corder.resize(ccorr_i.size());

      unsigned ni = 0;
      for (multimap<double, unsigned>::iterator i = ccorr_i.begin(); i != ccorr_i.end(); ++ i, ++ ni)
      {
        corder[ni] = ccomponent[i->second];
        cdisorder[ccomponent[i->second]] = ni;
      }
    }

    /**
        @note
        The screen size is: 640. x 480.
    */

    const double ny = 480.;
    const double nph = 1. / ny; // pixel height
    const double nys = 1. - 10. * font_height * nph / nzoom - 2. * font_height * nph / nzoom; // y axis size

    // subcontig
    unsigned npb = cdata.cabinfo[nc[0]].first.first.size(); // number of spectra per subcontig
    unsigned nbs = 1; // number of subcontigs

    if (cdata.cabinfo[nc[0]].first.first.size() > ny * nys / font_height)
      if (cdata.cabinfo[nc[0]].first.first.size() > 100)
      {
        nbs = cdata.cabinfo[nc[0]].first.first.size() / 100 + 1;
        npb = 100;
      }

    // iterate subsets
    for (unsigned nb = nba; nb < nbs; ++ nb)
    {
      // number of spectra in this subcontig
      const unsigned nps = min<unsigned>(cdata.cabinfo[nc[0]].first.first.size(), npb * (nb + 1)) - npb * nb;

      double nx = 640.;
      double ny = 480.;
      double npw = 1. / nx; // pixel width
      double nph = 1. / ny; // pixel height
      double nxs = 1. - 8. * font_width * npw / nzoom - 3. * font_width * npw / nzoom; // x axis size
      double nys = 1. - 12. * font_height * nph / nzoom - 2. * font_height * nph / nzoom; // y axis size
      double nxp = 8. * font_width * npw / nzoom; // x axis position
      double nyp = 2. * font_height * nph / nzoom; // y axis position
      double ntp = 1. - 12. * font_height * nph / nzoom; // top margin position

      // stretch graph
      if (nps * font_height / ny > nys)
      {
        ny = nps * font_height + nyp / nph + (1. - ntp) / nph;
        nph = 1. / ny; // pixel height
        nys = 1. - 12. * font_height * nph / nzoom - 2. * font_height * nph / nzoom; // y axis size
        nyp = 2. * font_height * nph / nzoom; // y axis position
        ntp = 1. - 12. * font_height * nph / nzoom; // top margin position
      }

      // consensus, homolog, reference & user defined sequences
      stringstream ssequence[cacid.back().m0 + 1];

      stringstream sns[2];
      sns[0] << nc[0] + 1;

      // implicit contig
      if (! bn[0])
        sns[1] << '.' << nc[0] + 1;

      // subcontig presence
      if (nbs > 1)
      {
        sns[0] << '-' << nb;

        // implicit subcontig
        if (! bn[1])
          sns[1] << '-' << nb;
      }

      const string ssns[] = {sns[0].str(), sns[1].str()};

      if (sopt[eformat] == "png")
        splot << "set term png enhanced crop font \"" << font_name << "\" " << font_size << " size " << nx * nzoom << ", " << ny * nzoom << '\n';
      else
        splot << "set term " << (sopt[eformat] == "eps" ? "postscript eps" : sopt[eformat] + " crop") << " enhanced font \"" << font_name << "\" " << font_size << " size " << 5. * nzoom << ", " << 3.5 * nzoom << '\n';

      splot << "set output \"" << slongfile << ssns[1] << "." << sopt[eformat] << "\"\n";
      splot << "unset xlabel\n";
      splot << "unset ylabel\n";
      splot << "set multiplot\n";

      /**
          @note
          The screen size is: 411x287.7
          The graph left margin / width / right margin is: 43-354-14
          The graph bottom margin / height / top margin is: 16-222-50
      */

      /// main graph
      if (nzoom >= 1.)
        splot << "set title \"Protein Contig " << ssns[0] << "\" 0," << 3 + (cacid.size() - 1) * 2 << " font \"" << font_name << "," << font_size + 2 << "\"\n";

      splot << "unset origin\n";
      splot << "unset size\n";
      splot << "set border 3\n";
      splot << "set lmargin 9\n";
      splot << "set rmargin 3\n";
      splot << "set tmargin 12\n";
      splot << "set bmargin 2\n";

      // intervals
      double xtemp = nbound[1] / nzoom / 7.;
      double xmag = floor(log10(xtemp));
      double xpow = pow(10., xmag);
      double xmsd = unsigned(xtemp / xpow + 0.5);

      if (xmsd > 5.)
        xmsd = 10.;
      else if (xmsd > 2.)
        xmsd = 5.;
      else if (xmsd > 1.)
        xmsd = 2.;

      splot << "set xtics border nomirror tc rgb \"dark-gray\" " << xmsd * xpow << '\n';
      splot << "set ytics border nomirror tc rgb \"dark-gray\"\n";

      // iterate spectra
      splot << "set xtics out\n";
      splot << "set ytics (";
      for (unsigned i = npb * nb; i < min<unsigned>(cdata.cabinfo[nc[0]].first.first.size(), npb * (nb + 1)); ++ i)
      {
        // global spec index
        //const unsigned ns = cdata.cabinfo[nc[0]].first.first[i];
        const unsigned ns = cdata.cabinfo[nc[0]].first.first[corder[i]];

        // tic renaming
        if (i != npb * nb)
          splot << ", ";

        splot << "\"" << ns + 1 << "\" " << i + 1;
      }
      splot << ") out\n";


      /// vertices
      double nmmin = 0;

      if (! cvmvi.empty())
      {
        const unsigned nv = cvmvi.rbegin()->second;
        const double nmvert = cdata.cabinfo[nc[0]].second[nv].second[0] - ccorr[cdata.cabinfo[nc[0]].second[nv].first[0]];
        const double nmcons = cdata.cconsset.specs[nc[0]].peakList[nv][0];
        nmmin = nmvert - nmcons;

        // absolute consensus
        if (nmmin < 0.)
          nbound[0] += nmmin;

        // iterate vertices
        for (unsigned nv = 0; nv < cdata.cabinfo[nc[0]].second.size(); ++ nv)
        {
          // iterate a vertex
          for (unsigned nk = 1; nk < cdata.cabinfo[nc[0]].second[nv].second.size(); ++ nk)
          {
            // global spec index
            const unsigned ns[] = {cdata.cabinfo[nc[0]].second[nv].first[nk - 1], cdata.cabinfo[nc[0]].second[nv].first[nk]};

            // component spec index
            unsigned ni[] = {cdisorder[ccomponent[ns[0]]], cdisorder[ccomponent[ns[1]]]};

            if (ni[0] < npb * nb || ni[0] >= npb * nb + nps && ni[1] < npb * nb || ni[1] >= npb * nb + nps)
              continue;

            // vertex mass
            const double nm[] =
            {
              cdata.cabinfo[nc[0]].second[nv].second[nk - 1] - ccorr[ns[0]],
              cdata.cabinfo[nc[0]].second[nv].second[nk] - ccorr[ns[1]]
            };

            // position
            double np[] =
            {
              ((nm[0] - nbound[0]) / (nbound[1] - nbound[0])) * nxs + nxp,
              (1. / nps * (ni[0] - npb * nb)) * nys + nyp,
              ((nm[1] - nbound[0]) / (nbound[1] - nbound[0])) * nxs + nxp,
              (1. / nps * (ni[1] - npb * nb)) * nys + nyp
            };

            // straighten
            if (np[1] > np[3])
            {
              swap(ni[0], ni[1]);
              swap(np[0], np[2]);
              swap(np[1], np[3]);
            }

            // clip lines
            struct
            {
              void operator () (double (& np)[4], double ny)
              {
                const double nm = (np[3] - np[1]) / (np[2] - np[0]);
                const double nb = np[1] - (nm * np[0]);

                if (ny > .5)
                {
                  np[3] = ny;

                  if (std::isfinite(nm))
                    np[2] = (ny - nb) / nm;
                }
                else
                {
                  np[1] = ny;

                  if (std::isfinite(nm))
                    np[0] = (ny - nb) / nm;
                }
              }
            } cclip;

            if (np[1] < nyp - 10. * nph)
              cclip(np, nyp - 10. * nph);
            if (np[3] >= nyp + nys)
              cclip(np, nyp + nys);

            if (sopt[eformat] == "png")
              splot << "set arrow from screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " to screen " << np[2] + 8. * npw << ", screen " << np[3] + 10. * nph << " nohead back lt 0 lc rgb \"#C00000\"\n";
            else
              splot << "set arrow from screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " to screen " << np[2] + 8. * npw << ", screen " << np[3] + 10. * nph << " nohead back lt 2 lc rgb \"#C00000\"\n";

            if (np[1] >= nyp - 10. * nph)
              splot << "set arrow from screen " << np[0] << ", screen " << np[1] + 10. * nph << " to screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " nohead back lt 1 lc rgb \"#C00000\"\n";
            if (np[3] < nyp + nys)
              splot << "set arrow from screen " << np[2] << ", screen " << np[3] + 10. * nph << " to screen " << np[2] + 8. * npw << ", screen " << np[3] + 10. * nph << " nohead back lt 1 lc rgb \"#C00000\"\n";
          }

          if (cdata.cabinfo[nc[0]].second[nv].second.size())
          {
            // header link
            const unsigned nk = cdata.cabinfo[nc[0]].second[nv].second.size() - 1;

            // global spec index
            const unsigned ns = cdata.cabinfo[nc[0]].second[nv].first[nk];

            // component spec index
            const unsigned ni = cdisorder[ccomponent[ns]];

            if (ni < npb * nb)
              continue;
            if (ni >= npb * nb + nps)
              continue;

            // vertex mass
            const double nm[] =
            {
              cdata.cabinfo[nc[0]].second[nv].second[nk] - ccorr[ns],
              cdata.cconsset.specs[nc[0]].peakList[nv][0] + nmmin
            };

            // position
            const double np[] =
            {
              ((nm[0] - nbound[0]) / (nbound[1] - nbound[0])) * nxs + nxp,
              (1. / nps * (ni - npb * nb)) * nys + nyp,
              ((nm[1] - nbound[0]) / (nbound[1] - nbound[0])) * nxs + nxp,
              ntp + font_height / ny / nzoom * 1.5
            };

            if (sopt[eformat] == "png")
              splot << "set arrow from screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " to screen " << np[2] + 8. * npw << ", screen " << np[3] - 10. * nph << " nohead back lt 0 lc rgb \"#C00000\"\n";
            else
              splot << "set arrow from screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " to screen " << np[2] + 8. * npw << ", screen " << np[3] - 10. * nph << " nohead back lt 2 lc rgb \"#C00000\"\n";

            splot << "set arrow from screen " << np[0] + 8. * npw << ", screen " << np[1] + 10. * nph << " to screen " << np[0] << ", screen " << np[1] + 10. * nph << " nohead back lt 1 lc rgb \"#C00000\"\n";
            splot << "set arrow from screen " << np[2] + 8. * npw << ", screen " << np[3] - 10. * nph << " to screen " << np[2] << ", screen " << np[3] << " nohead back lt 1 lc rgb \"#C00000\"\n";
          }
        }

        /// header
        set<uf> cpeak[cacid.back().m0 + 1];

        for (unsigned i = 1, k; k = i, i < cdata.cconsset.specs[nc[0]].peakList.size(); ++ i)
        {
          // peptide fragment
          stringstream sfragment[cacid.back().m0 + 1];
          sps::vector<unsigned> nfragment(cacid.back().m0 + 1);

          // header position
          const double np[] = {ntp + font_height / ny / nzoom * 1.5};

          // segment mass
          const double nm[] = {cdata.cconsset.specs[nc[0]].peakList[i - 1][0] + nmmin, cdata.cconsset.specs[nc[0]].peakList[i][0] + nmmin};

          cpeak[0].insert(uf(nm[0], npeakmasstol));
          cpeak[0].insert(uf(nm[1], npeakmasstol));

          if (nzoom >= 1. && i == 1)
            splot << "set label " << "\"Consensus:\" at graph 0, screen " << np[0] << " offset 0,1 right\n";

          splot << "set arrow from " << nm[0] << ", screen " << np[0] << " to " << nm[1] << ", screen " << np[0] << " head front lw 2 lc rgb \"#C00000\"\n";

          // text offset
          hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(nm[1] - nm[0], npeakmasstol));
          if (l != cdata.caa.end())
            splot << "set label " << "\"" << l->second << "\" at " << (nm[1] - nm[0]) / 2 + nm[0] << ", screen " << np[0] << " offset 0,1 center\n";
          else if (nzoom >= 1.)
            splot << "set label " << "\"[" << fixed << setprecision(nuo) << nm[1] - nm[0] << resetiosflags(ios::floatfield) << setprecision(6) << "]\" at " << (nm[1] - nm[0]) / 2 + nm[0] << ", screen " << np[0] << " offset 0,1 center\n";

          // de novo
          if (l != cdata.caa.end())
            ssequence[0] << l->second;
          else
            ssequence[0] << "[" << fixed << setprecision(nuo) << nm[1] - nm[0] << resetiosflags(ios::floatfield) << setprecision(6) << "]";

          // homolog / reference
          unsigned nrc = 0;
          for (cacid_t::iterator ir = ++ cacid.begin(); ir != cacid.end() && i - 1 < ir->m2.size(); ++ ir, ++ nrc)
          {
            // header position
            const double np[] = {ntp + font_height / ny / nzoom * (1.5 + 2. * nrc), ntp + font_height / ny / nzoom * (3.5 + 2. * nrc)};

            if (nzoom >= 1. && i == 1)
              splot << "set label " << "\"" << ir->m1 << ":\" at graph 0, screen " << np[1] << " offset 0,1 right\n";

            // skip override
            if (ir->m1 == "User")
              continue;

            // skip space fillers
            if (! ir->m2[i - 1].size())
              continue;

            // bounds
            unsigned u;
            for (u = k; u < ir->m2.size() - 1 && ir->m2[u].empty(); ++ u)
              ;

            double nm2[] = {nm[0], nm[0], cdata.cconsset.specs[nc[0]].peakList[u][0] + nmmin};

            if (sopt[eformat] == "png")
              splot << "set arrow from " << nm2[0] << ", screen " << np[1] << " to " << nm2[0] << ", screen " << np[0] << " nohead back lt 0 lc rgb \"#C00000\"\n";
            else
              splot << "set arrow from " << nm2[0] << ", screen " << np[1] << " to " << nm2[0] << ", screen " << np[0] << " nohead back lt 0 lc rgb \"#C00000\"\n";

            for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
            {
              map<char, double>::const_iterator l = cdata.ciaa.find(ir->m2[i - 1][j]);
              if (l != cdata.ciaa.end())
                nm2[1] += l->second;
            }

            if (uf(nm2[1], npeakmasstol) != uf(nm2[2], npeakmasstol))
            {
              splot << "set arrow from " << nm2[0] << ", screen " << np[1] << " to " << nm2[2] << ", screen " << np[1] << " head back lw 2 lc rgb \"#C00000\"\n";
              splot << "set label \"(";

              nfragment[ir->m0] += ir->m2[i - 1].size();
              sfragment[ir->m0] << "(";

              for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
              {
                splot << ir->m2[i - 1][j];
                sfragment[ir->m0] << ir->m2[i - 1][j];
              }

              if (nzoom >= 1.)
                splot << "," << fixed << setprecision(nuo) << nm2[2] - nm2[1] << noshowpos << resetiosflags(ios::floatfield) << setprecision(6);

              sfragment[ir->m0] << "," << fixed << setprecision(nuo) << nm2[2] - nm2[1] << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << ")";

              splot << ")\"" << " at " << (nm2[2] - nm2[0]) / 2. + nm2[0] << ", screen " << np[1] << " offset 0,1 center\n";

              // homolog
              ssequence[ir->m0] << "(";

              for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
                ssequence[ir->m0] << ir->m2[i - 1][j];

              ssequence[ir->m0] << "," << fixed << setprecision(nuo) << nm2[2] - nm2[1] << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << ")";
            }
            else
            {
              nm2[1] = nm[0];

              for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
              {
                map<char, double>::const_iterator l = cdata.ciaa.find(ir->m2[i - 1][j]);
                if (l != cdata.ciaa.end())
                  nm2[1] += l->second;

                splot << "set arrow from " << nm2[0] << ", screen " << np[1] << " to " << nm2[1] << ", screen " << np[1] << " head back lw 2 lc rgb \"#C00000\"\n";

                splot << "set label \"";

                map<char, double>::const_iterator m = cdata.ciaa.find(ir->m2[i - 1][j]);
                if (m == cdata.ciaa.end())
                {
                  splot << "[" << fixed << setprecision(nuo) << l->second << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]";
                  ++ nfragment[ir->m0];
                  sfragment[ir->m0] << "[" << fixed << setprecision(nuo) << l->second << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]";
                }
                else
                {
                  splot << ir->m2[i - 1][j];
                  ++ nfragment[ir->m0];
                  sfragment[ir->m0] << ir->m2[i - 1][j];
                }

                splot << "\" at " << (nm2[1] - nm2[0]) / 2. + nm2[0] << ", screen " << np[1] << " offset 0,1 center\n";

                // homolog
                if (m == cdata.ciaa.end())
                  ssequence[ir->m0] << "[" << fixed << setprecision(nuo) << l->second << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]";
                else
                  ssequence[ir->m0] << ir->m2[i - 1][j];

                nm2[0] = nm2[1];
              }
            }

            cpeak[ir->m0].insert(uf(nm2[2], npeakmasstol));

            if (sopt[eformat] == "png")
              splot << "set arrow from " << nm2[2] << ", screen " << np[1] << " to " << nm2[2] << ", screen " << np[0] << " nohead back lt 0 lc rgb \"#C00000\"\n";
            else
              splot << "set arrow from " << nm2[2] << ", screen " << np[1] << " to " << nm2[2] << ", screen " << np[0] << " nohead back lt 2 lc rgb \"#C00000\"\n";
          }

          // save top sequence
          csequence[nc[0]][i - 1] = make_pair(nfragment[cacid.rbegin()->m0], sfragment[cacid.rbegin()->m0].str());
        }

        // override
        if (cacid.rbegin()->m1 == "User")
        {
          // header position
          const double np[] = {ntp + font_height / ny / nzoom * (1.5 + 2. * (cacid.size() - 2)), ntp + font_height / ny / nzoom * (3.5 + 2. * (cacid.size() - 2))};

          // segment mass
          double nm = 0.;

          for (unsigned i = 0; i < cmass.size() - 1; nm += cmass[i], ++ i)
          {
            if (nzoom >= 1. && i == 0)
              splot << "set label " << "\"User:\" at graph 0, screen " << np[1] << " offset 0,1 right\n";

            splot << "set arrow from " << nm << ", screen " << np[1] << " to " << nm + cmass[i] << ", screen " << np[1] << " head front lw 2 lc rgb \"#C00000\"\n";

            // text offset
            hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(cmass[i], npeakmasstol));
            if (l != cdata.caa.end())
              splot << "set label " << "\"" << l->second << "\" at " << cmass[i] / 2 + nm << ", screen " << np[1] << " offset 0,1 center\n";
            else if (nzoom >= 1.)
              splot << "set label " << "\"[" << fixed << setprecision(nuo) << cmass[i] << resetiosflags(ios::floatfield) << setprecision(6) << "]\" at " << cmass[i] / 2 + nm << ", screen " << np[1] << " offset 0,1 center\n";

            // de novo
            if (l != cdata.caa.end())
              ssequence[cacid.rbegin()->m0] << l->second;
            else
              ssequence[cacid.rbegin()->m0] << "[" << fixed << setprecision(nuo) << cmass[i] << resetiosflags(ios::floatfield) << setprecision(6) << "]";

            if (cpeak[cacid.rbegin()->m0].find(uf(nm, npeakmasstol)) != cpeak[cacid.rbegin()->m0].end())
              if (sopt[eformat] == "png")
                splot << "set arrow from " << nm << ", screen " << np[1] << " to " << nm << ", screen " << np[0] << " nohead back lt 0 lc rgb \"#C00000\"\n";
              else
                splot << "set arrow from " << nm << ", screen " << np[1] << " to " << nm << ", screen " << np[0] << " nohead back lt 2 lc rgb \"#C00000\"\n";
          }
        }
      }

      splot << "plot [" << nbound[0] << ":" << nbound[1] << "] [" << npb * nb << " + 0.5:" << min<unsigned>(cdata.cabinfo[nc[0]].first.first.size(), npb * (nb + 1)) + 0.5 << "] NaN notitle\n";
      splot << "unset arrow\n";
      splot << "unset label\n";
      splot << "unset title\n";


      /// spectra graphs
      splot << "set lmargin 0\n";
      splot << "set rmargin 0\n";
      splot << "set bmargin 0\n";
      splot << "set tmargin 0\n";
      splot << "unset ytics\n";
      splot << "set autoscale fix\n";
      splot << "set yrange [0:" << nbound[2] << "]\n";
      splot << "set tics out\n";
      splot << "set noxtics\n";
      splot << "set noytics\n";
      splot << "set border 1\n";

      // iterate spectra
      for (unsigned i = npb * nb; i < npb * nb + nps; ++ i)
      {
        // global spec index
        //const unsigned ns = cdata.cabinfo[nc[0]].first.first[i];
        const unsigned ns = cdata.cabinfo[nc[0]].first.first[corder[i]];

        // position / size
        const double np[] =
        {
          (((* cdata.cspecset.specs[ns].peakList.begin())[0] - ccorr[ns] - nbound[0]) / (nbound[1] - nbound[0])) * nxs + nxp,
          (1. / nps * (i - npb * nb)) * nys + nyp,
          (std::abs((* cdata.cspecset.specs[ns].peakList.rbegin())[0] - (* cdata.cspecset.specs[ns].peakList.begin())[0]) / (nbound[1] - nbound[0])) * nxs,
          (1. / nps) * nys
        };

        splot << "set origin " << np[0] << ", " << np[1] << '\n';
        splot << "set size " << np[2] << ", " << np[3] << '\n';

        // variables
        unsigned nplotr = 0;
        unsigned nplotb = 0;
        stringstream splotr;
        stringstream splotb;

        // process line
        for (unsigned j = 0; j < cdata.cspecset.specs[ns].peakList.size(); ++ j)
        {
          const double m = cdata.cspecset.specs[ns].peakList[j][0], h = cdata.cspecset.specs[ns].peakList[j][1];

          if (cvsknot_hash[ns].find(uf(m, npeakmasstol)) != cvsknot_hash[ns].end())
          {
            splotr.write((const char *)(const void *)(& m), sizeof(m));
            splotr.write((const char *)(const void *)(& h), sizeof(h));
            nplotr ++;
          }
          else
          {
            splotb.write((const char *)(const void *)(& m), sizeof(m));
            splotb.write((const char *)(const void *)(& h), sizeof(h));
            nplotb ++;
          }
        }
        splot << "set arrow from " << (* cdata.cspecset.specs[ns].peakList.begin())[0] << ", graph 0 to " << (* cdata.cspecset.specs[ns].peakList.begin())[0] << ", graph 1 nohead back lw 4 lc rgb \"#0000FF\"\n";
        splot << "set arrow from " << (* cdata.cspecset.specs[ns].peakList.rbegin())[0] << ", graph 0 to " << (* cdata.cspecset.specs[ns].peakList.rbegin())[0] << ", graph 1 nohead back lw 4 lc rgb \"#0000FF\"\n";

        // speptide[0]
        sps::vector<double> peakList(cdata.cconsset.specs[nc[0]].peakList.size() + 2);

        (* peakList.begin()) = - numeric_limits<double>::max();

        for (unsigned j = 0; j < cdata.cconsset.specs[nc[0]].peakList.size(); ++ j)
          peakList[j + 1] = cdata.cconsset.specs[nc[0]].peakList[j][0] + nmmin + ccorr[ns] + nbound[0];

        (* peakList.rbegin()) = numeric_limits<double>::max();

        for (map<double, unsigned>::iterator j = cvsknot_map[ns].begin(); j != cvsknot_map[ns].end(); ++ j)
          peakList[j->second + 1] = j->first;

        ostringstream speptide;
        speptide.precision(nuo);
        speptide.setf(ios::fixed, ios::floatfield);

        for (unsigned j = 1; j < peakList.size(); ++ j)
        {
          // segment mass
          double * nm[] = {& peakList[j - 1], & peakList[j]};

          if (uf(* nm[1], npeakmasstol) < uf(0., npeakmasstol))
            continue;

          if (uf((* cdata.cspecset.specs[ns].peakList.rbegin())[0] - AAJumps::massH2O, npeakmasstol) < uf(* nm[0], npeakmasstol))
            break;

          if (uf(* nm[0], npeakmasstol) < uf(0., npeakmasstol))
            * nm[0] = 0;

          if (uf((* cdata.cspecset.specs[ns].peakList.rbegin())[0] - AAJumps::massH2O, npeakmasstol) < uf(* nm[1], npeakmasstol))
            * nm[1] = (* cdata.cspecset.specs[ns].peakList.rbegin())[0] - AAJumps::massH2O;

          if (uf(* nm[1] - * nm[0], npeakmasstol) == uf(0, npeakmasstol))
            continue;

          hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(* nm[1] - * nm[0], npeakmasstol));
          if (l == cdata.caa.end())
            speptide << "[" << * nm[1] - * nm[0] << "]";
          else
            speptide << l->second;
        }

        cpeptide[0][ns] = make_pair(ssns[0], speptide.str());

        // peptide db
        for (cacid_t::iterator ir = ++ cacid.begin(); ir != cacid.end(); ++ ir)
        {
          // skip override
          if (ir->m1 == "User")
            continue;

          ostringstream speptide;
          speptide.precision(nuo);
          speptide.setf(ios::fixed, ios::floatfield);

          // header
          for (unsigned i = 0, k; k = i, i < peakList.size() - 1; i = k + 1)
          {
            // segment mass
            double nm[] = {peakList[i], peakList[i + 1]};

            if (uf(nm[1], npeakmasstol) < uf(0., npeakmasstol))
              continue;

            if (uf(nm[0], npeakmasstol) < uf(0., npeakmasstol))
              nm[0] = 0;

            if (uf(nm[1] - nm[0], npeakmasstol) == uf(0., npeakmasstol))
              continue;

            // peptide db
            if (! i || i >= ir->m2.size() || ir->m2[i - 1].empty())
            {
              hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(nm[1] - nm[0], npeakmasstol));
              if (l != cdata.caa.end())
                speptide << l->second;
              else
                speptide << "[" << fixed << setprecision(nuo) << nm[1] - nm[0] << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]";
            }
            else
            {
              // bounds
              for (; k < ir->m2.size() - 1 && ir->m2[k].empty(); ++ k)
                ;

              double nm2[] = {nm[0], nm[0], peakList[k + 1]};

              for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
              {
                map<char, double>::const_iterator l = cdata.ciaa.find(ir->m2[i - 1][j]);
                if (l != cdata.ciaa.end())
                  nm2[1] += l->second;
              }

              if (uf(nm2[1], npeakmasstol) != uf(nm2[2], npeakmasstol))
              {
                // homolog
                speptide << "(";

                for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
                  speptide << ir->m2[i - 1][j];

                speptide << "," << fixed << setprecision(nuo) << nm2[2] - nm2[1] << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << ")";
              }
              else
              {
                nm2[1] = nm[0];

                for (unsigned j = 0; j < ir->m2[i - 1].size(); ++ j)
                {
                  map<char, double>::const_iterator l = cdata.ciaa.find(ir->m2[i - 1][j]);
                  if (l != cdata.ciaa.end())
                    nm2[1] += l->second;

                  map<char, double>::const_iterator m = cdata.ciaa.find(ir->m2[i - 1][j]);

                  // homolog
                  if (m == cdata.ciaa.end())
                    speptide << "[" << fixed << setprecision(nuo) << l->second << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]";
                  else
                    speptide << ir->m2[i - 1][j];

                  nm2[0] = nm2[1];
                }
              }
            }
          }

          cpeptide[ir->m0][ns] = make_pair(ssns[0], speptide.str());
        }

        // override
        if (cacid.rbegin()->m1 == "User")
        {
          ostringstream speptide;
          speptide.precision(nuo);
          speptide.setf(ios::fixed, ios::floatfield);

          if (uf(ccorr[ns], npeakmasstol) != uf(0., npeakmasstol))
            speptide << "[" << ccorr[ns] << "]";

          // segment mass
          double nm = ccorr[ns];

          for (unsigned i = 0; i < cmass.size() - 1; nm += cmass[i], ++ i)
          {
            // text offset
            hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(cmass[i], npeakmasstol));

            // acid
            if (l != cdata.caa.end())
              speptide << l->second;
            else
              speptide << "[" << fixed << setprecision(nuo) << cmass[i] << resetiosflags(ios::floatfield) << setprecision(6) << "]";
          }

          cpeptide[3U][ns] = make_pair(ssns[0], speptide.str());
        }

/*
              // cancellable acids
              if (i == 1)
              {
                map<char, double>::iterator l = cdata.ciaa.find(acid[na[0]]);
                if (l->second == uf(nm2[1] - nm2[2], npeakmasstol))
                  na[0] = 1;
              }
              else if (i == peakList.size() - 2)
              {
                map<char, double>::iterator l = cdata.ciaa.find(acid[na[1] - 1]);
                if (l->second == uf(nm2[1] - nm2[2], npeakmasstol))
                  na[1] = acid.size() - 1;
              }
*/

        // iterate vertex/spectrum knots
        for (map<double, unsigned>::iterator j = cvsknot_map[ns].begin(), k = cvsknot_map[ns].end(); j != cvsknot_map[ns].end(); k = j, ++ j)
        {
          if (k == cvsknot_map[ns].end())
            continue;

          // segment mass
          const double nm[] = {k->first, j->first};

          // arrow color
          const double nd = cdata.cconsset.specs[nc[0]].peakList[j->second][0] - cdata.cconsset.specs[nc[0]].peakList[k->second][0];
          if (uf(nm[1] - nm[0], npeakmasstol) == uf(nd, npeakmasstol))
            splot << "set arrow from " << nm[0] << ", graph 0 to " << nm[1] << ", graph 0 head front lw 2 lc rgb \"#00A000\"\n";
          else
          {
            splot << "set arrow from " << nm[0] << ", graph 0 to " << nm[1] << ", graph 0 head front lw 2 lc rgb \"#C00000\"\n";

            hash_map<uf, char>::const_iterator l = cdata.caa.find(uf(nm[1] - nm[0], npeakmasstol));
            if (l != cdata.caa.end())
              splot << "set label " << "\"" << l->second << "\" at " << (nm[1] - nm[0]) / 2 + nm[0] << ", graph 0 offset 0, 1 center tc rgb \"#C00000\"\n";
            else
              splot << "set label " << "\"[" << showpos << fixed << setprecision(nuo) << nm[1] - nm[0] - nd << noshowpos << resetiosflags(ios::floatfield) << setprecision(6) << "]\" at " << (nm[1] - nm[0]) / 2 + nm[0] << ", graph 0 offset 0, 1 center tc rgb \"#C00000\"\n";
          }
        }

        // ...
        if (splotb.rdbuf()->in_avail() && splotr.rdbuf()->in_avail())
        {
          splot << "plot \"-\" binary record=" << nplotr << " format=\"%double%double\" with impulses lt 1 lw 4 lc rgb \"#C00000\" notitle, \"-\" binary record=" << nplotb << " format=\"%double%double\" with impulses lt 1 lw 4 lc rgb \"black\" notitle\n";
          splot << splotr.rdbuf() << sps::clear;
          splot << splotb.rdbuf() << sps::clear;
        }
        else if (splotr.rdbuf()->in_avail())
        {
          splot << "plot \"-\" binary record=" << nplotr << " format=\"%double%double\" with impulses lt 1 lw 4 lc rgb \"#C00000\" notitle\n";
          splot << splotr.rdbuf() << sps::clear;
        }
        else if (splotb.rdbuf()->in_avail())
        {
          splot << "plot \"-\" binary record=" << nplotb << " format=\"%double%double\" with impulses lt 1 lw 4 lc rgb \"black\" notitle\n";
          splot << splotb.rdbuf() << sps::clear;
        }

        splot << "unset arrow\n";
        splot << "unset label\n";
      }

      splot << "set nomultiplot\n";
      splot << "unset output\n";

      shtml << "<tr align=\"center\" valign=\"middle\">\n";
      shtml << "<td>" << ssns[0] << "</td>\n";
      shtml << "<td>" << nps << "</td>\n";
      shtml << "<td>\n";

      if (sopt[efollow].empty())
        shtml << "<h3 align=\"center\"><img src=\"" << sfilename << ssns[1] << "." << sopt[eformat] << "\"/></h3>\n";
      else
        shtml << "<h3 align=\"center\"><a href=\"" << sopt[efollow] << ssns[1] << ".html\"><img src=\"" << sfilename << ssns[1] << "." << sopt[eformat] << "\" style=\"border-style: none\"/></a></h3>\n";

      shtml << "</td>\n";
      shtml << "<td>\n";
      shtml << "  <table>\n";

      for (cacid_t::reverse_iterator ir = cacid.rbegin(); ir != cacid.rend(); ++ ir)
        shtml << "    <tr align=\"center\"><td><b>" << ir->m1 << "</b></td></tr><tr align=\"center\"><td>" << ssequence[ir->m0].rdbuf() << sps::clear << "</td></tr>\n";

      shtml << "    <tr align=\"center\">\n";
      shtml << "      <td><br>\n";
//      shtml << "        <form method=\"POST\" action=\"/cgi-bin/spsplot.fcgi\">";
//      shtml << "          <script type=\"text/javascript\">common();</script>\n";
//      shtml << "          <input type=\"hidden\" name=\"--contig\" value=\"" << ssns[0] << "\" />\n";
//      shtml << "          <input type=\"text\" name=\"--peptide\" style=\"text-transform: uppercase; width:100%\" /><br>\n";
//      shtml << "          <input type=submit value=\"Update\"/>\n";
//      shtml << "        </form>\n";
      shtml << "      </td>\n";
      shtml << "    </tr>\n";
      shtml << "  </table>\n";
      shtml << "</td>\n";

      if (! sopt[efasta].empty())
      {
        unsigned nidx = 0;

        shtml << "<td>";

        if (! sopt[emp].empty())
        {
          if (sopt[erefmp].empty() || ! cdata.reference(nc[0]))
            nidx = cdata.cproteinidx[nc[0]][0];
          else
            nidx = cdata.creferenceidx[cdata.crefindex.find(nc[0])->second][0];

          if (cdata.cfasta.getID(nidx))
            shtml << "<i>" << cdata.cfasta.getID(nidx) << "</i>";

          if (cdata.cfasta.getDesc(nidx))
            shtml << "<br><br>" << cdata.cfasta.getDesc(nidx);
        }

        shtml << "</td>\n";
      }

      shtml << "</tr>\n";

      // explicit subcontig
      if (bn[1])
        break;
    } // iterate subsets

    // explicit contig
    if (bn[0])
      break;
  } // iterate contigs

  return 0;
}


}; // namespece specnets
