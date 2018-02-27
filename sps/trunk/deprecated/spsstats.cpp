#include <errno.h>

#include <cmath>
#include <numeric>

#include "spsplot.h"
#include "range.h"


using namespace std;
using namespace specnets;

namespace specnets
{

struct type_t
{
  const double & nparentmasstol;
  map<unsigned, sps::vector<double> > & cmass;
  map<unsigned, sps::vector<string> > & cname;

  type_t(map<unsigned, sps::vector<double> > & cmass, map<unsigned, sps::vector<string> > & cname, double & nparentmasstol) : cmass(cmass), cname(cname), nparentmasstol(nparentmasstol) {}

  int operator () (const unsigned ns[2], const char * p[2], const char * q[2])
  {
    // PTM count
    struct ptm_t
    {
      typedef map<const char *, double> cptm_t;

      cptm_t operator () (const char * p, sps::vector<string> & cname, sps::vector<double> & cmass)
      {
        cptm_t cptm;

        for (unsigned i = 0, j = 0; i < cname.size() - 1; ++ i)
          if (cname[i].empty())
            cptm.insert(make_pair(p + j, cmass[i]));
          else
            ++ j;

        return cptm;
      }
    } ptm;

    ptm_t::cptm_t cptm[] = {ptm(p[0], cname[ns[0]], cmass[ns[0]]), ptm(p[1], cname[ns[1]], cmass[ns[1]])};

    // iterate acids
    for (ptm_t::cptm_t::iterator i = cptm[0].begin(), j; j = i, i != cptm[0].end(); i = j)
    {
      ++ j;

      if (cptm[1].find(i->first) != cptm[1].end())
        if (uf(cptm[1][i->first], nparentmasstol) == uf(i->second, nparentmasstol))
        {
          cptm[1].erase(i->first);
          cptm[0].erase(i->first);
        }
    }

    int ntype = int();
    const unsigned ndiff = cptm[0].size() + cptm[1].size();

    if (p[0] == p[1] && q[0] == q[1])
      switch (ndiff)
      {
      case 0: ntype = 1; break;
      case 1: ntype = 2; break;
      }
    else if (p[0] == p[1] || q[0] == q[1])
      switch (ndiff)
      {
      case 0: ntype = 3; break;
      default: ntype = 4; break;
      }
    else if (range1d<const char *>(p[0], q[0]) == range1d<const char *>(p[1], q[1]))
      switch (ndiff)
      {
      case 0: ntype = 6; break;
      default: ntype = 4; break;
      }

    return ntype;
  }
};

struct shift_t
{
  const double & nparentmasstol;
  map<unsigned, sps::vector<double> > & cmass;
  map<unsigned, sps::vector<string> > & cname;

  shift_t(map<unsigned, sps::vector<double> > & cmass, map<unsigned, sps::vector<string> > & cname, double & nparentmasstol) : cmass(cmass), cname(cname), nparentmasstol(nparentmasstol) {}

  pair<double, double> operator () (const unsigned ns[2], const char * p[2], const char * q[2])
  {
    // shift count
    struct shift_t
    {
      typedef map<const char *, double> cshift_t;

      cshift_t operator () (const char * p, sps::vector<string> & cname, sps::vector<double> & cmass)
      {
        cshift_t cshift;

        for (unsigned i = 0, j = 0; i < cname.size() - 1; ++ i)
        {
          cshift.insert(make_pair(p + j, cmass[i]));

          if (! cname[i].empty())
            ++ j;
        }

        return cshift;
      }
    } shift;

    pair<double, double> nshift;
    shift_t::cshift_t cshift[] = {shift(p[0], cname[ns[0]], cmass[ns[0]]), shift(p[1], cname[ns[1]], cmass[ns[1]])};

    // iterate acids
    for (shift_t::cshift_t::iterator i = cshift[0].begin(), j = cshift[1].begin(); i != cshift[0].end() || j != cshift[1].end(); )
      if (i == cshift[0].end())
      {
        nshift.second -= j->second;
        ++ j;
      }
      else if (j == cshift[1].end())
      {
        nshift.second += i->second;
        ++ i;
      }
      else if (i->first < j->first)
      {
        nshift.first += i->second;
        ++ i;
      }
      else if (j->first < i->first)
      {
        nshift.first -= j->second;
        ++ j;
      }
      else
      {
        ++ i;
        ++ j;
      }

    return nshift;
  }
};

} // namespace


int SpsPlot::stats()
{
  double nparentmasstol;
  istringstream(sopt[eparentmasstol]) >> nparentmasstol;

  const string sname[] = {"ms", "prm", "stars", "align", "all", "prm_pairs_hist", "stars_pairs_hist", "prm_pairs", "stars_pairs", "stars_seq", "contig_seq"};
  ofstream stxt[sizeof(sname) / sizeof(* sname)];
  SpecSet * pset[] = {& cdata.cmsset, & cdata.cprm, & cdata.cspecset};

  for (unsigned nt = 0; nt < sizeof(stxt) / sizeof(* stxt); ++ nt)
  {
    const string sfilename = "stats-" + sname[nt] + ".txt";
    const string slongfile = sopt[eoutdir] + "/" + sfilename;

    if (stxt[nt].open(slongfile.c_str()), ! stxt[nt])
    {
      cerr << "error: cannot open '" << sfilename << "': " << strerror(errno) << endl;
      return -2;
    }

    const unsigned nuo = unsigned(fabs(floor(log10(nparentmasstol))));

    stxt[nt].precision(nuo);
    stxt[nt].setf(ios::fixed, ios::floatfield);
  }

  map<unsigned, sps::vector<double> > cmass;
  map<unsigned, sps::vector<string> > cname;

  for (Data::cinspect_t::iterator ii = cdata.cinspect.begin(); ii != cdata.cinspect.end(); ++ ii)
    getMasses(ii->second.m0, cmass[ii->first], cname[ii->first], cdata.ciaa, true);

  type_t ctype(cmass, cname, nparentmasstol);
  shift_t cshift(cmass, cname, nparentmasstol);

  /// spectrum statistics
  for (unsigned nt = 0; nt < 3; ++ nt)
  {
    const string sfilename = "stats-" + sname[nt] + ".txt";
    const string slongfile = sopt[eoutdir] + "/" + sfilename;

    cout << slongfile << endl;
  }

  unsigned nline[3][2] = {0, 0};
  sps::vector<double> nave[3][2] = {{sps::vector<double>(14), sps::vector<double>(14)}, {sps::vector<double>(14), sps::vector<double>(14)}, {sps::vector<double>(14), sps::vector<double>(14)}};

  // stats
  Data::cstat_t cstat[3];

  // iterate lines
  for (Data::cinspect_t::iterator ii = cdata.cinspect.begin(); ii != cdata.cinspect.end(); ++ ii)
  {
    const unsigned ns = ii->first;

    ostringstream sns;
    sns << ns + 1;

    // iterate input types
    for (unsigned nt = 0; nt < 3; ++ nt)
    {
      // skip
      if (! pset[nt]->size())
        continue;

      switch (nt)
      {
      case 2:
        sps::vector< sps::vector<unsigned> >::iterator j;

        // iterate indices
        for (j = cdata.cstarsindex.begin(); j != cdata.cstarsindex.end(); ++ j)
          if (find(j->begin(), j->end(), ns) != j->end())
            break;

        if (j == cdata.cstarsindex.end())
          continue;

        break;
      }

      double nstat[14] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

      Data::cpeptide_t cpeptide;

      ofstream splot;
      ofstream shtml;
      SpecPlot cplot(cdata, * pset[nt], splot, shtml, cpeptide, cstat[nt]);

      cplot.sopt[SpecPlot::ezoom] = sopt[ezoom];
      cplot.sopt[SpecPlot::efont] = sopt[efont];
      cplot.sopt[SpecPlot::espectrum] = sns.str();
      cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];

      switch (nt)
      {
      case 1:
      case 2:
        cplot.bopt[SpecPlot::eshift] = true;
        break;
      }

      cpeptide[0][ns] = make_pair(string(), ii->second.m0);

      if (cplot())
      {
        cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;
        return -4;
      }

      nstat[0] = cstat[nt][ns].second[10] * 100. / cstat[nt][ns].second[8]; // % explained intensity
      nstat[1] = cstat[nt][ns].second[11] * 100. / cstat[nt][ns].second[9]; // % explained peaks
      nstat[2] = cstat[nt][ns].second[16] * 100. / cstat[nt][ns].second[7]; // % observed b ions
      nstat[3] = cstat[nt][ns].second[17] * 100. / cstat[nt][ns].second[7]; // % observed y ions
      nstat[4] = cstat[nt][ns].second[19] * 100. / cstat[nt][ns].second[7]; // % observed neutral-loss
      nstat[5] = cstat[nt][ns].second[18] * 100. / cstat[nt][ns].second[7]; // % observed charge 2
      nstat[6] = cstat[nt][ns].second[12] * 100. / cstat[nt][ns].second[20]; // % intensity explained by b-ions
      nstat[7] = cstat[nt][ns].second[13] * 100. / cstat[nt][ns].second[20]; // % intensity explained by y-ions
      nstat[8] = cstat[nt][ns].second[14] * 100. / cstat[nt][ns].second[20]; // % intensity explained by neutral-loss ions
      nstat[9] = cstat[nt][ns].second[15] * 100. / cstat[nt][ns].second[20]; // % intensity explained by doubly-charged ions
      nstat[10] = cstat[nt][ns].second[9]; // # spectrum peaks
      nstat[11] = cstat[nt][ns].second[3]; // mass error
      nstat[12] = (cstat[nt][ns].second[23] + cstat[nt][ns].second[24]) * 100. / cstat[nt][ns].second[7]; // % observed b, y, b^2 and y^2 ions
      nstat[13] = cstat[nt][ns].second[28]; // # consecutive breaks

      if (! std::isfinite(nstat[0])) nstat[0] = 0.;
      if (! std::isfinite(nstat[1])) nstat[1] = 0.;
      if (! std::isfinite(nstat[2])) nstat[2] = 0.;
      if (! std::isfinite(nstat[3])) nstat[3] = 0.;
      if (! std::isfinite(nstat[4])) nstat[4] = 0.;
      if (! std::isfinite(nstat[5])) nstat[5] = 0.;
      if (! std::isfinite(nstat[6])) nstat[6] = 0.;
      if (! std::isfinite(nstat[7])) nstat[7] = 0.;
      if (! std::isfinite(nstat[8])) nstat[8] = 0.;
      if (! std::isfinite(nstat[9])) nstat[9] = 0.;
      if (! std::isfinite(nstat[10])) nstat[10] = 0.;
      if (! std::isfinite(nstat[11])) nstat[11] = 0.;
      if (! std::isfinite(nstat[12])) nstat[12] = 0.;
      if (! std::isfinite(nstat[13])) nstat[13] = 0.;

      stxt[nt] << ns + 1 << '\t'; // index
      stxt[nt] << nstat[0] << '\t';
      stxt[nt] << nstat[1] << '\t';
      stxt[nt] << nstat[2] << '\t';
      stxt[nt] << nstat[3] << '\t';
      stxt[nt] << nstat[4] << '\t';
      stxt[nt] << nstat[5] << '\t';
      stxt[nt] << nstat[6] << '\t';
      stxt[nt] << nstat[7] << '\t';
      stxt[nt] << nstat[8] << '\t';
      stxt[nt] << nstat[9] << '\t';
      stxt[nt] << nstat[10] << '\t';
      stxt[nt] << nstat[11] << '\t';
      stxt[nt] << nstat[12] << '\t';
      stxt[nt] << nstat[13] << endl;

      const unsigned ng = nstat[7] <= nstat[6] ? 0 : 1;

      ++ nline[nt][ng];

      nave[nt][ng][0] += nstat[0];
      nave[nt][ng][1] += nstat[1];
      nave[nt][ng][2] += nstat[2];
      nave[nt][ng][3] += nstat[3];
      nave[nt][ng][4] += nstat[4];
      nave[nt][ng][5] += nstat[5];
      nave[nt][ng][6] += nstat[6];
      nave[nt][ng][7] += nstat[7];
      nave[nt][ng][8] += nstat[8];
      nave[nt][ng][9] += nstat[9];
      nave[nt][ng][10] += nstat[10];
      nave[nt][ng][11] += nstat[11];
      nave[nt][ng][12] += nstat[12];
      nave[nt][ng][13] += nstat[13];
    }
  }

  for (unsigned nt = 4; nt < 5; ++ nt)
  {
    const string sfilename = "stats-" + sname[nt] + ".txt";
    const string slongfile = sopt[eoutdir] + "/" + sfilename;

    cout << slongfile << endl;
  }

  // iterate input types
  for (unsigned nt = 0; nt < 2; ++ nt)
  {
    stxt[4] << '\t';
    stxt[4] << (nave[nt][0][0] + nave[nt][1][0]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][1] + nave[nt][1][1]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][2] + nave[nt][1][2]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][3] + nave[nt][1][3]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][4] + nave[nt][1][4]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][5] + nave[nt][1][5]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][6] + nave[nt][1][6]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][7] + nave[nt][1][7]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][8] + nave[nt][1][8]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][9] + nave[nt][1][9]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][10] + nave[nt][1][10]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][11] + nave[nt][1][11]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][12] + nave[nt][1][12]) / (nline[nt][0] + nline[nt][1]) << '\t';
    stxt[4] << (nave[nt][0][13] + nave[nt][1][13]) / (nline[nt][0] + nline[nt][1]) << '\t' << endl;
  }

  stxt[4] << '\t';
  stxt[4] << (nave[2][0][0]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][1]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][2]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][3]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][4]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][5]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][6]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][7]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][8]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][9]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][10]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][11]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][12]) / (nline[2][0]) << '\t';
  stxt[4] << (nave[2][0][13]) / (nline[2][0]) << '\t' << endl;

  stxt[4] << '\t';
  stxt[4] << (nave[2][1][0]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][1]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][2]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][3]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][4]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][5]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][6]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][7]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][8]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][9]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][10]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][11]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][12]) / (nline[2][1]) << '\t';
  stxt[4] << (nave[2][1][13]) / (nline[2][1]) << '\t' << endl;

  /// alignment statistics
  const string sfilename = "stats-" + sname[3] + ".txt";
  const string slongfile = sopt[eoutdir] + "/" + sfilename;

  cout << slongfile << endl;

  // missed pairs
  sps::vector< sps::vector<unsigned> > cmatrix[] = {sps::vector< sps::vector<unsigned> >(10, sps::vector<unsigned>(10)), sps::vector< sps::vector<unsigned> >(10, sps::vector<unsigned>(10))};

  // iterate pairs
  for (unsigned np = 0; np < cdata.cpairs.size(); ++ np)
  {
    const unsigned ns[] = {cdata.cpairs[np].spec1, cdata.cpairs[np].spec2};
    const char * p[] = {0, 0};
    unsigned j[] = {0, 0};

    if (cdata.cinspect.find(ns[0]) != cdata.cinspect.end())
      for (; j[0] < cdata.cfasta.size(); ++ j[0])
        if (p[0] = strstr(cdata.cfasta[j[0]], cdata.cinspect[ns[0]].m1.c_str()))
          break;

    if (cdata.cinspect.find(ns[1]) != cdata.cinspect.end())
      for (; j[1] < cdata.cfasta.size(); ++ j[1])
        if (p[1] = strstr(cdata.cfasta[j[1]], cdata.cinspect[ns[1]].m1.c_str()))
          break;

    if (cdata.cinspect.find(ns[0]) == cdata.cinspect.end() || cdata.cinspect.find(ns[1]) == cdata.cinspect.end())
      stxt[3] << -1 << '\t' << -1;
    else if (j[0] != j[1] || ! p[0] || ! p[1])
      stxt[3] << -1 << '\t' << -1;
    else
    {
      const char * q[] = {p[0] + cdata.cinspect[ns[0]].m1.size(), p[1] + cdata.cinspect[ns[1]].m1.size()};

      const int ntype = ctype(ns, p, q);
      const pair<double, double> nshift = cshift(ns, p, q);

      stxt[3] << ntype << '\t';

      if (ntype != 0)
        if (uf(cdata.cpairs[np].shift1, nparentmasstol) == uf(nshift.first, nparentmasstol) || uf(cdata.cpairs[np].shift1, nparentmasstol) == uf(nshift.second, nparentmasstol) || uf(cdata.cpairs[np].shift2, nparentmasstol) == uf(nshift.first, nparentmasstol) || uf(cdata.cpairs[np].shift2, nparentmasstol) == uf(nshift.second, nparentmasstol))
          stxt[3] << 1;
        else
          stxt[3] << 0;
      else
        stxt[3] << -1;

      switch (ntype)
      {
      case 0:
        stxt[3] << '\t' << -1;
        stxt[3] << '\t' << -1 << '\t' << -1;
        stxt[3] << '\t' << -1 << '\t' << -1;
        break;

      default:
        double nintensity[] = {0., 0., 0., 0.};

        for (unsigned i = 0; i < cdata.cprm[cdata.cpairs[np].spec1].size(); ++ i)
          nintensity[0] += cdata.cprm[cdata.cpairs[np].spec1][i][1];

        for (unsigned i = 0; i < cdata.cprm[cdata.cpairs[np].spec2].size(); ++ i)
          nintensity[1] += cdata.cprm[cdata.cpairs[np].spec2][i][1];

        for (unsigned i = 0; i < cdata.cmatchpairs[np * 2].size(); ++ i)
          nintensity[2] += cdata.cmatchpairs[np * 2][i][1];

        for (unsigned i = 0; i < cdata.cmatchpairs[np * 2 + 1].size(); ++ i)
          nintensity[3] += cdata.cmatchpairs[np * 2 + 1][i][1];

        stxt[3] << '\t' << min(cdata.cmatchpairs[np * 2].size(), cdata.cmatchpairs[np * 2 + 1].size());
        stxt[3] << '\t' << nintensity[2] * 100. / nintensity[0] << '\t' << nintensity[3] * 100. / nintensity[1];

        double nmass[] = {0., 0.};

        if (cdata.cpairs[np].shift1 < 0.)
          nmass[0] = min((* cdata.cmatchpairs[np * 2 + 1].peakList.rbegin())[0] + cdata.cpairs[np].shift1, (* cdata.cmatchpairs[np * 2].peakList.rbegin())[0]);
        else
          nmass[0] = min((* cdata.cmatchpairs[np * 2].peakList.rbegin())[0] - cdata.cpairs[np].shift1, (* cdata.cmatchpairs[np * 2 + 1].peakList.rbegin())[0]);

        if (cdata.cpairs[np].shift1 < 0.)
          nmass[1] = min((* cdata.cmatchpairs[np * 2 + 1].peakList.rbegin())[0] - cdata.cpairs[np].shift1, (* cdata.cmatchpairs[np * 2].peakList.rbegin())[0]);
        else
          nmass[1] = min((* cdata.cmatchpairs[np * 2].peakList.rbegin())[0] + cdata.cpairs[np].shift1, (* cdata.cmatchpairs[np * 2 + 1].peakList.rbegin())[0]);

        stxt[3] << '\t' << nmass[0] * 100. / (* cdata.cmatchpairs[np * 2].peakList.rbegin())[0] << '\t' << nmass[1] * 100. / (* cdata.cmatchpairs[np * 2 + 1].peakList.rbegin())[0];
        //stxt[3] << '\t' << ns[0] << '\t' << ns[1] << '\t' << cdata.cpairs[np].shift1 << '\t' << cdata.cpairs[np].shift2 << '\t' << min(nshift[0], nshift[1]) << '\t' << max(nshift[0], nshift[1]) << '\t' << cdata.cinspect[ns[0]].m0 << '\t' << cdata.cinspect[ns[1]].m0;
        break;
      }
    }

    stxt[3] << endl;
  }

  if (sopt[epairs].empty())
    return 0;

  /// missed pairs
  for (unsigned nt = 5; nt < 9; ++ nt)
  {
    const string sfilename = "stats-" + sname[nt] + ".txt";
    const string slongfile = sopt[eoutdir] + "/" + sfilename;

    cout << slongfile << endl;
  }

  unsigned nc = 0;

  // cached pairs
  set< pair<unsigned, unsigned> > ccache[2];

  for (sps::vector<Results_PA>::iterator ip = cdata.cpairs.begin(); ip != cdata.cpairs.end(); ++ ip)
    ccache[0].insert(make_pair(ip->spec1, ip->spec2));

  for (sps::vector<Results_PA>::iterator ip = cdata.cstarspairs.begin(); ip != cdata.cstarspairs.end(); ++ ip)
    ccache[1].insert(make_pair(ip->spec1, ip->spec2));

  // iterate lines
  for (Data::cinspect_t::iterator ii = cdata.cinspect.begin(); ii != cdata.cinspect.end(); ++ ii)
    // iterate lines
    for (Data::cinspect_t::iterator ij = cdata.cinspect.begin(); ij != cdata.cinspect.end(); ++ ij)
    {
      const unsigned ns[] = {ii->first, ij->first};
      const char * p[] = {0, 0};
      unsigned j[] = {0, 0};

      if (ns[0] == ns[1])
        continue;

      // protein mapping
      if (cdata.cinspect.find(ns[0]) != cdata.cinspect.end())
        for (; j[0] < cdata.cfasta.size(); ++ j[0])
          if (p[0] = strstr(cdata.cfasta[j[0]], cdata.cinspect[ns[0]].m1.c_str()))
            break;

      if (cdata.cinspect.find(ns[0]) != cdata.cinspect.end())
        for (; j[1] < cdata.cfasta.size(); ++ j[1])
          if (p[1] = strstr(cdata.cfasta[j[1]], cdata.cinspect[ns[1]].m1.c_str()))
            break;

      if (cdata.cinspect.find(ns[0]) == cdata.cinspect.end() || cdata.cinspect.find(ns[1]) == cdata.cinspect.end())
        continue;
      else if (j[0] != j[1] || ! p[0] || ! p[1])
        continue;

      const char * q[] = {p[0] + cdata.cinspect[ns[0]].m1.size(), p[1] + cdata.cinspect[ns[1]].m1.size()};

      if (range1d<const char *>(p[0], q[0]) != range1d<const char *>(p[1], q[1]))
        continue;

      // overlap percentage
      const char * nk[] = {p[0], p[1], q[0], q[1]};
      double nm[] = {max(cstat[0][ns[0]].second[2], cstat[0][ns[1]].second[2]), 0.};

      sort(nk, nk + sizeof(nk) / sizeof(* nk));

      // mass intersection
      for (const char * pc = nk[1]; pc != nk[2]; ++ pc)
        if (cdata.ciaa.find(* pc) != cdata.ciaa.end())
          nm[1] += cdata.ciaa.find(* pc)->second;

      double no = nm[1] * 100. / nm[0];

      if (! std::isfinite(no)) no = 0.;

      if (no < 20.)
        continue;

      // iterate input types
      for (unsigned nt = 5; nt < 7; ++ nt)
      {
        // skip if detected
        switch (nt)
        {
        case 5:
          if (ccache[0].find(make_pair(ns[0], ns[1])) != ccache[0].end())
            continue;
          if (ccache[0].find(make_pair(ns[1], ns[0])) != ccache[0].end())
            continue;
          break;

        case 6:
          if (ccache[1].find(make_pair(ns[0], ns[1])) != ccache[1].end())
            continue;
          if (ccache[1].find(make_pair(ns[1], ns[0])) != ccache[1].end())
            continue;
          break;
        }

        double ni = min(cstat[nt - 4][ns[0]].second[10] * 100. / cstat[nt - 4][ns[0]].second[8], cstat[nt - 4][ns[1]].second[10] * 100. / cstat[nt - 4][ns[1]].second[8]);

        if (! std::isfinite(ni)) ni = 0.;

        ++ cmatrix[nt - 5][min(unsigned(ni) / 10u, 9u)][min(unsigned(no) / 10u, 9u)];

        if (no < 45.)
          continue;

        const int ntype = ctype(ns, p, q);
        const pair<double, double> nshift = cshift(ns, p, q);

        stxt[nt + 2] << ns[0] + 1 << '\t' << ns[1] + 1 << '\t' << ntype << '\t' << nshift.first << '\t' << ni << '\t' << no << '\t' << string(nk[0], nk[1] - nk[0]) << '[' << string(nk[1], nk[2] - nk[1]) << ']' << string(nk[2], nk[3] - nk[2]) << '\n';
      }
    }

  // missed pairs
  for (unsigned nt = 5; nt < 7; ++ nt)
  {
    for (unsigned nx = 0; nx < cmatrix[nt - 5].size(); ++ nx)
    {
      stxt[nt] << '\t';

      if (nx < 9)
        stxt[nt] << nx * 10 << "-" << (nx + 1) * 10 - 1 << "%";
      else
        stxt[nt] << nx * 10 << "-100%";
    }
    stxt[nt] << '\n';

    unsigned ny = 0;
    for (sps::vector< sps::vector<unsigned> >::iterator iy = cmatrix[nt - 5].begin(); iy != cmatrix[nt - 5].end(); ++ iy, ++ ny)
    {
      if (ny < 9)
        stxt[nt] << ny * 10 << "-" << (ny + 1) * 10 - 1 << "%";
      else
        stxt[nt] << ny * 10 << "-100%";

      for (sps::vector<unsigned>::iterator ix = iy->begin(); ix != iy->end(); ++ ix)
      {
        stxt[nt] << '\t';

        if (* ix)
          stxt[nt] << * ix;
      }
      stxt[nt] << '\n';
    }
  }

  /// sequencing statistics
  for (unsigned nt = 9; nt < 11; ++ nt)
  {
    const string sfilename = "stats-" + sname[nt] + ".txt";
    const string slongfile = sopt[eoutdir] + "/" + sfilename;

    cout << slongfile << endl;
  }

  // iterate contigs
  for (abinfo_t::iterator ic = cdata.cabinfo.begin(); ic != cdata.cabinfo.end(); ++ ic)
  {
    const unsigned nc = ic->first;

    // skip
    if (cdata.cconsset.specs[nc].peakList.empty())
      continue;

    Data::cpeptide_t cpeptide;

    stringstream snc;
    snc << nc + 1;

    const string ssnc = snc.str();

    {
      ofstream splot;
      ofstream shtml;
      ContPlot cplot(cdata, splot, shtml, cpeptide);

      cplot.sopt[ContPlot::econtig] = ssnc;
      cplot.sopt[ContPlot::emp] = sopt[emp];
      cplot.sopt[ContPlot::emidx] = sopt[emidx];
      cplot.sopt[ContPlot::efasta] = sopt[efasta];
      cplot.sopt[ContPlot::erefmp] = sopt[erefmp];
      cplot.sopt[ContPlot::erefmidx] = sopt[erefmidx];
      cplot.sopt[ContPlot::epeakmasstol] = sopt[epeakmasstol];

      if (cplot())
      {
        cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;
        return -4;
      }
    }

    // stats
    unsigned nline = 0;
    double nave[14] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // iterate spectra
    for (unsigned i = 0; i < cdata.cabinfo[nc].first.first.size(); ++ i)
    {
      const unsigned ns = cdata.cabinfo[nc].first.first[i];

      ostringstream sns;
      sns << ns + 1;

      double nstat[14] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

      {
        ofstream splot;
        ofstream shtml;
        SpecPlot cplot(cdata, * pset[2], splot, shtml, cpeptide);

        cplot.sopt[SpecPlot::espectrum] = sns.str();
        cplot.sopt[SpecPlot::epeakmasstol] = sopt[epeakmasstol];
        cplot.bopt[SpecPlot::eshift] = true;
        cplot.bopt[SpecPlot::ebounds] = true;

        if (cplot())
        {
          cout << "error: see '" << sopt[eoutdir] + "/spsplot.log'" << endl;
          return -4;
        }
      }

      // unidentified
      if (cdata.cinspect.find(ns) == cdata.cinspect.end())
        continue;

      nstat[0] = cstat[2][ns].second[10] * 100. / cstat[2][ns].second[8]; // % explained intensity
      nstat[1] = cstat[2][ns].second[11] * 100. / cstat[2][ns].second[9]; // % explained peaks
      nstat[2] = cstat[2][ns].second[16] * 100. / cstat[2][ns].second[7]; // % observed b ions
      nstat[3] = cstat[2][ns].second[17] * 100. / cstat[2][ns].second[7]; // % observed y ions
      nstat[4] = cstat[2][ns].second[19] * 100. / cstat[2][ns].second[7]; // % observed neutral-loss
      nstat[5] = cstat[2][ns].second[18] * 100. / cstat[2][ns].second[7]; // % observed charge 2
      nstat[6] = cstat[2][ns].second[12] * 100. / cstat[2][ns].second[20]; // % intensity explained by b-ions
      nstat[7] = cstat[2][ns].second[13] * 100. / cstat[2][ns].second[20]; // % intensity explained by y-ions
      nstat[8] = cstat[2][ns].second[14] * 100. / cstat[2][ns].second[20]; // % intensity explained by neutral-loss ions
      nstat[9] = cstat[2][ns].second[15] * 100. / cstat[2][ns].second[20]; // % intensity explained by doubly-charged ions
      nstat[10] = cstat[2][ns].second[9]; // number of spectrum peaks
      nstat[11] = cstat[2][ns].second[3]; // mass error
      nstat[12] = (cstat[2][ns].second[23] + cstat[2][ns].second[24]) * 100. / cstat[2][ns].second[7]; // % observed b, y, b^2 and y^2 ions
      nstat[13] = cstat[2][ns].second[28]; // # consecutive breaks

      if (! std::isfinite(nstat[0])) nstat[0] = 0.;
      if (! std::isfinite(nstat[1])) nstat[1] = 0.;
      if (! std::isfinite(nstat[2])) nstat[2] = 0.;
      if (! std::isfinite(nstat[3])) nstat[3] = 0.;
      if (! std::isfinite(nstat[4])) nstat[4] = 0.;
      if (! std::isfinite(nstat[5])) nstat[5] = 0.;
      if (! std::isfinite(nstat[6])) nstat[6] = 0.;
      if (! std::isfinite(nstat[7])) nstat[7] = 0.;
      if (! std::isfinite(nstat[8])) nstat[8] = 0.;
      if (! std::isfinite(nstat[9])) nstat[9] = 0.;
      if (! std::isfinite(nstat[10])) nstat[10] = 0.;
      if (! std::isfinite(nstat[11])) nstat[11] = 0.;
      if (! std::isfinite(nstat[12])) nstat[12] = 0.;
      if (! std::isfinite(nstat[13])) nstat[13] = 0.;

      stxt[9] << ns + 1 << '\t'; // index
      stxt[9] << nstat[0] << '\t';
      stxt[9] << nstat[1] << '\t';
      stxt[9] << nstat[2] << '\t';
      stxt[9] << nstat[3] << '\t';
      stxt[9] << nstat[4] << '\t';
      stxt[9] << nstat[5] << '\t';
      stxt[9] << nstat[6] << '\t';
      stxt[9] << nstat[7] << '\t';
      stxt[9] << nstat[8] << '\t';
      stxt[9] << nstat[9] << '\t';
      stxt[9] << nstat[10] << '\t';
      stxt[9] << nstat[11] << '\t';
      stxt[9] << nstat[12] << '\t';
      stxt[9] << nstat[13] << endl;

      ++ nline;

      nave[0] += nstat[0];
      nave[1] += nstat[1];
      nave[2] += nstat[2];
      nave[3] += nstat[3];
      nave[4] += nstat[4];
      nave[5] += nstat[5];
      nave[6] += nstat[6];
      nave[7] += nstat[7];
      nave[8] += nstat[8];
      nave[9] += nstat[9];
      nave[10] += nstat[10];
      nave[11] += nstat[11];
      nave[12] += nstat[12];
      nave[13] += nstat[13];
    }

    sps::vector<double> ctype[] = {sps::vector<double>(cdata.cabinfo[nc].second.size()), sps::vector<double>(cdata.cabinfo[nc].second.size()), sps::vector<double>(cdata.cabinfo[nc].second.size())};

    // iterate vertices
    for (unsigned nv = 0; nv < cdata.cabinfo[nc].second.size(); ++ nv)
      // iterate a vertex
      for (unsigned nk = 0; nk < cdata.cabinfo[nc].second[nv].second.size(); ++ nk)
      {
        // global spec index
        const unsigned ns = cdata.cabinfo[nc].second[nv].first[nk];

        // vertex mass
        const double nm = cdata.cabinfo[nc].second[nv].second[nk];

        // undefined
        if (cdata.cinspect.find(ns) == cdata.cinspect.end())
          ctype[2][nv] = cdata.cconsset.specs[nc].peakList[nv][1];
        else
        {
          // incorrect
          ctype[1][nv] = cdata.cconsset.specs[nc].peakList[nv][1];

          for (unsigned j = 0; j < (* pset[2])[ns].size(); ++ j)
            if (uf(nm, nparentmasstol) == uf((* pset[2])[ns][j][0], nparentmasstol))
            {
              ctype[0][nv] = cdata.cconsset.specs[nc].peakList[nv][1];
              break;
            }
        }
      }

    unsigned ni[] = {0, 0, 0};
    double nd[] = {0., 0., 0., 0.};

    for (unsigned nv = 0; nv < cdata.cabinfo[nc].second.size(); ++ nv)
      if (ctype[0][nv] != 0.)
      {
        ++ ni[0];
        nd[0] += ctype[0][nv];
      }
      else if (ctype[1][nv] != 0.)
      {
        ++ ni[1];
        nd[1] += ctype[1][nv];
      }
      else
      {
        ++ ni[2];
        nd[2] += ctype[2][nv];
      }

    nd[3] = nd[0] + nd[1] + nd[2];

    nave[0] = nave[0] / nline;
    nave[1] = nave[1] / nline;
    nave[2] = nave[2] / nline;
    nave[3] = nave[3] / nline;
    nave[4] = nave[4] / nline;
    nave[5] = nave[5] / nline;
    nave[6] = nave[6] / nline;
    nave[7] = nave[7] / nline;
    nave[8] = nave[8] / nline;
    nave[9] = nave[9] / nline;
    nave[10] = nave[10] / nline;
    nave[11] = nave[11] / nline;
    nave[12] = nave[12] / nline;
    nave[13] = nave[13] / nline;

    if (! std::isfinite(nave[0])) nave[0] = 0.;
    if (! std::isfinite(nave[1])) nave[1] = 0.;
    if (! std::isfinite(nave[2])) nave[2] = 0.;
    if (! std::isfinite(nave[3])) nave[3] = 0.;
    if (! std::isfinite(nave[4])) nave[4] = 0.;
    if (! std::isfinite(nave[5])) nave[5] = 0.;
    if (! std::isfinite(nave[6])) nave[6] = 0.;
    if (! std::isfinite(nave[7])) nave[7] = 0.;
    if (! std::isfinite(nave[8])) nave[8] = 0.;
    if (! std::isfinite(nave[9])) nave[9] = 0.;
    if (! std::isfinite(nave[10])) nave[10] = 0.;
    if (! std::isfinite(nave[11])) nave[11] = 0.;
    if (! std::isfinite(nave[12])) nave[12] = 0.;
    if (! std::isfinite(nave[13])) nave[13] = 0.;

    stxt[10] << nc + 1 << '\t' << ni[0] << '\t' << ni[1] << '\t' << ni[2] << '\t' << nd[0] * 100. / nd[3] << '\t' << nd[1] * 100. / nd[3] << '\t' << nd[2] * 100. / nd[3] << '\t';
    stxt[10] << nave[0] << '\t';
    stxt[10] << nave[1] << '\t';
    stxt[10] << nave[2] << '\t';
    stxt[10] << nave[3] << '\t';
    stxt[10] << nave[4] << '\t';
    stxt[10] << nave[5] << '\t';
    stxt[10] << nave[6] << '\t';
    stxt[10] << nave[7] << '\t';
    stxt[10] << nave[8] << '\t';
    stxt[10] << nave[9] << '\t';
    stxt[10] << nave[10] << '\t';
    stxt[10] << nave[11] << '\t';
    stxt[10] << nave[12] << '\t';
    stxt[10] << nave[13] << '\n';
  }

  return 0;
}
