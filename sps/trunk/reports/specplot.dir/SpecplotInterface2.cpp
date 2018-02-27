///////////////////////////////////////////////////////////////////////////////
#include "SpecplotInterface2.h"
#include "mzxml.h"
#include "PWizInterface.h"
#include "Timer.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecplotInterface2::SpecplotInterface2()
{
}

SpecplotInterface2::~SpecplotInterface2()
{
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface2::load(string &fn, string &ext, bool spectrumScanPresent, bool spectrumIdPresent, bool usePwizFirst)
{
  int fileLoaded = 0;
  // load te file
  if(!specSet) {
    // create object
    specSet = new SpecSet();
    // start the timer
    Timer_c timer;

    // load using pwiz
    if(usePwizFirst && fileLoaded < 1 && (ext.compare("mzxml") || ext.compare("mzml"))) {
      if(m_verbose)
        cout << "Loading file using piwz....";
      // start the timer
      timer.start();
      PWizInterface pwiz;
      fileLoaded = pwiz.loadDataUsingPWiz(fn, *specSet, 2);
      if(m_verbose)
        if(fileLoaded < 1)
          cout << "exited in error." << endl;
        else
          cout << "loaded " << fileLoaded << " spectra; took " << timer.restart() << endl;
    }

    // load MGF, given and id or scan
    if(fileLoaded < 1 && !m_loadFull) {
      if(ext.compare("mgf") == 0) {
        if(m_verbose)
          cout << "Loading MGF file using id or scan....";
        int idx = (m_spectrumIndex > 0 ? m_spectrumIndex : -1);
        unsigned scan = getInt(m_spectrumScan.c_str());
        scan = (scan > 0 ? scan : 0);
        fileLoaded = specSet->LoadSpectrum_mgf(fn.c_str(), scan, idx, m_useindex);
        m_origspectrumIndex = m_spectrumIndex;
        m_spectrumIndex = (m_spectrumIndex > 0 ? 1 : 0);
        if(m_verbose)
          if(fileLoaded < 1)
            cout << "exited in error." << endl;
          else
            cout << "loaded " << fileLoaded << " spectra; took " << timer.restart() << endl;
      }
    }

    //int sz = specSet->size();
    //cout << "Read specset has " << sz << endl;
    //if(sz) {
    //  cout << "  Spectrum scan is = " << (*specSet)[0].scan << endl;
    //}


    // load mzxml file. Load using first the SPS mzxml loader. Use the appropriate method in case scan # is specified
    if(fileLoaded < 1) {
      if(ext.compare("mzxml") == 0) {
        if(m_verbose)
          cout << "Loading MZXML file....";
        // start the timer
        timer.start();
        // auxilizary vector needed
        vector<short> msLevel;
        // hold the return value
        int ret;
        // load method depends on the presence of spectrumscan specification
        try {
          if(m_spectrumScan.empty()) {
            fileLoaded = LoadMzxml( (char * const)(fn.c_str()), fn, *specSet, & msLevel, 2);
          } else {
            if(m_verbose)
              cout << "using spectrumscan fast parsing....";
            fileLoaded = LoadMzxml(fn.c_str(), fn, *specSet, m_spectrumScan.c_str(), & msLevel, 2);
          }
        } catch (...) {
          if(m_verbose)
             cout << "exited in error." << endl;
          fileLoaded = 0;
        }
        if(m_verbose)
          if(fileLoaded < 1)
            cout << "loaded " << fileLoaded << " spectra; took " << timer.restart() << endl;
      }
    }

    // If not file was loaded at this point, load it using the specsset load methods.
    if(fileLoaded < 1) {
      if(m_verbose)
        cout << "Loading specset file ....";
      timer.start();
      // load other formats
      fileLoaded = specSet->Load(fn.c_str(), NULL);
      if(m_verbose)
        if(fileLoaded < 1)
          cout << "exited in error." << endl;
        else
          cout << "loaded " << fileLoaded << " spectra; took " << timer.restart() << endl;
    }

    // load using pwiz - final attempt
    if(fileLoaded < 1) {
      timer.start();
      if(m_verbose)
        cout << "Loading using pwiz ....";
      PWizInterface pwiz;
      fileLoaded = pwiz.loadDataUsingPWiz(fn, *specSet, 2);
      if(m_verbose)
        if(fileLoaded < 1)
          cout << "exited in error." << endl;
        else
          cout << "loaded " << fileLoaded << " spectra; took " << timer.restart() << endl;
    }

  }

  // if no file was loaded at this point, fail
  if(fileLoaded < 1) {
    stringstream err;
    err << "Error loading file: " << fn;
    return error(err.str());
  }

  // Check if spectrum ID is present. If it is, find the corresponding spectrum index
  if(spectrumIdPresent) {
    bool found = false;
    for(int i = 0 ; i < specSet->size() ; i++){
      if(  (*specSet)[i].psmList.size() == 1 ){
        if ( (*specSet)[i].psmList.front()->m_spectrumID == m_spectrumID) {
          m_spectrumIndex = i+1;
          plotSpectrum.setSpectrumIndex(m_spectrumIndex);
          found = true;
        }
      }
    }
    // in case we didn't find it, exit in error
    if(!found) {
        stringstream err; err << "ERROR: Spectrum ID " << m_spectrumID << " not found.";
        return error(err.str());
    }
  }

  // If spectrumscan was specified,
  if(spectrumScanPresent) {
    // get spectrumscan from index
    bool found = false;
    Timer_c timer;
    if(m_verbose)
      cout << "Looking for scan # " << m_spectrumScan << " in " << specSet->size() << " spectra" << endl;
    for(int i = 0 ; i < specSet->size() ; i++){
      if(m_verbose)
        cout << "  Scan # at index " << i << " : " << (*specSet)[i].scan << endl;
      if((*specSet)[i].scan == getInt(m_spectrumScan.c_str())) {
        m_spectrumIndex = i+1;
        plotSpectrum.setSpectrumIndex(m_spectrumIndex);
        found = true;
        if(m_verbose)
          cout << "Found at index (0-based) : " << i << " -- Took " << timer.restart() << endl;
        break;
      }
    }
    // in case we didn't find it, exit in error
    if(!found) {
      stringstream err; err << "ERROR: Spectrumscan not found.";
      return error(err.str());
    }
  }

  if(fileLoaded)
    m_inputSpectraFilename = fn;

  return fileLoaded;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
