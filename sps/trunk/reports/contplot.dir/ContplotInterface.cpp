///////////////////////////////////////////////////////////////////////////////
#include "ContplotInterface.h"
#include "PlotContig.h"

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "aminoacid.h"
#include "Logger.h"
//#include "mzxml.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
ContplotInterface::ContplotInterface() :
  abinfo(NULL), abinfo_own(true),
  star(NULL), star_own(true),
  seqs(NULL), seqs_own(true)
{
  // Sets draw object in plotContig as PlotGnu
  plotContig.setDrawMethod(RENDERER_TYPE_GNU);
}

ContplotInterface::~ContplotInterface()
{
  if(abinfo_own && abinfo)
    delete abinfo;
  if(star_own && star)
    delete star;
  if(seqs_own && seqs)
    delete seqs;
}
///////////////////////////////////////////////////////////////////////////////
void ContplotInterface::setData(int type, void *data)
{
  switch(type) {
  case 1:
    if(abinfo) return;
    abinfo = (abinfo_t *)data;
    abinfo_own = false;
    plotContig.setAbinfo(abinfo);
    break;
  case 2:
    if(star) return;
    star = new SpecSet();
    *star = (*(SpecSet *)data);
    star_own = true;
    plotContig.setStar(star);
    break;
  case 3:
    if(seqs) return;
    seqs = (SpecSet *)data;
    seqs_own = false;
    plotContig.setSeqs(seqs);
    break;
  }
}
///////////////////////////////////////////////////////////////////////////////
int ContplotInterface::plot(void)
{
  // set the spectrum to plot
  //plotContig.setSpectrum(&specSet[m_spectrumIndex - 1]);

  // Draws the image
  plotContig.draw();

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ContplotInterface::processOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------
  // Get the execultable directory. Unless specified, renderer directory and font directory is expected
  string str = argv[0];

  //convert path
  size_t found;
  found = str.find_last_of("/\\");
  string aux;
  aux = '.';
  if(found != string::npos)
    aux = str.substr(0, found);
  string exePath = aux;

//#if  defined(__MINGW32__) || defined(__CYGWIN__)
  exePath = replaceAll(exePath, "\\", "/");
//#endif

  // set renderer location
  plotContig.setRendererLocation(aux);
  // Font location also defaults to executale directory
  plotContig.setFontLocation(exePath);
  // Annotation file directory also defaults to executable directory
  plotContig.setAnnotationModelDirectory(aux);
  // Amino acid file directory also defaults to executable directory
  plotContig.setAminoAcidsFileDirectory(aux);

  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
//  listOptions.push_back(CommandLineParser::Option("-verbose",         "verbose",          false));

  listOptions.push_back(CommandLineParser::Option("-notitle",         "NOTITLE",          false));
  listOptions.push_back(CommandLineParser::Option("-title",           "TITLE",            true));
  listOptions.push_back(CommandLineParser::Option("-shift-value",     "SHIFT_VALUE",      true));

  listOptions.push_back(CommandLineParser::Option("-aa",              "AMINO_ACID_LIST",  true));
  listOptions.push_back(CommandLineParser::Option("-annotation-model","ANNOTATION_MODEL", true));
  listOptions.push_back(CommandLineParser::Option("-renderer-dir",    "RENDERER_LOCATION",true));
  listOptions.push_back(CommandLineParser::Option("-font-dir",        "FONT_LOCATION",    true));

  listOptions.push_back(CommandLineParser::Option("-star",            "STAR",             true));
  listOptions.push_back(CommandLineParser::Option("-abruijn",         "ABRUIJN",          true));
  listOptions.push_back(CommandLineParser::Option("-seqs",            "SEQS",             true));

  listOptions.push_back(CommandLineParser::Option("-contig",          "CONTIG",           true));
  listOptions.push_back(CommandLineParser::Option("-reference",       "REFERENCE",        true));
  listOptions.push_back(CommandLineParser::Option("-homolog",         "HOMOLOG",          true));
  listOptions.push_back(CommandLineParser::Option("-denovo",          "DENOVO",           true));
  listOptions.push_back(CommandLineParser::Option("-user",            "USER",             true));

  listOptions.push_back(CommandLineParser::Option("-mass-reference",  "MASS_REFERENCE",   true));
  listOptions.push_back(CommandLineParser::Option("-mass-homolog",    "MASS_HOMOLOG",     true));
  listOptions.push_back(CommandLineParser::Option("-mass-user",       "MASS_USER",        true));

  listOptions.push_back(CommandLineParser::Option("-offset-reference","OFFSET_REFERENCE", true));
  listOptions.push_back(CommandLineParser::Option("-offset-homolog",  "OFFSET_HOMOLOG",   true));
  listOptions.push_back(CommandLineParser::Option("-offset-user",     "OFFSET_USER",      true));

  listOptions.push_back(CommandLineParser::Option("-reverse",         "REVERSE",          false));


  listOptions.push_back(CommandLineParser::Option("-annotation-style-inspect","ANNOT_STYLE_INSPECT",       false));

  listOptions.push_back(CommandLineParser::Option("-peakmasstol",     "TOLERANCE_PEAK",   true));
  listOptions.push_back(CommandLineParser::Option("-individual-peakmasstol","INDIVIDUAL_TOLERANCE_PEAK",   false));
  listOptions.push_back(CommandLineParser::Option("-parentmasstol",   "TOLERANCE_PM",     true));

  listOptions.push_back(CommandLineParser::Option("-outdir",          "OUTDIR",           true));
  listOptions.push_back(CommandLineParser::Option("-prefix",          "PREFIX",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",         "OUTFILE",          true));
  listOptions.push_back(CommandLineParser::Option("-format",          "FORMAT",           true));
  listOptions.push_back(CommandLineParser::Option("-encoding",        "ENCODING",         true));
  listOptions.push_back(CommandLineParser::Option("-target",          "TARGET",           true));

  listOptions.push_back(CommandLineParser::Option("-zoom",            "ZOOM",             true));
  listOptions.push_back(CommandLineParser::Option("-image-width",     "IMAGE_WIDTH",      true));
  listOptions.push_back(CommandLineParser::Option("-image-height",    "IMAGE_HEIGHT",     true));
  listOptions.push_back(CommandLineParser::Option("-image-stretch-height", "IMAGE_STRETCH_HEIGHT",     true));
  listOptions.push_back(CommandLineParser::Option("-image-stretch-width",  "IMAGE_STRETCH_WIDTH",      true));

  listOptions.push_back(CommandLineParser::Option("-range",           "RANGE",            true));

  listOptions.push_back(CommandLineParser::Option("-request-id",      "REQ_ID",           true));

  listOptions.push_back(CommandLineParser::Option("-gnuplot-file",    "GNUPLOT_FILE",     false));

  listOptions.push_back(CommandLineParser::Option("-eps-correction",    "EPS_CORRECTION",     false));
  listOptions.push_back(CommandLineParser::Option("-minimum-intensity", "MINIMUM_INTENSITY",  false));
  listOptions.push_back(CommandLineParser::Option("-peak-thickness",    "PEAK_THICKNESS",     true));



  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",               "PARAMETER_FILE",   true));


  ////////////////////////////////////////////////////////////////////////////////
  // Execute the command line parser
  CommandLineParser clp(argc, argv, 0, listOptions);

  ////////////////////////////////////////////////////////////////////////////////
  // Checking for errors
  string parser_error;
  if (!clp.validate(parser_error)) {
    ERROR_MSG(parser_error);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // The parameters' values
  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  return processOptions(commandLineParams);
}
////////////////////////////////////////////////////////////////////////////////
int ContplotInterface::processOptions(ParameterList &ip)
{
  ParameterList commandLineParams;

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (ip.exists("PARAMETER_FILE")) {

    string parameterFilename = ip.getValue("PARAMETER_FILE");

    commandLineParams.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag set)
  }

  commandLineParams.addList(ip, true);

  ////////////////////////////////////////////////////////////////////////////////
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // Boolean parameters

  if (commandLineParams.exists("shift"))
    plotContig.setMassShift(0.0);

  if (commandLineParams.exists("SHIFT_VALUE"))
    plotContig.setMassShift(getFloat(commandLineParams.getValue("SHIFT_VALUE").c_str()));

  if (commandLineParams.exists("NOTITLE"))
    plotContig.setTitlePresence(0);

  if (commandLineParams.exists("ANNOT_STYLE_INSPECT"))
    plotContig.setAnnotatinStyle(ANNOTATION_STYLE_INSPECT);

  if (commandLineParams.exists("REVERSE"))
    plotContig.setReverseFlag();


  ////////////////////////////////////////////////////////////////////////////////
  // File locations

  // Renderer location
  if (commandLineParams.exists("RENDERER_LOCATION")) {
    string aux = commandLineParams.getValue("RENDERER_LOCATION");
    aux = replaceAll(aux, "\\", "/");
    plotContig.setRendererLocation(aux);
  }

  // Font file location
  if (commandLineParams.exists("FONT_LOCATION")) {
    string aux = commandLineParams.getValue("FONT_LOCATION");
    aux = replaceAll(aux, "\\", "/");
    plotContig.setFontLocation(aux);
  }


  // Amino-acid masses file
  if (commandLineParams.exists("AMINO_ACID_LIST")) {
    string aux = commandLineParams.getValue("AMINO_ACID_LIST");
    aux = replaceAll(aux, "\\", "/");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      plotContig.setAminoAcidsFile(aux.substr(found + 1));
      plotContig.setAminoAcidsFileDirectory(aux.substr(0, found));
    } else {
      plotContig.setAminoAcidsFile(aux);
      plotContig.setAminoAcidsFileDirectory(".");
    }
  }

  // Annotation model
  if (commandLineParams.exists("ANNOTATION_MODEL")) {
    string aux = commandLineParams.getValue("ANNOTATION_MODEL");
    aux = replaceAll(aux, "\\", "/");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      plotContig.setAnnotationModel(aux.substr(found + 1));
      plotContig.setAnnotationModelDirectory(aux.substr(0, found));
    } else {
      plotContig.setAnnotationModel(aux);
      plotContig.setAnnotationModelDirectory(".");
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Peptide and amino acid parameters

  // contig index
  if (commandLineParams.exists("CONTIG")) {
    int ctig = getInt(commandLineParams.getValue("CONTIG").c_str());
    plotContig.setContigIndex(ctig - 1); // contig is 1-based
    stringstream title;
    title << "Contig " << ctig;
    plotContig.setTitle(title.str());
    //cout << "--------------------- Contig " << ctig << "--------------------- " << endl;
  } else {
    return error("contig index not specified.");
  }

  // reference sequence
  if (commandLineParams.exists("REFERENCE")) {
    string reference = commandLineParams.getValue("REFERENCE");
    transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
    plotContig.setReference(reference);
    //cout << "Reference: " << reference << endl;
  }

  // homolog sequence
  if (commandLineParams.exists("HOMOLOG")) {
    string homolog = commandLineParams.getValue("HOMOLOG");
    transform(homolog.begin(), homolog.end(), homolog.begin(), ::toupper);
    plotContig.setHomolog(homolog);
    //cout << "Homolog: " << homolog << endl;
  }

  // denovo sequence
  if (commandLineParams.exists("DENOVO")) {
    string denovo = commandLineParams.getValue("DENOVO");
    transform(denovo.begin(), denovo.end(), denovo.begin(), ::toupper);
    plotContig.setDenovo(denovo);
    //cout << "deNovo: " << denovo << endl;
  }

  // user sequence
  if (commandLineParams.exists("USER")) {
    string user = commandLineParams.getValue("USER");
    transform(user.begin(), user.end(), user.begin(), ::toupper);
    plotContig.setUser(user);
    //cout << "User: " << user << endl;
  }



  // reference mass intervals
  if (commandLineParams.exists("MASS_REFERENCE")) {
    string reference = commandLineParams.getValue("MASS_REFERENCE");
    plotContig.setReferenceMass(reference);
    //cout << "Reference masses: " << reference << endl;
  }

  // homolog mass intervals
  if (commandLineParams.exists("MASS_HOMOLOG")) {
    string homolog = commandLineParams.getValue("MASS_HOMOLOG");
    plotContig.setHomologMass(homolog);
    //cout << "Homolog masses: " << homolog << endl;
  }

  // user mass intervals
  if (commandLineParams.exists("MASS_USER")) {
    string user = commandLineParams.getValue("MASS_USER");
    plotContig.setUserMass(user);
    //cout << "User masses: " << user << endl;
  }



  // reference offset
  if (commandLineParams.exists("OFFSET_REFERENCE")) {
    string reference = commandLineParams.getValue("OFFSET_REFERENCE");
    plotContig.setReferenceOffset(getFloat(reference.c_str()));
  }

  // homolog offset
  if (commandLineParams.exists("OFFSET_HOMOLOG")) {
    string homolog = commandLineParams.getValue("OFFSET_HOMOLOG");
    plotContig.setHomologOffset(getFloat(homolog.c_str()));
  }

  // user offset
  if (commandLineParams.exists("OFFSET_USER")) {
    string user = commandLineParams.getValue("OFFSET_USER");
    plotContig.setUserOffset(getFloat(user.c_str()));
  }


  ////////////////////////////////////////////////////////////////////////////////
  // File load section

  // Abruijn graph
  if(!abinfo)
    if (commandLineParams.exists("ABRUIJN")) {
      abinfo = new abinfo_t();
      string fn = commandLineParams.getValue("ABRUIJN");
      if(Load_abinfo(fn.c_str(), *abinfo) == 0) {
        cerr << "Error loading Abruijn graph." << endl;
        return -1;
      }
      plotContig.setAbinfo(abinfo);
    } else {
      return error("Abruijn graph not specified.");
    }

  // Star spectra file
  if(!star)
    if (commandLineParams.exists("STAR")) {
      star = new SpecSet();
      string starFilename = commandLineParams.getValue("STAR");
      star->loadPklBin(starFilename.c_str());
      plotContig.setStar(star);
    } else {
      return error("Star spectra file not specified.");
    }

  // seqs file
  if(!seqs)
    if (commandLineParams.exists("SEQS")) {
      seqs = new SpecSet();
      string starFilename = commandLineParams.getValue("SEQS");
      seqs->loadPklBin(starFilename.c_str());
      plotContig.setSeqs(seqs);
    } else {
      return error("Sps seqs file not specified.");
    }

  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("OUTDIR")) {
    string aux = commandLineParams.getValue("OUTDIR");
    aux = replaceAll(aux, "\\", "/");
    plotContig.setFileOutDir(aux);
  }

  // Output file name.
  if (commandLineParams.exists("OUTFILE")) {
    string aux = commandLineParams.getValue("OUTFILE");
    aux = replaceAll(aux, "\\", "/");
    plotContig.setFileName(aux);
  }

  // Filename prefix
  if (commandLineParams.exists("PREFIX"))
   plotContig.setFilePrefix(commandLineParams.getValue("PREFIX"));

  // Output file format: png, ...
  if (commandLineParams.exists("FORMAT"))
    plotContig.setFileFormat(commandLineParams.getValue("FORMAT"));

  // Output file encoding: uu64, ...
  if (commandLineParams.exists("ENCODING"))
    plotContig.setEncoding(commandLineParams.getValue("ENCODING"));

  // Output target: file, cout, cout, internal
  if (commandLineParams.exists("TARGET"))
    plotContig.setTarget(commandLineParams.getValue("TARGET"));


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters

  // peak mass tolerance
  if (commandLineParams.exists("TOLERANCE_PEAK"))
    plotContig.setPeakMassTol(getFloat(commandLineParams.getValue("TOLERANCE_PEAK").c_str()));

  if (commandLineParams.exists("INDIVIDUAL_TOLERANCE_PEAK"))
    plotContig.useIndividualPeakMassTol(true);

  ////////////////////////////////////////////////////////////////////////////////
  //  Image related parameters

  // Title
  if (commandLineParams.exists("TITLE"))
    plotContig.setTitle(commandLineParams.getValue("TITLE"));

  // Output image zoom factor
  if (commandLineParams.exists("ZOOM")) {
    plotContig.setZoom(getFloat(commandLineParams.getValue("ZOOM").c_str()));
    //cout << "Zoom: " << commandLineParams.getValue("ZOOM") << endl;
  }

  // Output image width, in pixels
  if (commandLineParams.exists("IMAGE_WIDTH"))
    plotContig.setImageWidth(getInt(commandLineParams.getValue("IMAGE_WIDTH").c_str()));

  // Output image height, in pixels
  if (commandLineParams.exists("IMAGE_HEIGHT"))
    plotContig.setImageHeight(getInt(commandLineParams.getValue("IMAGE_HEIGHT").c_str()));

  // Output image stretch width, in pixels
  if (commandLineParams.exists("IMAGE_STRETCH_HEIGHT"))
    plotContig.setImageStretchHeight(getInt(commandLineParams.getValue("IMAGE_STRETCH_HEIGHT").c_str()));

  // Output image stretch height, in pixels
  if (commandLineParams.exists("IMAGE_STRETCH_WIDTH"))
    plotContig.setImageStretchWidth(getInt(commandLineParams.getValue("IMAGE_STRETCH_WIDTH").c_str()));


/*
  // m/z  axis range.
  if (commandLineParams.exists("RANGE")) {
    string range = commandLineParams.getValue("RANGE");
    vector<string> aux;
    stringSplit(range, aux, ":");
    float faux;
    switch(aux.size()) {
    case 1:
      faux = getFloat(aux[0].c_str());
      if(range[0] ==':')
        plotContig.setRangeMax(faux);
      else
        plotContig.setRangeMin(faux);
      break;
    case 2:
      float rangeMin = getFloat(aux[0].c_str());
      float rangeMax = getFloat(aux[1].c_str());
      plotContig.setRangeMin(rangeMin);
      plotContig.setRangeMax(rangeMax);
      if(rangeMin >= rangeMax)
        return error("range minimun value greater or egual than range maximun value.");
      break;
    }
  }
*/

  //////////////////////////////////////////////////////////////////////////////
  // Other

  // keep gnuplot file
  if (commandLineParams.exists("GNUPLOT_FILE"))
    plotContig.setDebug(1);

  // eps file correction
  if (commandLineParams.exists("EPS_CORRECTION"))
    plotContig.setEpsCorrection(true);

  // eps file correction
  if (commandLineParams.exists("MINIMUM_INTENSITY"))
    plotContig.useMinimunIntentity(true);

  // peak line thickness
  if (commandLineParams.exists("PEAK_THICKNESS"))
    plotContig.setLineThickness(getInt(commandLineParams.getValue("PEAK_THICKNESS").c_str()));



  //////////////////////////////////////////////////////////////////////////////
  // Execute

  plotContig.draw();

  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int ContplotInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: specplot [OPTION]\n";
  sout << "Options:\n";

  sout << "  --abruijn FILE                 Abruijn graph file\n";
  sout << "  --star FILE                    Stars file\n";
  sout << "  --seqs FILE                    Seqs file\n";
  sout << '\n';
  sout << "  --font-dir DIRECTORY           Font files location in the file system.\n";
  sout << "  --renderer-dir DIRECTORY       Renderer location in the file system.\n";
  sout << '\n';
  sout << "  --contig INDEX                 Contig (X) 1-based\n";
  sout << "  --reverse                      Reversed contig\n";
  sout << "  --reference SEQUENCE           Reference sequence\n";
  sout << "  --homolog SEQUENCE             Homolog sequence\n";
  sout << "  --denovo SEQUENCE              Denovo sequence\n";
  sout << "  --user SEQUENCE                User sequence\n";
  sout << '\n';
  sout << "  --mass-reference VALUE LIST    Reference sequence mass intervals, separated by the '|' character\n";
  sout << "  --mass-homolog VALUE LIST      Homolog sequence mass intervals, separated by the '|' character\n";
//  sout << "  --mass-user VALUE LIST       User sequence mass intervals, separated by the '|' character\n";
  sout << '\n';
  sout << "  --offset-reference VALUE       Reference sequence mass offset\n";
  sout << "  --offset-homolog VALUE         Homolog sequence mass offset\n";
  sout << "  --offset-user VALUE            User sequence mass offset\n";
  sout << '\n';
  sout << "  --peakmasstol NUMBER           Amino acid peak mass tolerance.\n";
  sout << "  --individual-peakmasstol       If specified, individual peak mass tolerance is used.\n";
  sout << "  --parentmasstol NUMBER         Parent mass tolerance\n";
  sout << "  --annotation-style-inspect     Set annotation style to INSPECT. Default is SPECNETS\n";
  sout << "  --notitle                      Disables title generation\n";
  sout << "  --title STRING                 Specifies title\n";
  sout << "  --zoom NUMBER                  Zoom factor\n";
  sout << "  --image-width NUMBER           Image width, in pixels\n";
  sout << "  --image-height NUMBER          Image height, in pixels\n";
  sout << "  --image-stretch-width NUMBER   Image width limit, in pixels. -1 for unlimited stretch.\n";
  sout << "  --image-stretch-height NUMBER  Image height limit, in pixels. -1 for unlimited stretch.\n";
//  sout << "  --range X:Y                  Zoom spectrum to m/z range between m/z=X and m/z=Y\n";
  sout << '\n';
  sout << "  --format STRING                Image output format (png or eps)\n";
  sout << "  --outdir PATH                  Output directory\n";
  sout << "  --outfile FILE                 Output filename\n";
  sout << "  --target STRING                Output target (file, cout, internal)\n";
  sout << "  --encoding STRING              Output encoding (uu64, default is no encoding)\n";
  sout << '\n';
  sout << "  --eps-correction               Use shift correction for eps files\n";
  sout << "  --minimum-intensity            Use minimum intensity values to draw the peaks\n";
  sout << "  --peak-thickness NUMBER        Specify peak line thickness (Default is 4)\n";
  sout << '\n';
  sout << "  --verbose                      Output additional info\n";
  sout << '\n';
  sout << "  --p FILE                       Read parameters from file\n";
  sout << '\n';
  sout << "  --help                         Display this help and exit\n";
  sout << "  --version                      Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ContplotInterface::version(ostream & sout)
{
  sout << "contplot 2.0." << XSTR(SPS_VERSION) << '\n';
  sout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
  sout << "SH1: " <<  XSTR(GIT_SH1) << endl;
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ContplotInterface::error(const string & a)
{
  cerr << "contplot: invalid or missing option: " << a << endl;
  cerr << "Type 'contplot --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
