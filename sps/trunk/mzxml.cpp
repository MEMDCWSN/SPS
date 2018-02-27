#include "mzxml.h"
//#include "regex2.h"
#include "Base64.h"

#if defined(__MINGW32__)
#include <Winsock2.h>
#include <Winsock.h>
#else
#include <netinet/in.h>
#endif

#ifndef HAVE_U_INT32_T
typedef uint32_t u_int32_t;
#define HAVE_U_INT32_T
#endif

#include <cmath>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

using namespace specnets;

XERCES_CPP_NAMESPACE_USE

short mzxmlSAX2Handler::STATE_SKIPPING = 0, mzxmlSAX2Handler::STATE_SCAN = 1,
    mzxmlSAX2Handler::STATE_PRECURSORMZ = 2, mzxmlSAX2Handler::STATE_PEAKS = 3;

string mzxmlSAX2Handler::DEFAULT_INSTRUMENT_ID = "default_mass_analyzer";

typedef union
{
  uint32_t u32;
  float flt;
} U32;

using namespace std;

//unsigned int LoadMzxml(char *filename, vector<Spectrum> &specs, vector<unsigned int> &scanNums,
//                  vector<short> *msLevels, short minMsLevel) {
template<typename T>
unsigned int LoadMzxml(const T & xmlFile,
                       const string& filename,
                       SpecSet &specs,
                       vector<short> *msLevels,
                       short minMsLevel)
{
  specs.resize(0);

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Error during initialization! :\n";
    cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return 1;
  }

  SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
  parser->setFeature(XMLUni::fgSAX2CoreValidation, false);   // optional
  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);   // optional
  parser->setFeature(XMLUni::fgXercesSchema, false);   // optional

//    DefaultHandler* defaultHandler = new DefaultHandler();
  mzxmlSAX2Handler* defaultHandler = new mzxmlSAX2Handler(specs,
                                                          filename,
                                                          msLevels,
                                                          minMsLevel);
  parser->setContentHandler(defaultHandler);
  parser->setErrorHandler(defaultHandler);

  try
  {
    parser->parse(xmlFile);
  }
  catch (const XMLException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return 0;
  }
  catch (const SAXParseException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return 0;
  }
  catch (const mzxmlEx e)
  {
    cout << "Exception: " << e.description << "\n";
    return 0;
  }
  catch (...)
  {
    cout << "Unexpected Exception \n";
    return 0;
  }

  if (defaultHandler->numSpecs < specs.size())
  {
    specs.resize(defaultHandler->numSpecs);
    if (msLevels)
      msLevels->resize(defaultHandler->numSpecs);
  }

  delete parser;
  delete defaultHandler;

  return specs.size();
}

/**
 Explicit template instanciations
 */

template
unsigned int LoadMzxml(char * const & xmlFile,
                       const string& filename,
                       SpecSet &specs,
                       vector<short> *msLevels,
                       short minMsLevel);

template
unsigned int LoadMzxml(const char * const & xmlFile,
                       const string& filename,
                       SpecSet &specs,
                       vector<short> *msLevels,
                       short minMsLevel);

template
unsigned int LoadMzxml(const MemBufInputSource & xmlFile,
                       const string& filename,
                       SpecSet &specs,
                       vector<short> *msLevels,
                       short minMsLevel);

namespace
{

  /**
   Grep utility.

   @note The text must be shorter than \code sizeof(sbuffer) \endcode
   */

  streampos grep(ifstream & sin, const char * pc)
  {
    for (char sbuffer[80]; sin.peek() != EOF;)
    {
      const streampos np = sin.tellg();

      sin.getline(sbuffer, sizeof(sbuffer));

      if (sin.fail())
      {
        sin.clear();
        sin.ignore(numeric_limits<int>::max(), '\n');
      }

      if (strstr(sbuffer, pc))
        return np;
    }

    return sin.tellg();
  }

  /**
   Egrep utility.

   @note The text must be shorter than \code sizeof(sbuffer) \endcode
   */

  streampos egrep2(ifstream & sin)
  {
    for (char sbuffer[80]; sin.peek() != EOF;)
    {
      const streampos np = sin.tellg();
      sin.getline(sbuffer, sizeof(sbuffer));
      if (sin.fail())
      {
        sin.clear();
        sin.ignore(numeric_limits<int>::max(), '\n');
      }

      if (strstr(sbuffer, "</ scan"))
        return np;
      if (strstr(sbuffer, "</scan"))
        return np;
      if (strstr(sbuffer, "<scan"))
        return (streampos)((int)np);
    }

    return sin.tellg();
  }

/*
 streampos egrep(ifstream & sin, const sps::regex & pc)
 {
 for (char sbuffer[80]; sin.peek() != EOF; )
 {
 const streampos np = sin.tellg();

 sin.getline(sbuffer, sizeof(sbuffer));

 if (sin.fail())
 {
 sin.clear();
 sin.ignore(numeric_limits<int>::max(), '\n');
 }

 if (pc == sbuffer)
 return np;
 }

 return sin.tellg();
 } */

} // namespace

/**
 Fast spectrum load.

 Used to retrieve only 1 spectrum out of a large mzXML file.

 @param  filename    Name of file to fetch from
 @param  specs       Loaded spectrum
 @param  sscan       Scan number (string format for convenience)
 */

unsigned int LoadMzxml(const char *filename,
                       const string& filenameStr,
                       SpecSet &specs,
                       const string & sscan,
                       vector<short> *msLevels,
                       short minMsLevel)
{
  ifstream smzxml(filename);

  //static const sps::regex ctag("</?scan");
  //streampos np[] = {grep(smzxml, ("<scan num=\"" + sscan).c_str()), egrep(smzxml, ctag)};
  streampos np[] = { grep(smzxml, ("<scan num=\"" + sscan).c_str()),
                     egrep2(smzxml) };

  char sbuffer[np[1] - np[0]];

  //cout << np[0] << endl;
  //cout << np[1] << endl;

  smzxml.seekg(np[0]);
  smzxml.read(sbuffer, sizeof(sbuffer));
  smzxml.close();

  ostringstream sxml;
  sxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
  sxml << "<mzXML>\n";
  sxml << " <msRun scanCount=\"1\">\n";
  sxml.write(sbuffer, sizeof(sbuffer));
  sxml << "  </scan>\n";
  sxml << " </msRun>\n";
  sxml << "</mzXML>\n";

  const string ssxml = sxml.str();

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    cout << "Error during initialization! :\n";
    cout << "Exception message is: \n" << message << "\n";
    XMLString::release(&message);
    return 1;
  }

  MemBufInputSource cmzxml((const XMLByte *)(ssxml.c_str()),
                           ssxml.size(),
                           __PRETTY_FUNCTION__);

  unsigned int res = LoadMzxml(cmzxml,
                               filenameStr,
                               specs,
                               msLevels,
                               minMsLevel);
  if (res > 0)
  {
    FilenameManager mngr(filenameStr);
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      if (specs[i].fileName.length() == 0)
      {
        specs[i].fileName = mngr.getFilenameWithExtension();
      }
    }
  }
  return res;
}

void mzxmlSAX2Handler::startElement(const XMLCh* const uri,
                                    const XMLCh* const localname,
                                    const XMLCh* const qname,
                                    const Attributes& attrs)
{

  char* name = XMLString::transcode(localname);

  char *str;
  XMLCh *xmlStr;
  int intValue;
  float floatValue;
  if (!strcmp(name, "msRun"))
  {
    xmlStr = XMLString::transcode("scanCount");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    intValue = atoi(str);
    if (intValue >= 0)
    {
      specs->resize(intValue);
      if (msLevels)
        msLevels->resize(intValue);
    }
    XMLString::release(&str);
    XMLString::release(&xmlStr);
  }

  if (!strcmp(name, "msInstrument"))
  {
    xmlStr = XMLString::transcode("msInstrumentID");
    if (attrs.getValue(xmlStr) != NULL)
    {
      str = XMLString::transcode(attrs.getValue(xmlStr));
      curInstrumentId = str;
      XMLString::release(&str);
    }
    else
    {
      curInstrumentId = DEFAULT_INSTRUMENT_ID;
    }
    XMLString::release(&xmlStr);
  }

  if (!strcmp(name, "msMassAnalyzer") && curInstrumentId.length() > 0)
  {
    xmlStr = XMLString::transcode("value");
    if (attrs.getValue(xmlStr) != NULL)
    {
      str = XMLString::transcode(attrs.getValue(xmlStr));
      instrumentIdMap[curInstrumentId] = str;

      // THIS LINE WAS BREAKING DOWN SPECPLOT AS A SERVER. IN THE FUTURE, USE cerr INSTEAD OF cout FOR CODE DEBUGGING.
      //cout << "Instrument map: " << curInstrumentId << " -> " << str << endl;
      curInstrumentId = "";
      XMLString::release(&str);
    }
    XMLString::release(&xmlStr);
  }

  if (!strcmp(name, "scan"))
  {
    if (specs->size() == 0)
    {
      stringstream aux;
      aux
          << "ERROR: Can't read scans before knowing how many scans in the file (from msRun.scanCount)!\n";
      throw mzxmlEx(1, aux.str().c_str());
    }
    if (numSpecs >= specs->size())
    {
      stringstream aux;
      aux << "ERROR: Can't read more scans than indicated in msRun.scanCount! "
          << numSpecs << " : " << specs->size() << "\n";
      throw mzxmlEx(1, aux.str().c_str());
    }

    xmlStr = XMLString::transcode("msLevel");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    intValue = atoi(str);
    if (intValue >= minMsLevel)
    {
      curSpec.push_front(numSpecs++);
      if (msLevels)
        (*msLevels).at(curSpec.front()) = max((unsigned int)intValue,
                                              (unsigned int)curSpec.size());
      (*specs)[curSpec.front()].msLevel = max((unsigned int)intValue,
                                              (unsigned int)curSpec.size());
      ;
      state = STATE_SCAN;
    }
    else
      state = STATE_SKIPPING;
    XMLString::release(&str);
    XMLString::release(&xmlStr);

    if (state == STATE_SCAN)
    {
      xmlStr = XMLString::transcode("num");
      str = XMLString::transcode(attrs.getValue(xmlStr));
      intValue = atoi(str);
      (*specs)[curSpec.front()].scan = intValue;
      XMLString::release(&str);
      XMLString::release(&xmlStr);

      (*specs)[curSpec.front()].fileName = filename;

      xmlStr = XMLString::transcode("peaksCount");
      str = XMLString::transcode(attrs.getValue(xmlStr));
      intValue = atoi(str);
      (*specs)[curSpec.front()].resize(intValue);
      XMLString::release(&str);
      XMLString::release(&xmlStr);

      xmlStr = XMLString::transcode("msInstrumentID");
      if (attrs.getValue(xmlStr) != NULL)
      {
        str = XMLString::transcode(attrs.getValue(xmlStr));
        if (instrumentIdMap.count(str) > 0)
        {
          (*specs)[curSpec.front()].msMassAnalyzerType =
              Spectrum::parseMassAnalyzer(instrumentIdMap[str]);
        }
        else
        {
          WARN_MSG("Failed to locate instrument ID \'" << str << "\' in mzXML header");
        }
        XMLString::release(&str);
      }
      else if (instrumentIdMap.count(DEFAULT_INSTRUMENT_ID) > 0)
      {
        (*specs)[curSpec.front()].msMassAnalyzerType =
            Spectrum::parseMassAnalyzer(instrumentIdMap[DEFAULT_INSTRUMENT_ID]);
      }
      XMLString::release(&xmlStr);

      xmlStr = XMLString::transcode("collisionEnergy");
      if (attrs.getValue(xmlStr) != NULL)
      {
        str = XMLString::transcode(attrs.getValue(xmlStr));
        floatValue = atof(str);
        (*specs)[curSpec.front()].collision_energy = floatValue;
        XMLString::release(&str);
      }
      XMLString::release(&xmlStr);

      xmlStr = XMLString::transcode("retentionTime");
      str = XMLString::transcode(attrs.getValue(xmlStr));
      string str_string(str);
      size_t found = str_string.find("PT");
      if (found != string::npos)
        str_string.erase(found, 2);
      found = str_string.find("S");
      if (found != string::npos)
        str_string.erase(found, 1);

      (*specs)[curSpec.front()].retention_time = atof(str_string.c_str());

      XMLString::release(&str);
      XMLString::release(&xmlStr);
    }
  }
  if (state == STATE_SKIPPING)
  {
    XMLString::release(&name);
    return;
  }
  if (!strcmp(name, "precursorMz"))
  {
    state = STATE_PRECURSORMZ;
    xmlStr = XMLString::transcode("precursorCharge");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    (*specs)[curSpec.front()].parentCharge = str ? atoi(str) : 0;

    xmlStr = XMLString::transcode("precursorIntensity");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    (*specs)[curSpec.front()].precursor_intensity = str ? atof(str) : 0.0;

    xmlStr = XMLString::transcode("activationMethod");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    (*specs)[curSpec.front()].msFragType = Spectrum::FragType_CID; // Default value
    if (str and strncmp(str, "ETD", 3) == 0)
    {
      (*specs)[curSpec.front()].msFragType = Spectrum::FragType_ETD;
    }
    if (str and strcmp(str, "HCD") == 0)
      (*specs)[curSpec.front()].msFragType = Spectrum::FragType_HCD;
    XMLString::release(&str);
    XMLString::release(&xmlStr);
  }
  if (!strcmp(name, "peaks"))
  {
    state = STATE_PEAKS;
    xmlStr = XMLString::transcode("precision");
    str = XMLString::transcode(attrs.getValue(xmlStr));
    if ((str ? atoi(str) : 0) != 32)
    {
      stringstream aux;
      aux << "ERROR: " << atoi(str)
          << "-bit precision not supported, only 32-bit.\n";
      throw mzxmlEx(1, aux.str().c_str());
    }
    XMLString::release(&str);
    XMLString::release(&xmlStr);
  }
  XMLString::release(&name);

}

void mzxmlSAX2Handler::characters(const XMLCh* const chars,
                                  const XMLSize_t length)
{
  char *str = XMLString::transcode(chars);
  if (curSpec.empty()
      or str != (char *)0
          and (str[0] == (char)0 or str[0] == (char)13 or str[0] == (char)10))
  {
    XMLString::release(&str);
    return;
  }

  Spectrum *spec = &(*specs)[curSpec.front()];
  if (state == STATE_PRECURSORMZ)
  {
    //spec->parentMass = (float)atof(str);
    spec->parentMass = (float)atof(str);
    spec->parentMZ = (float)atof(str);
    if (spec->parentCharge > 0)
      spec->parentMass = spec->parentMass * spec->parentCharge
          - 1.0072763 * (spec->parentCharge - 1);  // Convert to standard MH+
  }

  Base64 b64;
  float fValue;
  if (state == STATE_PEAKS and spec->size() > 0)
  {
    unsigned char *decoded;
    int bytesRead = b64.b64_decode(str, strlen(str));
    //cout<<str<<endl;
    if (((unsigned int)bytesRead) != spec->size() * 8)
    {
      stringstream aux;
      aux << "ERROR: Could not decode peak list for scan " << spec->scan
          << " - needed " << spec->size() * 8 + 8 << " bytes, got " << bytesRead
          << " bytes (" << spec->size() << " peaks).\n";
      throw mzxmlEx(1, aux.str().c_str());
    }
    decoded = b64.getOutputBuffer();
    if (decoded == (unsigned char *)0)
    {
      stringstream aux;
      aux << "ERROR: Could not decode peak list for scan " << spec->scan
          << ".\n";
      throw mzxmlEx(1, aux.str().c_str());
    }

    unsigned int valIdx = 0;
    U32 tmp;
    for (unsigned int peakIdx = 0; peakIdx < spec->size(); peakIdx++)
    {
      for (unsigned int colIdx = 0; colIdx < 2; colIdx++, valIdx++)
      {
        tmp.u32 = (u_int32_t)ntohl((u_int32_t)((u_int32_t *)decoded)[valIdx]);
        try
        {
          (*spec)[peakIdx][colIdx] = tmp.flt;
        }
        catch (...)
        {
          stringstream aux;
          aux << "FNB2: " << peakIdx << "/" << spec->size() << ", col "
              << colIdx << "\n";
          throw mzxmlEx(1, aux.str().c_str());
        }
      }
    }
    b64.freeOutputBuffer();
  }

  XMLString::release(&str);
}

void mzxmlSAX2Handler::endElement(const XMLCh* const uri,
                                  const XMLCh* const localname,
                                  const XMLCh* const qname)
{

  char* name = XMLString::transcode(localname);
  if (!strcmp(name, "scan") and state != STATE_SKIPPING)
  {
    curSpec.pop_front();
    if (curSpec.empty())
      state = STATE_SKIPPING;
    else
      state = STATE_SCAN;
  }
  XMLString::release(&name);
}
;

void mzxmlSAX2Handler::endDocument()
{
}

void mzxmlSAX2Handler::fatalError(const SAXParseException& exception)
{
  char* message = XMLString::transcode(exception.getMessage());
  cout << "Fatal Error: " << message << " at line: "
      << exception.getLineNumber() << endl;
}
