////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_COLUMN_TYPES_H__
#define __REPORT_COLUMN_TYPES_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <list>

#include "spsFiles.h"
////////////////////////////////////////////////////////////////////////////////
using namespace std;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportParamsOption

   Defines parameters options for renderes. Used when defined report table views.
   This is a simple triplet of string containing the parameter, the option and
   the validator.
   The parser uses data from each option to build the command line. For each
   option, the validator is evaluated first, and if it evaulates to an emtpy
   string, the option is excluded.

   */
class ReportParamsOption {

 public:

  /*! \brief Parameter idetifier
  Specifies the parameter contents
   */
  string param;
  /*! \brief Specified option
  Specifies the option (if available)
   */
  string option;
  /*! \brief Validator.
    If validator evaluates to an empty string, then the parameter is excluded
   */
  string validator;

  /*! \brief Default contructor
   */
  ReportParamsOption(string p, string o, string v) :
    param(p), option(o), validator(v)
    {};

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportParamsFiles

   Defines parameters input files for renderes. Used when defined report table
   views.
   Each instance of the record below contains information about a file, including
   the file contents.
   This is used to simplify the process of sharing files between spsReports and
   specplot/contplot. Initially, the drawing tools would load each file when
   generating every single image in the reports. To avoid this extremelly
   time consuming and pointless process, files started to be loaded by spsReports,
   and then shared with specplot and contplot. In order to do so, and as they would
   reside in the same process, data is shared by sending a pointer to a vector of
   file containers, being a file container the structure defined below.

   */
class ReportParamsFiles {

 public:

  /*! \brief File type
   */
  int         type;
  /*! \brief File ID.
    Used to identify files in SPS context. The allows a process to
    search for a file, given it's id.
   */
  SpsFileID   file;
  /*! \brief Parameter
   */
  string      param;
  /*! \brief Validator.
    If validator evaluates to an empty string, then the parameter is excluded
   */
  string      validator;
  /*! \brief The file contents.
   */
  void       *data;
  /*! \brief If the file is loaded and valid.
   */
  bool        valid;

  /*! \brief Default constructor
   */
  ReportParamsFiles(int t, SpsFileID f, string p, string v) :
    type(t), file(f), param(p), validator(v), data(NULL), valid(false)
    {};

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportColumnTypeBase

   Base class for column types.
   This structure defines the base classs for structures defining the possible
   column types that can be used in a report table cell.
   Common data types are defined here, which include the links, column labels,
   etc..

   */
class ReportColumnTypeBase {

 public:

  /*! \brief CSS class name for HTML formating
   */
  string  cssClass;       // CSS class name for HTML formating
  /*! \brief If a divifder is to be renderer to the left
   */
  bool    leftDivider;    //
  /*! \brief specifies a cell with content dynamically used (content sent to server)
   */
  bool    dynamic;        // specifies a cell with content dynamically used (content sent to server)
  /*! \brief label on the table column
   */
  string  columnLabel;    // label on the table column
  /*! \brief URL template for the link
   */
  string  link;           // URL template for the link
  /*! \brief template for onClick HTML method
   */
  string  onClick;        // template for onClick HTML method
  /*! \brief template for field ID, needed to read or write data to
   */
  string  id;             // template for field ID, needed to read or write data to
  /*! \brief validator must be not null in order for the cell to be displayed
   */
  string  validator;      // validator must be not null in order for the cell to be displayed
  /*! \brief display level. If input level < display level, then the item is not displayed
   */
  int     displayLevel;   // display level. If input level < display level, then the item is not displayed

  // constructor used to initialize methods
  /*! \brief constructor used to initialize methods
   */
  ReportColumnTypeBase() : leftDivider(false), dynamic(false), displayLevel(0) {};

  // virtual destructor to make class polymorfic
  /*! \brief virtual destructor to make class polymorfic
   */
  virtual ~ReportColumnTypeBase() {};

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportColumnTypeImageOnDemand

   Column type for Image On Demand (uses AJAX for requests to the server).

   When using a CGI call, the command is constructed in the following way:
   /cgi-bin/<renderer> <params>

   when rendering local static pages, the renderer name is used to generate/request a render object by name (using a object factory model)
   and <params> are passed to build the image
   */
class ReportColumnTypeImageOnDemand : public ReportColumnTypeBase {

 public:

  /*! \brief renderer used to generate the icon. If empty, iconParams treated as image/URL
   */
  string                      iconRenderer;       // renderer used to generate the icon. If empty, iconParams treated as image/URL
  //string                      iconParams;         // Icon path/image/URL
  /*! \brief Icon path/image/URL
   */
  vector<ReportParamsOption>  iconParams;         // Icon path/image/URL
  /*! \brief display level for the icon. If input level < display level, then the item is not displayed
   */
  int                         iconDisplayLevel;   // display level for the icon. If input level < display level, then the item is not displayed
  /*! \brief alternative to display in place of icon
   */
  string                      alt;              // alternative to display in place of icon


  /*! \brief label to show for the link (defined by renderer and params)
   */
  string                      label;            // label to show for the link (defined by renderer and params)
  /*! \brief Object name used for rendering the Image On Demand
   */
  string                      renderer;         // Object name used for rendering the Image On Demand
  /*! \brief parameters passed to the renderer object on Image On Demand
   */
  vector<ReportParamsOption>  params;           // parameters passed to the renderer object on Image On Demand
  /*! \brief display level for the link. If input level < display level, then the item is not displayed
   */
  int                         linkDisplayLevel; // display level for the link. If input level < display level, then the item is not displayed
  /*! \brief used to specify if label should be split into chunks (i.e. if it is a sequence)
   */
  bool                        splitLabel;       // used to specify if label should be split into chunks (i.e. if it is a sequence)

  /*! \brief Files passed to the renderer object on Image On Demand
   */
  vector<ReportParamsFiles>   files ;           // Files passed to the renderer object on Image On Demand


  // When using a CGI call, the command is constructed in the following way:
  // /cgi-bin/<renderer> <params>
  //
  // when rendering local static pages, the renderer name is used to generate/request a render object by name (using a object factory model)
  // and <params> are passed to build the image

  /*! \brief Default constructor
   */
  ReportColumnTypeImageOnDemand() : iconDisplayLevel(0), linkDisplayLevel(0), splitLabel(false) {}

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportColumnTypeString

   Column type for strings.

   */
class ReportColumnTypeString : public ReportColumnTypeBase {

 public:

  /*! \brief Text template for cell contents, button, input box.
   */
  string  text;       // Text template for cell contents, button, input box.
  /*! \brief If true, a button is drawn with the text in the "text" field
   */
  bool    isButton;   // If true, a button is drawn with the text in the "text" field
  /*! \brief if True, an input box is drawn
   */
  bool    isInput;    // if True, an input box is drawn
  /*! \brief used to specify if text should be split into chunks (i.e. if it is a sequence)
   */
  bool    splitText;  // used to specify if text should be split into chunks (i.e. if it is a sequence)


  // constructor used to initialize methods
  /*! \brief constructor used to initialize
   */
  ReportColumnTypeString() : isButton(false), isInput(false), splitText(false) {};

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportColumnTypeStringMultiple

   Column type for multiple strings. Used to generate input spectra file list in initial page.

   */
class ReportColumnTypeStringMultiple : public ReportColumnTypeBase {

 public:

  /*! \brief link filename prefix
   */
  string  linkPrefix; // link filename prefix
  /*! \brief link filename suffix
   */
  string  linkSuffix; // link filename suffix
  /*! \brief Text template for cell contents.
   */
  string  text;       // Text template for cell contents.

  // constructor used to initialize methods
  /*! \brief constructor used to initialize
   */
  ReportColumnTypeStringMultiple() {};

};
////////////////////////////////////////////////////////////////////////////////
 /*! \brief ReportColumnTypeBox

   Column type for containers. Containers may have an arbitrary number of ReportColumnTypeBase elements

   */
class ReportColumnTypeBox : public ReportColumnTypeBase {

 public:

  /*! \brief Vector of several column types.
   */
  vector<ReportColumnTypeBase *> sequences; // Vector of several column types.

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
