////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>

#include "ReportBase64.h"

////////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////////
//
// This executable implements a simple uuEncode64 tool that reads any given file,
// encodes it, and dumps it to the standard output. It relies on the function stored
// in the reports library, named base64_encode().
//

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // test for empty parameters. Only one is expected, which is the input filename
  if(argc != 2) {
    cout << "uuEncode64 - a uu64 file encoding tool." << endl;
    cout << "Usage: uuEncode64 <input file>" << endl;
    return 0;
  }

  // get the filename from the command line (as the 1st parameter)
  string fn = argv[1];

  // declare variables to store pointers to the binary data
  int length;
  char * buffer;

  // open the file
  ifstream is;
  is.open (fn.c_str(), ios::binary );

  // check if the file was open, and exits in error in case it is not
  if(!is.is_open()) {
    cout << "Unable to open file " << fn << endl;
    return false;
  }

  // get length of file
  is.seekg (0, ios::end);
  length = is.tellg();
  is.seekg (0, ios::beg);

  // allocate the memory needed
  buffer = new char [length+1];

  // read data as a block, and close the file
  is.read(buffer, length);
  is.close();
  // Add a terminator, just to make sure.
  buffer[length] = 0;

  // uuEncode the file and store it in a string. As uuEncoded data is a string with no extra characters, a string object may be used
  string m_image = base64_encode(reinterpret_cast<const unsigned char*>(buffer), length);

  // dump encoded file to stdout
  cout <<  m_image;

  // return ok
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
