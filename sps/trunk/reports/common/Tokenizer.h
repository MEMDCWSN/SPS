#ifndef __TOKENIZER_H__
#define __TOKENIZER_H__

#include <string>
#include <vector>

using namespace std;
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////

/*! \brief Split a string

  Splits a string into sub strings given a set of separators

  @param str The string to split
  @param results The vector of resulting substrings
  @param removeDelimiters If delimeters are to be removed
  @param separators A string containing all the 1 char separators
*/
void StringExplode(string str, vector<string> &results, bool removeDelimiters = false, string separators = " \t\f\v\n\r", string delimiters = "\"'");

/*! \brief Compares two strings

  Compares two strings, and consideres string length (length has priority over order)

  @param s1 The 1st string to compare
  @param s2 The 2nd string to compare
  @return -1, 0 or 1
*/
bool stringSortCompare(const string &s1, const string &s2);

/*! \brief Compares two strings

  Compares two strings

  @param s1 The 1st string to compare
  @param s2 The 2nd string to compare
  @return -1, 0 or 1
*/
bool stringUniqueCompare(const string &s1, const string &s2);

////////////////////////////////////////////////////////////////////////////////
};
////////////////////////////////////////////////////////////////////////////////
#endif // __TOKENIZER_H__
////////////////////////////////////////////////////////////////////////////////
