#include "ion.h"

using namespace std;


const Ion Ion::ion[15] =
{
  Ion(true,  Ion::a, -44.998203008, "-H_20",       0.113),
  Ion(true,  Ion::a, -44.014187423, "-NH_3",       0.150),
  Ion(true,  Ion::b, -35.013853072, "-H_20-H_20",  0.145),
  Ion(true,  Ion::b, -34.029837487, "-H_20-NH_3",  0.163),
  Ion(true,  Ion::a, -26.987638322, "",            0.174),
  Ion(true,  Ion::b, -17.003288386, "-H_20",       0.300),
  Ion(true,  Ion::b, -16.019272801, "-NH_3",       0.298),
  Ion(true,  Ion::b, 1.007276300,   "",            0.502),
  Ion(true,  Ion::b, 2.010631137,   "(iso)",       0.445),
  Ion(false,  Ion::y, -17.003288386, "-H_20-H_20",  0.101),
  Ion(false,  Ion::y, -16.019272801, "-H_20-NH_3",  0.122),
  Ion(false,  Ion::y, 1.007276300,   "-H_20",       0.210),
  Ion(false,  Ion::y, 1.991291885,   "-NH_3",       0.220),
  Ion(false,  Ion::y, 19.017840986,  "",            0.515),
  Ion(false,  Ion::y, 20.021195823,  "(iso)",       0.465),
};

static pair<Ion::Type, std::string> name_value[] =
{
  pair<Ion::Type, std::string>(Ion::a, "a"),
  pair<Ion::Type, std::string>(Ion::b, "b"),
  pair<Ion::Type, std::string>(Ion::c, "c"),
  pair<Ion::Type, std::string>(Ion::x, "x"),
  pair<Ion::Type, std::string>(Ion::y, "y"),
  pair<Ion::Type, std::string>(Ion::z, "z"),
  pair<Ion::Type, std::string>(Ion::phos, "phos"),
};

map<Ion::Type, string> Ion::cname(name_value, name_value + sizeof(name_value) / sizeof(* name_value));

static pair<std::string, Ion::Type> type_value[] =
{
  pair<std::string, Ion::Type>("a", Ion::a),
  pair<std::string, Ion::Type>("b", Ion::b),
  pair<std::string, Ion::Type>("c", Ion::c),
  pair<std::string, Ion::Type>("x", Ion::x),
  pair<std::string, Ion::Type>("y", Ion::y),
  pair<std::string, Ion::Type>("z", Ion::z),
  pair<std::string, Ion::Type>("phos", Ion::phos),
};

map<string, Ion::Type> Ion::ctype(type_value, type_value + sizeof(type_value) / sizeof(* type_value));
