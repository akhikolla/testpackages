#ifndef HIPREAD_STRINGUTILS_H_
#define HIPREAD_STRINGUTILS_H_

#include "boost.h"

// Adapted from readr's QiParsers.h
template <typename Iterator, typename Attr>
inline bool parseDouble(Iterator& first, Iterator& last, Attr& res) {
  return boost::spirit::qi::parse(first, last, boost::spirit::qi::double_, res);
}

template <typename Iterator, typename Attr>
inline bool parseInteger(Iterator& first, Iterator& last, Attr& res) {
  return boost::spirit::qi::parse(first, last, boost::spirit::qi::long_, res);
}


class IpStringUtils {
public:
  // Adapted from https://stackoverflow.com/questions/216823
  static void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
      return !std::isspace(ch);
    }));
  }

  static void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
      return !std::isspace(ch);
    }).base(), s.end());
  }

  static void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
  }

  static void newtrim(const char* &start, const char* &end) {
    start = std::find_if(start, end, [](int ch) {
      return !std::isspace(ch);
    });

    while(end > start) {
      if (!std::isspace(end[-1])) break;
      --end;
    }
  }

};

#endif
