#ifndef kvh_h
#define kvh_h

#include <string>
#include <Rcpp.h>

//using namespace std;
using namespace Rcpp;

// type declarations
typedef struct {
    std::string key;
    RObject val;
    bool tab_found;
} keyval;

typedef struct {
    List res;
    std::string line;
} list_line;

// function declarations
std::string unescape(std::string s);
bool indent_lacking(std::string& buf, size_t& lev);
bool escaped_eol(std::string& buf);
std::string kvh_get_line(std::ifstream& fin, const std::string& comment_str);
inline keyval kvh_parse_kv(std::string& line, size_t& lev, const bool strip_white, const std::string& split_str);
list_line kvh_read(FILE *fin, size_t lev, const std::string& comment_str, const bool strip_white, const bool skip_blank, const std::string& split_str, const bool follow_url);
RObject kvh_read(std::string fn, const std::string& comment_str, const bool strip_white, const bool skip_blank, const std::string& split_str, const bool follow_url);

#endif
