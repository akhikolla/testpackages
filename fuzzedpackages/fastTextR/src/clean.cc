#include <Rcpp.h>

// [[Rcpp::plugins("cpp11")]]


int is_break(std::string s) {
    if ( s.length() < 4 ) return 0;
    if ( s.compare(0, 4, "<br>") == 0 ) return 3;
    if ( s.compare(0, 5, "<br >") == 0 ) return 4;
    if ( s.compare(0, 6, "<br />") == 0 ) return 5;
    return 0;
}


// s1 == " " | s1 == "\a" | s1 == "\b" | s1 == "\f" 
// | s1 == "\n" | s1 == "\r" | s1 == "\t" | s1 == "\v"
int is_control_char(std::string s1) {
    if ( (s1 == " ") | (s1 == "\a") | (s1 == "\b") | (s1 == "\f") 
       | (s1 == "\n") | (s1 == "\r") | (s1 == "\t") | (s1 == "\v") ) {
        return 1;
    }
    return 0;
}


// s1 == "\'" | s1 == "\"" | s1 == "." | s1 == "," | s1 == "("
//   | s1 == ")" | s1 == "!" | s1 == "\?"
int is_punctation(std::string s1) {
    if ( (s1 == "\'") | (s1 == "\"") | (s1 == ".") | (s1 == ",") | (s1 == "(")
       | (s1 == ")") | (s1 == "!") | (s1 == "\?") ) {
        return 1;
    }
    return 0;
}


// ## add space
// ## "'", '"', ".", ",", "(", ")", "!", "?", ""
// ## replace with " "
// ## "<br />", ";", ":", 
// ## "<br />"
// [[Rcpp::export]]
std::vector<std::string> clean_text(std::vector<std::string> x) {
    std::vector<std::string> y(x.size());
    std::string s0, s1, s2;
    for (size_t i = 0; i < x.size(); i++) {
        s1 = ""; s2 = ""; y[i] = "";
        size_t n = x[i].length();
        for (size_t j = 0; j < n; j++) {
            s1 = x[i].at(j);
            if ( j + 1 <  n ) {
                s2 = x[i].at(j + 1);
            }

            if ( is_punctation(s1) ) {
                if ( s0 != " " ) {
                    y[i] += " ";
                }
                y[i] += s1;
                if ( j + 1 <  n ) {
                    if ( s2 != " ") {
                        y[i] += " ";       
                    }
                }
            } else if ( is_control_char(s1) ) {
                s1 = " ";
                if ( s0 != " " ) {
                    y[i] += " ";
                }
            } else if ( (s1 == ";") | (s1 == ":") ) {
                if ( s0 != " " ) {
                    y[i] += " ";
                }
            } else if ( (s1 == "<") & (s2 == "b") & (n - j > 4) ) {
                size_t len = ( 6 < (n - j) ) ? 6 : (n-j);
                j += is_break(x[i].substr(j, len));
                if ( j > 0 ) {
                    if ( s0 != " " ) {
                        y[i] += " ";       
                    }
                } else {
                    y[i] += s1;    
                }
            } else {
                y[i] += s1;
            }
            s0 = y[i].back();
        }
    }
    return y;
}
