#ifndef Utility_H
#define Utility_H


#include <iostream>

#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif

/* All code is based on
 * http://www.codeproject.com/articles/514443/debug-print-in-variadic-template-style
*/

/* my_print single (RCPP EXPOSED)
 * Arguments: one C-style string
 * This is the base-case for recursively printing a variable
 * number of string arguments between C++ and R to stdout.
 * Returns 0 if no errors in formatting detected.
*/
inline int my_print(const char *s)
{
    int rv = 0; // By default, assume success

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
                //throw std::runtime_error("invalid format string: missing arguments");
                rv = 1;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }

#ifndef STANDALONE
    Rcpp::Rcout.flush();
#else
    std::cout.flush();
#endif

    return rv;
}


/* my_print multiple (RCPP EXPOSED)
 * Arguments: a variable number of C-style strings
 * This is the recursive function to print a variable number of
 * string arguments between C++ and R to stdout.
 * Returns 0 if no errors in formatting detected.
*/
template<typename T, typename... Args>
inline int my_print(const char *s, T value, Args... args)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
#ifndef STANDALONE
                Rcpp::Rcout << value;
#else
                std::cout << value;
#endif
                rv = my_print(s + 1, args...); // call even when *s == 0 to detect extra arguments

                // TODO note: this may be an extraneous flush.
#ifndef STANDALONE
                Rcpp::Rcout.flush();
#else
                std::cout.flush();
#endif
                return rv;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to my_print");
    return 1;
}


/* my_printError single (RCPP EXPOSED)
 * Arguments: one C-style string
 * This is the base-case for recursively printing a variable
 * number of string arguments between C++ and R to stderr.
 * Returns 0 if no errors in formatting detected.
*/
inline int my_printError(const char *s)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
                //throw std::runtime_error("invalid format string: missing arguments");
                rv = 1;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }

#ifndef STANDALONE
    Rcpp::Rcerr.flush();
#else
    std::cerr.flush();
#endif

    return rv;
}


/* my_printError multiple (RCPP EXPOSED)
 * Arguments: a variable number of C-style strings
 * This is the recursive function to print a variable number of
 * string arguments between C++ and R to stderr.
 * Returns 0 if no errors in formatting detected.
*/
template<typename T, typename... Args>
inline int my_printError(const char *s, T value, Args... args)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
#ifndef STANDALONE
                Rcpp::Rcerr << value;
#else
                std::cerr << value;
#endif
                rv = my_printError(s + 1, args...); // call even when *s == 0 to detect extra arguments

                // TODO note: this may be an extraneous flush.
#ifndef STANDALONE
                Rcpp::Rcerr.flush();
#else
                std::cerr.flush();
#endif
                return rv;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to my_print");
    return 1;
}

//Blank header
#endif // Utility_H
