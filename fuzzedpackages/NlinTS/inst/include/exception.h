/**
 * @authors Hmamouche Youssef
 * @date    15/10/2017
 **/

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

class Exception : public std::exception
{
    std::string _message;
    
public:
    
    Exception (const std::string & message) throw ():
    _message (message)
    {};
    virtual ~Exception (void) {};
    virtual const char* what()  const throw()
    {
        return _message.c_str();
    };
};

#endif
