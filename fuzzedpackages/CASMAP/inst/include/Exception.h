/*
 * Exception.h
 *
 *  Created on: Aug 19, 2016
 *      Author: mabaker
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <string>
#include <exception>
#include <stdexcept>

namespace SignificantPattern
{

	class Exception : public std::exception
	{
	protected:
		std::string text;

	public:
		explicit Exception (const char *text_);
		explicit Exception (const std::string& text_);
		virtual ~Exception() throw (){}

		 const char *what () const throw ()
		 {
			 return text.c_str();
		 }

	};

} /* namespace SignificantPattern */

#endif /* EXCEPTION_H_ */
