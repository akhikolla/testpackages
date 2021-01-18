/*
 * Exception.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: mabaker
 */

#include "Exception.h"

namespace SignificantPattern
{

	Exception::Exception (const char *text_)
		: text(text_)
	{
	}
	Exception::Exception (const std::string& text_)
		: text(text_)
	{
	}

} /* namespace SignificantPattern */
