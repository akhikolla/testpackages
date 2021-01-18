/*
	Copyright (c) 2000 Paul T. Miller

	LGPL DISCLAIMER
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Library General Public
	License as published by the Free Software Foundation; either
	version 2 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Library General Public License for more details.

	You should have received a copy of the GNU Library General Public
	License along with this library; if not, write to the
	Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	Boston, MA  02111-1307, USA.

	http://www.gnu.org/copyleft/lgpl.html
*/

#ifndef XMLCONFIG_H
#define XMLCONFIG_H

#ifdef __cplusplus
	#define XML_BEGIN_NAMESPACE	namespace XML {
	#define XML_END_NAMESPACE	}
#endif

#include <stdlib.h>

#ifdef XML_UNICODE
	#error UNICODE not implemented
	typedef unsigned short XML_Char;	// UTF-16
#else
	typedef char XML_Char;				// UTF-8
#endif

#endif	// XMLCONFIG_H

