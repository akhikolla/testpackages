/* 
 *  global.cpp : Global variables for the Oja Median library.
 *  Copyright (C) 2000 Tommi Ronkainen
 *  CVS $Id: global.cpp,v 1.1 2008/01/25 11:47:49 ruthe Exp $
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include <string.h>
#include <stdlib.h>
#include "global.h"

//using namespace std; //df

/* Ei geneerinen osa parseria */
int debug=0;
bool trace=false;
bool verbose=false;
bool quiet=false;
bool adaptive=false;
/* Geneerinen jatkuu */

static char* _programname;
static int _parameters;
static char**_parameter;
static int _options;
static char**_option;

bool parse_arguments(int argc,char** argv)
{
	bool options_done=false;
	
	_parameters=0;
	_parameter=new char*[argc];
	_options=0;
	_option=new char*[argc];

	const char* s=argv[0];
	for(const char* d=argv[0]; *d; d++)
		if(*d=='/')
			s=d+1;

	_programname=new char[strlen(s)+1];
	strcpy(_programname,s);

	if(argc==2)
	{
		if(strcmp(argv[1],"-?")==0 ||
		  strcmp(argv[1],"-h")==0 ||
		  strcmp(argv[1],"--help")==0)
		{
			return true;
		}
	}
		
	for(int i=1; i<argc; i++)
	{
		/* Ei geneerinen osa parseria */
		if(strcmp(argv[i],"-D")==0)
			debug=1;
		else if(strcmp(argv[i],"-t")==0)
			trace=true;
		else if(strcmp(argv[i],"-v")==0)
			verbose=true;
		else if(strcmp(argv[i],"-q")==0)
			quiet=true;
		/* Geneerinen jatkuu */
		else if(options_done)
			_parameter[_parameters++]=argv[i];
		else if(strcmp(argv[i],"--")==0)
			options_done=true;
		else
		{
			if(argv[i][0]=='-')
			{
				if(i+2 < argc && argv[i+2][0]=='-' && argv[i+1][0]!='-')
				{
					_option[_options]=new char[strlen(argv[i])+strlen(argv[i+1])+1];
					strcpy(_option[_options],argv[i]);
					strcpy(_option[_options]+strlen(argv[i]),argv[i+1]);
					_options++;
					i++;
				}
				else
					_option[_options++]=argv[i];
			}
			else
			{				
				options_done=true;
				_parameter[_parameters++]=argv[i];				
			}
		}
	}

	return false;
}

const char* program_name()
{
	return _programname;
}

int options()
{
	return _options;
}

const char* option(int idx)
{
	if(idx < 0 || idx >= _options)
		return "";

	return _option[idx];
}

int parameters()
{
	return _parameters;
}

const char* parameter(int idx)
{
	if(idx < 0)
		idx = _parameters+idx;
	
	if(idx < 0 || idx >= _parameters)
		return "";

	return _parameter[idx];
}

const char* find_option(const char* name)
{
	for(int i=0; i<_options; i++)
		if(strncmp(_option[i],name,strlen(name))==0)
			return _option[i]+strlen(name);

	return 0;
}

bool has_option(const char* opt)
{
	return find_option(opt) != 0;
}
	
int int_option(const char* opt,int def)
{
	const char* s;

	s=find_option(opt);
	if(s)
		return atoi(s);
	else
		return def;
}

double double_option(const char* opt,double def)
{
	const char* s;

	s=find_option(opt);
	if(s)
		return atof(s);
	else
		return def;
}

const char* string_option(const char* opt,const char* def)
{
	const char* s;

	s=find_option(opt);
	if(s)
		return s;
	else
		return def;
}

