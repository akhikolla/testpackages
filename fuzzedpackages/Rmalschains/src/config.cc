/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include <exception>

extern Config *getConfig(void);

string Config::extractType(void) {
   string defaultvalue;

   if (!fileconfig.readInto(defaultvalue, m_strategy+".default"))
     defaultvalue = "";

   return defaultvalue;
}

string Config::extractName(void) {
   string begin = m_strategy +"." +m_type;
   string name;

   if (m_type == "") {
      return "";
   }
   // Comprueba que exista "<typestrategy>.<type>.id = <name>
   if (!fileconfig.readInto(name, begin +".id"))
      name = m_type;

   return name;
}
