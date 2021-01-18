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

#ifndef _CONFIG
#define _CONFIG 1

#include "ConfigFile.h"
#include "debug.h"
#include <string>
#include <iostream>

using namespace std;


/**
 *
 * @class Config 
 *
 * Esta clase permite leer un fichero de configucación de la aplicación.
 *
 * El modelo que admite es el siguiente:
 *
 * En un fichero de configuración se especifique los distintos
 * parámetros para cada uno de las estrategias del AG. 
 *
 * Se considera estrategia a cada uno de los métodos de cruce, 
 * de selección o de reemplazo a aplicar.
 *
 * Se considera <type> de una estrategia no sólo a un tipo de estrategia
 * (como oeprador de cruce BLX o alpha) sino también al conjunto de parámetros.
 *
 * La idea es poder asociar el mismo identificador de estrategia con distintos
 * conjuntos de parámetros, cada uno identificado con un nombre distintivo. 
 *
 * De esa forma, simplemente indicando el tipo de la estrategia a aplicar
 * se seleccionaría el total de parámetros. 
 *
 * formato:
 * <typeestrategia>.<type>.id = <name>
 *
 * Para asociar a este <type> los distintos parámetros la sintaxis es la 
 * siguiente (se supone que el objeto factory asociado identificará
 * la distinta ). 
 * 
 * <typeestrategia>.<type>.<param1> = <value1>
 * <typeestrategia>.<type>.<param2> = <value2>
 *
 * o bien mediante un fichero externo
 * <typeestrategia>.<type>.params = <filename>
 *
 * <filename> en este caso contendría la sintaxis
 * <param1> = <value1>
 * <param2> = <value2>
 * .
 * .
 * .
 * 
 * Se puede especificar parámetros comunes a los distintos tipos asociados
 * a una estrategia con la sintaxis. 
 * <typeestrategia>.<name>.<param1> = <value1>
 *
 * Dado que la idea de esta sintaxis es especificar parámetros comunes
 * tendrá preferencia la sintaxis <typeestrategia>.<type>.<paran> sobre
 * ésta. 
 * 
 * La elección de la estrategia (el identificador <type>) 
 * a aplicar se puede especificar mediante paso de parámetros
 * o bien mediante un valor default:
 *
 * <typeestrategia>.<default> = <type>
 * 
 */

class Config {
private:
    typedef ConfigFile::key_not_found key_not_found;
    ConfigFile fileconfig;
    string m_strategy;
    string m_type;
    string m_name;

private:
    /**
     * Permite obtener el valor por defecto si existe
     *
     * @return el tipo asociado
     */
    string extractType(void);

    /**
     * Permite obtener el nombre de la estrategia elegida.
     * Si no existe lanza una excepción (type_not_valide)
     * 
     * @return El nombre de la estrategia elegia ('' si todavía
     * no existe ninguna elegida). 
     */
    string extractName(void);

public:

    /**
     * @class config_error
     *
     * Esta clase es lanzada cuando se produce un error en la lectura del fichero
     * de configuración
     */
    struct config_error : public std::exception{
       string m_msg;
       string m_strategy;
       config_error( const string& strategy, const string& msg)
	 : m_msg(msg), m_strategy(strategy) {} 
       virtual const char* what() const throw() {
	  string output = "Config Error: " +m_msg;
	  return output.c_str();
       }

       virtual ~config_error() throw() {}
    };

public:
    /**
     * Constructor.
     * @param strategy
     * @param name Nombre del fichero
     */
    Config(string strategy, string name) : fileconfig(name) {
       m_strategy = strategy;
       m_type = extractType();
       m_name = extractName();
    }

    /**
     * Constructor
     * @param config Fichero de configuración
     */
    Config(string strategy, ConfigFile &config) : fileconfig(config) {
       print_info("Creando fileconfig");
       m_strategy = strategy;
       m_type = extractType();
       m_name = extractName();
    }

    Config(string strategy, Config *config) : fileconfig(config->fileconfig) {
       m_strategy = strategy;
       m_type = extractType();
       m_name = extractName();
    }


    /**
     * Devuelve el nombre de la estrategia elegida
     */
    string getName() {
       if (m_name != "") {
	 return m_name;
       }
       else {
	 return m_type; 
       }
    }

    /**
     * Devuelve el nombre del grupo elegido de la estrategia 
     */
    string getType() {
       return m_type;
    }


    /**
     * Especifico el tipo
     */
    void setType(string name) {
        cerr <<"setType : '" <<name <<"'" <<endl;
        if (name != "") {
	    m_type = name;
	    m_name = extractName();
            cerr <<"Hecho asignacion de extractName" <<endl;
	}

        cerr <<"Realizado setType : '" <<name <<"'" <<endl;
    }

    /**
     * Especifico el identificador
     */
    void setName(string name) {
        if (name != "") {
	    m_name = name;
	}
    }


    /**
     * Obtengo el parámetro asociado
     *
     * @param param el nombre del parámetro
     * @return Valor en forma de cadena
     */
    template <class T>
    void getParam(const char *param, T &value) {
       string begin = m_strategy +".";
       string last = ".";
       last += param;
       string key = begin +m_type +last;
       string key_name_default = begin +m_name +last;
       string key_type_default = begin + param;

       if (fileconfig.readInto(value, key) ) {
       }
       else if (fileconfig.readInto(value, key_name_default)) {
       }
       else if (fileconfig.readInto(value, key_type_default)) {
       }
       else  {
	  throw ConfigFile::key_not_found(key_type_default);
       }
    }

    template <class T>
      void getParam(const string &param, T &value) {
	 getParam(param.c_str(), value);
      }

};

typedef Config* ConfigPtr;

#endif
