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

#ifndef _EA_STATISTICS_H

#define _EA_STATISTICS_H 1

/**
 * Este fichero permite mostrar información estadística sobre la ejecución de un AE. 
 *
 * Permite mostrar información por cada % de nueva generación, y tras cada evaluación
 *
 */
#include "problem.h"
#include "populationreal.h"

namespace realea {

class Statistics {
public:
    /**
     * Constructor, permite especificar cada cuantos genes debe de mostrarse la información
     */
    Statistics (int num_gen);

    void setProblem(Problem *problem) {
	m_problem = problem;
    }

    /*
     * Inicia las estadísticas
     */
    void reset(void);

    /**
     * Avisa de un evento. Por defecto sólo muestra los eventos activados con
     * activeEvent
     *
     * @param event Nombre del evento
     */
    void newEvent(string event);
    /**
     * Indica que se activa el evento deseado
     */
    void activeEvent(string event);
    /**
     * Notifica que empieza un nuevo experimento
     */
    void newExperiment(void);
    /**
     * Notifica que empieza una nueva generación
     */
    void newGeneration(void);
    /**
     * Avisa del fin de la generación, indicando el mejor fitness encontrado hasta el momento
     */
    void endGeneration(tFitness best);
    /**
     * Avisa del fin del experimento
     */
    void endExperiment(void);
private:
    Problem *m_problem;
    tFitness m_lastbest;
    unsigned m_experiment, m_generation;
    unsigned m_rate;
    map<string,bool> m_events;
};

}
#endif
