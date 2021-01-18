/*
 * Objective.h
 *
 *  Created on: Aug 21, 2015
 *      Author: Javon
 */

#ifndef SRC_OBJECTIVE_H_
#define SRC_OBJECTIVE_H_

class Objective {
public:
	virtual ~Objective() {};

	Objective()
	: m_overalObj(0), m_loss(0), m_penalty(0){}

	Objective(double overalObj, double loss, double penalty)
	: m_overalObj(overalObj), m_loss(loss), m_penalty(penalty) {}

	inline void setOveralObj(double overalObj) {
		m_overalObj = overalObj;
	}

	inline double getOveralObj() {
		return m_overalObj;
	}

	inline void setLoss(double loss) {
		m_loss = loss;
	}

	inline double getLoss() {
		return m_loss;
	}

	inline void setPenalty(double penalty) {
		m_penalty = penalty;
	}

	inline double getPenalty() {
		return m_penalty;
	}

	/** initialize all components to be zero */
	void init() {
		m_overalObj = 0;
		m_loss = 0;
		m_penalty = 0;
	}

private:
	double m_overalObj;
	double m_loss;
	double m_penalty;
};

#endif /* SRC_OBJECTIVE_H_ */
