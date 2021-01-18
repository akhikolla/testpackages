#ifndef MANIFOLD_OPTIM_EXCEPTION_H
#define MANIFOLD_OPTIM_EXCEPTION_H

#include <exception>
#include <string>

class ManifoldOptimException : public std::exception
{
public:
	ManifoldOptimException(const std::string& msg)
		: m_msg(msg)
	{
	}

	~ManifoldOptimException() throw() {}

	virtual const char* what() const throw()
	{
		return m_msg.c_str();
	}

private:
	std::string m_msg;
};

#endif
