//
//  Logger.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.10.2013.
//
//

#include "config.h"

#include "Logger.h"
#include <RcppArmadillo.h>
#include <exception>

#ifdef HAVE_PTHREAD_H

#ifdef ENABLE_DEBUG_VERBOSITY
#define CHECK_PTHREAD_RETURN_CODE(expr) {int rc = expr; if((rc) != 0) { Rcpp::Rcout << "Warning: Call to pthread function failed with error code " << (rc) << " in " << __FILE__ << ":" << __LINE__ << std::endl; }}

#else
#define CHECK_PTHREAD_RETURN_CODE(expr) {expr;}
#endif

template <>
inline std::streamsize LoggerStreamBuffer<false>::xsputn(const char *s, std::streamsize n) {
	if(this->threadSafe) {
		this->tsBuffer.append(s, n);
	} else {
		Rprintf("%.*s", n, s);
	}
	return n;
}

template <>
inline std::streamsize LoggerStreamBuffer<true>::xsputn(const char *s, std::streamsize n) {
	if(this->threadSafe) {
		this->tsBuffer.append(s, n);
	} else {
		REprintf("%.*s", n, s);
	}
	return n;
}

template <>
inline int LoggerStreamBuffer<false>::overflow(int c) {
	if(c != traits_type::eof()) {
		if(this->threadSafe) {
			this->tsBuffer.append(1, (char) c);
		} else {
			Rprintf("%.1s", &c);
		}
	}
	return c;
}

template <>
inline int LoggerStreamBuffer<true>::overflow(int c) {
	if(c != traits_type::eof()) {
		if(this->threadSafe) {
			this->tsBuffer.append(1, (char) c);
		} else {
			Rprintf("%.1s", &c);
		}
	}
	return c;
}

template <>
inline int LoggerStreamBuffer<false>::sync() {
	if(!this->threadSafe) {
		R_FlushConsole();
	}
	return 0;
}

template <>
inline int LoggerStreamBuffer<true>::sync() {
	if(!this->threadSafe) {
		R_FlushConsole();
	}
	return 0;
}

template <>
void LoggerStreamBuffer<false>::flushThreadSafeBuffer() {
	if(this->tsBuffer.length() > 0) {
		Rprintf("%.*s", this->tsBuffer.length(), this->tsBuffer.c_str());
		R_FlushConsole();
		this->tsBuffer.clear();
	}
}

template <>
void LoggerStreamBuffer<true>::flushThreadSafeBuffer() {
	if(this->tsBuffer.length() > 0) {
		Rprintf("%.*s", this->tsBuffer.length(), this->tsBuffer.c_str());
		R_FlushConsole();
		this->tsBuffer.clear();
	}
}

/*
 * Logger
 * if pthreads are available
 */
template <>
Logger<false>::Logger() : std::ostream(new Buffer()), buf(static_cast<Buffer*>(rdbuf())), threadSafe(false) {
	int pthreadRC = pthread_mutex_init(&this->printMutex, NULL);
	if(pthreadRC != 0) {
		throw std::runtime_error("Mutex to synchronize printing could not be initialized");
	}
}

template <>
Logger<true>::Logger() : std::ostream(new Buffer()), buf(static_cast<Buffer*>(rdbuf())), threadSafe(false) {
	int pthreadRC = pthread_mutex_init(&this->printMutex, NULL);
	if(pthreadRC != 0) {
		throw std::runtime_error("Mutex to synchronize printing could not be initialized");
	}
}

template <>
Logger<false>::~Logger()  {
	if(this->buf != NULL) {
		delete this->buf;
		this->buf = NULL;
	}
	pthread_mutex_destroy(&this->printMutex);
}

template <>
Logger<true>::~Logger()  {
	if(this->buf != NULL) {
		delete this->buf;
		this->buf = NULL;
	}
	pthread_mutex_destroy(&this->printMutex);
}

template <>
void Logger<false>::flushThreadSafeBuffer() {
	if(this->buf != NULL) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		this->buf->flushThreadSafeBuffer();
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	}
}

template <>
void Logger<true>::flushThreadSafeBuffer() {
	if(this->buf != NULL) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		this->buf->flushThreadSafeBuffer();
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	}
}

template <>
void Logger<false>::placeMutexLock(bool lock) {
	if(this->threadSafe) {
		if(lock) {
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		} else {
			this->flush();
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
		}
	}
}


template <>
void Logger<true>::placeMutexLock(bool lock) {
	if(this->threadSafe) {
		if(lock) {
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		} else {
			this->flush();
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
		}
	}
}

#else
/*
 * LoggerStreamBuffer
 * if pthreads are NOT available
 */
template <>
inline std::streamsize LoggerStreamBuffer<false>::xsputn(const char *s, std::streamsize n) {
	Rprintf("%.*s", n, s);
	return n;
}

template <>
inline std::streamsize LoggerStreamBuffer<true>::xsputn(const char *s, std::streamsize n) {
	REprintf("%.*s", n, s);
	return n;
}

template <>
inline int LoggerStreamBuffer<false>::overflow(int c) {
	if(c != traits_type::eof()) {
		Rprintf("%.1s", &c);
	}
	return c;
}

template <>
inline int LoggerStreamBuffer<true>::overflow(int c) {
	if(c != traits_type::eof()) {
		REprintf("%.1s", &c);
	}
	return c;
}

template <>
inline int LoggerStreamBuffer<false>::sync() {
	R_FlushConsole();
	return 0;
}

template <>
inline int LoggerStreamBuffer<true>::sync() {
	R_FlushConsole();
	return 0;
}

template <>
void LoggerStreamBuffer<false>::flushThreadSafeBuffer() {
}

template <>
void LoggerStreamBuffer<true>::flushThreadSafeBuffer() {
}

/*
 * Logger
 * if pthreads are NOT available
 */

template <>
Logger<false>::Logger() : std::ostream(new Buffer()), buf(static_cast<Buffer*>(rdbuf())), threadSafe(false) {}

template <>
Logger<true>::Logger() : std::ostream(new Buffer()), buf(static_cast<Buffer*>(rdbuf())), threadSafe(false) {}

template <>
Logger<false>::~Logger()  {
	if(this->buf != NULL) {
		delete this->buf;
		this->buf = NULL;
	}
}

template <>
Logger<true>::~Logger()  {
	if(this->buf != NULL) {
		delete this->buf;
		this->buf = NULL;
	}
}

template <>
void Logger<false>::flushThreadSafeBuffer() {
	if(this->buf != NULL) {
		this->buf->flushThreadSafeBuffer();
	}
}

template <>
void Logger<true>::flushThreadSafeBuffer() {
	if(this->buf != NULL) {
		this->buf->flushThreadSafeBuffer();
	}
}

template <>
void Logger<false>::placeMutexLock(bool) {}


template <>
void Logger<true>::placeMutexLock(bool) {}

#endif
Logger<false> GAout;
Logger<true> GAerr;
