//
//  Logger.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.10.2013.
//
//

#ifndef GenAlgPLS_Logger_h
#define GenAlgPLS_Logger_h

#include <iostream>
#include <streambuf>
#include <string>

template <bool ERROR_STREAM>
class LoggerStreamBuffer : public std::streambuf
{
public:
	LoggerStreamBuffer() : threadSafe(false) {};
	virtual ~LoggerStreamBuffer() {};

	void flushThreadSafeBuffer();

	void enableThreadSafety(bool threadSafe) {
		this->flushThreadSafeBuffer();
		this->threadSafe = threadSafe;
	}
protected:
	virtual std::streamsize xsputn(const char *s, std::streamsize n);
	virtual int overflow(int c = traits_type::eof());
	virtual int sync();

private:
	bool threadSafe;
	std::string tsBuffer;
};

template <bool ERROR_STREAM>
class Logger : public std::ostream
{
private:
	typedef LoggerStreamBuffer<ERROR_STREAM> Buffer;
	Buffer* buf;
	bool threadSafe;
#ifdef HAVE_PTHREAD_H
	pthread_mutex_t printMutex;
#endif
public:
	Logger();
	~Logger();

	void flushThreadSafeBuffer();
	
	void placeMutexLock(bool lock);

	void enableThreadSafety(bool threadSafe = true) {
		this->threadSafe = threadSafe;
		if(this->buf != NULL) {
			this->buf->enableThreadSafety(threadSafe);
		}
	}

	class LogLocker {
	public:
		LogLocker(Logger& logger, bool lock) : logger(logger), lock(lock) {}

		friend std::ostream& operator<<(std::ostream& os, const LogLocker& locker) {
			locker.logger.placeMutexLock(locker.lock);
			return os;
		}

	private:
		Logger& logger;
		bool lock;
	};

	LogLocker lock() {
		return LogLocker(*this, true);
	}

	LogLocker unlock() {
		return LogLocker(*this, false);
	}

};

extern Logger<false> GAout;
extern Logger<true> GAerr;

#endif
