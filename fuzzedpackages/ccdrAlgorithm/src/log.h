#ifndef __LOG_H__
#define __LOG_H__

//
// This is a simple log header class
//
// To change reporting level:
//	 FILELog::ReportingLevel() = logDEBUG1;
// To log:
//   FILE_LOG(logINFO) << "text goes here";
// To make it output to a file:
//   FILE* pFile = fopen("log.txt", "w");
//   Output2FILE::Stream() = pFile;
//

#include <sstream>
#include <string>
#include <stdio.h>
#include <iomanip>

inline std::string NowTime();

enum TLogLevel {logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};

template <typename T>
class Log
{
public:
    Log();
    virtual ~Log();
    std::ostringstream& Get(TLogLevel level = logINFO, std::string file = "");
public:
    static TLogLevel& ReportingLevel();
    static std::string ToString(TLogLevel level);
    static TLogLevel FromString(const std::string& level);
protected:
    std::ostringstream os;
private:
	TLogLevel level;
	std::string file;
    Log(const Log&);
    Log& operator =(const Log&);
};

template <typename T>
Log<T>::Log()
{
}

template <typename T>
std::ostringstream& Log<T>::Get(TLogLevel level, std::string file)
{
	this->level = level;
    os << "- " << NowTime();
    os << " " << std::setw(6) << ToString(level) << ": ";
    os << std::string(level > logDEBUG ? level - logDEBUG : 0, '\t');
	this->file = file;
    return os;
}

template <typename T>
Log<T>::~Log()
{
	//add appropriate tabs to newlines
	std::string output = os.str();
	std::size_t i = 0;
	int numTabs = 2 + (level > logDEBUG ? level - logDEBUG : 0);
	while ( i < output.size() && i != std::string::npos) {
		i = output.find('\n',i+numTabs+1);
		if (i != std::string::npos) {
			std::string tabs = std::string(numTabs, '\t');
			output.insert(i+1,tabs);
		}
	}
	if (output.back() == '\n') output.erase(output.size()-1);

//    output += " @"  + file + "\n";
    output += " "  + file + "\n";
    T::Output(output);
}

template <typename T>
TLogLevel& Log<T>::ReportingLevel()
{
    static TLogLevel reportingLevel = logDEBUG4;
    return reportingLevel;
}

template <typename T>
std::string Log<T>::ToString(TLogLevel level)
{
	static const char* const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
    return buffer[level];
}

template <typename T>
TLogLevel Log<T>::FromString(const std::string& level)
{
    if (level == "DEBUG4")
        return logDEBUG4;
    if (level == "DEBUG3")
        return logDEBUG3;
    if (level == "DEBUG2")
        return logDEBUG2;
    if (level == "DEBUG1")
        return logDEBUG1;
    if (level == "DEBUG")
        return logDEBUG;
    if (level == "INFO")
        return logINFO;
    if (level == "WARNING")
        return logWARNING;
    if (level == "ERROR")
        return logERROR;
    Log<T>().Get(logWARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return logINFO;
}

class Output2FILE
{
public:
    static FILE*& Stream();
    static void Output(const std::string& msg);
};

inline FILE*& Output2FILE::Stream()
{
    static FILE* pStream = stderr;
    return pStream;
}

inline void Output2FILE::Output(const std::string& msg)
{
    FILE* pStream = Stream();
    if (!pStream)
        return;
    fprintf(pStream, "%s", msg.c_str());
    fflush(pStream);
}
class FILELog : public Log<Output2FILE> {};
//typedef Log<Output2FILE> FILELog;

#ifndef FILELOG_MAX_LEVEL
#define FILELOG_MAX_LEVEL logDEBUG4
#endif

#define FILE_LOG(level) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else if (level > FILELog::ReportingLevel() || !Output2FILE::Stream()) ; \
    else FILELog().Get(level)

#include <sys/time.h>

inline std::string NowTime()
{
    char buffer[11];
    time_t t;
    time(&t);
    tm r = {0};
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = {0};
    sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
}
#endif //__LOG_H__
