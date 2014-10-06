#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnException_h_
#define _gnException_h_

#include "libGenome/gnClone.h"
#include "libGenome/gnExceptionCode.h"
#include <string>
#include <list>


namespace genome {

class GNDLLEXPORT gnException
{
public:
	gnException( const char* const srcFile,
				 const unsigned int srcLine,
				 const char* const function,
				 gnExceptionCode& code, 
				 const char* const message );

	friend std::ostream& operator<<(std::ostream& os, const gnException& gne); //write to source.
	gnExceptionCode& GetCode(){return m_code;}
	std::string GetMessage(){return m_message;}
	void AddCaller(	const char* const function );
protected:
	gnExceptionCode& m_code;
	std::string m_message;
	const char* const m_file;
	unsigned int m_line;
	std::list<std::string> function_trace;
};

inline
gnException::gnException(const char* const srcFile, const unsigned int srcLine, const char* const function, gnExceptionCode& code, const char* const message ) : 
m_code(code), m_message(message), m_file(srcFile), m_line(srcLine){
	AddCaller(function);
}


GNDLLEXPORT 
std::ostream& operator<<(std::ostream& os, const gnException& gne); //write to source.


#ifdef __GNDEBUG__
#define STACK_TRACE_START try{
#define STACK_TRACE_END \
}catch(genome::gnException& e){\
	e.AddCaller(__PRETTY_FUNCTION__);\
	throw e;\
}
#else
#define STACK_TRACE_START
#define STACK_TRACE_END
#endif // ifdef __DEBUG__

#define Throw_gnEx(code) throw genome::gnException(__FILE__, __LINE__, __PRETTY_FUNCTION__, code, "")
#define Throw_gnExMsg(code,msg) throw genome::gnException(__FILE__, __LINE__, __PRETTY_FUNCTION__, code, msg)


}	// end namespace genome

#endif //_gnException_h_
