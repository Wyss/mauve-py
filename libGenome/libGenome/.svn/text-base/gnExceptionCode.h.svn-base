#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnExceptionCode_h_
#define _gnExceptionCode_h_

#include "libGenome/gnDefs.h"
#include <string>


namespace genome {

class GNDLLEXPORT gnExceptionCode{
public:
	gnExceptionCode(uint32 code, const char* name);
	boolean operator==(gnExceptionCode& gnec);
	boolean operator!=(gnExceptionCode& gnec);
	uint32 GetInt(){ return m_code; }
	std::string GetName(){ return m_name; }
private:
	gnExceptionCode();
	//prevent instances from being copied
	gnExceptionCode(const gnExceptionCode& gnec);
	gnExceptionCode& operator=(gnExceptionCode& gnec);
	uint32 m_code;
	std::string m_name;
};

inline
gnExceptionCode::gnExceptionCode(uint32 code, const char* name) :
m_code(code), m_name(name)
{}

inline
boolean gnExceptionCode::operator==(gnExceptionCode& gnec){
	return m_code == gnec.m_code;
}

inline
boolean gnExceptionCode::operator!=(gnExceptionCode& gnec){
	return m_code != gnec.m_code;
}

GNDLLEXPORT
uint32& GetNewExceptionCode();

inline
uint32& GetNewExceptionCode(){
	//static initializer is called only once
	static uint32 new_code = 0;
	//increment it each time the function is called
	new_code++;
	return new_code;
};

//Creates an exception code with the given name
//currently it chooses a unique id for each exception
//this may have to be changed in the future if the integer
//associated with each exception must be the same across compiles
#define CREATE_EXCEPTION(E_NAME) \
inline \
static genome::gnExceptionCode& E_NAME(){ \
	static genome::gnExceptionCode* m_excp = new genome::gnExceptionCode(genome::GetNewExceptionCode(), #E_NAME); \
	return *m_excp; \
}

//define a bunch of exception codes
//this must be done in a header file to work correctly

/**
 * Thrown when a generic array index is too large or otherwise invalid
 */
CREATE_EXCEPTION(IndexOutOfBounds)
/**
 * Thrown when a sequence index references an invalid coordinate
 */
CREATE_EXCEPTION(SeqIndexOutOfBounds)
/**
 * Thrown when a fragment index references an invalid fragment
 */
CREATE_EXCEPTION(FragmentIndexOutOfBounds)
/**
 * Thrown when a contig index references an invalid contig
 */
CREATE_EXCEPTION(ContigIndexOutOfBounds)
/**
 * Thrown when a header index references an invalid header
 */
CREATE_EXCEPTION(HeaderIndexOutOfBounds)
/**
 * Thrown when a spec index references an invalid spec
 */
CREATE_EXCEPTION(SpecIndexOutOfBounds)
/**
 * Thrown when a feature index references an invalid feature
 */
CREATE_EXCEPTION(FeatureIndexOutOfBounds)
/**
 * Thrown when a file can't be opened
 */
CREATE_EXCEPTION(FileNotOpened)
/**
 * Thrown when a URL can't be opened
 */
CREATE_EXCEPTION(URLNotFound)
/**
 * Thrown when a file's data is corrupt or unreadable
 */
CREATE_EXCEPTION(FileUnreadable)
/**
 * Thrown when an operation on a stream fails
 */
CREATE_EXCEPTION(IOStreamFailed)
/**
 * Thrown when an invalid pointer is given
 */
CREATE_EXCEPTION(NullPointer)



}	// end namespace genome

#endif  //_gnExceptionCode_h_
