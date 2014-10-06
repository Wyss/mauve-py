#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnException.h"
#include <iostream>


using namespace std;
namespace genome {


void gnException::AddCaller(const char* const function)
{
	string func(function);
	string pretty_func = func.substr(0, func.find('(')+1);
	pretty_func += func.substr(func.rfind(')'), func.length());
	function_trace.push_back(pretty_func);
}

ostream& operator<<(ostream& os, const gnException& gne){ 
	//write exception to stream.
	os << "Exception " << gne.m_code.GetName() << " thrown ";
	list<string>::const_iterator func_iter = gne.function_trace.begin();

	//print the original function
	if(func_iter != gne.function_trace.end()){
		os << "from\n" << *func_iter << " in " << gne.m_file << " " << gne.m_line;
		func_iter++;
	}

	//print the call stack
	while(func_iter != gne.function_trace.end()){
		os << "\nCalled by " << *func_iter;
		func_iter++;
	}
	if(gne.m_message.length() > 0)
		os <<"\n" << gne.m_message;
	os << "\n";
	return os;
}

}	// end namespace genome

