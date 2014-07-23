/////////////////////////////////////////////////////////////////////////////
// File:            gnFastTranslator.h
// Purpose:         Filter for all Sequences
// Description:     translates, converts sequence
// Changes:        
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFastTranslator.h"
#include <iostream>
#include <cstring>


using namespace std;
namespace genome {


//	static data access, avoids static initialization order fiasco
const gnFastTranslator *gnFastTranslator::ProteinDNATranslator(){
	const static gnFastTranslator* t_trans = new gnFastTranslator(ProteinDNATranslatorType);
	return t_trans;
}
const gnFastTranslator *gnFastTranslator::DNAProteinTranslator(){
	const static gnFastTranslator* t_trans = new gnFastTranslator(DNAProteinTranslatorType);
	return t_trans;
}

//	public:
gnFastTranslator::gnFastTranslator()
{
	use_default = false;
	m_defaultChar = 0;
}

gnFastTranslator::gnFastTranslator( const gnFastTranslator &sf )
{
	m_name = sf.m_name;
	use_default = sf.use_default;
	m_defaultChar = sf.m_defaultChar;
	m_transCache = sf.m_transCache;
}
gnFastTranslator::gnFastTranslator( gnTranslatorType t_type )
{
	use_default = false;
	m_defaultChar = 0;
	switch(t_type){
		case ProteinDNATranslatorType:
			CacheTranslator(gnTranslator::ProteinDNATranslator(), "FLIMVPTAY.HQNKDECGSR", 1);
			break;
		case DNAProteinTranslatorType:
			CacheTranslator(gnTranslator::DNAProteinTranslator(), "ACGTRYKMBVDHSWNX", 3);
			break;
	}
}

	// gnSeqC 
gnSeqC gnFastTranslator::Filter( const gnSeqC ch ) const{
/*	for(uint32 i=0; i < m_inputTable.size(); i++){
		if(m_inputTable[i].length() == 1)
			if(compare->Contains(m_inputTable[i][0], ch))
				return m_outputTable[i][0];
	}
*/	return m_defaultChar;
}

void gnFastTranslator::Filter( gnSeqC** seq, gnSeqI& len ) const{
/*	uint32 curpos = 0;
	string output;
	while(curpos < len){
		uint32 i=0;
		for(; i < m_inputTable.size(); i++){
			//don't compare if there aren't enough chars
			uint32 curlen = m_inputTable[i].length();
			if(len - curpos < curlen)
				continue;
			if(compare->Contains(m_inputTable[i].data(), *seq + curpos, curlen)){
				output += m_outputTable[i];
				curpos += curlen;
				break;
			}
		}
		if(i == m_inputTable.size()){
			//no match was found.  
			if(use_default)  //fill with the default char?
				output += m_defaultChar;
			curpos++;
		}
	}
	if(output.length() > len){
		delete[] *seq;
		*seq = new gnSeqC[output.length()];
	}
	len = output.length();
	memcpy(*seq, output.data(), len);
*/}
	// string
void gnFastTranslator::Filter( string &seq ) const{
	uint32 curpos = 0, outpos = 0;
	uint32 len = seq.length();
	uint32 width = m_transCache.begin()->first.length();
	uint32 out_width = m_transCache.begin()->second.length();
	uint32 out_size = (seq.length() / width) * out_width + seq.length() % width + 1;
	gnSeqC* output_array = new gnSeqC[out_size];
	output_array[out_size-1] = 0;
	string seq_upper;
	while(curpos < len){
		//transform to upper case
		seq_upper = seq.substr(curpos, width);
		for(uint32 i=0; i < seq_upper.size(); i++)
			seq_upper[i] = toupper(seq_upper[i]);
		
		map<string, string>::const_iterator iter = m_transCache.find(seq_upper);
		
		if(iter == m_transCache.end()){
			//no match was found.  
			if(use_default)  //fill with the default char?
				output_array[curpos] = m_defaultChar;
			curpos++;
		}else{
			strcpy( output_array + outpos, iter->second.c_str() );
			curpos += width;
			outpos += out_width;
		}
	}
	seq = output_array;
}

void gnFastTranslator::CacheTranslator(const gnTranslator* tranny, string inputs, const gnSeqI input_width){
	string cur_input;
	string cur_trans;
	vector<gnSeqI> index;
	gnSeqI cur_index = input_width;
	
	//fill the index array with input_width 0's
	for(gnSeqI curI = 0; curI < input_width; curI++)
		index.push_back(0);

	while(true){
		//ensure the validity of our indices
		cur_index = input_width - 1;
		while(index[cur_index] == inputs.length()){
			if(cur_index == 0){
				return;
			}
			index[cur_index] = 0;
			cur_index--;
			index[cur_index]++;
			continue;
		}
		
		//create a sequence to cache.
		for(gnSeqI i = 0; i < input_width; i++){
			cur_input += inputs[index[i]];
		}
		cur_trans = cur_input;
		tranny->Filter(cur_trans);
		m_transCache[cur_input] = cur_trans;
//		m_transCache.insert(map<string, string>::value_type(cur_input, cur_trans));
		// prepare for next time thru the loop
		cur_input = "";
		index[input_width - 1]++;
	}
}

}	// end namespace genome

