#include "libMems/IntervalList.h"
#include "libMems/Backbone.h"

using namespace mems;
using namespace genome;
using namespace std;

int main (int ARGC,char ** ARGV) {

  IntervalList input_alignment;
  ifstream align_file;
  if(ARGC != 4){
    cout <<"Usage:\nbackbone_global_to_local <xmfa file> <backbone file> <output file>\n";
    return 0;
  }
  align_file.open(ARGV[1]);
  if(!align_file.is_open()){
    cerr <<"Couldn't read xmfa file: "<<ARGV[1]<<"\n";
  }
  input_alignment.ReadStandardAlignment(align_file);
  LoadSequences(input_alignment,&cout);
  ifstream backbone_file;
  backbone_file.open(ARGV[2]);
  if(!backbone_file.is_open()){
    cerr <<"Couldn't read backbone file: "<<ARGV[2]<<"\n";
  }
  
  ofstream new_backbone(ARGV[3]);
  if(!align_file.is_open()){
    cerr <<"Couldn't write to output file: "<<ARGV[3]<<"\n";
  }

  vector< bb_seqentry_t > backbone_struct;
  readBackboneSeqFile(backbone_file, backbone_struct);

  for(int i=0; i < backbone_struct.size(); i++){
    for(int j=0; j < backbone_struct[i].size(); j++){
      uint64 start = absolut(backbone_struct[i][j].first);
      uint64 end = absolut(backbone_struct[i][j].second);
      uint32 contig_num1;
      uint32 contig_num2;
      if(start == 0){
	contig_num1=0;
	contig_num2=0;
      }else{
	input_alignment.seq_table[j]->globalToLocal(contig_num1,start);
        input_alignment.seq_table[j]->globalToLocal(contig_num2,end);
      }

      if(contig_num1 != contig_num2){
	//cerr <<"Not the same contig!" <<contig_num1 <<" "<<contig_num2;
      }
      if(j>0){
	new_backbone<<"\t";
      }
      new_backbone<<contig_num1 <<":"<<start<<"\t"<<contig_num2 <<":"<<end;
    }
    new_backbone <<"\n";
  }
}
