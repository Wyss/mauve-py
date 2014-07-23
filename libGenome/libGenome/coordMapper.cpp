#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include <iostream>
#include <fstream>
#include <algorithm>

struct ExMem{
	gnSeqI length;
	int64 ex_start;
	int64 in_start;

};

static boolean LocationLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetStart() < b.GetStart();
}

static boolean LocationEndLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationEndLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetEnd() < b.GetStart();
}

static boolean LocationSizeLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationSizeLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetEnd() - a.GetStart() < b.GetEnd() - b.GetStart();
}

static boolean ExMemLessthan(const ExMem& a, const ExMem& b){
	if(a.ex_start == b.ex_start){
		if(a.in_start == b.in_start){
			return a.length < b.length;
		}
		return a.in_start < b.in_start;
	}
	return a.ex_start < b.ex_start;
}

void print_usage(char* pname){
	cout << "Usage: " << pname << " <genbank file> <exon_list> <intron_list> <mem_list> <regulation_network> <minimum match size>\n";
}

void load_map_file(string& filename, vector<gnLocation>& loc_list){
	ifstream coord_file;
	coord_file.open(filename.c_str());
	if(!coord_file.is_open()){
		cout << "couldn't open file\n";
		return;
	}

	//read all the exon coordinates
	while(coord_file.good()){
		gnSeqI start_coord, end_coord;
		coord_file >> start_coord;
		coord_file >> end_coord;
		
		gnLocation new_location(start_coord, end_coord);
		loc_list.push_back(new_location);
	}
}

void map_coordinates(vector<gnLocation>& loc_list, vector<gnLocation>& map_list){
	gnSeqI curStart = 1;
	for(uint32 locationI = 0; locationI < loc_list.size(); locationI++){
		gnSeqI cur_len = loc_list[locationI].GetEnd() - loc_list[locationI].GetStart() + 1;
		gnLocation map_location(curStart, curStart+cur_len-1);
		map_list.push_back(map_location);
		curStart += cur_len;
	}
}

void map_list_coordinates(list<gnLocation>& loc_list, list<gnLocation>& map_list){
	gnSeqI curStart = 1;
	list<gnLocation>::iterator loc_iter = loc_list.begin();
	for(; loc_iter != loc_list.end(); loc_iter++){
		gnSeqI cur_len = loc_iter->GetEnd() - loc_iter->GetStart() + 1;
		gnLocation map_location(curStart, curStart+cur_len-1);
		map_list.push_back(map_location);
		curStart += cur_len;
	}
}

void print_feature(ostream& os, gnBaseFeature* cur_feat){
	os << cur_feat->GetName() << ": \n";
	for(uint32 i=0; i < cur_feat->GetQualifierListLength(); i++)
		os << cur_feat->GetQualifierName(i) << "\t" << cur_feat->GetQualifierValue(i) << "\n";
	os << "\n";

}

uint32 loc_binary_search(vector<gnLocation>& loc_list, uint32 startI, uint32 endI, gnLocation& query_loc){
	uint32 midI = ((endI - startI) / 2) + startI;
	if(startI == endI)
		return endI;

	if(loc_list[midI].Intersects(query_loc)){
		return midI;
	}else if(loc_list[midI].GetStart() < query_loc.GetStart())
		return loc_binary_search(loc_list, midI + 1, endI, query_loc);
	else
		return loc_binary_search(loc_list, startI, midI, query_loc);
}

int main(int argc, char* argv[]){
	
	boolean run_interactive = false;
	string seq_filename;
	string exon_list_filename;
	string intron_list_filename;
	string mem_list_filename;
	string reg_net_filename;
	vector<gnLocation> exon_list;
	vector<gnLocation> intron_list;
	vector<gnLocation> exon_map_list;
	vector<gnLocation> intron_map_list;
	vector<ExMem> mem_list;
	uint32 minimum_match_size;

	// check for correct calling semantics
	if(argc != 7){
		print_usage(argv[0]);
		return -1;
	}

	seq_filename = argv[1];
	exon_list_filename = argv[2];
	intron_list_filename = argv[3];
	mem_list_filename = argv[4];
	reg_net_filename = argv[5];
	minimum_match_size = atoi(argv[6]);
	
	ifstream mem_file(mem_list_filename.c_str());
	if(!mem_file.is_open()){
		cout << "Error opening " << mem_list_filename << "\n";
		return -1;
	}

	if(run_interactive){
		cout << "Give the name of the exon list to search\n";
		cin >> exon_list_filename;
		cout << "Give the name of the intron list to search\n";
		cin >> intron_list_filename;
		cout << "Give the name of the regulatory network to output\n";
		cin >> reg_net_filename;

	}
	ofstream net_file(reg_net_filename.c_str());
	
	if(!net_file.is_open()){
		cout << "Error opening regulatory network file: " << reg_net_filename << "\n";
		return -2;
	}
	
	load_map_file(exon_list_filename, exon_list);
	load_map_file(intron_list_filename, intron_list);
	
	cout << exon_list.size() << " unique exons loaded from file\n";
	cout << intron_list.size() << " unique introns loaded from file\n";
	
	//now load the genbank file
	gnSequence seq_file;
	if(run_interactive){
		cout << "Enter the name of the genbank sequence file you are using\n";
		cin >> seq_filename;
	}
	if(!seq_file.LoadSource(seq_filename)){
		cout << "Error loading file\n";
		return -1;
	}
	cout << "Sequence loaded successfully, " << seq_file.length() << " base pairs.\n";
	
	//construct a mapping between coordinates...
	map_coordinates(exon_list, exon_map_list);
	map_coordinates(intron_list, intron_map_list);
	
	//now read the mem file
	while(mem_file.good()){
		ExMem m;
		mem_file >> m.length;
		mem_file >> m.ex_start;
		mem_file >> m.in_start;
		if(m.length >= minimum_match_size)
			mem_list.push_back(m);
	}
	cout << mem_list.size() << " matches loaded.\n";
	//sort the mem list
	sort(mem_list.begin(), mem_list.end(), &ExMemLessthan);
	
	//now get the intersection for each mem in the list...
	uint32 exonmapI = 0;
	uint32 notify_percent = 10;
	uint32 notify_interval = mem_list.size() / notify_percent;
	uint32 percent_count = 0;
	cout << "Searching for complementary matches:\n";
	for(uint32 memI = 0; memI < mem_list.size(); memI++){
		if(memI % notify_interval == 0){
			cout << percent_count << "%.. ";
			percent_count += notify_percent;
		}
		//simple linear search for intersecting exon mapping
		gnLocation ex_map_loc(mem_list[memI].ex_start, mem_list[memI].ex_start + mem_list[memI].length - 1);
		for(; exonmapI < exon_map_list.size(); exonmapI++){
			if(exon_map_list[exonmapI].Intersects(ex_map_loc))
				break;
		}
		
		//continue to search for any other mappings that intersect
		uint32 mapEnd = exonmapI;
		for(; mapEnd < exon_map_list.size(); mapEnd++){
			if(!exon_map_list[exonmapI].Intersects(ex_map_loc))
				break;
		}
		mapEnd--;
		
		uint32 intronmapI = 0; // intronmapI will contain the index of the first intersecting intron
		//do a binary search for intersecting intron mappings
		int64 cur_in_start = mem_list[memI].in_start;
		if(cur_in_start < 0)
			cur_in_start = -cur_in_start;
		gnLocation in_map_loc(cur_in_start, cur_in_start + mem_list[memI].length - 1);
		uint32 search_mapI = loc_binary_search(intron_map_list, 0, intron_map_list.size()-1, in_map_loc);
		intronmapI = search_mapI;

		//search backwards for previous intersections
		for(; intronmapI >= 0; intronmapI--){
			if(!intron_map_list[intronmapI].Intersects(in_map_loc)){
				intronmapI++;
				break;
			}
			if(intronmapI == 0)
				break;
		}
		//continue to search for any other mappings that intersect
		uint32 intron_mapEnd = search_mapI;
		for(; intron_mapEnd < intron_map_list.size(); intron_mapEnd++){
			if(!intron_map_list[intronmapI].Intersects(in_map_loc))
				break;
		}
		intron_mapEnd--;
		
		//we have the mappings, now map the coordinates
		vector<uint32> ex_feat_index, in_feat_index;
		vector<gnBaseFeature*> ex_feat_list;
		vector<gnBaseFeature*> in_feat_list;
		gnSeqI cur_match_len = mem_list[memI].length;
		
		//find out how much of the first exon was matched
		//extra exon start has the number of unmatched bases at the beginning of the exon
		gnSeqI extra_exon_start = mem_list[memI].ex_start - exon_map_list[exonmapI].GetStart();
		gnSeqI cur_exon_len = exon_map_list[exonmapI].GetEnd() - exon_map_list[exonmapI].GetStart() + 1;
		gnSeqI max_exon_chunk = cur_exon_len - extra_exon_start;
		gnSeqI cur_exon_chunk = max_exon_chunk < cur_match_len ? max_exon_chunk : cur_match_len;
		
		//find out how much of the first intron was matched
		gnSeqI extra_intron_start, cur_intron_len, max_intron_chunk, cur_intron_chunk;
		boolean complement = false;
		if(mem_list[memI].in_start > 0){
			extra_intron_start = mem_list[memI].in_start - intron_map_list[intronmapI].GetStart();
			cur_intron_len = intron_map_list[intronmapI].GetEnd() - intron_map_list[intronmapI].GetStart() + 1;
			max_intron_chunk = cur_intron_len - extra_intron_start;
			cur_intron_chunk = max_intron_chunk < mem_list[memI].length ? max_intron_chunk : mem_list[memI].length;
		}else{
			//reverse complement, start at the end.
			if(cur_in_start >= intron_map_list[intron_mapEnd].GetStart())
				cur_intron_chunk = cur_match_len;
			else
				cur_intron_chunk = cur_in_start + cur_match_len - intron_map_list[intron_mapEnd].GetStart();
			complement = true;
			seq_file.getIntersectingFeatures(intron_list[intronmapI], in_feat_list, in_feat_index);
		}
		
		//the current chunk will be the smaller of the two mappings
		gnSeqI cur_chunk = cur_intron_chunk < cur_exon_chunk ? cur_intron_chunk : cur_exon_chunk;

		gnLocation cur_exon_loc(exon_list[exonmapI].GetStart() + extra_exon_start, exon_list[exonmapI].GetStart() + cur_chunk);
		seq_file.getIntersectingFeatures(cur_exon_loc, ex_feat_list, ex_feat_index);
		if(mem_list[memI].in_start > 0){
			gnLocation cur_intron_loc(intron_list[intronmapI].GetStart() - 1, intron_list[intronmapI].GetEnd() + 1);
			seq_file.getIntersectingFeatures(cur_intron_loc, in_feat_list, in_feat_index);
		}else{
			gnLocation cur_intron_loc(intron_list[intron_mapEnd].GetStart() - 1, intron_list[intron_mapEnd].GetEnd() + 1);
			seq_file.getIntersectingFeatures(cur_intron_loc, in_feat_list, in_feat_index);
		}
		vector<gnBaseFeature*> ex_forward;
		vector<gnBaseFeature*> ex_reverse;
		vector<gnBaseFeature*> in_forward;
		vector<gnBaseFeature*> in_reverse;

		for(uint32 featI = 0; featI < ex_feat_list.size(); featI++){
			string featName = ex_feat_list[featI]->GetName();
			if(featName == "mRNA" || featName == "CDS" || featName == "gene" )
				if(ex_feat_list[featI]->GetLocationType() == gnLocation::LT_Complement)
					ex_reverse.push_back(ex_feat_list[featI]);
				else
					ex_forward.push_back(ex_feat_list[featI]);
		}
		for(uint32 featI = 0; featI < in_feat_list.size(); featI++){
			string featName = in_feat_list[featI]->GetName();
			if(featName == "mRNA" || featName == "CDS" || featName == "gene" )
				if(in_feat_list[featI]->GetLocationType() == gnLocation::LT_Complement)
					in_reverse.push_back(in_feat_list[featI]);
				else
					in_forward.push_back(in_feat_list[featI]);
		}

		if(complement){
			vector<gnBaseFeature*> tmp_in = in_forward;
			in_forward = in_reverse;
			in_reverse = tmp_in;
		}

		//now print out all the complementary features
		if((ex_forward.size() > 0 && in_reverse.size() > 0) || (in_forward.size() > 0 && ex_reverse.size() > 0)){
			net_file << "================================\n";
			net_file << "Mem: " << mem_list[memI].length << "\n";
			net_file << "This exon/intron matching size: " << cur_chunk << "\n";
		}
		if(ex_forward.size() > 0 && in_reverse.size() > 0){
			net_file << "Forward Exons:\n";
			for(uint32 featI = 0; featI < ex_forward.size(); featI++)
				print_feature(net_file, ex_forward[featI]);
			net_file << "Matching introns:\n";
			for(uint32 featI = 0; featI < in_reverse.size(); featI++)
				print_feature(net_file, in_reverse[featI]);
		}
		if(in_forward.size() > 0 && ex_reverse.size() > 0){
			net_file << "Reverse Exons:\n";
			for(uint32 featI = 0; featI < ex_reverse.size(); featI++)
				print_feature(net_file, ex_reverse[featI]);
			net_file << "Matching introns:\n";
			for(uint32 featI = 0; featI < in_forward.size(); featI++)
				print_feature(net_file, in_forward[featI]);
		}
		
		//release memory
		for(uint32 featI = 0; featI < ex_feat_list.size(); featI++)
			delete ex_feat_list[featI];
		for(uint32 featI = 0; featI < in_feat_list.size(); featI++)
			delete in_feat_list[featI];
		
		//loop while there is stuff to match
//		while(cur_match_len > 0){
		
//		}
	}
}

