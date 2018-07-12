#ifndef PET_H
#define PET_H

#include "operation.h"

class PET
{
public:
	string chr1;
	int start1;
	int end1;
	char strand1;
	string chr2;
	int start2;
	int end2;
	int len;
	char strand2;
	int q;
};

class PET_bank
{
public:
	vector< PET > pet_ve;
	map<string, map<pair<int, int >, size_t > > chr_frg_pet;
	map<string, map<int, set<pair<int, int > > > > chr_index_pos;
	
	void readinPET( string infile, int down, int up, map<string, int > &chrlen );
	void readinPET( string infile, int down, int up );
	void readinPET( string infile ); 
	
	void generate_pos_map();
	
	void filter_HQNR_pet(int qlim);
	void filter_NR_pet();
	
	void getsubgroupbylen( vector<PET> &subpet, int down, int up );
	
	void getoverlappedPET_region( vector<size_t > &pet_ovl, map<string, set<pair<int, int > > > &regions );
	// at least one end of tag overlap with region
	void getoverlappedPET_region2( vector<size_t > &pet_ovl, map<string, set<pair<int, int > > > &regions, int down, int up );
	// at least one end of tag overlap with region, further require that intra-chr and length limits
	
	void outputPET( string outfile, vector<size_t > &pet );
	void outputPET( string outfile );
	
	void generateindex();
	
};

class Loop_bank
{
public:
	map<string, set<pair<int, int > > > peaks;
	map<string, map<int, set<pair<int, int > > > > chr_index_peaks;
	PET_bank petbank;
	
	map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > looped_peak;
	map<string, map<pair<int, int >, int > > peak_height;
	
	void readinpeak( string infile );
	void getloopedpeak( int minP );
	void generateindex();
	
	
	void getcluster_loop_by_pet( int win, int gap );
	void getcluster_height_from_loopedpeak();
	 
	void cal_overlap( vector< PET > &pet_ve_o );
};

void cal_overlap( map<string, map<pair<int, int >, int > > &peak1,
	map<string, map<pair<int, int >, int > > &peak2 );

void cal_overlap( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > &looped_peak1,
	map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > &looped_peak2 );

#endif
