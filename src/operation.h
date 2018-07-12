#ifndef OPERATION_H
#define OPERATION_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>

using namespace std;

typedef size_t Region_id;
typedef size_t Sample_id;

vector<string > parse_string( string & instr, char spl );
vector<string > parse_string( string & instr );

class Position
{
public:
	int start;
	int end;
	string chr;
	Position()
	{
	}
	Position(int inst, int inend, string inchr)
	{
		start = inst;
		end = inend;
		chr = inchr;
	}
	Position(string str)
	{
		vector<string > splstr1 = parse_string(str, ':');
		if ( splstr1.size() != 2 )
		{
			cout<<"error parse position: str "<< str<<endl;
			exit(1);
		}
		chr = splstr1[0];
		vector<string > splstr2 = parse_string(splstr1[1], '+');
		if ( splstr2.size() != 2 )
		{
			cout<<"error parse position: str "<< str<<endl;
			exit(1);
		}
		start = atoi( splstr2[0].c_str() );
		end = atoi( splstr2[1].c_str() );
	}

	void addpos(int inst, int inend, string inchr);
	void addpos(string str);


};


double log_2(double r );

void callendensity( vector<int > &len_ve, map<int, double > &len_den, int bin );

bool overlaptest(Position &p1, Position &p2);

bool overlaptest(Position &p1, Position &p2, int extension);

bool overlaptest(pair<int, int> p1, pair<int, int> p2);

bool overlaptest(pair<int, int> p1, pair<int, int> p2, int extension);

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap );

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve );

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve, int extension );

bool overlaptest( int site1, map<int, Region_id > &posmap, int extension, vector<int > &dis );

int matrixtrans( int stage, int elem, vector<int > &matrix );
vector<int > revtransmatrix( int stage, int elem, int digit );

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, string > > &region_map );

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, Region_id > > &region_map );

void filter_region_remain( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void filter_region_minus( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void filter_region_minus_remain( map<string, vector< pair<int, int> > > &overlapped, 
				   map<string, vector< pair<int, int> > > &nonoverlapped,
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void stitch_region( map<string, vector<vector<pair<int, int> > > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance );
	

void region_merge_naive( vector< map<string, vector< pair<int, int> > > > &reg_ve,
	map<string, vector< pair<int, int> > > &merged );

void transformregion( map<string, vector< pair<int, int> > > &newregion,
	map<string, set<pair<int, int > > > &oldregion );

int calsize( map<string, vector<pair<int, int > > >& regions );
					 
string inttostr(int i );

void outputtable( ofstream &outf, vector<vector<int> > &mat, vector<string > &rowname, vector<string> &colname );

void outputtable( ofstream &outf, vector<vector<double> > &mat, vector<string > &rowname, vector<string> &colname );

vector<vector<double> > transposmat( vector<vector<double> > &mat );

void normalizemat( vector<vector<double > > &mat, vector<double > &thr );

int findnearTSSinchr( set<int > &tssset, pair<int, int > region );

int findnearTSSinchr( set<int > &tssset, pair<int, int > region, vector< int > &containedTSS );


set<Region_id > vetoset(vector<Region_id > &ve );

vector<pair<int, int > > mergeragions( set<pair<int, int > > &r );

int shift_tagpos( int start, int end, char strand, int seg_len );

void cal_column_average( vector<vector<double> > &mat, vector<double > &ave );  // rows are in the same size;
void cal_column_sum( vector<vector<double> > &mat, vector<double > &addve ); // rows are in the same size;

void generate_wig( string outfile, map<string, map<int, int > > &chr_pos_sc, string name, int span );
void generate_wig( string outfile, map<string, map<int, int > > &chr_pos_sc, string name, int span, string tchr );
void generate_wig( string outfile, map<string, map<int, double > > &chr_pos_sc, string name, int span );
void generate_wig( string outfile, map<string, map<int, double > > &chr_pos_sc, string name, int span, string tchr );

void readwigfile( string infile, map<string, map<int, double > > &score_map );
void readwigfile_bunch( string infile, map<string, map<int, double > > &score_map );
void readwigandchangeres( string infile, map<string, map<int, double > > &score_map, int span );
void readwigandchangeres_bunch( string infile, map<string, map<int, double > > &score_map, int span );



void cal_summits_from_ve( vector<double > &ve, map<int, double > &summits );
void cal_summits_from_map( map<int, int> &m, map<int, int > &summits );
void cal_summits_from_map( map<int, double> &m, map<int, double > &summits );
void cal_summits_from_map_win( map<int, int> &m, map<int, int > &summits, int win );
void cal_summits_from_map_win( map<int, double> &m, map<int, double > &summits, int win );
void cal_bottoms_from_ve( vector<double > &ve, map<int, double > &bottoms );
void cal_boundaries_for_peaks( map<int, double > &summits, map<int, double > &bottoms, int leftsite, int rightsight, map<int, pair<int, int > > &peak_boundary );
void merge_tiny_peaks( map<int, pair<int, int > > &peak_boundary, map<int, double > &summit_score, map<int, double > &bottom_score, int win_str ) ;
	
void get_boundary_summit_plus1( map<int, pair<int, int > > &peak_boundary, int oricenter, int &upbound, int &downbound, int &summit );

void cal_density_and_win_smooth( map<int, int >& dis_count, map<int, double > &den, int halfw );

void smooth_ve( vector<double > &ve, int halfw );
void cal_deviation_peak( map<int, pair<int, int > > &peak_boundary, vector< vector<double > > &score_mat, map<int, double > &peak_deviation, int length, int extend );

map<int, double > find_max( vector<map<int, double > > close_map, int gap );
void remove_nearby_summit( map<string, map<int, double > > &chr_summit_map, map<string, map<int, double > > &rm_chr_summit_map );
void remove_nearby_summit( map<int, double > &summit_map, map<int, double > &rm_summit_map, int distance );

int getmean(vector<int> & );
double getmean(vector<double > & );
int getn50(vector<int> & );
int getmedium(vector<int > & );
int getmax(vector<int > & );
int getsum(vector<int > & ); 
pair<double, double> getmeanstd(vector<int > &);
pair<double, double> getmeanstd(vector<double > &ve);
double getstd( vector<double > &ve );
double getCV( vector<double > &ve );

double entropy( vector<double > &Pi );

#endif

