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

vector<string > parse_string( string & instr, char spl )
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == spl )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

vector<string > parse_string( string & instr)
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == '\t' || instr[i] == ' ' )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

void readintracloop( string infile, 
	map<string, map<pair<int, int >, set<pair<int, int > > > > &loops,
	map<string, set<pair<int, int > > > &anchors )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	int i = 0;
//	map<string, map<string, string > > r1_r2_line;
	getline(inf, line);
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string > ps = parse_string( line );
		
		string chr = ps[0];
		int start1 = atoi(ps[1].c_str() );
		int end1 = atoi(ps[2].c_str() );
		int start2 = atoi(ps[3].c_str() );
		int end2 = atoi(ps[4].c_str() );
		
		
		loops[chr][make_pair(start1, end1)].insert( make_pair(start2, end2 ) );
		anchors[chr].insert( make_pair(start1, end1 ) );
		anchors[chr].insert( make_pair(start2, end2 ) );
	}
	inf.close();
		
}

void merge_anchors( map<string, set<pair<int, int > > > &anchors,
	map<string, set<pair<int, int > > > &merged_anchors,
	map<string, map<pair<int, int >, pair<int, int > > > &anchor_map )
{
	
	for ( map<string, set<pair<int, int > > >::iterator ite = anchors.begin();
		ite != anchors.end(); ++ite )
	{
		string chr = ite->first;
		set<pair<int, int > > usedsets;
		set<pair<int, int > > ms;
		map<pair<int, int >, pair<int, int > > ms_map;
		set<pair<int, int > >::iterator si = ite->second.begin();
		int start = si->first;
		int end = si->second;
		usedsets.insert( *si );
		++si;
		for ( ; si != ite->second.end(); ++si )
		{
			if ( si->first <= end )
			{
				usedsets.insert( *si );
				end = si->second;
			} else
			{
				ms.insert( make_pair(start, end ) );
				for ( set<pair<int, int > >::iterator ui = usedsets.begin(); 
					ui != usedsets.end(); ++ui )
				{
					ms_map[*ui] = make_pair( start, end );
				}
				usedsets.clear();
				
				start = si->first;
				end = si->second;
				usedsets.insert( *si );
			}
		}
		ms.insert( make_pair(start, end ) );
		for ( set<pair<int, int > >::iterator ui = usedsets.begin(); 
			ui != usedsets.end(); ++ui )
		{
			ms_map[*ui] = make_pair( start, end );
		}
		merged_anchors[chr] = ms;
		anchor_map[chr] = ms_map;
	}
}

void match_merged_loop( map<string, map<pair<int, int >, set<pair<int, int > > > > &loops,
	map<string, map<pair<int, int >, pair<int, int > > > &anchor_map,
	map<string, map<pair<int, int >, set<pair<int, int > > > > &merged_loops )
{
	for ( map<string, map<pair<int, int >, set<pair<int, int > > > >::iterator ite = loops.begin();
		ite != loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, set<pair<int, int > > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > a1 = si->first;
			if ( anchor_map[chr].find( a1 ) == anchor_map[chr].end() )
			{
				cout<<"error cannot find map anchor "<<chr<<" "<<a1.first<<" "<<a1.second<<endl; exit(1);
			}
			pair<int, int > ma1 = anchor_map[chr][a1];
			for ( set<pair<int, int > >::iterator ti = si->second.begin();
				ti != si->second.end(); ++ti )
			{
				pair<int, int > a2 = *ti;
				if ( anchor_map[chr].find( a2 ) == anchor_map[chr].end() )
				{
					cout<<"error cannot find map anchor "<<chr<<" "<<a2.first<<" "<<a2.second<<endl; exit(1);
				}
				pair<int, int > ma2 = anchor_map[chr][a2];
				
				if ( ma2.first > ma1.second+1000 )
				{
					merged_loops[chr][ma1].insert( ma2);
				}
			}
		}
	}
}

void output_loops( string outfile, map<string, map<pair<int, int >, set<pair<int, int > > > > &loops )
{
	ofstream outf( outfile.data() );
	
	for ( map<string, map<pair<int, int >, set<pair<int, int > > > >::iterator ite = loops.begin(); 
		ite != loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, set<pair<int, int > > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			for (set<pair<int, int > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				
				outf<<chr<<"\t"<<si->first.first<<"\t"<<si->first.second<<"\t"<<ti->first<<"\t"<<ti->second<<endl;
			}
		}
	}
	outf.close();
}

void output_anchors( string outfile, map<string, set<pair<int, int > > > &anchors )
{
	ofstream outf( outfile.data() );
	
	for ( map<string, set<pair<int, int > > >::iterator ite = anchors.begin();
		ite != anchors.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			outf<<chr<<"\t"<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outf.close();
}

int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"merge traloop"<<endl;
		cout<<"Usage: prog outmergedloop outmergedanchor infile1 infile2 ..."<<endl;
		exit(1);
	}
	
	string outfile1 = argv[1];
	string outfile2 = argv[2];
	
	map<string, map<pair<int, int >, set<pair<int, int > > > > loops;
	map<string, set<pair<int, int > > > anchors;
	
	for ( int i = 3; i <= argc-1; ++i )
	{
		string infile = argv[i];
		readintracloop( infile, loops, anchors );
		
	}
	map<string, set<pair<int, int > > > merged_anchors;
	map<string, map<pair<int, int >, pair<int, int > > > anchor_map;
	merge_anchors( anchors, merged_anchors, anchor_map );
	
	map<string, map<pair<int, int >, set<pair<int, int > > > > merged_loops;
	match_merged_loop( loops, anchor_map, merged_loops );
	
	output_loops( outfile1, merged_loops );
	output_anchors( outfile2, merged_anchors );
	
	return 0;
}




