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

void readinregion( string infile, map<string, map<pair<int, int >, int > > &region)
{
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	cout<<"read file "<<infile<<endl;
	string line;
	

	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		
		region[chr][make_pair(start, end)] = 0;
	}
	
	inf.close();
}

void readinbed( string infile, map<string, map<pair<int, int >, int > > &region, int &totalsize )
{
	int L = 10000;
	map<string, map<int, set<pair<int, int > > > > index_region;
	for ( map<string, map<pair<int, int >, int > >::iterator ite = region.begin();
		ite != region.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int idx1 = start / L;
			int idx2 = end / L;
			index_region[chr][idx1].insert( make_pair(start, end ) );
			if (idx2 > idx1 ) 
				index_region[chr][idx2].insert( make_pair(start, end) ); 
		}
	}
	
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	cout<<"read file "<<infile<<endl;
	string line;
	
	totalsize = 0;
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		int c = start + (end-start)/2;
		
		int idx = c / L;
		totalsize += 1;
		if ( totalsize % 100000 == 0 )
			cout<<totalsize<<endl;
		
		if ( index_region[chr].find( idx ) != index_region[chr].end() )
		{
			
			for ( set<pair<int, int > >::iterator ci = index_region[chr][idx].begin();
				ci != index_region[chr][idx].end(); ++ci )
			{
				if ( ci->second < c )
					continue;
				if ( ci->first > c )
					break;
				
				region[chr][*ci] += 1;
				break;
			}
		}
	}
	inf.close();
	
	
}

void output( string outfile, map<string, map<pair<int, int >, int > > &region, int &totalsize )
{
	if ( totalsize == 0 )
	{
		cout<<"error total size == 0 "<<endl;
		exit(1);
	}
	ofstream outf( outfile.data() );
	for ( map<string, map<pair<int, int >, int > >::iterator ite = region.begin();
		ite != region.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			double nc = si->second * 1000000.0 / totalsize;
			outf<<chr<<"_"<<si->first.first<<"_"<<si->first.second<<"\t"<<nc<<endl;
		}
	}
	outf.close();
	
}

int main( int argc, char * argv[] )
{
	if ( argc == 1 )
	{
		cout<<"count tag density on regions"<<endl;
		cout<<"Usage: prog peakfile bedfile outfile"<<endl;
		exit(1);
	}
	
	string regionfile = argv[1];
	string bedfile = argv[2];
	string outfile = argv[3];
	
	cout<<"read in region"<<endl;
	map<string, map<pair<int, int >, int > > region;
	readinregion( regionfile, region);
	
	cout<<"count tag"<<endl;
	int totalsize = 0;
	readinbed( bedfile, region, totalsize );
	
	cout<<"output"<<endl;
	output( outfile, region, totalsize );
	
	return 0;
}






