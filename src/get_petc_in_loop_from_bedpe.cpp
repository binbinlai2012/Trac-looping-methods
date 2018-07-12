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

void readinloop( string infile, map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loops )
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
		
		loops[chr][make_pair(start1, end1)][make_pair(start2, end2)] = 0;
		
	}
	inf.close();
		
}


void filter_pets_in_loops( string infile, map< string, map<pair<int, int >, map<pair<int, int >, int >  > > &loops )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	int c = 0;
	
	
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		
		
		string chr1 = ps[0];
		int start1 = atoi(ps[1].c_str() );
		int end1 = atoi(ps[2].c_str() );
		
		string chr2 = ps[3];
		int start2 = atoi(ps[4].c_str() );
		int end2 = atoi(ps[5].c_str() );
		
		c += 1;
		if ( chr1 != chr2 )
			continue;
		int c1 = start1 + (end1 - start1)/2;
		int c2 = start2 + (end2 - start2)/2;
		if ( ( c2 - c1 ) < 1000)
			continue;
		if ( c % 100000 == 0 )
			cout<<c<<endl;
			
	//	if ( chr1=="chr3" && c1 < 37026876 && c1 > 37022876 && c2 > 37054748 && c2 < 37058748)
		//	cout<<"+1"<<endl;
		for ( map<pair<int, int >, map<pair<int, int >, int >  >::iterator ite = loops[chr1].begin();
			ite != loops[chr1].end(); ++ite )
		{
			pair<int, int > s = ite->first;
			if ( s.second < start1 )
				continue;
			if ( s.first > end1 )
				break;
			
			for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				pair<int, int > e = si->first;
				if ( e.second < start2 )
					continue;
				if ( e.first > end2 )
					break;
				
			//	outf<<line<<endl;
				si->second += 1;
				break;
			}
		//	break;
		}
	}
	
	inf.close();
}

void output( string outfile, map< string, map<pair<int, int >, map<pair<int, int >, int >  > > &loops )
{
	ofstream outf( outfile.data() );
	for ( map< string, map<pair<int, int >, map<pair<int, int >, int >  > >::iterator ite = loops.begin(); 
		ite != loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, int >  >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > s = si->first;
			for ( map<pair<int, int >, int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				pair<int, int > e = ti->first;
				int c = ti->second;
				outf<<chr<<"\t"<<s.first<<"\t"<<s.second<<"\t"<<e.first<<"\t"<<e.second<<"\t"<<c<<endl;
			}
		}
	}
}

int main( int argc,char* argv[] )
{
	if ( argc == 1)
	{
		cout<<"call pet count in loop from bedpe"<<endl;
		cout<<"Usage: prog loopfile(header=TRUE) bedpefile outfile"<<endl;
		exit(1);
	}
	
	string loopfile = argv[1];
	
	string bedpefile = argv[2];
	string outfile = argv[3];
	
	cout<<"read loop file"<<endl;
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > loops;
	readinloop( loopfile, loops );
	
	cout<<"read bedpe and get count"<<endl;
	filter_pets_in_loops( bedpefile, loops );
	
	cout<<"output"<<endl;
	output( outfile, loops );
	
	return 0;
}
	









