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

void rm_RD_pets( string infile, string outfile )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
//	map<string, map<string, string > > r1_r2_line;
	ofstream outf( outfile.data() );
	map<string, map<int, map<string, set<string> > > > chr1_index1_r1_r2;
	int L = 100000;
	int i = 0;
	int rdc = 0;
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		i += 1;
		if ( i % 1000000 == 0)
			cout<<i<<endl;
			
		vector<string> ps = parse_string( line );
		
		string chr1 = ps[0];
		string chr2 = ps[3];
		int start1 = atoi(ps[1].c_str() );
		int end1 = atoi(ps[2].c_str() );
		int start2 = atoi(ps[4].c_str() );
		int end2 = atoi(ps[5].c_str() );
		int c1 = start1;
		if ( ps[7] == "-" )
			c1 = end1;
		int c2 = start2;
		if ( ps[8] == "-" )
			c2 = end2;
			
		
		
		int index1 = c1 / L;
		int index2 = c2 / L;
		
		
		string r1 = ps[0]+"_"+ps[1]+"_"+ps[7];
		if (ps[7] == "-" )
			r1 = ps[0]+"_"+ps[2]+"_"+ps[7];
		string r2 = ps[3]+"_"+ps[4]+"_"+ps[8];
		if ( ps[8] == "-")
			r2 = ps[3]+"_"+ps[5]+"_"+ps[8];
		
		if ( chr1_index1_r1_r2[chr1][index1][r1].find( r2 ) 
			!= chr1_index1_r1_r2[chr1][index1][r1].end() )
		{
			rdc += 1;
			continue;
		}
	/*	if ( r1_r2_line[r1].find( r2 ) != r1_r2_line[r1].end() )
		{
			rdc += 1;
			continue;
		}  */
	//	r1_r2_line[r1][r2] = line;
	//	r1_r2_line[r2][r1] = line;
		chr1_index1_r1_r2[chr1][index1][r1].insert(r2);
		chr1_index1_r1_r2[chr2][index2][r2].insert(r1);
		outf<<line<<endl;
		
	}
	outf.close();
	inf.close();
	cout<<"redundant reads "<<rdc<<endl;
}

int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"filter redundant bedpe"<<endl;
		cout<<"Usage: prog inbedpe output"<<endl;
		exit(1);
	}
	string infile = argv[1];
	string outfile = argv[2];
	rm_RD_pets( infile, outfile );
	
	return 0;
}
