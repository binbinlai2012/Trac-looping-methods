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

void readinpetc( string infile, map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc, int &total_petc)
{
	ifstream inf(infile.data() );
	total_petc = 0;

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
		int petc = atoi(ps[5].c_str() );
		loop_petc[chr][make_pair(start1, end1)][make_pair(start2, end2)] = petc;
		total_petc += petc;
	}
	inf.close();
		
}

void readin_acc( string infile, map<string, map<pair<int, int>, double> >  &anchor_acc )
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
//	getline(inf, line);
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string > ps = parse_string( line );
		vector<string > pos_ve = parse_string( ps[0], '_');
		string chr = pos_ve[0];
		int start = atoi(pos_ve[1].c_str() );
		int end = atoi(pos_ve[2].c_str() );
		double acc = atof( ps[1].c_str() );
		
		anchor_acc[chr][make_pair(start, end)] = acc;
	}
	inf.close();
}

void assign_acc( map<string, map<pair<int, int>, double>  > &anchor_acc, 
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc,
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > &loop_acc )
{
	for ( map<string, map<pair<int, int >, map<pair<int, int >, int > > >::iterator ite = loop_petc.begin();
		ite != loop_petc.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, int > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > s = si->first;
			if ( anchor_acc[chr].find(s) == anchor_acc[chr].end() )
			{
				cout<<"error cannot find anchor acc "<<chr<<"\t"<<s.first<<"\t"<<s.second<<endl;
				exit(1);
			}
			double acc1 = anchor_acc[chr][s];
			for ( map<pair<int, int >, int >::iterator ti = si->second.begin();
				ti != si->second.end(); ++ti )
			{
				pair<int, int > e = ti->first;
				if ( anchor_acc[chr].find(e) == anchor_acc[chr].end() )
				{
					cout<<"error cannot find anchor acc "<<chr<<"\t"<<e.first<<"\t"<<e.second<<endl;
					exit(1);
				}
				double acc2 = anchor_acc[chr][e];
				loop_acc[chr][s][e] = make_pair(acc1, acc2 );
			}
		}
	}
	
}

double BC( int n, int k )
{
	if ( k == 0 )
		return 1;
	double a = 1;
	for ( int i = 0; i <= (k-1); ++i )
	{
		a *= ( (n-i)*1.0 / (i+1) );
	}
	
	return a;
}

double cdf( int k, int n, double p )
{
	double v = 0;
	for ( int i = 0; i <= k; ++i )
	{
		double c = BC(n, i);
		double a = pow(p, i*1.0);
		double b = pow(1-p, (n-i)*1.0);
		v += c*a*b;
	//	cout<<v<<" "<<c<<" "<<a<<" "<<b<<" "<<c*a*b<<endl;	
	}	
	if ( v > 1 )
		v = 1;
	return 1-v;
}

/*
pair<double, double> cal_binomial_pvalue( double acc1, double acc2, int pet, double sec_acc1, double sec_acc2, int sec_pet, int pet_lib1, int pet_lib2 )
{
	if ( pet == 0 )
		pet = 1;
	if ( sec_pet == 0 )
		sec_pet = 1;
	double exp_secpet = pet * pow( ((sec_acc1+0.01)*(sec_acc2+0.01))/((acc1+0.01)*(acc2+0.01)), 0.5) * (pet_lib2*1.0/pet_lib1);
	
	double parp = exp_secpet / pet_lib2;
	
	double fd = sec_pet*1.0/exp_secpet;
	
	double pv = 1;
	
	if ( fd > 1 )
		pv = cdf( sec_pet, pet_lib2, parp );
	
	return make_pair( fd, pv); 
	
}

void comp( map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc_A, 
	int total_petc_A,
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc_B, 
	int total_petc_B,
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > &loop_acc_A,
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > &loop_acc_B,
	string outfile )
{
	
	ofstream outf( outfile.data() );
	
	outf<<"chr\tAnchor1_s\tAnchor1_e\tAnchor2_s\tAnchor2_e\tAnchors_Acc_A\tAnchors_Acc_B\tAnchors_AccFC_BvsA\tAnchors_AccFC_AvsB";
	outf<<"\tPet_A\tPet_B\tNormed_PetFC_BvsA\tPvalue_BvsA\tNormed_PetFC_AvsB\tPvalue_AvsB"<<endl;
	for ( map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > >::iterator ite = loop_acc_A.begin();
		ite != loop_acc_A.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, pair<double, double> > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > s = si->first;
			for ( map<pair<int, int >, pair<double, double> >::iterator ti = si->second.begin();	
				ti != si->second.end(); ++ti )
			{
				pair<int, int > e = ti->first;
				pair<double, double > acc_A = ti->second;
				pair<double, double > acc_B = loop_acc_B[chr][s][e];
				int petc_A = loop_petc_A[chr][s][e];
				int petc_B = loop_petc_B[chr][s][e];
				pair<double, double> p_B_vs_A = cal_binomial_pvalue( acc_A.first, acc_A.second, petc_A, acc_B.first, acc_B.second, petc_B, total_petc_A, total_petc_B );
				pair<double, double> p_A_vs_B = cal_binomial_pvalue( acc_B.first, acc_B.second, petc_B, acc_A.first, acc_A.second, petc_A, total_petc_B, total_petc_A );
				double c_acc_A = pow((acc_A.first+0.01)*(acc_A.second+0.01), 0.5);
				double c_acc_B = pow((acc_B.first+0.01)*(acc_B.second+0.01), 0.5);
				double acc_fc_BvsA = c_acc_B/c_acc_A;
				double acc_fc_AvsB = c_acc_A/c_acc_B;
				outf<<chr<<"\t"<<s.first<<"\t"<<s.second<<"\t"<<e.first<<"\t"<<e.second;
				outf<<"\t"<<c_acc_A<<"\t"<<c_acc_B<<"\t"<<acc_fc_BvsA<<"\t"<<acc_fc_AvsB;
				outf<<"\t"<<petc_A<<"\t"<<petc_B<<"\t"<<p_B_vs_A.first<<"\t"<<p_B_vs_A.second<<"\t"<<p_A_vs_B.first<<"\t"<<p_A_vs_B.second<<endl;
			}
		}
	}
}
*/

pair<double, double> cal_binomial_pvalue( double acc1, double acc2, int pet, double sec_acc1, double sec_acc2, int sec_pet, int pet_lib1, int pet_lib2 )
{
	if ( pet == 0 )
		pet = 1;
	if ( sec_pet == 0 )
		sec_pet = 1;
	double exp_secpet = pet * ((sec_acc1+0.01)*(sec_acc2+0.01))/((acc1+0.01)*(acc2+0.01)) * (pet_lib2*1.0/pet_lib1);
	
	double parp = exp_secpet / pet_lib2;
	
	double fd = sec_pet*1.0/exp_secpet;
	
	double pv = 1;
	
	if ( fd > 1 )
		pv = cdf( sec_pet, pet_lib2, parp );
	
	return make_pair( fd, pv); 
	
}

void comp( map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc_A, 
	int total_petc_A,
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > &loop_petc_B, 
	int total_petc_B,
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > &loop_acc_A,
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > &loop_acc_B,
	string outfile )
{
	
	ofstream outf( outfile.data() );
	
	outf<<"chr\tAnchor1_s\tAnchor1_e\tAnchor2_s\tAnchor2_e\tAnchors_Acc_A\tAnchors_Acc_B\tAnchors_AccFC_BvsA\tAnchors_AccFC_AvsB";
	outf<<"\tPet_A\tPet_B\tNormed_PetFC_BvsA\tPvalue_BvsA\tNormed_PetFC_AvsB\tPvalue_AvsB"<<endl;
	for ( map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > >::iterator ite = loop_acc_A.begin();
		ite != loop_acc_A.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, pair<double, double> > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > s = si->first;
			for ( map<pair<int, int >, pair<double, double> >::iterator ti = si->second.begin();	
				ti != si->second.end(); ++ti )
			{
				pair<int, int > e = ti->first;
				pair<double, double > acc_A = ti->second;
				pair<double, double > acc_B = loop_acc_B[chr][s][e];
				int petc_A = loop_petc_A[chr][s][e];
				int petc_B = loop_petc_B[chr][s][e];
				pair<double, double> p_B_vs_A = cal_binomial_pvalue( acc_A.first, acc_A.second, petc_A, acc_B.first, acc_B.second, petc_B, total_petc_A, total_petc_B );
				pair<double, double> p_A_vs_B = cal_binomial_pvalue( acc_B.first, acc_B.second, petc_B, acc_A.first, acc_A.second, petc_A, total_petc_B, total_petc_A );
				double c_acc_A = (acc_A.first+0.01)*(acc_A.second+0.01);
				double c_acc_B = (acc_B.first+0.01)*(acc_B.second+0.01);
				double acc_fc_BvsA = c_acc_B/c_acc_A;
				double acc_fc_AvsB = c_acc_A/c_acc_B;
				outf<<chr<<"\t"<<s.first<<"\t"<<s.second<<"\t"<<e.first<<"\t"<<e.second;
				outf<<"\t"<<c_acc_A<<"\t"<<c_acc_B<<"\t"<<acc_fc_BvsA<<"\t"<<acc_fc_AvsB;
				outf<<"\t"<<petc_A<<"\t"<<petc_B<<"\t"<<p_B_vs_A.first<<"\t"<<p_B_vs_A.second<<"\t"<<p_A_vs_B.first<<"\t"<<p_A_vs_B.second<<endl;
			}
		}
	}
}


int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"call differential interactions"<<endl;
		cout<<"Usage: prog pets_in_loop_A pets_in_loop_B Acc_in_loop_A Acc_in_loop_B outfile"<<endl;
		exit(1);
	}
	
	string Petsfile1 = argv[1];
	string Petsfile2 = argv[2];
	string Accfile1 = argv[3];
	string Accfile2 = argv[4];
	string outfile = argv[5];
	
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > loop_petc_A;
	int total_petc_A = 0;
	readinpetc( Petsfile1, loop_petc_A, total_petc_A);
	
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > loop_petc_B;
	int total_petc_B = 0;
	readinpetc( Petsfile2, loop_petc_B, total_petc_B);
	
	map<string, map<pair<int, int>, double> > anchor_acc_A;
	readin_acc( Accfile1, anchor_acc_A );
	map<string, map<pair<int, int>, double> > anchor_acc_B;
	readin_acc( Accfile2, anchor_acc_B );
	
	
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > loop_acc_A;
	assign_acc( anchor_acc_A,  loop_petc_A, loop_acc_A );
	
	map<string, map<pair<int, int >, map<pair<int, int >, pair<double, double> > > > loop_acc_B;
	assign_acc( anchor_acc_B,  loop_petc_B, loop_acc_B );
	
	
	comp( loop_petc_A, total_petc_A, loop_petc_B, total_petc_B, loop_acc_A, loop_acc_B, outfile );
	
	return 0;
}










