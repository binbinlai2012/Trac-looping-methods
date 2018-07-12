#include "pet.h"

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

void readin_distance_dist0( string infile, map<int, double > &d_p )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	
	vector<string > ps = parse_string( line );
	int d = atoi(ps[0].c_str());
	double p = atof( ps[1].c_str() );
	d_p.insert( make_pair(d, p ) );
	inf.close();
	
}

void readin_distance_dist( string infile, map<int, double > &d_p )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		int d = atoi( ps[0].c_str() );
		double p = atof( ps[1].c_str() );
		d_p.insert( make_pair(d, p ) );
	}
	inf.close();
	
}

void readin_accm_dist( string infile, map<int, double > &d_p )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		int d = atoi( ps[0].c_str() );
		double p = atof( ps[1].c_str() );
		d_p.insert( make_pair(d, p ) );
	}
	inf.close();
}

void readin_accm_PETdist( string infile, map<int, double > &d_p )
{
	ifstream inf(infile.data() );

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
		int d = atoi( ps[0].c_str() );
		double p = atof( ps[1].c_str() );
		d_p.insert( make_pair(d, p ) );
	}
	inf.close();
}

void readin_total_space_size( string infile, long &space_size )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	space_size = atoi(line.c_str() );
	inf.close();
	cout<<"space size "<<space_size<<endl;
}

void readin_total_pets( string infile, int &pets_size )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	pets_size = atoi(line.c_str() );
	inf.close();
	cout<<"pets size "<<pets_size<<endl;
}

void readin_region( string infile, map<string, set<pair<int, int > > > &peaks, int win )
{
	ifstream inf(infile.data() );

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
		int start = atoi( ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		int c = start + (end-start+1)/2;
		start = c - win/2;
		end = c + win/2;
		peaks[chr].insert( make_pair( start, end ) );
	}
	inf.close();
}

void assign_pet_to_peak_pair( PET_bank & pet_bank, 
	map<string, set<pair<int, int > > > &peaks,
	map<string, map<pair<int, int >, double > > &peak_acc,
	map<string, vector<pair<pair<int,int>, pair<int, int > > > > &peak_pair, int minl, int maxl, int win, int total_c )
{
	map<string, map<int, set<pair<int, int > > > > chr_index_peak;
	for ( map<string, set<pair<int, int > > >::iterator ite = peaks.begin(); ite != peaks.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int index1 = si->first / win;
			int index2 = si->second / win;
			for ( int i = index1; i <= index2; ++i )
			{
				chr_index_peak[chr][i].insert( *si );
			}
		}
	}
	
	cout<<"total ini pets "<<pet_bank.pet_ve.size()<<endl;
	pair<int, int > tag1 = make_pair( 6055900, 6056900 );
	pair<int, int > tag2 = make_pair( 6093500, 6094500 );
	int ck = 0;
	map<string, map<pair<int, int >, int > > peak_petc;
	for ( map<string, set<pair<int, int > > >::iterator ite = peaks.begin();
		ite != peaks.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{	
			peak_petc[chr][*si] = 0;
		}
	}
	
	for ( size_t i = 0; i < pet_bank.pet_ve.size(); ++i )
	{
	/*	if ( i % 100000 == 0 )
		{
			cout<<"i"<<i<<endl;
		} */
		if ( pet_bank.pet_ve[i].chr1 == "chrM" || pet_bank.pet_ve[i].chr2 == "chrM" )
			continue;
		if ( pet_bank.pet_ve[i].chr1 != pet_bank.pet_ve[i].chr2 )
			continue;
		if ( pet_bank.pet_ve[i].len <= minl || pet_bank.pet_ve[i].len > maxl )
			continue;
		
		string chr = pet_bank.pet_ve[i].chr1;
		
		vector<pair<int, int > > left_hit_peak;
		vector<pair<int, int > > right_hit_peak;
		
		
		int pos1 = pet_bank.pet_ve[i].start1;
		if ( pet_bank.pet_ve[i].strand1 == '-' )	
			pos1 = pet_bank.pet_ve[i].end1;
		string chr1 = pet_bank.pet_ve[i].chr1;
		int index1 = pos1 / win;
		
		if ( chr_index_peak[chr1].find(index1) != chr_index_peak[chr1].end())
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peak[chr1][index1].begin();
				ci != chr_index_peak[chr1][index1].end(); ++ci )
			{
				if ( ci->first > pos1 )
					break;
				if ( ci->first <= pos1 && ci->second >= pos1 )
				{
					left_hit_peak.push_back( *ci );
				}
			}
		}
		
	//	int p_start = (pos1 / chr_Par[chr1].win) * chr_Par[chr1].win;
	//	int id1 = p_start / chr_Par[chr1].win;
		
		int pos2 = pet_bank.pet_ve[i].start2;
		if ( pet_bank.pet_ve[i].strand2 == '-' )	
			pos2 = pet_bank.pet_ve[i].end2;
		string chr2 = pet_bank.pet_ve[i].chr2;
		int index2 = pos2 / win;
		
		if ( pos2 < pos1 )
		{
			cout<<"error pos2 < pos1 "<<endl; exit(1);
		}
		
		if ( chr_index_peak[chr2].find(index2) != chr_index_peak[chr2].end())
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peak[chr2][index2].begin();
				ci != chr_index_peak[chr2][index2].end(); ++ci )
			{
				if ( ci->first > pos2 )
					break;
				if ( ci->first <= pos2 && ci->second >= pos2 )
				{
					right_hit_peak.push_back( *ci );
				}
			}
		}
		
		if ( !left_hit_peak.empty() )
		{
			for ( vector<pair<int, int > >::iterator ai = left_hit_peak.begin();
				ai != left_hit_peak.end(); ++ai )
			{
				peak_petc[chr][*ai] += 1;
			}
		}
		if ( !right_hit_peak.empty() )
		{
			for ( vector<pair<int, int > >::iterator bi = right_hit_peak.begin();
				bi != right_hit_peak.end(); ++bi )
			{
				peak_petc[chr][*bi] += 1;
			}
		}
		
		if ( !left_hit_peak.empty() && !right_hit_peak.empty() )
		{
			for ( vector<pair<int, int > >::iterator ai = left_hit_peak.begin();
				ai != left_hit_peak.end(); ++ai )
			{
				for ( vector<pair<int, int > >::iterator bi = right_hit_peak.begin();
					bi != right_hit_peak.end(); ++bi )
				{
					peak_pair[chr].push_back( make_pair( *ai, *bi) );
					if ( chr == "chr10" && *ai == tag1 && *bi == tag2 )
					{
						cout<<"Hit 1"<<endl;
					}
				}
			}
			ck += 1;
		}
		

		
		
	}
	
	cout<<ck<<endl;
	
	
	for ( map<string, map<pair<int, int >, int > >::iterator ite = peak_petc.begin(); ite != peak_petc.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = si->second;
			double r = c * 1000000.0 / total_c;
			peak_acc[chr][si->first] = r;
		}
	}
	
	
	
}

void assign_pet_to_peak_pair2( PET_bank & pet_bank, 
	map<string, set<pair<int, int > > > &peaks,
	map<string, map<pair<int, int >, double > > &peak_acc,
	map<string, map<pair<pair<int,int>, pair<int,int> >, int > > &peak_pair_c,
	int minl, int maxl, int win, int total_c )
{
	map<string, vector<pair<pair<int,int>, pair<int, int > > > > peak_pair;

	assign_pet_to_peak_pair( pet_bank, peaks, peak_acc, peak_pair, minl, maxl, win, total_c );
	
	
	for ( map<string, vector<pair<pair<int,int>, pair<int, int > > > >::iterator ite = peak_pair.begin();
		ite != peak_pair.end(); ++ite )
	{
		string chr = ite->first;
		
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( peak_pair_c[chr].find( ite->second[i] ) != peak_pair_c[chr].end() )
			{
				peak_pair_c[chr][ite->second[i]] += 1;
			} else
				peak_pair_c[chr][ite->second[i]] = 1;
			
		}  
	}
	
	map<int, int > rpm_c;
	for ( map<string, map<pair<pair<int,int>, pair<int,int> >, int > >::iterator ite = peak_pair_c.begin(); 
		ite != peak_pair_c.end(); ++ite )
	{	
		for ( map<pair<pair<int,int>, pair<int,int> >, int >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			double rpm1 = peak_acc[ite->first][si->first.first];
			double rpm2 = peak_acc[ite->first][si->first.second];
			double rpm = rpm1*rpm2;
			int logrpm = (int)(log2(rpm*1.0+1) * 10);
			if ( rpm_c.find(logrpm) != rpm_c.end() )
				rpm_c[logrpm] += 1;
			else
				rpm_c[logrpm] = 1;
		}
	}
	ofstream outf( "tmp.rpm_dist.txt" );
	for ( map<int, int >::iterator ite= rpm_c.begin(); ite != rpm_c.end(); ++ite )
		outf<<ite->first<<"\t"<<ite->second<<endl;
	outf.close();
}

void complete_dist_p( map<int, double > &dist_p, int skip )
{
	map<int, double > comp_dist_p;
	int minv = dist_p.begin()->first;
	if ( skip == 1 )
	{
		map<int, double >::iterator ite = dist_p.begin();
		comp_dist_p.insert(*ite);
		++ite;
		minv = ite->first;
	}
	int maxv = dist_p.rbegin()->first;
	
	for ( int v = minv; v <= maxv; ++v )
	{
		if ( dist_p.find( v ) != dist_p.end() )
			comp_dist_p[v] = dist_p[v];
		else
		{
			int v_f = v;
			do {
				--v_f;
				if ( v_f < minv )
					break;
			} while ( dist_p.find(v_f) == dist_p.end() );
			if ( v_f < minv )
			{
				cout<<"error v_f < minv "<<v_f<<" "<<minv<<endl;
				exit(1);
			}
			int v_b = v;
			do {
				++v_b;
				if ( v_b > maxv )
					break;
			} while ( dist_p.find(v_b) == dist_p.end() );
			if ( v_b > maxv )
			{
				cout<<"error v_b > maxv "<<v_b<<" "<<maxv<<endl; exit(1);
			}
			double a = min(dist_p[v_f], dist_p[v_b]);
			double b = max(dist_p[v_f], dist_p[v_b]);
			double y = a + ( v - v_f )*1.0*(b - a) / (v_b - v_f);
			if ( y == 0 )
			{
				cout<<a<<" "<<b<<" "<<y<<" "<<v<<" "<<v_f<<" "<<v_b<<"\t"<<minv<<"\t"<<maxv<<endl;
			}
			comp_dist_p[v] = y; 
		}
	}
	dist_p = comp_dist_p;
}

void output_d_p( string outfile, map<int, double > &d_p )
{
	ofstream outf( outfile.data() );
	for ( map<int, double >::iterator ite= d_p.begin(); ite != d_p.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
}

void complete_alldist_p( map<int, double > &dist_p,
	map<int, double > &dist_PETs_p,
	map<int, double > &rpm_p,
	map<int, double > &rpm_PETs_p )
{
	cout<<"comp distance p"<<endl;
	complete_dist_p( dist_p, 1 );
	
	output_d_p( "comp_distance_p.txt", dist_p );
	cout<<"comp dist PEts p"<<endl;
	complete_dist_p( dist_PETs_p, 1 );
	
	output_d_p("comp_distance_PETs_p.txt", dist_PETs_p );
	cout<<"comp rpm p"<<endl;
	complete_dist_p( rpm_p, 0 );
	output_d_p("comp_rpm_p.txt", rpm_p );
	cout<<"comp rpm PETs p"<<endl;
	complete_dist_p( rpm_PETs_p, 0 );
	output_d_p("comp_rpm_PETs_p.txt", rpm_PETs_p );
}

void cal_pv( map<string, map<pair<pair<int,int>, pair<int,int> >, int > > &peak_pair_c,
	map<string, map<pair<int, int >, double > > &peak_acc,
	map<int, double > &dist_p,
	map<int, double > &dist_PETs_p,
	map<int, double > &rpm_p,
	map<int, double > &rpm_PETs_p,
	int minl, int maxl, int win, long total_space, int valid_PETs_in_anchor, int min_PETs, 
	string outfile )
{
	ofstream outf( outfile.data() );
	int lowrpm = 0;
	int lowrpmv = rpm_p.begin()->first;
	outf<<"chr\tpeak1_s\tpeak1_e\tpeak2_s\tpeak2_e\trpm1\trpm2\trpm_product\tdistance\tPets\tl_20\tlogl\tlogrpm\tlp1\tlp2\tdp1\tdp2\tprobability\tpvalue\n";
	for ( map<string, map<pair<pair<int,int>, pair<int,int> >, int > >::iterator ite = peak_pair_c.begin();
		ite != peak_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<pair<int,int>, pair<int,int> >, int >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int c = si->second;
			if ( c < min_PETs )
				continue;
			int pos1 = si->first.first.first+(si->first.first.second - si->first.first.first) / 2;
			int pos2 = si->first.second.first + (si->first.second.second - si->first.second.first ) / 2;
			int length = pos2 - pos1;
			if ( length < minl || length >= maxl )
				continue;
			int l = length / win;
			int minstep = minl / win;
			int maxstep = maxl / win;
			if ( l < minstep || l >= maxstep )
				continue;
			int l_20 = l / 20;
			int loglength = (int)(log2(l*1.0)*10);
			
			double p1 = 0;
			if ( dist_p.find( l_20 ) != dist_p.end() )
			{
				p1 = dist_p[l_20];
			} else
			{
				cout<<"warning cannot find dist_p for "<<length<<" "<<l_20<<endl; exit(1);
			}
			double p2 = 0;
			if ( dist_PETs_p.find( loglength ) != dist_PETs_p.end() )
			{
				p2 = dist_PETs_p[loglength];
			} else
			{
				cout<<"warning cannot find dist_PETs_p for "<<length<<" "<<loglength<<endl; exit(1);
			}
			
			if ( peak_acc[chr].find(si->first.first ) == peak_acc[chr].end() )
			{
				cout<<"error peak_acc cannot find "<<chr<<" "<<si->first.first.first<<" "<<si->first.first.second<<endl; exit(1);
			}
			if ( peak_acc[chr].find(si->first.second ) == peak_acc[chr].end() )
			{
				cout<<"error peak_acc cannot find "<<chr<<" "<<si->first.second.first<<" "<<si->first.second.second<<endl; exit(1);
			}
			
			double rpm1 = peak_acc[chr][si->first.first];
			double rpm2 = peak_acc[chr][si->first.second];
			double rpm = rpm1*rpm2;
			int logrpm = (int)(log2(rpm*1.0+1) * 10);
			if ( logrpm > 50 )
				logrpm = 50;  
			if ( logrpm < lowrpmv )
			{
				lowrpm += 1;
				continue;	
			} 
			double rp1 = 0;
			if ( rpm_p.find(logrpm) != rpm_p.end() )
			{
				rp1 = rpm_p[logrpm];
			} else
			{
				if ( logrpm > rpm_p.rbegin()->first )
				{
					rp1 = rpm_p.rbegin()->second;
				//	cout<< "continue"<<endl;
				} else
				{
					cout<<"warning cannot find rpm_p for "<<rpm1<<" "<<rpm2<<" "<<rpm<<" "<<logrpm<<endl; 
					exit(1);
				}
			}
			double rp2 = 0;
			if ( rpm_PETs_p.find(logrpm) != rpm_PETs_p.end() )
			{
				rp2 = rpm_PETs_p[logrpm];
			} else
			{
				if ( logrpm > rpm_PETs_p.rbegin()->first )
				{
					rp2 = rpm_PETs_p.rbegin()->second;
				//	cout<<"continue"<<endl;
				} else
				{
					cout<<"warning cannot find rpm_PETs_p for "<<rpm1<<" "<<rpm2<<" "<<rpm<<" "<<logrpm<<endl; 
					exit(1);
				}
				
			}
			
			double p = (p2 * rp2) / (p1 * rp1 * total_space );
			
			double pv = cdf( c, valid_PETs_in_anchor, p );
			
			outf<<chr<<"\t"<<si->first.first.first<<"\t"<<si->first.first.second<<"\t"<<si->first.second.first<<"\t"<<si->first.second.second
				<<"\t"<<rpm1<<"\t"<<rpm2<<"\t"<<rpm<<"\t"<<length<<"\t"<<c<<"\t"<<l_20<<"\t"<<loglength<<"\t"<<logrpm<<"\t"<<p1<<"\t"<<p2<<"\t"<<rp1<<"\t"<<rp2<<"\t"<<p<<"\t"<<pv<<endl;
			
		}
	}
	outf.close();
	cout<<"lowrpm "<<lowrpm<<endl;
}

void exit_with_help()
{
	cerr <<"call significant interaction"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-R		region file"<<endl;
	cerr <<"-b		bedpe file" <<endl;
	cerr <<"-w		bin size (default: 1K)"<<endl;
	cerr <<"-p		output prefix file input" <<endl;
	cerr <<"-l		lower bound pet length (default: 2K)"<<endl;
	cerr <<"-u		upper bound pet length (default: 2M)"<<endl;
	cerr <<"-c		minimal pet count in loop for consideration (default: 3)"<<endl;
	cerr <<"-0		distance_dist0 file"<<endl;
	cerr <<"-1		distance dist file"<<endl;
	cerr <<"-2		distance PETs dist0 file"<<endl;
	cerr <<"-3		distance PETs dist file"<<endl;
	cerr <<"-4		accm dist file"<<endl;
	cerr <<"-5		accm PETs dist file"<<endl;
	cerr <<"-6		total_space_file"<<endl;
	cerr <<"-7		total_PETs_file"<<endl;
	cerr <<"-8		total PETs_in_anchor file"<<endl;
	exit(1);
}

void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}


int main(int argc, char* argv[] )
{
	string regionfile = "";
	string bedpefile = "";
	string outprefix = "";
	int win = 1000;
	int minl = 2000;
	int maxl = 2000000;
	int min_count = 3;
	string l_dist0file = "";
	string l_distfile = "";
	string l_petdist0file = "";
	string l_petdistfile = "";
	string d_distfile = "";
	string d_petdistfile = "";
	string spacefile = "";
	string totalpetsfile = "";
	string totalpetsanchorfile = "";
	if (argc == 1)
	{
		exit_with_help();
	}
	
	for(int i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			exit_with_help("Options must start with \'-\'.");

		if(argv[i][2] != '\0')
			exit_with_help("The option should be exactly one letter.");
		int option = argv[i][1];

		i++;

		switch(option)
		{
		
		case 'R':
			regionfile = argv[i];
			break;
		case 'b':
			bedpefile = argv[i];
			break;
		case 'p':
			outprefix = argv[i];
			break;
		case 'w':
			win = atoi(argv[i]);
			break;
		case 'l':
			minl = atoi(argv[i]);
			break;
		case 'u':
			maxl = atoi(argv[i]);
			break;
		case 'c':
			min_count = atoi(argv[i]);
			break;
		case '0':
			l_dist0file = argv[i];
			break;
		case '1':
			l_distfile = argv[i];
			break;
		case '2':
			l_petdist0file = argv[i];
			break;
		case '3':
			l_petdistfile = argv[i];
			break;
		case '4':
			d_distfile = argv[i];
			break;
		case '5':
			d_petdistfile = argv[i];
			break;
		case '6':
			spacefile = argv[i];
			break;
		case '7':
			totalpetsfile = argv[i];
			break;
		case '8':
			totalpetsanchorfile = argv[i];
			break;
		default:
			exit_with_help();
		
		}
	}
	
	cout<<"read in length p and density p and space size"<<endl;
	map<int, double > l_p;
	readin_distance_dist0( l_dist0file, l_p );
	readin_distance_dist( l_distfile, l_p );
	map<int, double > l_pet_p;
	readin_distance_dist0( l_petdist0file, l_pet_p );
	readin_distance_dist( l_petdistfile, l_pet_p );
	map<int, double > d_p;
	readin_accm_dist( d_distfile, d_p );
	map<int, double > d_pet_p;
	readin_accm_PETdist( d_petdistfile, d_pet_p );
	long space_size= 1;
	readin_total_space_size( spacefile, space_size );
	int total_pet_count = 0;
	readin_total_pets(totalpetsfile, total_pet_count );
//	int total_pet_in_anchor = 0;
//	readin_total_pets(totalpetsanchorfile, total_pet_in_anchor );
	
	cout<<"read in region"<<endl;
	map<string, set<pair<int, int > > > peaks;
	readin_region( regionfile, peaks, win );
	
	cout<<"read in pets"<<endl;
	PET_bank pet_bank;
	pet_bank.readinPET( bedpefile, minl, maxl );
	
	cout<<"assign pet to peak pair"<<endl;
	map<string, map<pair<int, int >, double > > peak_acc;
	map<string, map<pair<pair<int,int>, pair<int,int> >, int > > peak_pair_c;

	assign_pet_to_peak_pair2( pet_bank, peaks, peak_acc, peak_pair_c, minl, maxl, win, total_pet_count );
	
	cout<<"complete probability"<<endl;
	complete_alldist_p( l_p, l_pet_p, d_p, d_pet_p );
	
	cout<<"call p value"<<endl;
	string outfile = outprefix+".interaction.summary";
//	cal_pv( peak_pair_c, peak_acc, l_p, l_pet_p, d_p, d_pet_p, minl, maxl, win, space_size, total_pet_in_anchor, min_count, outfile );
	cal_pv( peak_pair_c, peak_acc, l_p, l_pet_p, d_p, d_pet_p, minl, maxl, win, space_size, total_pet_count, min_count, outfile );
	
	return 0;
}





