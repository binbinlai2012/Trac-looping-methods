#include "pet.h"

void readinloop( string infile, map<string, map<pair<int, int >, map<pair<int, int >,  double  > > > &loops )
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
		
		loops[chr][make_pair(start1, end1)][make_pair(start2, end2)] = atof(ps[19].c_str() );
		
	}
	inf.close();
		
}

void merge_overlaploops( map<string, map<pair<int, int >, map<pair<int, int >, double  > > > &loops,
	map<string, map<pair<int, int >, map<pair<int, int >, double > > > &m_loops )
{
	map<string, set<pair<int, int > > > partitions;
	for ( map<string, map<pair<int, int >, map<pair<int, int >, double  > > >::iterator ite = loops.begin();
		ite != loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >,  double  > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			partitions[chr].insert( si->first );
			for ( map<pair<int, int >,  double  >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				partitions[chr].insert( ti->first );
			}
		}
	}
	
	map<string, map<pair<int, int >, int > > ori_par_id;
	map<string, vector<pair<int, int > > > merged_par;
	for ( map<string, set<pair<int, int > > >::iterator ki = partitions.begin(); ki != partitions.end(); ++ki )
	{
		string chr = ki->first;
		set<pair<int, int > >::iterator ite = ki->second.begin();
		
		int start = ite->first;
		int end = ite->second;
		int i = 0;
		ori_par_id[chr][*ite] = i;
		++ite;
	
		for ( ; ite != ki->second.end(); ++ite )
		{
			if ( ite->first <= end )
			{
				end = ite->second;
				ori_par_id[chr][*ite] = i;
			} else
			{
				merged_par[chr].push_back( make_pair( start, end ) );
				start = ite->first;
				end = ite->second;
				i+=1;
				ori_par_id[chr][*ite] = i;
			}
		}
		merged_par[chr].push_back( make_pair(start, end ) );
	}
	
	
	for ( map<string, map<pair<int, int >, map<pair<int, int >, double > > >::iterator ite = loops.begin();
		ite != loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >,  double > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( ori_par_id[chr].find( si->first) == ori_par_id[chr].end() )
			{
				cout<<"error cannot find ori_par_id: "<<si->first.first<<" "<<si->first.second<<endl;
				exit(1);
			}
			int mid = ori_par_id[chr][si->first];
			pair<int,int > mp = merged_par[chr][mid];
			
			for ( map<pair<int, int >, double >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				if ( ori_par_id[chr].find( ti->first) == ori_par_id[chr].end() )
				{
					cout<<"error cannot find ori_par_id: "<<ti->first.first<<" "<<ti->first.second<<endl;
					exit(1);
				}
				int smid = ori_par_id[chr][ti->first];
				if ( smid == mid )
					continue; 
				pair<int, int > smp = merged_par[chr][smid];
				if ( m_loops[chr].find( mp ) != m_loops[chr].end() )
				{
					if ( m_loops[chr][mp].find( smp ) != m_loops[chr][mp].end() )
					{
						if ( ti->second < m_loops[chr][mp][smp] )
						{
							m_loops[chr][mp][smp] = ti->second;
						} 
					} else
					{
						m_loops[chr][mp][smp] = ti->second;
					}
				} else
				{
					m_loops[chr][mp][smp] = ti->second;
				}
			}
		}
	}
	
	
}

void assignpettoloop( vector< PET > &pet_ve, map<string, map<pair<int, int >, map<pair<int, int >, double > > > &m_loops,
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > &m_loops_c, int minl, int maxl )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_firstanchor;
	for ( map<string, map<pair<int, int >, map<pair<int, int >, double > > >::iterator ite = m_loops.begin();
		ite != m_loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, double > > ::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			pair<int, int > anc = si->first;
			int index1 = anc.first / L;
			int index2 = anc.second / L;
			chr_index_firstanchor[chr][index1].insert( anc );
			if ( index2 > index1 )
				chr_index_firstanchor[chr][index2].insert( anc );
			
		}
	}
	
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].chr1 == "chrM" || pet_ve[i].chr2 == "chrM" )
			continue;
		if ( pet_ve[i].chr1 != pet_ve[i].chr2 )
			continue;
		if ( pet_ve[i].len <= minl || pet_ve[i].len > maxl )
			continue;
		if ( i % 100000 ==  0)
		{
			cout<<i<<endl;
		}
		string chr = pet_ve[i].chr1;
		
		int pos1 = pet_ve[i].start1;
		if ( pet_ve[i].strand1 == '-' )	
			pos1 = pet_ve[i].end1;
		
		int pos2 = pet_ve[i].start2;
		if ( pet_ve[i].strand2 == '-' )	
			pos2 = pet_ve[i].end2;
		
		int startpos = pos1;
		int endpos = pos2;
		if ( pos1 > pos2 )
		{
			startpos = pos2;
			endpos = pos1;
		}
		
		int index = startpos / L;
		if ( chr_index_firstanchor[chr].find( index ) != chr_index_firstanchor[chr].end() )
		{
			bool fd = false;
			for ( set< pair<int, int > >::iterator ci = chr_index_firstanchor[chr][index].begin();
				ci != chr_index_firstanchor[chr][index].end(); ++ci )
			{
				if ( ci->first > startpos )
					break;
				if ( ci->first <= startpos && ci->second >= startpos )
				{
					for ( map<pair<int, int >, double >::iterator ti = m_loops[chr][*ci].begin();
						ti != m_loops[chr][*ci].end(); ++ti )
					{
						if ( ti->first.first > endpos )
							break;
						if ( ti->first.first <= endpos && ti->first.second >= endpos )
						{
							fd = true;
							
							if ( m_loops_c[chr][*ci].find( ti->first )  != m_loops_c[chr][*ci].end() )
							{
								m_loops_c[chr][*ci][ti->first] += 1;
							} else
							{
								m_loops_c[chr][*ci][ti->first] = 1;
							}
							
							break;
						}
					}
					if ( fd )
						break;
				} 
			}
		}
		
	}
	
}

void output( string prefix,
	map<string, map<pair<int, int >, map<pair<int, int >, double > > > &m_loops,
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > &m_loops_c,
	int libsize )
{
	string outfile1 = prefix+".merged_sig_interaction.table.txt";
	string outfile2 = prefix+".merged_sig_interaction.hammock";
	
	double fac = 1000000.0 / libsize;
	
	ofstream outf1( outfile1.data() );
	ofstream outf2( outfile2.data() );
	
	outf1<<"chr\tregion1_start\tregion1_end\tregion2_start\tregion2_end\tpets\tnorm_pets\tFDR"<<endl;
	
	int ID = 0;
	for ( map<string, map<pair<int, int >, map<pair<int, int >, double > > >::iterator ite = m_loops.begin(); 
		ite != m_loops.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map<pair<int, int >, double > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			for ( map<pair<int, int >, double >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int c = 0;
				if ( m_loops_c[chr][si->first].find( ti->first ) != m_loops_c[chr][si->first].end() )
					c = m_loops_c[chr][si->first][ti->first];
				if ( c == 0 )
				{
					cout<<"warning pet = 0 "<<chr<<" "<<si->first.first<<" "<<si->first.second<<" "<<ti->first.first<<" "<<ti->first.second<<endl;
				}
				double nc = c * fac;
				outf1<<chr<<"\t"<<si->first.first<<" "<<si->first.second<<" "<<ti->first.first<<" "<<ti->first.second;
				outf1<<"\t"<<c<<"\t"<<nc<<"\t"<<ti->second<<endl;
				
				
				outf2<<chr<<"\t"<<si->first.first<<"\t"<<ti->first.second;
				outf2<<"\tstruct:{thick:[["<<si->first.first<<","<<si->first.second<<"],["<<ti->first.first<<","<<ti->first.second<<"],],}"<<endl;
			}
		}
	}
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

int main( int argc, char* argv[] )
{
	cout<<"merge significant loops and re calculate pet counts"<<endl;
	if ( argc == 1 )
	{
		cout<<"Usage: prog ori_loops bedpefile minl maxL totalpetsizefile outprefix"<<endl;
		exit(1);
	}
	
	string orifile = argv[1];
	string bedpefile = argv[2];
	int minl = atoi(argv[3]);
	int maxl = atoi(argv[4]);
	string totalpetsfile = argv[5];
	string prefix = argv[6];
	
	cout<<"read in pets"<<endl;
	PET_bank pet_bank;
	pet_bank.readinPET( bedpefile, minl, maxl );
	
	cout<<"read in loops"<<endl;
	map<string, map<pair<int, int >, map<pair<int, int >,  double  > > > loops;
	readinloop( orifile, loops );
	
	cout<<"merge loops "<<endl;
	map<string, map<pair<int, int >, map<pair<int, int >, double > > > m_loops;
	merge_overlaploops( loops, m_loops );
	
	cout<<"assign pet "<<endl;
	map<string, map<pair<int, int >, map<pair<int, int >, int > > > m_loops_c;
	assignpettoloop( pet_bank.pet_ve, m_loops, m_loops_c, minl, maxl );
	
	int total_pet_count = 0;
	readin_total_pets(totalpetsfile, total_pet_count );
	
	cout<<"output"<<endl;
	output( prefix, m_loops, m_loops_c, total_pet_count );
	
	return 0;
}







