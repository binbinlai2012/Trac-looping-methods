#include "pet.h"

void PET_bank::readinPET( string infile, int down, int up )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		PET pet;
		
		pet.chr1 = ps[0];
		pet.start1 = atoi(ps[1].c_str());
		pet.end1 = atoi(ps[2].c_str());
		pet.strand1 = ps[7][0];
		pet.chr2 = ps[3];
		pet.start2 = atoi(ps[4].c_str());
		pet.end2 = atoi(ps[5].c_str());
		pet.len = pet.end2 - pet.start1;
		pet.strand2 = ps[8][0];
		pet.q = atoi(ps[6].c_str() );
		if ( pet.chr1 == pet.chr2 && pet.len >= down && pet.len <= up )
		{
			pet_ve.push_back( pet );
			
		}
	}
	inf.close();
	
}

void PET_bank::readinPET( string infile, int down, int up, map<string, int > &chrlen )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		PET pet;
		
		pet.chr1 = ps[0];
		pet.start1 = atoi(ps[1].c_str());
		pet.end1 = atoi(ps[2].c_str());
		pet.strand1 = ps[7][0];
		pet.chr2 = ps[3];
		pet.start2 = atoi(ps[4].c_str());
		pet.end2 = atoi(ps[5].c_str());
		pet.len = pet.end2 - pet.start1;
		pet.strand2 = ps[8][0];
		pet.q = atoi(ps[6].c_str() );
		if ( chrlen.find( pet.chr1 ) == chrlen.end() || chrlen.find(pet.chr2) == chrlen.end() )
			continue;
		if ( pet.chr1 == pet.chr2 && pet.len >= down && pet.len <= up )
		{
			pet_ve.push_back( pet );
			
		}
	}
	inf.close();
	
}

void PET_bank::readinPET( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		PET pet;
		pet.chr1 = ps[0];
		pet.start1 = atoi(ps[1].c_str());
		pet.end1 = atoi(ps[2].c_str());
		pet.strand1 = ps[7][0];
		pet.chr2 = ps[3];
		pet.start2 = atoi(ps[4].c_str());
		pet.end2 = atoi(ps[5].c_str());
		pet.len = pet.end2 - pet.start1;
		pet.strand2 = ps[8][0];
		pet.q = atoi(ps[6].c_str() );
		pet_ve.push_back( pet );
			
		
	}
	inf.close();
	
}

void PET_bank::generate_pos_map( )
{
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		chr_frg_pet[pet_ve[i].chr1][make_pair(pet_ve[i].start1, pet_ve[i].end1) ] = i;
		chr_frg_pet[pet_ve[i].chr2][make_pair(pet_ve[i].start2, pet_ve[i].end2) ] = i;
		
	} 
}

void PET_bank::filter_HQNR_pet(int qlim)
{
	map <string, map<pair<int, int >, vector<size_t > > > hqnr;
	vector<size_t > hqnr_id;
	int rdc = 0;
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].q < qlim )
			continue;
		bool rd = false;
		if ( hqnr[pet_ve[i].chr1].find( make_pair( pet_ve[i].start1, pet_ve[i].end1 ) ) != hqnr[pet_ve[i].chr1].end() )
		{
			
			for ( size_t j = 0; j < hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )].size(); ++j )
			{
				size_t tagid = hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )][j];
				if ( pet_ve[i].chr1 == pet_ve[tagid].chr1 
					&& pet_ve[i].start1 == pet_ve[tagid].start1
					&& pet_ve[i].end1 == pet_ve[tagid].end1 
					&& pet_ve[i].chr2 == pet_ve[tagid].chr2 
					&& pet_ve[i].start2 == pet_ve[tagid].start2 
					&& pet_ve[i].end2 == pet_ve[tagid].end2 )
				{
					rd = true;
					break;
				} 
				if ( pet_ve[i].chr1 == pet_ve[tagid].chr2 
					&& pet_ve[i].start1 == pet_ve[tagid].start2
					&& pet_ve[i].end1 == pet_ve[tagid].end2 
					&& pet_ve[i].chr2 == pet_ve[tagid].chr1 
					&& pet_ve[i].start2 == pet_ve[tagid].start1 
					&& pet_ve[i].end2 == pet_ve[tagid].end1 )
				{
					rd = true;
					break;
				} 
			}
			
		} 
		
		if ( !rd )
		{
			hqnr_id.push_back( i );
			hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )].push_back( i );
			hqnr[pet_ve[i].chr2][make_pair( pet_ve[i].start2, pet_ve[i].end2 )].push_back( i );
		} else
		{
			rdc += 1;
		}
	}
	
	cout<<"RD count = "<<rdc<<endl;
	
	vector<PET> newpet;
	for ( size_t i = 0; i < hqnr_id.size(); ++i )
	{
		size_t id = hqnr_id[i];
		newpet.push_back( pet_ve[id] );
		
	}
	pet_ve = newpet;
	cout<<pet_ve.size()<<endl;
}

/*
void PET_bank::filter_NR_pet()
{
	map <string, map<pair<int, int >, vector<size_t > > > hqnr;
	vector<size_t > hqnr_id;
	int rdc = 0;
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].q < qlim )
	//		continue;
		bool rd = false;
		if ( hqnr[pet_ve[i].chr1].find( make_pair( pet_ve[i].start1, pet_ve[i].end1 ) ) != hqnr[pet_ve[i].chr1].end() )
		{
			
			for ( size_t j = 0; j < hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )].size(); ++j )
			{
				size_t tagid = hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )][j];
				if ( pet_ve[i].chr1 == pet_ve[tagid].chr1 
					&& pet_ve[i].start1 == pet_ve[tagid].start1
					&& pet_ve[i].end1 == pet_ve[tagid].end1 
					&& pet_ve[i].chr2 == pet_ve[tagid].chr2 
					&& pet_ve[i].start2 == pet_ve[tagid].start2 
					&& pet_ve[i].end2 == pet_ve[tagid].end2 )
				{
					rd = true;
					break;
				} 
				if ( pet_ve[i].chr1 == pet_ve[tagid].chr2 
					&& pet_ve[i].start1 == pet_ve[tagid].start2
					&& pet_ve[i].end1 == pet_ve[tagid].end2 
					&& pet_ve[i].chr2 == pet_ve[tagid].chr1 
					&& pet_ve[i].start2 == pet_ve[tagid].start1 
					&& pet_ve[i].end2 == pet_ve[tagid].end1 )
				{
					rd = true;
					break;
				} 
			}
			
		} 
		
		if ( !rd )
		{
			hqnr_id.push_back( i );
			hqnr[pet_ve[i].chr1][make_pair( pet_ve[i].start1, pet_ve[i].end1 )].push_back( i );
			hqnr[pet_ve[i].chr2][make_pair( pet_ve[i].start2, pet_ve[i].end2 )].push_back( i );
		} else
		{
			rdc += 1;
		}
	}
	
	cout<<"RD count = "<<rdc<<endl;
	
	vector<PET> newpet;
	for ( size_t i = 0; i < hqnr_id.size(); ++i )
	{
		size_t id = hqnr_id[i];
		newpet.push_back( pet_ve[id] );
		
	}
	pet_ve = newpet;
	cout<<pet_ve.size()<<endl;
}
*/

void PET_bank::getsubgroupbylen( vector<PET> &subpet, int down, int up )
{
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].chr1 == pet_ve[i].chr2 )
		{
			if ( pet_ve[i].len >= down && pet_ve[i].len <= up )
			{
				subpet.push_back( pet_ve[i] );
			}
		}
	}
	cout<<"size of subpet: "<<subpet.size()<<endl;
}

void PET_bank::getoverlappedPET_region( vector<size_t > &pet_ovl, map<string, set<pair<int, int > > > &regions )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_region;
	for ( map<string, set<pair<int, int > > >::iterator ite = regions.begin(); 
		ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int s = si->first;
			int e = si->second;
			int indexs = s / L;
			int indexe = e / L;
			chr_index_region[chr][indexs].insert( *si );
			if ( indexe > indexs )
			{
				for ( int k = indexs+1; k <= indexe; ++k )
					chr_index_region[chr][k].insert( *si );
			}
		}
	}
	
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		bool end1_ovl = false;
		string chr1 = pet_ve[i].chr1;
		int c1 = pet_ve[i].start1+(pet_ve[i].end1-pet_ve[i].start1)/2;
		int index1 = c1 / L;
	//	cout<<"PET"<<i<<endl;
	//	cout<<chr1<<" "<<pet_ve[i].start1<<" "<<pet_ve[i].end1<<" "<<c1<<" "<<index1<<endl;
		if ( chr_index_region[chr1].find( index1 ) != chr_index_region[chr1].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_region[chr1][index1].begin();
				ci != chr_index_region[chr1][index1].end(); ++ci )
			{
				if ( ci->second < c1 )
					continue;
				if ( ci->first > c1 )
					break;
			//	cout<<"Ovl at "<<chr1<<" "<<ci->first<<" "<<ci->second<<endl;
				end1_ovl = true;
				break;
			}
			
		} 
		
		if ( end1_ovl )
		{
			pet_ovl.push_back( i );
			continue;
		}
		
		bool end2_ovl = false;
		string chr2 = pet_ve[i].chr2;
		int c2 = pet_ve[i].start2+(pet_ve[i].end2-pet_ve[i].start2)/2;
		int index2 = c2 / L;
	//	cout<<chr2<<" "<<pet_ve[i].start2<<" "<<pet_ve[i].end2<<" "<<c2<<" "<<index2<<endl;
		if ( chr_index_region[chr2].find( index2 ) != chr_index_region[chr2].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_region[chr2][index2].begin();
				ci != chr_index_region[chr2][index2].end(); ++ci )
			{
				if ( ci->second < c2 )
					continue;
				if ( ci->first > c2 )
					break;
			//	cout<<"Ovl at "<<chr2<<" "<<ci->first<<" "<<ci->second<<endl;
				end2_ovl = true;
				break;
			}
			
		} 
		if ( end2_ovl )
		{
			pet_ovl.push_back(i );
		}
	//	exit(1);
	}
}

void PET_bank::getoverlappedPET_region2( vector<size_t > &pet_ovl, map<string, set<pair<int, int > > > &regions, int down, int up )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_region;
	for ( map<string, set<pair<int, int > > >::iterator ite = regions.begin(); 
		ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int s = si->first;
			int e = si->second;
			int indexs = s / L;
			int indexe = e / L;
			chr_index_region[chr][indexs].insert( *si );
			if ( indexe > indexs )
			{
				for ( int k = indexs+1; k <= indexe; ++k )
					chr_index_region[chr][k].insert( *si );
			}
		}
	}
	
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].chr1 != pet_ve[i].chr2 )
			continue;
		if ( pet_ve[i].len < down || pet_ve[i].len > up )
			continue;
			
		bool end1_ovl = false;
		string chr1 = pet_ve[i].chr1;
		int c1 = pet_ve[i].start1+(pet_ve[i].end1-pet_ve[i].start1)/2;
		int index1 = c1 / L;
		
		if ( chr_index_region[chr1].find( index1 ) != chr_index_region[chr1].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_region[chr1][index1].begin();
				ci != chr_index_region[chr1][index1].end(); ++ci )
			{
				if ( ci->second < c1 )
					continue;
				if ( ci->first > c1 )
					break;
				end1_ovl = true;
				break;
			}
			
		} 
		
		if ( end1_ovl )
		{
			pet_ovl.push_back( i );
			continue;
		}
		
		bool end2_ovl = false;
		string chr2 = pet_ve[i].chr2;
		int c2 = pet_ve[i].start2+(pet_ve[i].end2-pet_ve[i].start2)/2;
		int index2 = c2 / L;
		if ( chr_index_region[chr2].find( index2 ) != chr_index_region[chr2].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_region[chr2][index2].begin();
				ci != chr_index_region[chr2][index2].end(); ++ci )
			{
				if ( ci->second < c2 )
					continue;
				if ( ci->first > c2 )
					break;
				end2_ovl = true;
				break;
			}
			
		} 
		if ( end2_ovl )
		{
			pet_ovl.push_back(i );
		}
	}
}

void PET_bank::outputPET( string outfile, vector<size_t > & pet )
{
	ofstream outf( outfile.data() );
	for ( size_t i = 0; i < pet.size(); ++i )
	{
		size_t id = pet[i];
		outf<<pet_ve[id].chr1<<"\t"<<pet_ve[id].start1<<"\t"<<pet_ve[id].end1;
		outf<<"\t"<<pet_ve[id].chr2<<"\t"<<pet_ve[id].start2<<"\t"<<pet_ve[id].end2;
		outf<<"\t"<<pet_ve[id].q<<"\t"<<pet_ve[id].strand1<<"\t"<<pet_ve[id].strand2<<endl;
	}
}

void PET_bank::outputPET( string outfile )
{
	ofstream outf( outfile.data() );
	for ( size_t id = 0; id < pet_ve.size(); ++id )
	{
		
		outf<<pet_ve[id].chr1<<"\t"<<pet_ve[id].start1<<"\t"<<pet_ve[id].end1;
		outf<<"\t"<<pet_ve[id].chr2<<"\t"<<pet_ve[id].start2<<"\t"<<pet_ve[id].end2;
		outf<<"\t"<<pet_ve[id].q<<"\t"<<pet_ve[id].strand1<<"\t"<<pet_ve[id].strand2<<endl;
	}
}

void PET_bank::generateindex()
{
	int L = 100000;
	for ( map<string, map<pair<int, int >, size_t > >::iterator ite = chr_frg_pet.begin();
		ite != chr_frg_pet.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, size_t >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int c = si->first.first+ (si->first.second-si->first.first)/2;
			int index = c / L;
			chr_index_pos[chr][index].insert( si->first );	
		}
		
	}
}

void Loop_bank::readinpeak( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	string line;
	
	
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> parsed_items = parse_string( line );
		
		string chr = parsed_items[0];
		int start = atoi(parsed_items[1].c_str());
		int end = atoi(parsed_items[2].c_str());
		peaks[chr].insert( make_pair(start, end ) );
	}
	inf.close();
}

void Loop_bank::generateindex()
{
	int L = 100000;
	for ( map<string, set<pair<int, int > > >::iterator ite = peaks.begin(); ite != peaks.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int start = si->first;
			int end = si->second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_peaks[chr][index1].insert( *si );
			if ( index2 > index1 )
			{
				chr_index_peaks[chr][index2].insert( *si );
			
			}
		}
	}
}


void Loop_bank::getloopedpeak( int minP )
{
	
	int L = 100000;
	
	map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > looped_peak_all;
	
	for ( size_t i = 0; i < petbank.pet_ve.size(); ++i) 
	{
		string chr1 = petbank.pet_ve[i].chr1;
		int start1 = petbank.pet_ve[i].start1;
		int end1 = petbank.pet_ve[i].end1;
		int c1 = start1 + (end1-start1)/2;
		string chr2 = petbank.pet_ve[i].chr2;
		int start2 = petbank.pet_ve[i].start2;
		int end2 = petbank.pet_ve[i].end2;
		int c2 = start2 + (end2-start2)/2;
		if ( chr1 != chr2 )
			continue;
			
		int index1 = c1 / L;
		bool fd1 = false;
		int peak1_start = 0;
		int peak1_end = 0;
		if ( chr_index_peaks[chr1].find( index1 ) != chr_index_peaks[chr1].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peaks[chr1][index1].begin();
				ci != chr_index_peaks[chr1][index1].end(); ++ci )
			{
				if ( ci->second+200 < c1 )
					continue;
				if ( ci->first-200 > c1 )
					break;
				fd1 = true;
				peak1_start = ci->first;
				peak1_end = ci->second;
				break;
			}
		}
		if ( !fd1 )
			continue;
			
		int index2 = c2/ L;
		bool fd2 = false;
		int peak2_start = 0;
		int peak2_end = 0;
		if ( chr_index_peaks[chr2].find( index2 ) != chr_index_peaks[chr2].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peaks[chr2][index2].begin();
				ci != chr_index_peaks[chr2][index2].end(); ++ci )
			{
				if ( ci->second+200 < c2 )
					continue;
				if ( ci->first-200 > c2 )
					break;
				fd2 = true;
				peak2_start = ci->first;
				peak2_end = ci->second;
				break;
			}
		}
		
		if ( fd2 )
		{
			if ( peak1_start != peak2_start && peak1_end != peak2_end )
			{
				
				looped_peak_all[chr1][make_pair(peak1_start, peak1_end) ][make_pair(peak2_start, peak2_end)].push_back( i );
				looped_peak_all[chr1][make_pair(peak2_start, peak2_end) ][make_pair(peak1_start, peak1_end)].push_back( i );
				 
			}
		}
		
	}
		
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak_all.begin();
		ite != looped_peak_all.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			for ( map<pair<int, int >, vector<size_t > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				if ( (int)ti->second.size() >= minP )
				{
					looped_peak[chr][si->first][ti->first] = ti->second;
				}
			}
		}
	}
}

void Loop_bank::getcluster_loop_by_pet( int win, int gap )  // gap unit: times 0, 1, 2, ... of win
{
	map<string, map<int, map<int, vector<size_t > > > > chr_bin_bin_pet;
	for ( size_t i = 0; i < petbank.pet_ve.size(); ++i )
	{
		if ( petbank.pet_ve[i].chr1 != petbank.pet_ve[i].chr2 )
			continue;
		string chr = petbank.pet_ve[i].chr1;
		int b1 = petbank.pet_ve[i].start1 / win;
		int b2 = petbank.pet_ve[i].end2 / win;
		chr_bin_bin_pet[chr][b1][b2].push_back( i );
		chr_bin_bin_pet[chr][b2][b1].push_back( i );
	}
	
	map<string, set<pair<int, int > > > clustered_peak;
	map<string, map<int, pair<int, int > > > chr_bin_peak;
	for ( map<string, map<int, map<int, vector<size_t > > > >::iterator ite = chr_bin_bin_pet.begin();
		ite != chr_bin_bin_pet.end(); ++ite )
	{
		string chr = ite->first;
		vector<vector<int > > clustered_bin;
		vector<int > running_bin;
		for ( map<int, map<int, vector<size_t > > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int b = si->first;
			if ( running_bin.empty() )
				running_bin.push_back( b );
			else
			{
			/*	if ( b > running_bin.back()+gap+1 )
				{
					clustered_bin.push_back( running_bin );
					running_bin.clear();
					running_bin.push_back( b );
				} else
				{
					running_bin.push_back( b );
				}  */
				
				// not cluster
				clustered_bin.push_back( running_bin );
				running_bin.clear();
				running_bin.push_back( b );
			}
		}
		if ( !running_bin.empty() )
		{
			clustered_bin.push_back( running_bin );
		}
		
		for ( size_t i = 0; i < clustered_bin.size(); ++i )
		{
			int fb = clustered_bin[i].front();
			int bb = clustered_bin[i].back();
			int start = fb * win;
			int end = (bb+1)*win - 1;
			clustered_peak[chr].insert(  make_pair(start, end ) );
			for ( size_t j = 0; j < clustered_bin[i].size(); ++j )
			{
				chr_bin_peak[chr][clustered_bin[i][j]] = make_pair( start, end ); 
			}
		}
	}
	
	for ( map<string, map<int, map<int, vector<size_t > > > >::iterator ite = chr_bin_bin_pet.begin();
		ite != chr_bin_bin_pet.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, map<int, vector<size_t > > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int b1 = si->first;
			if ( chr_bin_peak[chr].find( b1 ) == chr_bin_peak[chr].end() )
			{
				cout<<"error chr_bin_peak cannot find bin1 "<<chr<<" "<<b1<<endl;
				exit(1);
			}
			pair<int, int > peak1 = chr_bin_peak[chr][b1];
			for ( map<int, vector<size_t > >::iterator ti = si->second.begin();
				ti != si->second.end(); ++ti )
			{
				int b2 = ti->first;
				if ( chr_bin_peak[chr].find( b2 ) == chr_bin_peak[chr].end() )
				{
					cout<<"error chr_bin_peak cannot find bin2 "<<chr<<" "<<b2<<endl;
					exit(1);
				}
				pair<int, int > peak2 = chr_bin_peak[chr][b2];
				
				// if peak1 == peak2 it is self looping. it is possible when pet length is short and clustered peak is long.
				for ( size_t k = 0; k < ti->second.size(); ++k )
				{
					looped_peak[chr][peak1][peak2].push_back( ti->second[k] );
				}
			}
		}
	}
	
}

void Loop_bank::getcluster_height_from_loopedpeak()
{
	int n = 0;
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak.begin();
		ite != looped_peak.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			pair<int, int > peak1 = si->first;
			int num = 0;
			for ( map<pair<int, int >, vector<size_t > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				num += (int)ti->second.size();
			}
			peak_height[chr][peak1] = num;
			n+= 1;
		}
	}
	cout<<"peak num "<<n<<endl;
}

void Loop_bank::cal_overlap( vector< PET > &pet_ve_o )
{
	int L = 100000;
	int olp = 0;
	int in = 0;
	for ( size_t i = 0; i < pet_ve_o.size(); ++i) 
	{
		string chr1 = pet_ve_o[i].chr1;
		int start1 = pet_ve_o[i].start1;
		int end1 = pet_ve_o[i].end1;
		int c1 = start1 + (end1-start1)/2;
		string chr2 = pet_ve_o[i].chr2;
		int start2 = pet_ve_o[i].start2;
		int end2 = pet_ve_o[i].end2;
		int c2 = start2 + (end2-start2)/2;
		if ( chr1 != chr2 )
			continue;
			
		int index1 = c1 / L;
		bool fd1 = false;
		int peak1_start = 0;
		int peak1_end = 0;
		if ( chr_index_peaks[chr1].find( index1 ) != chr_index_peaks[chr1].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peaks[chr1][index1].begin();
				ci != chr_index_peaks[chr1][index1].end(); ++ci )
			{
				if ( ci->second+200 < c1 )
					continue;
				if ( ci->first-200 > c1 )
					break;
				fd1 = true;
				peak1_start = ci->first;
				peak1_end = ci->second;
				break;
			}
		}
		
		int index2 = c2/ L;
		bool fd2 = false;
		int peak2_start = 0;
		int peak2_end = 0;
		if ( chr_index_peaks[chr2].find( index2 ) != chr_index_peaks[chr2].end() )
		{
			for ( set<pair<int, int > >::iterator ci = chr_index_peaks[chr2][index2].begin();
				ci != chr_index_peaks[chr2][index2].end(); ++ci )
			{
				if ( ci->second+200 < c2 )
					continue;
				if ( ci->first-200 > c2 )
					break;
				fd2 = true;
				peak2_start = ci->first;
				peak2_end = ci->second;
				break;
			}
		}
		
		if ( fd1 && fd2 )
			in += 1;
		
		if ( !fd1 )
			continue;
		
		if ( looped_peak[chr1].find( make_pair(peak1_start, peak1_end ) ) == looped_peak[chr1].end() )
			continue;
		
		if ( !fd2 )
			continue;
		
		if ( looped_peak[chr1][make_pair(peak1_start, peak1_end )].find(make_pair(peak2_start, peak2_end ) ) 
			!= looped_peak[chr1][make_pair(peak1_start, peak1_end )].end() )
		{
			olp += 1;
		}
	}
	
	double r = olp*1.0/(int)pet_ve_o.size();
	cout<<"total pets "<<pet_ve_o.size()<<"\tolp: "<<olp<<"\trate: "<<r<<endl;
	cout<<"in pets "<<in<<"\t"<<in*1.0/(int)pet_ve_o.size()<<endl;
}

void cal_overlap( map<string, map<pair<int, int >, int > > &peak1,
	map<string, map<pair<int, int >, int > > &peak2 )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_peak1;
	for ( map<string, map<pair<int, int >, int > >::iterator ite = peak1.begin();
		ite != peak1.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_peak1[chr][index1].insert( make_pair(start, end ) );
			if ( index2 > index1 )
			{
				chr_index_peak1[chr][index2].insert( make_pair(start, end ) );
			}
		}
		
	}
	
	map<string, map<int, set<pair<int, int > > > > chr_index_peak2;
	for ( map<string, map<pair<int, int >, int > >::iterator ite = peak2.begin();
		ite != peak2.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_peak2[chr][index1].insert( make_pair(start, end ) );
			if ( index2 > index1 )
			{
				chr_index_peak2[chr][index2].insert( make_pair(start, end ) );
			}
		}
	}
	
	int totalpeak1 = 0;
	int ovlpeak1 = 0;
	for ( map<string, map<pair<int, int >, int > >::iterator ite = peak1.begin();
		ite != peak1.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			int index = c / L;
			bool fd = false;
			if ( chr_index_peak2[chr].find( index ) != chr_index_peak2[chr].end() )
			{
				for ( set<pair<int, int > >::iterator ci = chr_index_peak2[chr][index].begin();
					ci != chr_index_peak2[chr][index].end(); ++ci )
				{
					if ( ci->second < start )
						continue;
					if ( ci->first > end )
						break;
					fd = true;
					break;
				}
				
			}
			
			if ( fd )
				ovlpeak1 += 1;
			totalpeak1 += 1;
		}
	}
	
	cout<<"total peak from lib1: "<<totalpeak1<<endl;
	cout<<"overlap peak from lib1: "<<ovlpeak1<<endl;
	
	int totalpeak2 = 0;
	int ovlpeak2 = 0;
	for ( map<string, map<pair<int, int >, int > >::iterator ite = peak2.begin();
		ite != peak2.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, int >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			int index = c / L;
			bool fd = false;
			if ( chr_index_peak1[chr].find( index ) != chr_index_peak1[chr].end() )
			{
				for ( set<pair<int, int > >::iterator ci = chr_index_peak1[chr][index].begin();
					ci != chr_index_peak1[chr][index].end(); ++ci )
				{
					if ( ci->second < start )
						continue;
					if ( ci->first > end )
						break;
					fd = true;
					break;
				}
				
			}
			
			if ( fd )
				ovlpeak2 += 1;
			totalpeak2 += 1;
		}
	}
	
	cout<<"total peak from lib2: "<<totalpeak2<<endl;
	cout<<"overlap peak from lib2: "<<totalpeak2<<endl;
}

void cal_overlap( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > &looped_peak1,
	map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > > &looped_peak2 )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_peak1;
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak1.begin();
		ite != looped_peak1.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_peak1[chr][index1].insert( make_pair(start, end ) );
			if ( index2 > index1 )
			{
				chr_index_peak1[chr][index2].insert( make_pair(start, end ) );
			}
		}
		
	}
	
	map<string, map<int, set<pair<int, int > > > > chr_index_peak2;
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak2.begin();
		ite != looped_peak2.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_peak2[chr][index1].insert( make_pair(start, end ) );
			if ( index2 > index1 )
			{
				chr_index_peak2[chr][index2].insert( make_pair(start, end ) );
			}
		}
	}
	
	int totalpet1 = 0;
	int ovlpet1 = 0;
	int totalloop1 = 0;
	int ovlloop1 = 0;
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak1.begin();
		ite != looped_peak1.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			int index = c / L;
			vector<pair<int, int > > ovl_peak2;
			if ( chr_index_peak2[chr].find( index ) != chr_index_peak2[chr].end() )
			{
				for ( set<pair<int, int > >::iterator ci = chr_index_peak2[chr][index].begin();
					ci != chr_index_peak2[chr][index].end(); ++ci )
				{
					if ( ci->second < start )
						continue;
					if ( ci->first > end )
						break;
					ovl_peak2.push_back( *ci );
				}
			} 
			
			for ( map< pair<int, int >, vector<size_t> >::iterator ti = si->second.begin();
				ti != si->second.end(); ++ti )
			{
				int pstart = ti->first.first;
				int pend = ti->first.second;
				bool ovl = false;
				for ( vector<pair<int, int > >::iterator ci = ovl_peak2.begin(); 
					ci != ovl_peak2.end(); ++ci )
				{
					for ( map< pair<int, int >, vector<size_t> >::iterator sci = looped_peak2[chr][*ci].begin();
						sci != looped_peak2[chr][*ci].end(); ++sci )
					{
						if ( sci->first.second < pstart )
							continue;
						if ( sci->first.first > pend )
							break;
						ovl = true;
						break;
					}
					if ( ovl )
						break;
				}
				
				
				if ( ovl )
				{
					ovlloop1 += 1;
					ovlpet1 += (int)ti->second.size();
				} 
				
				totalloop1 += 1;
				totalpet1 += (int)ti->second.size();
				
			}
		}
	}
	
	cout<<"total loop from lib1: "<<totalloop1<<endl;
	cout<<"overlap loop from lib1: "<<ovlloop1<<endl;
	cout<<"total pet from lib1: "<<totalpet1<<endl;
	cout<<"overlap pet from lib1: "<<ovlpet1<<endl;
	
	int totalpet2 = 0;
	int ovlpet2 = 0;
	int totalloop2 = 0;
	int ovlloop2 = 0;
	for ( map<string, map<pair<int, int >, map< pair<int, int >, vector<size_t> > > >::iterator ite = looped_peak2.begin();
		ite != looped_peak2.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, map< pair<int, int >, vector<size_t> > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			int index = c / L;
			vector<pair<int, int > > ovl_peak1;
			if ( chr_index_peak1[chr].find( index ) != chr_index_peak1[chr].end() )
			{
				for ( set<pair<int, int > >::iterator ci = chr_index_peak1[chr][index].begin();
					ci != chr_index_peak1[chr][index].end(); ++ci )
				{
					if ( ci->second < start )
						continue;
					if ( ci->first > end )
						break;
					ovl_peak1.push_back( *ci );
				}
			} 
			
			for ( map< pair<int, int >, vector<size_t> >::iterator ti = si->second.begin();
				ti != si->second.end(); ++ti )
			{
				int pstart = ti->first.first;
				int pend = ti->first.second;
				bool ovl = false;
				for ( vector<pair<int, int > >::iterator ci = ovl_peak1.begin(); 
					ci != ovl_peak1.end(); ++ci )
				{
					for ( map< pair<int, int >, vector<size_t> >::iterator sci = looped_peak1[chr][*ci].begin();
						sci != looped_peak1[chr][*ci].end(); ++sci )
					{
						if ( sci->first.second < pstart )
							continue;
						if ( sci->first.first > pend )
							break;
						ovl = true;
						break;
					}
					if ( ovl )
						break;
				}
				
				
				if ( ovl )
				{
					ovlloop2 += 1;
					ovlpet2 += (int)ti->second.size();
				} 
				
				totalloop2 += 1;
				totalpet2 += (int)ti->second.size();
				
			}
		}
	}
	cout<<"total loop from lib2: "<<totalloop2<<endl;
	cout<<"overlap loop from lib2: "<<ovlloop2<<endl;
	cout<<"total pet from lib2: "<<totalpet2<<endl;
	cout<<"overlap pet from lib2: "<<ovlpet2<<endl;
	
}





