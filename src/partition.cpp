#include "partition.h"

size_t Partition::get_bin_by_start( int start )
{
	int id = start / win;
	if ( id <0 || id >= par_ve.size() )
	{
		cout<<"error cannot locate Par_Bin by start "<<start<<" "<<id<<"\t"<<chr<<endl;
		exit(1);
	}
	if ( par_ve[id].start != start )
	{
		cout<<"error in get bin by start "<<endl;
	}
	size_t i = id;
	return i;
}

void Partition_Bank::ini_Partition_from_genome( Genome &genome, int win )
{
	
	
	for ( map<string, string >::iterator ite = genome.genomeseq.begin(); ite != genome.genomeseq.end(); ++ite )
	{
		string chr = ite->first;
		if ( chr.find("_") != std::string::npos || chr == "chrM")
			continue;
		int l = (int)ite->second.size();
		int step = l / win;
		cout<<chr<<" "<<l<<endl;
		Partition P;
		P.chr = chr;
		P.win = win;
		for ( int i = 0; i < step; ++i )
		{
			Par_Bin pbin;
			pbin.start = i*win;
			P.par_ve.push_back( pbin );
		}
		chr_Par.insert( make_pair(chr, P ) );
	}
	
	
	
}

void Partition_Bank::ini_Partition_from_chrlen( string infile, int win )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	map<string, int > chr_len;
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		int l = atoi(ps[1].c_str() );
		chr_len.insert( make_pair(chr, l ) );
	}
	inf.close();
	
	for ( map<string, int >::iterator ite = chr_len.begin(); ite != chr_len.end(); ++ite )
	{
		string chr = ite->first;
	/*	if ( chr.find("_") != std::string::npos || chr == "chrM")
			continue;   */
		int l = ite->second;
		int step = l / win;
		cout<<chr<<" "<<l<<endl;
		Partition P;
		P.chr = chr;
		P.win = win;
		
		for ( int i = 0; i < step; ++i )
		{
			Par_Bin pbin;
			pbin.start = i*win;
			pbin.end_count = 0;
			pbin.fea_ACC = 0;
			pbin.peak_anchor = false;
			P.par_ve.push_back( pbin );
			
		}
		chr_Par.insert( make_pair(chr, P ) );
	}
	
	
}





void Partition_Bank::add_Partition_freq_GC( Genome &genome, string inkmer )
{
	kmer = inkmer;
	map<string, vector<bool > > plus;
	map<string, vector<bool > > minus;
	
	cout<<"get kmer distri"<<endl;
	genome.get_kmer_distri( kmer, plus, minus );
	
	cout<<"add kmer freq"<<endl;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite  )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			int c = 0;
			int start = ite->second.par_ve[i].start;
			int win = ite->second.win;
			for ( int pos = start; pos < start+win; ++pos )
			{
				if ( pos >= (int)plus[chr].size() || pos >= (int)minus[chr].size() )
				{
					cout<<"error in Partition_Bank::add_Partition_freq_GC index exceed "<<endl;
					cout<<pos<<" "<<plus[chr].size()<<" "<<minus[chr].size()<<endl;
					exit(1);
				}
				if ( plus[chr][pos] )
					c += 1;
				if ( minus[chr][pos] )
					c += 1;
					
			}
			
			ite->second.par_ve[i].freq_GC = c *1.0 / (win*2);
		}
	}
}

void Partition_Bank::add_Partition_mapa( string infile )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	map<string, vector<double > > mapa;
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		cout<<"chr"<<endl;
		vector<double > ve;
		for ( size_t i = 4; i < ps.size(); ++i )
		{
			if ( ps[i] == "NA" )
			{
				ve.push_back(0);
			} else
			{
				double s = atof(ps[i].c_str() );
				if ( s < 0 || s > 1 )
				{
					cout<<"unexpected score in Partition_Bank::add_Partition_mapa"<<endl;
					cout<<"chr"<<" "<<i<<" "<<ps[i]<<" "<<s<<endl;
					exit(1);
				}
				ve.push_back(s);
			}
		}
		
		mapa[chr] = ve;
	}
	
	cout<<"assign mapa"<<endl;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite  )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			int c = 0;
			int start = ite->second.par_ve[i].start;
			int win = ite->second.win;
			double s_t = 0;
			for ( int pos = start; pos < start+win; ++pos )
			{
				if ( pos >= (int)mapa[chr].size() )
				{
					cout<<"error in Partition_Bank::add_Partition_mapa index exceed "<<endl;
					cout<<pos<<" "<<mapa[chr].size()<<endl;
					exit(1);
				}
				s_t += mapa[chr][pos];
				c += 1;
			}
			double s_a = s_t / c;
			ite->second.par_ve[i].mapa = s_a;
		}
	}
	
}

void Partition_Bank::write_partition( string outfile )
{
	ofstream outf( outfile.data() );
	
	outf<<"# Partition_Bank, kmer="<<kmer<<endl;
	outf<<"# pos_start freq_kmer mappability_score"<<endl;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		outf<<">chrome="<<chr<<"\twin="<<ite->second.win<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			outf<<ite->second.par_ve[i].start<<"\t"<<ite->second.par_ve[i].freq_GC<<"\t"<<ite->second.par_ve[i].mapa<<endl;
		}
	}
	outf.close();
	
}

void Partition_Bank::read_partition_with_genomic_fea( string infile )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline( inf, line );
	getline( inf, line );
	
	getline(inf,line);
	while(!inf.eof())
	{
		
		vector<string > ps = parse_string( line );
		if ( ps.size() != 2 || line[0] != '>' )
		{
			cout<<"error wrong chr line "<<line<<endl;
			exit(1);
		}
		vector<string > chr_s = parse_string(ps[0], '=' );
		string chr = chr_s[1];
		vector<string > win_s = parse_string( ps[1], '=' );
		int win = atoi(win_s[1].c_str() );
		vector< Par_Bin > parve;
		Partition p;
		p.chr = chr;
		p.win = win;
		
		getline( inf, line );
		
		while ( !inf.eof() )
		{
			
			if ( line.empty() )
				break;
			
			vector<string > pss = parse_string( line );
			if ( pss.size() != 3 )
			{
				cout<<"error body line "<<line<<endl;
				exit(1);
			}
			Par_Bin pbin;
			
			pbin.start = atoi(pss[0].c_str() );
			pbin.freq_GC = atof( pss[1].c_str() );
			pbin.mapa = atof(pss[2].c_str() );
			parve.push_back(pbin);
			getline( inf, line );
			if ( line[0] == '>' )
				break;
		}
		p.par_ve = parve;
		chr_Par[chr] = p;
		cout<<chr<<"\t"<<p.win<<endl;
		if ( line.empty() )
			break;
	}
	inf.close();
}

void Partition_Bank::add_Partition_fea_ACC_from_PET_bank( PET_bank & pet_bank )  // all pets included
{
	// ini fea_ACC
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			ite->second.par_ve[i].fea_ACC = 0;
			ite->second.par_ve[i].end_count = 0;
		}
	}
/*	cout<<pet_bank.pet_ve.size()<<endl;
	cout<<"test "<<chr_Par.size()<<endl;
	for ( map<string, Partition>::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		cout<<ite->first<<"\t"<<ite->second.win<<endl;
	} */
	
	for ( size_t i = 0; i < pet_bank.pet_ve.size(); ++i )
	{
	/*	if ( i % 100000 == 0 )
		{
			cout<<"i"<<i<<endl;
		} */
		if ( pet_bank.pet_ve[i].chr1 == "chrM" || pet_bank.pet_ve[i].chr2 == "chrM" )
			continue;
			
		int pos1 = pet_bank.pet_ve[i].start1;
		if ( pet_bank.pet_ve[i].strand1 == '-' )	
			pos1 = pet_bank.pet_ve[i].end1;
		string chr1 = pet_bank.pet_ve[i].chr1;
		if ( chr_Par[chr1].win == 0 )
		{
			cout<<"win == 0 "<<chr1<<endl; exit(1);
		}
		int p_start = (pos1 / chr_Par[chr1].win) * chr_Par[chr1].win;
		int id1 = p_start / chr_Par[chr1].win;
		if ( id1 < chr_Par[chr1].par_ve.size() )
		{
			size_t pid = chr_Par[chr1].get_bin_by_start(p_start);
		
			chr_Par[chr1].par_ve[pid].fea_ACC += 1;
			chr_Par[chr1].par_ve[pid].end_count += 1;
		}
		
		int pos2 = pet_bank.pet_ve[i].start2;
		if ( pet_bank.pet_ve[i].strand2 == '-' )	
			pos2 = pet_bank.pet_ve[i].end2;
		string chr2 = pet_bank.pet_ve[i].chr2;
		if ( chr_Par[chr2].win == 0 )
		{
			cout<<"win == 0 "<<chr2<<endl; exit(1);
		}
		int p_start2 = (pos2 / chr_Par[chr2].win ) * chr_Par[chr2].win;
		int id2 = p_start2 / chr_Par[chr2].win;
		if ( id2 < chr_Par[chr2].par_ve.size() )
		{
			size_t pid2 = chr_Par[chr2].get_bin_by_start(p_start2 );
		
			chr_Par[chr2].par_ve[pid2].fea_ACC += 1;
			chr_Par[chr2].par_ve[pid2].end_count += 1;
		}
	}
	
//	cout<<"b"<<endl;
	// normalize
	double fac = 1000000.0 / ( (int)pet_bank.pet_ve.size() * 2 );
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			ite->second.par_ve[i].fea_ACC *= fac;
		}
	}
	
}

void Partition_Bank::write_partition_ext2( string outfile )
{
	ofstream outf( outfile.data() );
	
	outf<<"# Partition_Bank, kmer="<<kmer<<endl;
	outf<<"# pos_start freq_kmer mappability_score fea_ACC tag_end_count"<<endl;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		outf<<">chrome="<<chr<<"\twin="<<ite->second.win<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			outf<<ite->second.par_ve[i].start<<"\t"<<ite->second.par_ve[i].freq_GC;
			outf<<"\t"<<ite->second.par_ve[i].mapa<<"\t"<<ite->second.par_ve[i].fea_ACC;
			outf<<"\t"<<ite->second.par_ve[i].end_count<<endl;
		}
	}
	outf.close();
	
}

void Partition_Bank::read_partition_from_par2( string infile )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline( inf, line );
	getline( inf, line );
	
	getline(inf,line);
	while(!inf.eof())
	{
		
		vector<string > ps = parse_string( line );
		if ( ps.size() != 2 || line[0] != '>' )
		{
			cout<<"error wrong chr line "<<line<<endl;
			exit(1);
		}
		vector<string > chr_s = parse_string(ps[0], '=' );
		string chr = chr_s[1];
		vector<string > win_s = parse_string( ps[1], '=' );
		int win = atoi(win_s[1].c_str() );
		vector< Par_Bin > parve;
		Partition p;
		p.chr = chr;
		p.win = win;
		
		getline( inf, line );
		
		while ( !inf.eof() )
		{
			
			if ( line.empty() )
				break;
			
			vector<string > pss = parse_string( line );
			if ( pss.size() != 5 )
			{
				cout<<"error body line "<<line<<endl;
				exit(1);
			}
			Par_Bin pbin;
			
			pbin.start = atoi(pss[0].c_str() );
			pbin.freq_GC = atof( pss[1].c_str() );
			pbin.mapa = atof(pss[2].c_str() );
			pbin.fea_ACC = atof(pss[3].c_str() );
			pbin.end_count = atoi(pss[4].c_str() );
			pbin.loop = false;
			parve.push_back(pbin);
			getline( inf, line );
			if ( line[0] == '>' )
				break;
		}
		p.par_ve = parve;
		chr_Par[chr] = p;
		cout<<chr<<"\t"<<p.win<<endl;
		if ( line.empty() )
			break;
	}
	inf.close();
}

void Partition_Bank::add_partition_anchor_from_region( string infile )
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
		if ( (int)ps.size() < 3 )
		{
			cout<<"error wrong line"<<line<<endl;
			exit(1);
		}
		string chr = ps[0];
		int start = atoi( ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
	/*	if ( chr == "chrM" )
			continue;
		if ( chr_Par.find( chr ) == chr_Par.end() )
		{
			cout<<"Unknown chr "<<chr<<endl;
			continue;
		}  */
		int p_start1 = (start / chr_Par[chr].win ) * chr_Par[chr].win;
		int p_start2 = (end / chr_Par[chr].win ) * chr_Par[chr].win;
		int id1 = p_start1 / chr_Par[chr].win;
		int id2 = p_start2 / chr_Par[chr].win;
		if ( id1 >= chr_Par[chr].par_ve.size() || id2 >= chr_Par[chr].par_ve.size() )
		{
			continue;
		}
		size_t pid1 = chr_Par[chr].get_bin_by_start(p_start1 );
		size_t pid2 = chr_Par[chr].get_bin_by_start(p_start2 );
		
		chr_Par[chr].par_ve[pid1].peak_anchor = true;
		chr_Par[chr].par_ve[pid2].peak_anchor = true;
	}
	inf.close();
}

void Partition_Bank::assign_pet_to_par_pair( PET_bank & pet_bank, 
	map<string, vector<pair<size_t, size_t > > > &par_pair, int minl, int maxl, int &ck, bool anchor )
{
	cout<<"total ini pets "<<pet_bank.pet_ve.size()<<endl;
	ck = 0;
	
	ofstream outf("tmp_pets.txt");
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
			
		int pos1 = pet_bank.pet_ve[i].start1;
		if ( pet_bank.pet_ve[i].strand1 == '-' )	
			pos1 = pet_bank.pet_ve[i].end1;
		string chr1 = pet_bank.pet_ve[i].chr1;
		if ( chr_Par[chr1].win == 0 )
		{
			cout<<"win == 0 "<<chr1<<endl; exit(1);
		}
		int p_start = (pos1 / chr_Par[chr1].win) * chr_Par[chr1].win;
		int id1 = p_start / chr_Par[chr1].win;
		
		int pos2 = pet_bank.pet_ve[i].start2;
		if ( pet_bank.pet_ve[i].strand2 == '-' )	
			pos2 = pet_bank.pet_ve[i].end2;
		string chr2 = pet_bank.pet_ve[i].chr2;
		if ( chr_Par[chr2].win == 0 )
		{
			cout<<"win == 0 "<<chr2<<endl; exit(1);
		}
		
		int p_start2 = (pos2 / chr_Par[chr2].win ) * chr_Par[chr2].win;
		int id2 = p_start2 / chr_Par[chr2].win;
		
			
		if ( id1 >= chr_Par[chr1].par_ve.size() || id2 >= chr_Par[chr2].par_ve.size() )
		{
			continue;
		}
		
		
		size_t pid1 = chr_Par[chr1].get_bin_by_start(p_start );
		size_t pid2 = chr_Par[chr2].get_bin_by_start(p_start2 );
		
		bool cont = true;
		if ( anchor )
		{
			if ( !chr_Par[chr1].par_ve[pid1].peak_anchor 
				|| !chr_Par[chr2].par_ve[pid2].peak_anchor )
				cont = false;
		}
		
		if ( cont )
		{
			chr_Par[chr1].par_ve[pid1].loop = true;
			chr_Par[chr2].par_ve[pid2].loop = true;
		
			par_pair[chr1].push_back( make_pair( pid1, pid2 ) );
		/*	if ( (chr1 == "chr10" && pid1 == 6056 && pid2 == 6092) || 
				(chr1 == "chr10" && pid2 == 6056 && pid1 == 6092) )
			{
				outf<<pet_bank.pet_ve[i].chr1<<"\t"<<pet_bank.pet_ve[i].start1
					<<"\t"<<pet_bank.pet_ve[i].end1<<"\t"<<pet_bank.pet_ve[i].chr2
					<<"\t"<<pet_bank.pet_ve[i].start2<<"\t"<<pet_bank.pet_ve[i].end2
					<<"\t"<<42<<"\t"<<pet_bank.pet_ve[i].strand1<<"\t"<<pet_bank.pet_ve[i].strand2<<endl;
				cout<<pet_bank.pet_ve[i].chr1<<"\t"<<pet_bank.pet_ve[i].start1
					<<"\t"<<pet_bank.pet_ve[i].end1<<"\t"<<pet_bank.pet_ve[i].chr2
					<<"\t"<<pet_bank.pet_ve[i].start2<<"\t"<<pet_bank.pet_ve[i].end2
					<<"\t"<<42<<"\t"<<pet_bank.pet_ve[i].strand1<<"\t"<<pet_bank.pet_ve[i].strand2<<endl;
			}  */
			ck += 1;
		}
	}
	
	cout<<ck<<endl;
	
	
	
}

void Partition_Bank::assign_pet_to_par_pair( PET_bank & pet_bank, 
	map<string, map<pair<size_t, size_t >, int > > &par_pair_c, string prefix )
{
	map<string, vector<pair<size_t, size_t > > > par_pair;
	int ck = 0; 
	assign_pet_to_par_pair( pet_bank, par_pair, 2000, 2000000, ck, false );
	
	int n = 0;
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = par_pair.begin();
		ite != par_pair.end(); ++ite )
	{
		string chr = ite->first;
		n+=(int)ite->second.size();
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( par_pair_c[chr].find( ite->second[i] ) != par_pair_c[chr].end() )
			{
				par_pair_c[chr][ite->second[i]] += 1;
			} else
				par_pair_c[chr][ite->second[i]] = 1;
		}
	}
	cout<<"assigned "<<n<<" pets"<<endl;
	
	string outfile = prefix+".raw_pets.longrange.txt";
	double fac = 1000000.0 / n;
	int ID = 0;
	ofstream outf( outfile.data() );
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); 
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			double nc = c * fac;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 2000 || pid1 > pid2 + 2000 )
				continue;
				
			ID += 1;	
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[pid1].start
				<<"\t"<<chr_Par[chr].par_ve[pid1].start+chr_Par[chr].win-1
				<<"\t"<<chr<<":"<<chr_Par[chr].par_ve[pid2].start
				<<"-"<<chr_Par[chr].par_ve[pid2].start+chr_Par[chr].win-1
				<<","<<nc<<"\t"<<ID<<"\t."<<endl;
			
			ID += 1;
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[pid2].start
				<<"\t"<<chr_Par[chr].par_ve[pid2].start+chr_Par[chr].win-1
				<<"\t"<<chr<<":"<<chr_Par[chr].par_ve[pid1].start
				<<"-"<<chr_Par[chr].par_ve[pid1].start+chr_Par[chr].win-1
				<<","<<nc<<"\t"<<ID<<"\t."<<endl;
			
		}
	}
		
}

void Partition_Bank::assign_pet_to_par_pair2( PET_bank & pet_bank, 
	int minl, int maxl, string prefix )
{
	map<string, vector<pair<size_t, size_t > > > par_pair;
	int total_pets = 0;
	assign_pet_to_par_pair( pet_bank, par_pair, minl, maxl, total_pets, false );
	
	string outfile0= prefix+".total_valid_PETs.txt";
	ofstream outf0(outfile0.data() );
	outf0<<total_pets<<endl;
	outf0.close();
	
	map<string, map<pair<size_t, size_t >, int > > par_pair_c;
	
	int total_valid_pets = 0;
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = par_pair.begin();
		ite != par_pair.end(); ++ite )
	{
		string chr = ite->first;
		total_valid_pets+=(int)ite->second.size();
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( par_pair_c[chr].find( ite->second[i] ) != par_pair_c[chr].end() )
			{
				par_pair_c[chr][ite->second[i]] += 1;
			} else
				par_pair_c[chr][ite->second[i]] = 1;
			chr_Par[chr].par_ve[ite->second[i].first].end_count += 1;
			chr_Par[chr].par_ve[ite->second[i].second].end_count += 1;
		}  
	}
	
	cout<<"assigned "<<total_valid_pets<<" pets"<<endl;
	
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			double rpm = ite->second.par_ve[i].end_count * 1000000.0 / total_valid_pets;
			ite->second.par_ve[i].fea_ACC = rpm;
		}
	}
	
	string outfile = prefix+".raw_pets.longrange.txt";
	double fac = 1000000.0 / total_valid_pets;
	int ID = 0;
	ofstream outf( outfile.data() );
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); 
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			double nc = c * fac;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 2000 || pid1 > pid2 + 2000 )
				continue;
				
			ID += 1;	
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[pid1].start
				<<"\t"<<chr_Par[chr].par_ve[pid1].start+chr_Par[chr].win-1
				<<"\t"<<chr<<":"<<chr_Par[chr].par_ve[pid2].start
				<<"-"<<chr_Par[chr].par_ve[pid2].start+chr_Par[chr].win-1
				<<","<<nc<<"\t"<<ID<<"\t."<<endl;
			
			ID += 1;
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[pid2].start
				<<"\t"<<chr_Par[chr].par_ve[pid2].start+chr_Par[chr].win-1
				<<"\t"<<chr<<":"<<chr_Par[chr].par_ve[pid1].start
				<<"-"<<chr_Par[chr].par_ve[pid1].start+chr_Par[chr].win-1
				<<","<<nc<<"\t"<<ID<<"\t."<<endl;
			
		}
	}
		
}

void Partition_Bank::assign_pet_to_par_pair3( PET_bank & pet_bank, 
	map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
	int minl, int maxl )
{
	map<string, vector<pair<size_t, size_t > > > par_pair;
	int valid_PETs= 0;
	assign_pet_to_par_pair( pet_bank, par_pair, minl, maxl, valid_PETs, true );
	
	total_valid_pets_in_anchor = 0;
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = par_pair.begin();
		ite != par_pair.end(); ++ite )
	{
		string chr = ite->first;
		total_valid_pets_in_anchor+=(int)ite->second.size();
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( par_pair_c[chr].find( ite->second[i] ) != par_pair_c[chr].end() )
			{
				par_pair_c[chr][ite->second[i]] += 1;
			} else
				par_pair_c[chr][ite->second[i]] = 1;
		
		}  
	}
	
	cout<<"assigned "<<total_valid_pets_in_anchor<<" pets in anchor"<<endl;
	

		
}

void Partition_Bank::statistics_fea_ACC()
{
	multiset< double > accset;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			if ( ite->second.par_ve[i].fea_ACC > 0 )
				accset.insert( ite->second.par_ve[i].fea_ACC );
		}
	}
	
	cout<<"statistics fea_ACC "<<endl;
	
	double total = 0;
	
	for ( multiset< double >::iterator ite = accset.begin(); ite != accset.end(); ++ite )
	{
		total += *ite;
	}
	
	double ave = total / (int)accset.size();
	cout<<"ave="<<ave<<endl;
	cout<<"max="<<*accset.rbegin()<<endl;
	
	ofstream outf("ACC_density.dist.txt");
	map<int, int > acc_d_g_c;
	for ( int i = 0; i < 10000; ++i )
	{
		acc_d_g_c.insert( make_pair(i, 0) );
	}
	for ( multiset< double >::iterator ite = accset.begin(); ite != accset.end(); ++ite )
	{
		int g = (int)( (*ite) * 100);
		if ( g < 10000 )
			acc_d_g_c[g] += 1;
	}
	for ( map<int, int >::iterator ite = acc_d_g_c.begin(); ite != acc_d_g_c.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
	
}

void Partition_Bank::statistics_freq_GC()
{
	multiset< double > freqset;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			if ( ite->second.par_ve[i].freq_GC > 0 )
				freqset.insert( ite->second.par_ve[i].freq_GC );
		}
	}
	
	cout<<"statistics freq_GC "<<endl;
	
/*	double total = 0;
	
	for ( multiset< double >::iterator ite = freqset.begin(); ite != freqset.end(); ++ite )
	{
		total += *ite;
	}
	
	double ave = total / (int)accset.size();
	cout<<"ave="<<ave<<endl;
	cout<<"max="<<*accset.rbegin()<<endl;  */
	
	ofstream outf("freq_kmer.dist.txt");
	map<int, int > freq_d_g_c;
	for ( int i = 0; i < 10; ++i )
	{
		freq_d_g_c.insert( make_pair(i, 0) );
	}
	for ( multiset< double >::iterator ite = freqset.begin(); ite != freqset.end(); ++ite )
	{
		int g = (int)( (*ite) * 1000);
		if ( g < 10 )
			freq_d_g_c[g] += 1;
	}
	for ( map<int, int >::iterator ite = freq_d_g_c.begin(); ite != freq_d_g_c.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
	
}

void Partition_Bank::statistics_mapa()
{
	multiset< double > mapaset;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			
			mapaset.insert( ite->second.par_ve[i].mapa );
		}
	}
	
	cout<<"statistics mapa "<<endl;
	

	
	ofstream outf("mappability.dist.txt");
	map<int, int > mapa_d_g_c;
	for ( int i = 0; i < 10; ++i )
	{
		mapa_d_g_c.insert( make_pair(i, 0) );
	}
	for ( multiset< double >::iterator ite = mapaset.begin(); ite != mapaset.end(); ++ite )
	{
		int g = (int)( (*ite) * 10);
		if ( g < 10 )
			mapa_d_g_c[g] += 1;
	}
	for ( map<int, int >::iterator ite = mapa_d_g_c.begin(); ite != mapa_d_g_c.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
	
}

void Partition_Bank::generate_simulate_space( map<int, vector< pair<int, int > > > &length_accm_map,
	map<int, vector<int > > &space, int mings )
{
	map< int, map<int, int > > length_accm_count;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
	//	if ( chr != "chr1" )
		//	continue;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			if ( i % 10000 == 0 )
				cout<<i<<endl;
			double acc_1 = ite->second.par_ve[i].fea_ACC;
			if ( acc_1 == 0 )
				continue;
			if ( !ite->second.par_ve[i].loop )
				continue;
				
			for ( size_t j = 2; j <= 2000; ++j )
			{
				size_t k = i+j;
				if ( k >= ite->second.par_ve.size() )
					break;
				double acc_2 = ite->second.par_ve[k].fea_ACC;
				if ( acc_2 == 0 )
					continue; 
				if ( !ite->second.par_ve[k].loop )
					continue;
				
				int length = j;
				int loglength = (int)(log(length*1.0));
				
				double acc_m = acc_1 * acc_2;
				int ind = (int)(acc_m * 1000);
				if ( length_accm_count[loglength].find(ind ) == length_accm_count[loglength].end() )
					length_accm_count[loglength][ind]=1;
				else
					length_accm_count[loglength][ind]+=1;
			}
		}
	}
	
	
	for ( map< int, map<int, int > >::iterator ite = length_accm_count.begin(); ite != length_accm_count.end(); ++ite )
	{
		int l = ite->first;
		map<int, int >::iterator si = ite->second.begin();
		
		int s = si->first;
		int c = si->second;
		int e = s;
		++si;
		for ( ; si != ite->second.end(); ++si )
		{
			if ( si->first > 1000 )
			{
				mings = 500;
				
			}
			
			if ( c >= mings )
			{
				length_accm_map[l].push_back( make_pair(s, e ) );
				space[l].push_back(c);
				
				s = si->first;
				c = si->second;
				e = s;
			} else
			{
				e = si->first;
				c += si->second;
			}
			
			
		}
		length_accm_map[l].push_back( make_pair( s, e) );
		space[l].push_back(c);
	}
	
	ofstream outf( "sm_space.txt" );
	for ( map<int, vector< pair<int, int > > >::iterator ite = length_accm_map.begin(); ite != length_accm_map.end(); ++ite )
	{
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			outf<<ite->first<<"\t"<<ite->second[i].first<<"\t"<<ite->second[i].second<<"\t"<<space[ite->first][i]<<endl;
		}
	}
	outf.close();
	
}

void Partition_Bank::generate_observed_space( map<int, vector< pair<int, int > > > &length_accm_map,
		map<int, vector<int > > &sm_space,
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		double ex_rate,
		map<int, vector<int > > &space )
{
//	map<int, map<int, set<size_t> > > length_accm_index;
	map<int, vector<multiset<int, greater<int > > > > space_det;
	for ( map<int, vector< pair<int, int> > >::iterator ite = length_accm_map.begin(); 
		ite != length_accm_map.end(); ++ite )
	{
		vector<multiset<int, greater<int > > > sp;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
		/*	int ind1 = (int)(ite->second[i].first / 1000);
			int ind2 = (int)(ite->second[i].second / 1000 );
			for ( int k = ind1; k <= ind2; ++k )
			{
				length_accm_index[ite->first][k].insert( i );
			}  */
			multiset<int, greater<int > > nu;
			sp.push_back( nu );
		}
		space_det[ite->first] = sp;
	}    
	
	int i = 0;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin();
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 2000 || pid1 > pid2 + 2000 )
				continue;
			if ( pid2 == pid1 + 1 || pid1 == pid2 + 1 )
				continue;
			i+=1;
			if ( i % 100000 == 0 )
				cout<<"i"<<i<<endl;
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			
			if ( acc_1 == 0  )
				continue;		
			
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			
			if ( acc_2 == 0  )
				continue;
			
			double acc_m = acc_1 * acc_2;
			int l = (pid2 - pid1);
			int loglength = (int)(log(l*1.0));
			int accm_id = (int)(acc_m*1000);
			bool fd = false;
			for ( size_t gid = 0; gid < length_accm_map[loglength].size(); ++gid )
			{
				if (  length_accm_map[loglength][gid].second < accm_id )
					continue;
				if ( length_accm_map[loglength][gid].first > accm_id )
					break;
				fd = true;
				space_det[loglength][gid].insert( c );
				break;
			}
			
			
			if ( !fd )
			{
				cout<<"error cannot find accm "<<accm_id<<" for length "<<loglength<<endl; exit(1);
			}
		}
	}
	
	vector<int > groups;
	ofstream outf( "space.txt" );
	outf<<"#loglength\tacc_m_s\tacc_m_e\tsimu_space\tobs_loop_space\tobs_loop_pets\tleft_simu_space\tobs_loop_pets_left"<<endl;
	for ( map<int, vector<multiset<int, greater<int > > > >::iterator ite = space_det.begin(); ite != space_det.end(); ++ite )
	{
		int l = ite->first;
		if ( ite->second.size() != length_accm_map[l].size() )
		{
			cout<<"error space_det size != length_acc_map size "<<ite->second.size() <<" "<< length_accm_map[l].size()<<" for length "<<l<<endl;
			exit(1);
		}
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			int loopsize = (int)ite->second[i].size();
			if ( loopsize == 0 )
			{
			//	cout<<"warning looping size == 0 for l="<<l<<" accmid_range "<<length_accm_map[l][i].first<<" "<<length_accm_map[l][i].second<<endl;
				outf<<l<<"\t"<<length_accm_map[l][i].first<<" "<<length_accm_map[l][i].second<<"\t"<<sm_space[l][i]<<"\t0\t0\t0\t0"<<endl;
				continue;
			}
			
			int ex = (int)(loopsize*ex_rate);
			int left_s = sm_space[l][i] - ex;
			int total_c = 0;
			int left_c = 0;
			int c = 0;
			for ( multiset<int, greater<int > >::iterator si = ite->second[i].begin(); si != ite->second[i].end(); ++si )
			{
				total_c += *si;
				if ( c < ex )
				{
					
					c += 1;
					continue;
				}
				left_c += *si;
			}
			outf<<l<<"\t"<<length_accm_map[l][i].first<<" "<<length_accm_map[l][i].second<<"\t"<<sm_space[l][i]<<"\t"<<loopsize<<"\t"<<total_c<<"\t"<<left_s<<"\t"<<left_c<<endl;
			sm_space[l][i] = left_s;
			space[l].push_back( left_c );
			groups.push_back( left_s );
		}
	}
	
	cout<<"group ve: size:"<<groups.size()<<" sum:"<<getsum( groups )<<" max:"<<getmax(groups )<<" medium:"<<getmedium(groups)<<endl;
	
	
}


void Partition_Bank::simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		map<int, vector< pair<int, int > > > &length_accm_map,
		map<int, vector<int > > &sm_space,
		map<int, vector<int > > &obsspace,
		int times,
		string outfile )
{
/*	map<int, map<int, set<size_t> > > length_accm_index;
	map<int, vector<multiset<int, greater<int > > > > space_det;
	for ( map<int, vector< pair<double, double > > >::iterator ite = length_accm_map.begin(); 
		ite != length_accm_map.end(); ++ite )
	{
		int l = ite->first;
		vector<multiset<int, greater<int > > > sp;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			int ind1 = (int)(ite->second[i].first * 10);
			int ind2 = (int)(ite->second[i].second * 10 );
			for ( int k = ind1; k <= ind2; ++k )
			{
				length_accm_index[l][k].insert( i );
			}
			multiset<int, greater<int > > nu;
			sp.push_back( nu );
		}
		space_det[l] = sp;
	} */
	
	int tc = 0;
	map<int, map<size_t, vector<pair<string, pair<size_t, size_t > > > > > length_accmid_par_pair;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin();
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 2000 || pid1 > pid2 + 2000 )
				continue;
			if ( pid2 == pid1 + 1 || pid1 == pid2 + 1 )
				continue;
			if ( c < 3 )
				continue; 
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			
			if ( acc_1 == 0  )
				continue;		
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			
			if ( acc_2 == 0 )
				continue;
			
			double acc_m = acc_1 * acc_2;
			int l = (pid2 - pid1);
			int loglength = (int)(log(l*1.0));
			int accm_id = (int)(acc_m*1000);
			bool fd = false;
			for ( size_t gid = 0; gid < length_accm_map[loglength].size(); ++gid )
			{
				if (  length_accm_map[loglength][gid].second < accm_id )
					continue;
				if ( length_accm_map[loglength][gid].first > accm_id )
					break;
				fd = true;
				length_accmid_par_pair[loglength][gid].push_back( make_pair(ite->first, si->first ) );
				break;
			}
			
			
			if ( !fd )
			{
				cout<<"error cannot find accm "<<accm_id<<" for length "<<loglength<<endl; exit(1);
			}
			
			
			tc += 1;
			if ( tc % 100000 == 0 )
				cout<<"tc "<<tc<<endl;
			
			
			if ( !fd )
			{
				cout<<"IN simulation_and_cal_pv: error cannot find accm "<<accm_id<<" for length "<<loglength<<endl; exit(1);
			}
		}
	}
	
	cout<<"simulation"<<endl;
//	estimate by binomial distribution
	int o = 0;
	int o_tag = 0;
	srand(time(NULL));
	map<string, map<pair<size_t, size_t >, pair<int, double > > > par_pair_pv;
	for ( map<int, map<size_t, vector<pair<string, pair<size_t, size_t > > > > >::iterator ite = length_accmid_par_pair.begin();
		ite != length_accmid_par_pair.end(); ++ite )
	{
		int l = ite->first;
		for ( map<size_t, vector<pair<string, pair<size_t, size_t > > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t acc_id = si->first;
			
			o += si->second.size();
			if ( o / 10000 > o_tag )
			{
				cout<<o<<endl;
				o_tag = o/ 10000;
			}
			
			int ve_size = sm_space[l][acc_id];
			int total_pet_n = obsspace[l][acc_id];
		//	vector< vector<int > > mat;
		//	sim_func( ve_size, total_pet_n, times, mat );
					
			
			
		/*	int rd = std::rand();
			int rdn = rd % (ve_size - 100);
			int cc = 0;
			multiset<int, greater<int> > total_col_set;
			for ( int j = rdn; j < ve_size; ++j )
			{
				for ( int i = 0; i < times; ++i )
				{
					if ( cc >= 100000 )
						break;
					total_col_set.insert( mat[i][j]);
					cc += 1;
					
				}
				if ( cc >= 100000 )
					break;
			}	*/
			for ( size_t k = 0; k < si->second.size(); ++k )
			{
				int c = par_pair_c[si->second[k].first][si->second[k].second];
				if ( c == 0 )
				{
					cout<<"error c == 0 "<<endl; exit(1);
				}
			/*	int r2 = 0;
				
				
				for ( multiset<int, greater<int> >::iterator ii = total_col_set.begin(); 
					ii != total_col_set.end(); ++ii )
				{

					if ( c > *ii )
					{
						break;
					}
					r2 += 1;
				}
				double p2 = r2 * 1.0 / (int)total_col_set.size();
				*/
				
				double bp = 1.0 / ve_size;
				double p2 = cdf( c-1, total_pet_n, bp );
			
				par_pair_pv[si->second[k].first][si->second[k].second] = make_pair( ve_size, p2 );
			}
			
		}
	}
	
	ofstream outf(outfile.data() );
	outf<<"## simulation times "<<times<<endl;
	outf<<"#chrom\tpar1_start\tpart1_end\tpar2_start\tpart2_end\tdistance\tpart1_PE_count\tpart2_PE_count\tpart1_PE_density\tpart2_PE_density\tcongener_space\tpet_count\tp_value\n";
	for ( map<string, map<pair<size_t, size_t >, pair<int, double > > >::iterator ite = par_pair_pv.begin();
		ite != par_pair_pv.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, pair<int, double > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start - chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].fea_ACC
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].fea_ACC
				<<"\t"<<si->second.first
				<<"\t"<<par_pair_c[chr][si->first]
				<<"\t"<<si->second.second<<endl;
			
		}
	}
	outf.close();
}

void Partition_Bank::generate_simulate_space(map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > &space )
{
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
	//	if ( chr != "chr1" )
		//	continue;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			if ( i % 10000 == 0 )
				cout<<i<<endl;
			double acc_1 = ite->second.par_ve[i].fea_ACC;
			double mapa_1 = ite->second.par_ve[i].mapa;
			double freq_kmer_1 = ite->second.par_ve[i].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;
			int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;
				
			for ( size_t j = 1; j <= 1000; ++j )
			{
				size_t k = i+j;
				if ( k >= ite->second.par_ve.size() )
					break;
				
				double acc_2 = ite->second.par_ve[k].fea_ACC;
				double mapa_2 = ite->second.par_ve[k].mapa;
				double freq_kmer_2 = ite->second.par_ve[k].freq_GC;
				if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
					continue;
				
				int mapa_id_2 = (int)(mapa_2 * 10);
				int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000 );
				if ( freq_kmer_id_2 >= 10 )
					freq_kmer_id_2 = 10;
					
				double acc_m = acc_1 * acc_2;
			/*	int acc_m_id = (int)( acc_m * 10 );	
				if ( acc_m_id >= 1000 )
					acc_m_id = 1000;  */
				
				int log_acc_m_10_id = transfer_acc_m( acc_m );
				
				pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
				pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
				int length = j;
				
				if ( space[j][log_acc_m_10_id][mapa_pair].find( freq_pair ) ==  space[j][log_acc_m_10_id][mapa_pair].end() )
					space[j][log_acc_m_10_id][mapa_pair][freq_pair] = 1;
				else
					space[j][log_acc_m_10_id][mapa_pair][freq_pair] += 1;
					
			}
		}
	}
	
/*	vector<int > min_ve;
	vector<int > sec_ve;
	vector<int > thr_ve;
	map<int, vector<int > > acc_m_ve;
	map<int, vector<int > > length_ve;
	for ( map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > >::iterator ite = space.begin(); ite != space.end(); ++ite )
	{
		int t_c = 0;
		for ( map<int, map<pair<int, int >, map<pair<int, int >, int > > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int s_c = 0;
			for ( map<pair<int, int >, map<pair<int, int >, int  > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<pair<int, int >, int >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{	
					int c = fi->second;
					min_ve.push_back(c);
					s_c += c;
				}
			}
			sec_ve.push_back( s_c );
			acc_m_ve[si->first].push_back( s_c );
			length_ve[ite->first].push_back( s_c );
			t_c += s_c;
			
		}
		
		thr_ve.push_back(t_c );
	}
	
	cout<<"minimal ve: size:"<<min_ve.size()<<" sum:"<<getsum( min_ve )<<" max:"<<getmax(min_ve )<<" medium:"<<getmedium(min_ve)<<endl;
	cout<<"second ve: size:"<<sec_ve.size()<<" sum:"<<getsum( sec_ve )<<" max:"<<getmax(sec_ve )<<" medium:"<<getmedium(sec_ve)<<endl;
	cout<<"third ve: size:"<<thr_ve.size()<<" sum:"<<getsum( thr_ve )<<" max:"<<getmax(thr_ve )<<" medium:"<<getmedium(thr_ve)<<endl;
	
	ofstream outf( "acc_m_space.density.txt" );
	for ( map<int, vector<int > >::iterator ite = acc_m_ve.begin(); ite != acc_m_ve.end(); ++ite )
	{
		double acc_m = pow( 10, ((ite->first*1.0 / 10)) );
		
		outf<<acc_m<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf.close();
	
	ofstream outf2("length_space.density.txt");
	for ( map<int, vector<int > >::iterator ite = length_ve.begin(); ite != length_ve.end(); ++ite )
	{
		outf2<<ite->first<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf2.close();
	*/
}

void Partition_Bank::generate_simulate_space(map<int, map<int, map<int, map<int, int > > > > &space )
{
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
	//	if ( chr != "chr1" )
		//	continue;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
			if ( i % 10000 == 0 )
				cout<<i<<endl;
			double acc_1 = ite->second.par_ve[i].fea_ACC;
			double mapa_1 = ite->second.par_ve[i].mapa;
			double freq_kmer_1 = ite->second.par_ve[i].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;
		/*	int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;  */
				
			for ( size_t j = 1; j <= 1000; ++j )
			{
				size_t k = i+j;
				if ( k >= ite->second.par_ve.size() )
					break;
				
				double acc_2 = ite->second.par_ve[k].fea_ACC;
				double mapa_2 = ite->second.par_ve[k].mapa;
				double freq_kmer_2 = ite->second.par_ve[k].freq_GC;
				if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
					continue;
				
			/*	int mapa_id_2 = (int)(mapa_2 * 10);
				int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000 );
				if ( freq_kmer_id_2 >= 10 )
					freq_kmer_id_2 = 10;  */
					
				double acc_m = acc_1 * acc_2;
			/*	int acc_m_id = (int)( acc_m * 10 );	
				if ( acc_m_id >= 1000 )
					acc_m_id = 1000;  */
				
				int log_acc_m_10_id = transfer_acc_m( acc_m );
				
				double mapa_m = mapa_1 * mapa_2;
				int mapa_m_id = (int)( mapa_m * 1000);
				if ( mapa_m_id > 10 )
					mapa_m_id = 10;
				
				double freq_m = freq_kmer_1 * freq_kmer_2;
				int freq_m_id = (int)( freq_m * 10000000 );
				if ( freq_m_id > 10 )
					freq_m_id = 10;
				
			//	pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
			//	pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
				int length = j;
				
				if ( space[j][log_acc_m_10_id][mapa_m_id].find( freq_m_id ) ==  space[j][log_acc_m_10_id][mapa_m_id].end() )
					space[j][log_acc_m_10_id][mapa_m_id][freq_m_id] = 1;
				else
					space[j][log_acc_m_10_id][mapa_m_id][freq_m_id] += 1;
					
			}
		}
	}
	
	vector<int > min_ve;
	vector<int > sec_ve;
	vector<int > thr_ve;
	map<int, vector<int > > acc_m_ve;
	map<int, vector<int > > length_ve;
	for ( map<int, map<int, map<int, map<int, int > > > >::iterator ite = space.begin(); ite != space.end(); ++ite )
	{
		int t_c = 0;
		for ( map<int, map<int, map<int, int > > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int s_c = 0;
			for ( map<int, map<int, int  > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<int, int >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{	
					int c = fi->second;
					min_ve.push_back(c);
					s_c += c;
				}
			}
			sec_ve.push_back( s_c );
			acc_m_ve[si->first].push_back( s_c );
			length_ve[ite->first].push_back( s_c );
			t_c += s_c;
			
		}
		
		thr_ve.push_back(t_c );
	}
	
	cout<<"minimal ve: size:"<<min_ve.size()<<" sum:"<<getsum( min_ve )<<" max:"<<getmax(min_ve )<<" medium:"<<getmedium(min_ve)<<endl;
	cout<<"second ve: size:"<<sec_ve.size()<<" sum:"<<getsum( sec_ve )<<" max:"<<getmax(sec_ve )<<" medium:"<<getmedium(sec_ve)<<endl;
	cout<<"third ve: size:"<<thr_ve.size()<<" sum:"<<getsum( thr_ve )<<" max:"<<getmax(thr_ve )<<" medium:"<<getmedium(thr_ve)<<endl;
	
/*	ofstream outf( "acc_m_space.density.txt" );
	for ( map<int, vector<int > >::iterator ite = acc_m_ve.begin(); ite != acc_m_ve.end(); ++ite )
	{
		double acc_m = pow( 10, ((ite->first*1.0 / 10)) );
		
		outf<<acc_m<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf.close();
	
	ofstream outf2("length_space.density.txt");
	for ( map<int, vector<int > >::iterator ite = length_ve.begin(); ite != length_ve.end(); ++ite )
	{
		outf2<<ite->first<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf2.close();
	*/
}


void Partition_Bank::generate_observed_space( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
	map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > &space )
{
	
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
	//	cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 1000 || pid1 > pid2 + 1000 )
				continue;
			
		//	if ( i % 100000 == 0 )
			//	cout<<"i"<<i<<endl;
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			double mapa_1 = chr_Par[chr].par_ve[pid1].mapa;
			double freq_kmer_1 = chr_Par[chr].par_ve[pid1].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;		
			int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			double mapa_2 = chr_Par[chr].par_ve[pid2].mapa;
			double freq_kmer_2 = chr_Par[chr].par_ve[pid2].freq_GC;
			if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
				continue;
			int mapa_id_2 = (int)(mapa_2 * 10);
			int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000);
			if ( freq_kmer_id_2 >= 10 )
				freq_kmer_id_2 = 10;
			
			double acc_m = acc_1 * acc_2;
			int log_acc_m_10_id = transfer_acc_m( acc_m );
			
			pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
			pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
			int length = pid2 - pid1;
			if ( length < 0 )
				length *= (-1);
			
			if ( space[length][log_acc_m_10_id][mapa_pair].find( freq_pair ) ==  space[length][log_acc_m_10_id][mapa_pair].end() )
				space[length][log_acc_m_10_id][mapa_pair][freq_pair] = c;
			else
				space[length][log_acc_m_10_id][mapa_pair][freq_pair] += c;
		}
	}
	
	vector<int > min_ve;
	vector<int > sec_ve;
	vector<int > thr_ve;
	map<int, vector<int > > acc_m_ve;
	map<int, vector<int > > length_ve;
	for ( map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > >::iterator ite = space.begin(); ite != space.end(); ++ite )
	{
		int t_c = 0;
		for ( map<int, map<pair<int, int >, map<pair<int, int >, int > > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int s_c = 0;
			for ( map<pair<int, int >, map<pair<int, int >, int  > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<pair<int, int >, int >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{	
					int c = fi->second;
					min_ve.push_back(c);
					s_c += c;
				}
			}
			sec_ve.push_back( s_c );
			acc_m_ve[si->first].push_back( s_c );
			length_ve[ite->first].push_back( s_c );
			t_c += s_c;
			
		}
		
		thr_ve.push_back(t_c );
	}
	
	cout<<"minimal ve: size:"<<min_ve.size()<<" sum:"<<getsum( min_ve )<<" max:"<<getmax(min_ve )<<" medium:"<<getmedium(min_ve)<<endl;
	cout<<"second ve: size:"<<sec_ve.size()<<" sum:"<<getsum( sec_ve )<<" max:"<<getmax(sec_ve )<<" medium:"<<getmedium(sec_ve)<<endl;
	cout<<"third ve: size:"<<thr_ve.size()<<" sum:"<<getsum( thr_ve )<<" max:"<<getmax(thr_ve )<<" medium:"<<getmedium(thr_ve)<<endl;
	
/*	ofstream outf( "acc_m_OBSspace.density.txt" );
	for ( map<int, vector<int > >::iterator ite = acc_m_ve.begin(); ite != acc_m_ve.end(); ++ite )
	{
		double acc_m = pow( 10, ((ite->first*1.0 / 10)) );
		
		outf<<acc_m<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf.close();
	
	ofstream outf2("length_OBSspace.density.txt");
	for ( map<int, vector<int > >::iterator ite = length_ve.begin(); ite != length_ve.end(); ++ite )
	{
		outf2<<ite->first<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf2.close();  */
}

void Partition_Bank::generate_observed_space( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
	map<int, map<int, map<int, map<int, int > > > > &space )
{
	
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
	//	cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int c = si->second;
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 1000 || pid1 > pid2 + 1000 )
				continue;
			
		//	if ( i % 100000 == 0 )
			//	cout<<"i"<<i<<endl;
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			double mapa_1 = chr_Par[chr].par_ve[pid1].mapa;
			double freq_kmer_1 = chr_Par[chr].par_ve[pid1].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;		
		/*	int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;  */
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			double mapa_2 = chr_Par[chr].par_ve[pid2].mapa;
			double freq_kmer_2 = chr_Par[chr].par_ve[pid2].freq_GC;
			if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
				continue;
		/*	int mapa_id_2 = (int)(mapa_2 * 10);
			int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000);
			if ( freq_kmer_id_2 >= 10 )
				freq_kmer_id_2 = 10;  */
			
			double acc_m = acc_1 * acc_2;
			int log_acc_m_10_id = transfer_acc_m( acc_m );
			double mapa_m = mapa_1 * mapa_2;
			int mapa_m_id = (int)( mapa_m * 1000);
			if ( mapa_m_id > 10 )
				mapa_m_id = 10;
			
			double freq_m = freq_kmer_1 * freq_kmer_2;
			int freq_m_id = (int)( freq_m * 10000000 );
			if ( freq_m_id > 10 )
				freq_m_id = 10;
			
		//	pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
		//	pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
			int length = pid2 - pid1;
			if ( length < 0 )
				length *= (-1);
			
			if ( space[length][log_acc_m_10_id][mapa_m_id].find( freq_m_id ) ==  space[length][log_acc_m_10_id][mapa_m_id].end() )
				space[length][log_acc_m_10_id][mapa_m_id][freq_m_id] = c;
			else
				space[length][log_acc_m_10_id][mapa_m_id][freq_m_id] += c;
		}
	}
	
	vector<int > min_ve;
	vector<int > sec_ve;
	vector<int > thr_ve;
	map<int, vector<int > > acc_m_ve;
	map<int, vector<int > > length_ve;
	for ( map<int, map<int, map<int, map<int, int > > > >::iterator ite = space.begin(); ite != space.end(); ++ite )
	{
		int t_c = 0;
		for ( map<int, map<int, map<int, int > > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int s_c = 0;
			for ( map<int, map<int, int  > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<int, int >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{	
					int c = fi->second;
					min_ve.push_back(c);
					s_c += c;
				}
			}
			sec_ve.push_back( s_c );
			acc_m_ve[si->first].push_back( s_c );
			length_ve[ite->first].push_back( s_c );
			t_c += s_c;
			
		}
		
		thr_ve.push_back(t_c );
	}
	
	cout<<"minimal ve: size:"<<min_ve.size()<<" sum:"<<getsum( min_ve )<<" max:"<<getmax(min_ve )<<" medium:"<<getmedium(min_ve)<<endl;
	cout<<"second ve: size:"<<sec_ve.size()<<" sum:"<<getsum( sec_ve )<<" max:"<<getmax(sec_ve )<<" medium:"<<getmedium(sec_ve)<<endl;
	cout<<"third ve: size:"<<thr_ve.size()<<" sum:"<<getsum( thr_ve )<<" max:"<<getmax(thr_ve )<<" medium:"<<getmedium(thr_ve)<<endl;
	
/*	ofstream outf( "acc_m_OBSspace.density.txt" );
	for ( map<int, vector<int > >::iterator ite = acc_m_ve.begin(); ite != acc_m_ve.end(); ++ite )
	{
		double acc_m = pow( 10, ((ite->first*1.0 / 10)) );
		
		outf<<acc_m<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf.close();
	
	ofstream outf2("length_OBSspace.density.txt");
	for ( map<int, vector<int > >::iterator ite = length_ve.begin(); ite != length_ve.end(); ++ite )
	{
		outf2<<ite->first<<"\t"<<ite->second.size()<<"\t"<<getsum( ite->second )<<"\t"<<getmedium( ite->second )<<endl;
		
	}
	outf2.close();  */
}


void Partition_Bank::simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > & par_pair_c,
	map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > & space,
	map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > & obsspace,
	int times,
//	map<string, map<pair<size_t, size_t >, pair<double, double > > > &par_pair_pv,
	string outfile )
{
	int omit_1 = 0;
	int omit_2 = 0;
//	map< int, map<int, map<pair<int, int>, set<pair<int, int > > > > > targets;
	map<string, map< pair<size_t, size_t >, vector<int > > > par_pair_target;
	map<int, map<int, map<pair<int, int>, map<pair<int, int >, vector<pair<string, pair<size_t, size_t > > > > > > > targets;
	
	cout<<"get par pair target"<<endl;
	int o = 0;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
	//	cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 1000 || pid1 > pid2 + 1000 )
				continue;
			int c = si->second;
			if ( c < 3 )
				continue;
			
			o += 1;
			if ( o % 100000 == 0 )
				cout<<o<<endl;
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			double mapa_1 = chr_Par[chr].par_ve[pid1].mapa;
			double freq_kmer_1 = chr_Par[chr].par_ve[pid1].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;		
			int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			double mapa_2 = chr_Par[chr].par_ve[pid2].mapa;
			double freq_kmer_2 = chr_Par[chr].par_ve[pid2].freq_GC;
			if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
				continue;
			int mapa_id_2 = (int)(mapa_2 * 10);
			int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000);
			if ( freq_kmer_id_2 >= 10 )
				freq_kmer_id_2 = 10;
			
			double acc_m = acc_1 * acc_2;
			int log_acc_m_10_id = transfer_acc_m( acc_m );
			
			pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
			pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
			int length = pid2 - pid1;
			if ( length < 0 )
				length *= (-1);
				
			if ( space[length][log_acc_m_10_id][mapa_pair].find( freq_pair ) == space[length][log_acc_m_10_id][mapa_pair].end() )
			{
				omit_1 += 1;
				continue;
			}
			if ( obsspace[length][log_acc_m_10_id][mapa_pair].find( freq_pair ) == obsspace[length][log_acc_m_10_id][mapa_pair].end() )
			{
				omit_2 += 1;
				continue;
			}
			
			targets[length][log_acc_m_10_id][mapa_pair][freq_pair].push_back( make_pair(chr, si->first ) );
			
			
		/*	vector<int > t_ve;
			t_ve.push_back( length );
			t_ve.push_back( log_acc_m_10_id );
			t_ve.push_back( mapa_pair.first );
			t_ve.push_back( mapa_pair.second );
			t_ve.push_back( freq_pair.first );
			t_ve.push_back( freq_pair.second );
			par_pair_target[chr][si->first] = t_ve;  */
			
		}
	}
	cout<<"total "<<o<<" pairs"<<endl;
	cout<<"omit "<<omit_1<<" "<<omit_2<<endl;
	
	cout<<"simulation for each target"<<endl;
// 	estimate by binomial distribution
	o = 0;
	int o_tag = 0;
//	map< int, map<int, map<pair<int, int>, map<pair<int, int >, vector< vector<int > > > > > > target_sim;
//	map<string, map<pair<size_t, size_t >, pair<int, pair<double, double > > > > par_pair_pv;
	map<string, map<pair<size_t, size_t >, pair<int, double > > > par_pair_pv;
	for ( map< int, map<int, map<pair<int, int>, map<pair<int, int >, vector< pair<string, pair<size_t, size_t > > > > > > >::iterator ite = targets.begin();
		ite != targets.end(); ++ite )
	{
		for ( map<int, map<pair<int, int>, map<pair<int, int >, vector< pair<string, pair<size_t, size_t > > > > > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			for ( map< pair<int, int >, map<pair<int, int >, vector< pair<string, pair<size_t, size_t > > > > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<pair<int, int >, vector< pair<string, pair<size_t, size_t > > > >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{
					o += fi->second.size();
					if ( o / 10000 > o_tag )
					{
						cout<<o<<endl;
						o_tag = o/ 10000;
					}
					int ve_size = space[ite->first][si->first][ti->first][fi->first];
					int total_pet_n = obsspace[ite->first][si->first][ti->first][fi->first];
				/*	vector< vector<int > > mat;
					sim_func( ve_size, total_pet_n, times, mat );
					
					int size_ve = (int)mat[0].size();
				
					multiset<int, greater<int> > total_col_set;
				
					int cc = 0;
					for ( int j = 0; j < ve_size; ++j )
					{
						for ( int i = 0; i < times; ++i )
						{
							if ( cc >= 100000 )
								break;
							total_col_set.insert( mat[i][j]);
							cc += 1;
							
						}
						if ( cc >= 100000 )
							break;
					}	*/
					
					for ( size_t k = 0; k < fi->second.size(); ++k )
					{
						int c = par_pair_c[fi->second[k].first][fi->second[k].second];
					/*	int r2 = 0;
						for ( multiset<int, greater<int> >::iterator ii = total_col_set.begin(); 
							ii != total_col_set.end(); ++ii )
						{
		
							if ( c > *ii )
							{
								break;
							}
							r2 += 1;
						}  
						double p2 = r2 * 1.0 / (int)total_col_set.size();
	
						
						*/
						double bp = 1.0 / ve_size;
						
						double p2 = cdf( c, total_pet_n, bp );
						
						par_pair_pv[fi->second[k].first][fi->second[k].second] = make_pair( ve_size, p2 );
					}
				}
			}
		}
	}
	
	cout<<"output pv"<<endl;
	
	ofstream outf(outfile.data() );
	outf<<"## simulation times "<<times<<endl;
	outf<<"#chrom\tpar1_start\tpart1_end\tpar2_start\tpart2_end\tdistance\tpart1_PE_count\tpart2_PE_count\tpart1_PE_density\tpart2_PE_density\tpart1_mappability\tpart2_mappability\tpart1_CATG_freq\tpart2_CATG_freq\tcongener_space\tpet_count\tp_value\n";
	for ( map<string, map<pair<size_t, size_t >, pair<int, double > > >::iterator ite = par_pair_pv.begin();
		ite != par_pair_pv.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, pair<int, double > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start - chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].fea_ACC
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].fea_ACC
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].mapa
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].mapa
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].freq_GC
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].freq_GC
				<<"\t"<<si->second.first
				<<"\t"<<par_pair_c[chr][si->first]
				<<"\t"<<si->second.second<<endl;
			
		}
	}
	outf.close();
/*	for ( map<string, map< pair<size_t, size_t >, vector<int > > >::iterator ite = par_pair_target.begin(); 
		ite != par_pair_target.end(); ++ite )
	{
	
		string chr = ite->first;
		for ( map< pair<size_t, size_t >, vector<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
		
			int c = par_pair_c[chr][si->first];
			int length = si->second[0];
			int acc = si->second[1];
			pair<int, int > mapa = make_pair( si->second[2], si->second[3] );
			pair<int, int > freq = make_pair( si->second[4], si->second[5] );
			vector< vector<int > > mat = target_sim[length][acc][mapa][freq];
			
			int size_ve = (int)mat[0].size();
			
			pair<double, double > pv = cal_pv( c, mat );
			
		//	par_pair_pv[chr][si->first] = pv;
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start+chr_Par[chr].win-1
				<<"\t"<<size_ve<<"\t"<<pv.first<<"\t"<<pv.second<<endl;
			
		}
		
	}  */
	
}

void Partition_Bank::simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > & par_pair_c,
	map<int, map<int, map<int, map<int, int > > > > & space,
	map<int, map<int, map<int, map<int, int > > > > & obsspace,
	int times,
	string outfile )
{
	int omit_1 = 0;
	int omit_2 = 0;
//	map<string, map< pair<size_t, size_t >, vector<int > > > par_pair_target;
	map<int, map<int, map<int, map<int, vector<pair<string, pair<size_t, size_t > > > > > > > targets;
	
	cout<<"get par pair target"<<endl;
	int o = 0;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
	//	cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			
			if ( pid1 == pid2 )
				continue;
			if ( pid2 > pid1 + 1000 || pid1 > pid2 + 1000 )
				continue;
			int c = si->second;
			if ( c < 3 )
				continue;
			
			o += 1;
			if ( o % 100000 == 0 )
				cout<<o<<endl;
			
			double acc_1 = chr_Par[chr].par_ve[pid1].fea_ACC;
			double mapa_1 = chr_Par[chr].par_ve[pid1].mapa;
			double freq_kmer_1 = chr_Par[chr].par_ve[pid1].freq_GC;
			if ( acc_1 == 0 || mapa_1 == 0 || freq_kmer_1 == 0 )
				continue;		
		//	int mapa_id_1 = (int)(mapa_1 * 10);
		//	int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
		//	if ( freq_kmer_id_1 >= 10 )
		//		freq_kmer_id_1 = 10;
				
			double acc_2 = chr_Par[chr].par_ve[pid2].fea_ACC;
			double mapa_2 = chr_Par[chr].par_ve[pid2].mapa;
			double freq_kmer_2 = chr_Par[chr].par_ve[pid2].freq_GC;
			if ( acc_2 == 0 || mapa_2 == 0 || freq_kmer_2 == 0 )
				continue;
		/*	int mapa_id_2 = (int)(mapa_2 * 10);
			int freq_kmer_id_2 = (int)(freq_kmer_2 * 1000);
			if ( freq_kmer_id_2 >= 10 )
				freq_kmer_id_2 = 10;  */
			
			double acc_m = acc_1 * acc_2;
			int log_acc_m_10_id = transfer_acc_m( acc_m );
			double mapa_m = mapa_1 * mapa_2;
			int mapa_m_id = (int)( mapa_m * 1000);
			if ( mapa_m_id > 10 )
				mapa_m_id = 10;
			
			double freq_m = freq_kmer_1 * freq_kmer_2;
			int freq_m_id = (int)( freq_m * 10000000 );
			if ( freq_m_id > 10 )
				freq_m_id = 10;
			
			
		//	pair<int, int > mapa_pair = make_pair( mapa_id_1, mapa_id_2 );
		//	pair<int, int > freq_pair = make_pair( freq_kmer_1, freq_kmer_2 );
			int length = pid2 - pid1;
			if ( length < 0 )
				length *= (-1);
				
			if ( space[length][log_acc_m_10_id][mapa_m_id].find( freq_m_id ) == space[length][log_acc_m_10_id][mapa_m_id].end() )
			{
				omit_1 += 1;
				continue;
			}
			if ( obsspace[length][log_acc_m_10_id][mapa_m_id].find( freq_m_id ) == obsspace[length][log_acc_m_10_id][mapa_m_id].end() )
			{
				omit_2 += 1;
				continue;
			}
			
			targets[length][log_acc_m_10_id][mapa_m_id][freq_m_id].push_back( make_pair(chr, si->first ) );
			
			
		/*	vector<int > t_ve;
			t_ve.push_back( length );
			t_ve.push_back( log_acc_m_10_id );
			t_ve.push_back( mapa_pair.first );
			t_ve.push_back( mapa_pair.second );
			t_ve.push_back( freq_pair.first );
			t_ve.push_back( freq_pair.second );
			par_pair_target[chr][si->first] = t_ve;  */
			
		}
	}
	cout<<"total "<<o<<" pairs"<<endl;
	cout<<"omit "<<omit_1<<" "<<omit_2<<endl;
	
	cout<<"simulation for each target"<<endl;
	o = 0;
	int o_tag = 0;
//	map< int, map<int, map<pair<int, int>, map<pair<int, int >, vector< vector<int > > > > > > target_sim;
//	map<string, map<pair<size_t, size_t >, pair<int, pair<double, double > > > > par_pair_pv;
	map<string, map<pair<size_t, size_t >, pair<int, double > > > par_pair_pv;
	for ( map< int, map<int, map<int, map<int, vector< pair<string, pair<size_t, size_t > > > > > > >::iterator ite = targets.begin();
		ite != targets.end(); ++ite )
	{
		for ( map<int, map<int, map<int, vector< pair<string, pair<size_t, size_t > > > > > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			for ( map< int, map<int, vector< pair<string, pair<size_t, size_t > > > > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				for ( map<int, vector< pair<string, pair<size_t, size_t > > > >::iterator fi = ti->second.begin(); fi != ti->second.end(); ++fi )
				{
					o += fi->second.size();
					if ( o / 10000 > o_tag )
					{
						cout<<o<<endl;
						o_tag = o/ 10000;
					}
					int ve_size = space[ite->first][si->first][ti->first][fi->first];
					int total_pet_n = obsspace[ite->first][si->first][ti->first][fi->first];
					vector< vector<int > > mat;
					sim_func( ve_size, total_pet_n, times, mat );
				//	target_sim[ite->first][si->first][ti->first][*fi] = mat;
					
					int size_ve = (int)mat[0].size();
				//	pair<double, double > pv = cal_pv( c, mat );
					multiset<int, greater<int> > total_col_set;
				/*	for ( int i = 0; i < times; ++i )
					{
						for ( int j = 0; j < ve_size; ++j )
						{
							total_col_set.insert( mat[i][j]);
						}
					} */
					int cc = 0;
					for ( int j = 0; j < ve_size; ++j )
					{
						for ( int i = 0; i < times; ++i )
						{
							if ( cc >= 100000 )
								break;
							total_col_set.insert( mat[i][j]);
							cc += 1;
							
						}
						if ( cc >= 100000 )
							break;
					}	
					
					for ( size_t k = 0; k < fi->second.size(); ++k )
					{
						int c = par_pair_c[fi->second[k].first][fi->second[k].second];
						int r2 = 0;
						for ( multiset<int, greater<int> >::iterator ii = total_col_set.begin(); 
							ii != total_col_set.end(); ++ii )
						{
		
							if ( c > *ii )
							{
								break;
							}
							r2 += 1;
						}
						double p2 = r2 * 1.0 / (int)total_col_set.size();
	
						
						par_pair_pv[fi->second[k].first][fi->second[k].second] = make_pair( size_ve, p2 );
					}
				}
			}
		}
	}
	
	cout<<"output pv"<<endl;
	
	ofstream outf(outfile.data() );
	outf<<"## simulation times "<<times<<endl;
	outf<<"#chrom\tpar1_start\tpart1_end\tpar2_start\tpart2_end\tdistance\tpart1_PE_count\tpart2_PE_count\tpart1_PE_density\tpart2_PE_density\tpart1_mappability\tpart2_mappability\tpart1_CATG_freq\tpart2_CATG_freq\tcongener_space\tpet_count\tp_value\n";
	for ( map<string, map<pair<size_t, size_t >, pair<int, double > > >::iterator ite = par_pair_pv.begin();
		ite != par_pair_pv.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<size_t, size_t >, pair<int, double > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start - chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].end_count
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].fea_ACC
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].fea_ACC
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].mapa
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].mapa
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].freq_GC
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].freq_GC
				<<"\t"<<si->second.first
				<<"\t"<<par_pair_c[chr][si->first]
				<<"\t"<<si->second.second<<endl;
			
		}
	}
	outf.close();
/*	for ( map<string, map< pair<size_t, size_t >, vector<int > > >::iterator ite = par_pair_target.begin(); 
		ite != par_pair_target.end(); ++ite )
	{
	
		string chr = ite->first;
		for ( map< pair<size_t, size_t >, vector<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
		
			int c = par_pair_c[chr][si->first];
			int length = si->second[0];
			int acc = si->second[1];
			pair<int, int > mapa = make_pair( si->second[2], si->second[3] );
			pair<int, int > freq = make_pair( si->second[4], si->second[5] );
			vector< vector<int > > mat = target_sim[length][acc][mapa][freq];
			
			int size_ve = (int)mat[0].size();
			
			pair<double, double > pv = cal_pv( c, mat );
			
		//	par_pair_pv[chr][si->first] = pv;
			outf<<chr<<"\t"<<chr_Par[chr].par_ve[si->first.first].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.first].start+chr_Par[chr].win-1
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start
				<<"\t"<<chr_Par[chr].par_ve[si->first.second].start+chr_Par[chr].win-1
				<<"\t"<<size_ve<<"\t"<<pv.first<<"\t"<<pv.second<<endl;
			
		}
		
	}  */
	
}

void Partition_Bank::generate_potential_space( map<string, vector<pair<size_t, size_t > > > &space,
	size_t minstep, size_t maxstep, int min_count, double min_feaACC )
{	
	int c = 0;
	for ( map<string, Partition >::iterator ite = chr_Par.begin(); ite != chr_Par.end(); ++ite )
	{
		string chr = ite->first;
	//	if ( chr != "chr1" )
		//	continue;
		cout<<chr<<endl;
		for ( size_t i = 0; i < ite->second.par_ve.size(); ++i )
		{
		//	if ( i % 10000 == 0 )
			//	cout<<i<<endl;
			int end_count_1 = ite->second.par_ve[i].end_count;
			if ( end_count_1 < min_count )
				continue;
			double feaACC1 = ite->second.par_ve[i].fea_ACC;
			if ( feaACC1 < min_feaACC )
				continue;
			if ( !ite->second.par_ve[i].peak_anchor )
				continue;
		/*	int mapa_id_1 = (int)(mapa_1 * 10);
			int freq_kmer_id_1 = (int)(freq_kmer_1 * 1000);
			if ( freq_kmer_id_1 >= 10 )
				freq_kmer_id_1 = 10;  */
				
			for ( size_t j = minstep; j < maxstep; ++j )
			{
				size_t k = i+j;
				if ( k >= ite->second.par_ve.size() )
					break;
				if ( !ite->second.par_ve[k].peak_anchor )
					continue;
				int end_count_2 = ite->second.par_ve[k].end_count;
				if ( end_count_2 < min_count )
					continue;
				double feaACC2 = ite->second.par_ve[k].fea_ACC;
				if ( feaACC2 < min_feaACC )
					continue;
				
				space[chr].push_back( make_pair( i, k ) );
				c += 1;
			}
		}
	}
	cout<<c<<endl;
}

void Partition_Bank::space_dis_accm_distr( map<string, vector<pair<size_t, size_t > > > &space,
	string prefix )
{
	string outfile0 = prefix+".distance.dist0.txt";
	ofstream outf0( outfile0.data() );
	string outfile1 = prefix+".distance.dist.txt";
	ofstream outf1( outfile1.data() );
	string outfile2 = prefix+".accm.dist.txt";
	ofstream outf2( outfile2.data() );
	string outfile3 = prefix+".total_valid_space_size.txt";
	ofstream outf3(outfile3.data() );
	map<int, double > l_20_c;
	map<int, double > logrpm_c;
	long total = 0;
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = space.begin(); 
		ite != space.end(); ++ite )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		
		for ( vector<pair<size_t, size_t > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			double rpm_m = chr_Par[chr].par_ve[si->first].fea_ACC * chr_Par[chr].par_ve[si->second].fea_ACC;
			int logrpm = (int)(log2(rpm_m*1.0+1) * 10);
		//	if ( logrpm < 10 )
			//	continue;
			
			if ( logrpm_c.find( logrpm ) != logrpm_c.end() )
			{
				logrpm_c[logrpm] += 1;
			} else
			{
				logrpm_c[logrpm] = 1;
			}
			
			total += 1;
			int l = (int)(si->second - si->first );
			int l_20 = l / 20;
			if ( l_20_c.find( l_20 ) != l_20_c.end() )
			{
				l_20_c[l_20] += 1;
				
			} else
			{
				l_20_c[l_20] = 1;
			}
			
			
		}
	}
	
	for ( map<int, double >::iterator ite = l_20_c.begin(); ite != l_20_c.end(); ++ite )
	{
		double r = ite->second / total;
		if ( ite == l_20_c.begin() )
		{
			outf0<<ite->first<<"\t"<<r<<endl;
			continue;
		}
		
		outf1<<ite->first<<"\t"<<r<<endl;
	}
	outf0.close();
	outf1.close();
	for ( map<int, double >::iterator ite = logrpm_c.begin(); ite != logrpm_c.end(); ++ite  )
	{
	//	if ( ite->first < 10 )
		//	continue;
		double r = ite->second / total;
		outf2<<ite->first<<"\t"<<r<<endl;
	}
	outf2.close();
	outf3<<total<<endl;
	outf3.close();
}

void Partition_Bank::space_dis_accm_distr( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
	size_t minstep, size_t maxstep, string prefix )
{
	string outfile0 = prefix+".distance.dist0.txt";
	ofstream outf0( outfile0.data() );
	string outfile1 = prefix+".distance.dist.txt";
	ofstream outf1( outfile1.data() );
	string outfile2 = prefix+".accm.dist.txt";
	ofstream outf2( outfile2.data() );
	string outfile3 = prefix+".total_valid_space_size.txt";
	ofstream outf3(outfile3.data() );
	map<int, double > l_20_c;
	map<int, double > logrpm_c;
	long total = 0;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); 
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int p = si->second;
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			int l = (int)(pid2 - pid1 );
			
			if ( l < minstep || l >= maxstep )
				continue;
			
			double rpm_m = chr_Par[chr].par_ve[pid1].fea_ACC * chr_Par[chr].par_ve[pid2].fea_ACC;
			int logrpm = (int)(log2(rpm_m*1.0+1) * 10);
		//	if ( logrpm < 10 )
			//	continue;
			
			if ( logrpm_c.find( logrpm ) != logrpm_c.end() )
			{
				logrpm_c[logrpm] += 1;
			} else
			{
				logrpm_c[logrpm] = 1;
			}
			
			total += 1;
			
			int l_20 = l / 20;
			if ( l_20_c.find( l_20 ) != l_20_c.end() )
			{
				l_20_c[l_20] += 1;
				
			} else
			{
				l_20_c[l_20] = 1;
			}
			
			
		}
	}
	
	for ( map<int, double >::iterator ite = l_20_c.begin(); ite != l_20_c.end(); ++ite )
	{
		double r = ite->second / total;
		if ( ite == l_20_c.begin() )
		{
			outf0<<ite->first<<"\t"<<r<<endl;
			continue;
		}
		
		outf1<<ite->first<<"\t"<<r<<endl;
	}
	outf0.close();
	outf1.close();
	for ( map<int, double >::iterator ite = logrpm_c.begin(); ite != logrpm_c.end(); ++ite  )
	{
	//	if ( ite->first < 10 )
		//	continue;
		double r = ite->second / total;
		outf2<<ite->first<<"\t"<<r<<endl;
	}
	outf2.close();
	outf3<<total<<endl;
	outf3.close();
}


void Partition_Bank::pets_dis_accm_distr( map<string, vector<pair<size_t, size_t > > > &space,
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		string prefix)
{
	string outfile0 = prefix+".distance.PETs_dist0.txt";
	ofstream outf0( outfile0.data() );
	string outfile1 = prefix+".distance.PETs_dist.txt";
	ofstream outf1( outfile1.data() );
	string outfile2 = prefix+".accm.PETs_dist.txt";
	ofstream outf2( outfile2.data() );
	string outfile3 = prefix+".total_valid_PETs_in_anchor.txt";
	ofstream outf3(outfile3.data() );
	map<int, int > logl_c;
	map<int, int > logrpm_c;
	int total = 0;
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = space.begin(); 
		ite != space.end(); ++ite )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<pair<size_t, size_t > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			int p = 0;
			if ( par_pair_c[chr].find( *si ) != par_pair_c[chr].end() )
			{
				p = par_pair_c[chr][*si];
			} 
			
			double rpm_m = chr_Par[chr].par_ve[si->first].fea_ACC * chr_Par[chr].par_ve[si->second].fea_ACC;
			int logrpm = (int)(log2(rpm_m*1.0+1) * 10);
			
		//	if ( logrpm < 10 )
			//	continue;
				
			if ( logrpm_c.find( logrpm ) != logrpm_c.end() )
			{
				logrpm_c[logrpm] += p;
			} else
				logrpm_c[logrpm] = p;
				
			
			int l = (int)(si->second - si->first);
			int loglength = (int)(log2(l*1.0)*10);
			if ( logl_c.find(loglength) != logl_c.end() )
			{
				logl_c[loglength] += p;
			} else
				logl_c[loglength] = p;
			
			total += p;
		}
	}
	
	cout<<"total pets in anchor "<<total<<endl;
	
	for ( map<int, int >::iterator ite = logl_c.begin(); ite != logl_c.end(); ++ite )
	{
		
		double r = (ite->second+1) * 1.0 / total;
		if ( ite == logl_c.begin() )
		{
			outf0<<ite->first<<"\t"<<r<<endl;
			continue;
		}
		outf1<<ite->first<<"\t"<<r<<endl;
	}
	outf0.close();
	outf1.close();
	for ( map<int, int >::iterator ite = logrpm_c.begin(); ite != logrpm_c.end(); ++ite  )
	{
	//	if ( ite->first < 10 )
		//	continue;
		double r = (ite->second+1) * 1.0 / total;
		outf2<<ite->first<<"\t"<<r<<endl;
	}
	outf2.close();
	outf3<<total<<endl;
	outf3.close();
}

void Partition_Bank::pets_dis_accm_distr( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
	size_t minstep, size_t maxstep,	string prefix)
{
	string outfile0 = prefix+".distance.PETs_dist0.txt";
	ofstream outf0( outfile0.data() );
	string outfile1 = prefix+".distance.PETs_dist.txt";
	ofstream outf1( outfile1.data() );
	string outfile2 = prefix+".accm.PETs_dist.txt";
	ofstream outf2( outfile2.data() );
	string outfile3 = prefix+".total_valid_PETs_in_anchor.txt";
	ofstream outf3(outfile3.data() );
	map<int, int > logl_c;
	map<int, int > logrpm_c;
	int total = 0;
	for ( map<string, map<pair<size_t, size_t >, int > >::iterator ite = par_pair_c.begin(); 
		ite != par_pair_c.end(); ++ite )
	{
		string chr = ite->first;
		cout<<chr<<endl;
		for ( map<pair<size_t, size_t >, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int p = si->second;
			size_t pid1 = si->first.first;
			size_t pid2 = si->first.second;
			
			int l = (int)(pid2 - pid1 );
			
			if ( l < minstep || l >= maxstep )
				continue;
				
			double rpm_m = chr_Par[chr].par_ve[pid1].fea_ACC * chr_Par[chr].par_ve[pid2].fea_ACC;
			int logrpm = (int)(log2(rpm_m*1.0+1) * 10);
			
		//	if ( logrpm < 10 )
			//	continue;
				
			if ( logrpm_c.find( logrpm ) != logrpm_c.end() )
			{
				logrpm_c[logrpm] += p;
			} else
				logrpm_c[logrpm] = p;
				
			int loglength = (int)(log2(l*1.0)*10);
			if ( logl_c.find(loglength) != logl_c.end() )
			{
				logl_c[loglength] += p;
			} else
				logl_c[loglength] = p;
			
			total += p;
		}
	}
	
	cout<<"total pets in anchor "<<total<<endl;
	
	for ( map<int, int >::iterator ite = logl_c.begin(); ite != logl_c.end(); ++ite )
	{
		
		double r = (ite->second+1) * 1.0 / total;
		if ( ite == logl_c.begin() )
		{
			outf0<<ite->first<<"\t"<<r<<endl;
			continue;
		}
		outf1<<ite->first<<"\t"<<r<<endl;
	}
	outf0.close();
	outf1.close();
	for ( map<int, int >::iterator ite = logrpm_c.begin(); ite != logrpm_c.end(); ++ite  )
	{
		if ( ite->first < 10 )
			continue;
		double r = (ite->second+1) * 1.0 / total;
		outf2<<ite->first<<"\t"<<r<<endl;
	}
	outf2.close();
	outf3<<total<<endl;
	outf3.close();
}

void Partition_Bank::cal_space_dis_acc_stat( map<string, vector<pair<size_t, size_t > > > &space,
	string prefix )
{
	string outfile1 = prefix+".length.txt";
	ofstream outf1( outfile1.data() );
	string outfile2 = prefix+".accm.txt";
	ofstream outf2( outfile2.data() );
	map<int, int > logl_c;
	map<int, int > accm_c;
	map<int, double > accm_int_rpm;
	map<int, int > logl_int;
	
	for ( map<string, vector<pair<size_t, size_t > > >::iterator ite = space.begin(); 
		ite != space.end(); ++ite )
	{
		string chr = ite->first;
		int i = 0;
		cout<<chr<<endl;
		for ( vector<pair<size_t, size_t > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
		//	if ( i % 100000 == 0 )
			//	cout<<i<<endl;
			++i;
			int l = (int)(si->second - si->first);
			int loglength = (int)(log2(l*1.0)*10);
			if ( logl_c.find( loglength ) == logl_c.end() )
			{
				logl_c[loglength] = 1;
				logl_int.insert( make_pair(loglength, l) );
			} else
				logl_c[loglength] += 1;
				
			int acc_m = chr_Par[chr].par_ve[si->first].end_count * chr_Par[chr].par_ve[si->second].end_count;
			int acc_m_100 = acc_m / 100;
			if ( accm_c.find( acc_m_100 ) == accm_c.end() )
			{
				accm_c[acc_m_100] = 1;
				double acc_m_rpm = chr_Par[chr].par_ve[si->first].fea_ACC * chr_Par[chr].par_ve[si->second].fea_ACC;
				accm_int_rpm.insert( make_pair(acc_m_100, acc_m_rpm) );
			} else
				accm_c[acc_m_100] += 1;
		}
	}
	for ( map<int, int >::iterator ite = logl_c.begin(); ite != logl_c.end(); ++ite )
	{
		outf1<<ite->first<<"\t"<<logl_int[ite->first]<<"\t"<<ite->second<<endl;
	}
	outf1.close();
	for ( map<int, int >::iterator ite = accm_c.begin(); ite != accm_c.end(); ++ite )
	{
		outf2<<ite->first<<"\t"<<accm_int_rpm[ite->first]<<"\t"<<ite->second<<endl;
	}
	outf2.close();
}

int transfer_acc_m( double acc_m )
{
	double acc_m_10 = acc_m * 10;
	if ( acc_m_10 >= 1000 )
		acc_m_10 = 1000;
	if ( acc_m_10 < 1 )
		acc_m_10 = 1;
	double log_acc_m_10 = log10 ( acc_m_10 );
	int log_acc_m_10_id = (int)(log_acc_m_10 * 10);
	return log_acc_m_10_id;
}

void sim_func( int ve_size, int total_pets, int times,
	vector< vector<int > > &mat )
{
	srand(time(NULL));
 
	for ( int i = 0; i < times; ++i )
	{
		vector<int > ve;
		for ( int j = 0; j < ve_size; ++j )
		{
			ve.push_back(0);
		}
		for ( int m = 0; m < total_pets; ++m )
		{
			int rd = std::rand();
			int rdn = rd % ve_size;
			ve[rdn] += 1;
		}
		mat.push_back(ve );
	}
	
/*	int rd = std::rand();
	int rdn = rd % ve_size;
	
	for ( int i = 0; i < times; ++i )
	{
		single_col_set.insert( mat[i][rdn] );
	}
	
	for ( int i = 0; i < times; ++i )
	{
		for ( int j = 0; j < ve_size; ++j )
		{
			total_col_set.insert( mat[i][j]);
		}
	} */
	
/*	int r1 = 0;
	for ( multiset<int, greater<int> >::iterator ite = single_col_set.begin(); ite != single_col_set.end(); ++ite )
	{
		
		if ( t_pet_count > *ite )
		{
			break;
		}
		r1 += 1;
	}
	double p1 = r1 * 1.0 / (int)single_col_set.size();
	
	int r2 = 0;
	for ( multiset<int, greater<int> >::iterator ite = total_col_set.begin(); ite != total_col_set.end(); ++ite )
	{
		
		if ( t_pet_count > *ite )
		{
			break;
		}
		r2 += 1;
	}
	double p2 = r2 * 1.0 / (int)total_col_set.size();
	
	return make_pair( p1, p2 );  */
}

pair<double, double > cal_pv( int c, vector< vector<int > > &mat )
{	
	multiset<int, greater<int> > single_col_set;
	multiset<int, greater<int> > total_col_set;
//	srand(time(NULL));
	int ve_size = (int)mat[0].size();
	int rd = std::rand();
	int rdn = rd % ve_size;
	
	for ( int i = 0; i < (int)mat.size(); ++i )
	{
		single_col_set.insert( mat[i][rdn] );
	}
	
	for ( int i = 0; i < (int)mat.size(); ++i )
	{
		for ( int j = 0; j < ve_size; ++j )
		{
			total_col_set.insert( mat[i][j]);
		}
	}
	
	int r1 = 0;
	for ( multiset<int, greater<int> >::iterator ite = single_col_set.begin(); ite != single_col_set.end(); ++ite )
	{
		
		if ( c > *ite )
		{
			break;
		}
		r1 += 1;
	}
	double p1 = r1 * 1.0 / (int)single_col_set.size();
	
	int r2 = 0;
	for ( multiset<int, greater<int> >::iterator ite = total_col_set.begin(); ite != total_col_set.end(); ++ite )
	{
		
		if ( c > *ite )
		{
			break;
		}
		r2 += 1;
	}
	
	double p2 = r2 * 1.0 / (int)total_col_set.size();
	
	return make_pair( p1, p2 );
	
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









