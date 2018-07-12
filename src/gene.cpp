#include "gene.h"


void Transcript::addTranscript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName)
{
	start = instart;
	end = inend;
	chr = inchr;
	strand = instrand;
	name = inname;
	geneName = ingeneName;
	cdsStart = incdsStart;
	cdsEnd = incdsEnd;
	exonCount = inexonCount;
	exonStarts = inexonStarts;
	exonEnds = inexonEnds;
}

void Transcript::addTranscript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName)
{
	start = instart;
	end = inend;
	chr = inchr;
	strand = instrand;
	name = inname;
	geneName = ingeneName;
}

int Transcript::getTSS()
{
	if ( strand == '+' )
		return start;
	else if ( strand == '-' )
		return end;
	else
	{
		cout<<"error getTSS unrecognized strand: "<<strand<<" "<<name<<endl;
		exit(1);
	}
}

int Transcript::getTES()
{
	if ( strand == '+' )
		return end;
	else if ( strand == '-' )
		return start;
	else
	{
		cout<<"error getTES unrecognized strand: "<<strand<<" "<<name<<endl;
		exit(1);
	}
}

pair<int, int > Transcript::getPromoter( int extension )
{
	int tss = getTSS();
	return make_pair( max(1, tss-extension), tss+extension-1 );
}

pair<int, int > Transcript::getPromoter( int up, int down )
{
	int tss = getTSS();
	if ( strand == '+' )
		return make_pair( max(1, tss-up), tss+down );
	else if ( strand == '-' )
		return make_pair( max(1, tss-down), tss+up );
	else
	{
		cout<<"error unknown strand "<<strand<<endl; 
		exit(1);
	}
}

pair<int, int > Transcript::getTESregion( int up, int down )
{
	int tes = getTES();
	if ( strand == '+' )
		return make_pair( max(1, tes-up), tes+down );
	else if ( strand == '-' )
		return make_pair( max(1, tes-down), tes+up );
	else
	{
		cout<<"error unknown strand "<<strand<<endl; 
		exit(1);
	}
}

pair<int, int > Transcript::getgenebody( int up, int down )
{
	int tss = getTSS();
	int tes = getTES();
	if ( strand == '+' )
		return make_pair( max(1, tss-up), tes+down );
	else if ( strand == '-' )
		return make_pair( max(1, tes-down), tss+up );
	else
	{
		cout<<"error unknown strand "<<strand<<endl; 
		exit(1);
	}
}

pair<int, int> Transcript::getPromoter()
{
	return getPromoter( 1000 );

}

void Gene::addTranscript(std::string name)
{
	TranscriptName.insert( name );
}

string Transcript::getgenetype()
{
	return name.substr(0,2);
}

map<string, map<pair<int, int>, string > > Genome::get_chr_promoter_TranscriptName_map( int up, int down )
{
	map<string, map<pair<int, int>, string > > pt;
/*	for ( map<string, map<pair<int, int >, string > >::iterator ite = chr_pos_TranscriptName_map.begin();
		ite != chr_pos_TranscriptName_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			pair<int, int > tpromoter = name_Transcript_map[seci->second].getPromoter();
			pt[chr][tpromoter] = seci->second;
		}
	} */
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpromoter = name_Transcript_map[repren].getPromoter( up, down );
			pt[chr][tpromoter] = repren;
		}
	}
	return pt;
}

map<string, map<pair<int, int>, string > > Genome::get_chr_TES_TranscriptName_map(int up, int down)
{
	map<string, map<pair<int, int>, string > > pt;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tesregion = name_Transcript_map[repren].getTESregion( up, down );
			pt[chr][tesregion] = repren;
		}
	}
	return pt;
}

void Genome::getintergenicregion( )
{
	cout<<"pppp"<<endl;
//	intergenicregion.clear();
	cout<<"ddd"<<endl;
	map<string, vector< pair<int, int> > > coll;
	cout<<chr_pos_TranscriptNameS_map.size()<<endl;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin(); ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
	//	cout<<ite->first<<endl;
		for ( map<pair<int, int >, vector<string> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			coll[ite->first].push_back( si->first );
		}
	}
	cout<<"1"<<endl;
	vector< map<string, vector< pair<int, int> > > > reg_ve;
	reg_ve.push_back( coll );
	map<string, vector< pair<int, int> > > merged;
	cout<<"s"<<endl;
	region_merge_naive( reg_ve, merged );
	cout<<"2"<<endl;
	cout<<merged.size()<<endl;
	map<string, vector<pair<int, int > > > intergen;
	for ( map<string, vector< pair<int, int> > >::iterator ite = merged.begin(); ite != merged.end(); ++ite )
	{
	//	cout<<"m "<<ite->first<<"\t "<<ite->second.size()<<endl;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( i == ite->second.size()-1 )
				break;
			int s = ite->second[i].second+1;
			int e = ite->second[i+1].first-1;
			intergen[ite->first].push_back( make_pair(s, e) );
			if ( i < 5 )
			{
			//	cout<<ite->first<<"\t"<<s<<"\t"<<e<<endl;
			}
		}
	}
	cout<<"dd"<<endl;
	intergenicregion = intergen;
	cout<<"kk"<<endl;
//	exit(1);
}

map<string, set<pair<int, int > > > Genome::get_chr_promoter( int up, int down )
{
	map<string, set<pair<int, int > > > pt;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpromoter = name_Transcript_map[repren].getPromoter( up, down );
			pt[chr].insert(tpromoter);
		}
	}
	return pt;
}

map<string, set<pair<int, int > > > Genome::get_genebody( int up, int down )
{
	map<string, set<pair<int, int > > > pt;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpgenebody = name_Transcript_map[repren].getgenebody( up, down );
		/*	if ( tpgenebody.first == 102513726 )
			{
				cout<<chr<<" "<<repren<<" "<<seci->first.first<<" "<<seci->first.second<<endl;
			} */
			pt[chr].insert(tpgenebody);
		}
	}
	return pt;
}

map<string, set<int > > Genome::get_chr_TSS()
{
	map<string, set<int > > pt;
/*	for ( map<string, map<pair<int, int >, string > >::iterator ite = chr_pos_TranscriptName_map.begin();
		ite != chr_pos_TranscriptName_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			int tss = name_Transcript_map[seci->second].getTSS();
			pt[chr].insert(tss);
		}
	} */
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			int tss = name_Transcript_map[repren].getTSS();
			pt[chr].insert(tss);
		}
	}
	return pt;
}

map<string, set<int > > Genome::get_chr_TES()
{
	map<string, set<int > > pt;

	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			int tss = name_Transcript_map[repren].getTES();
			pt[chr].insert(tss);
		}
	}
	return pt;
}

vector<int > Genome::get_TSS_from_Gene( string gene)
{
	vector<int > res;
	if ( name_Gene_map.find(gene) == name_Gene_map.end() )
	{	
		cout<<"error name_Gene_map not find gene "<<gene<<endl; exit(1);
	}
	set<string > transname = name_Gene_map[gene].TranscriptName;
	for ( set<string >::iterator ite = name_Gene_map[gene].TranscriptName.begin(); ite != name_Gene_map[gene].TranscriptName.end(); ++ite )
	{
		if ( name_Transcript_map.find( *ite ) ==  name_Transcript_map.end() )
		{	
			cout<<"error name_Transcript_map not find trans "<<*ite<<endl; exit(1);
		}
		int tss = name_Transcript_map[*ite].getTSS();
		res.push_back(tss);
	}
	
	return res;
}

void Genome::transcripttogene()
{
	map< string, vector<string > > gene_transcript_map;
	for ( map<string, Transcript >::iterator ite = name_Transcript_map.begin(); ite != name_Transcript_map.end(); ++ite )
	{
		string genename = ite->second.geneName;
		gene_transcript_map[genename].push_back( ite->first );
	}
	for ( map< string, vector<string > >::iterator ite = gene_transcript_map.begin(); ite != gene_transcript_map.end(); ++ite )
	{
		string genename = ite->first;
		Gene gn;
		int start = -1;
		int end = -1;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			string trans = ite->second[i];
			gn.addTranscript( trans );
			gn.chr = name_Transcript_map[trans].chr;
			gn.strand = name_Transcript_map[trans].strand;
			gn.type = name_Transcript_map[trans].getgenetype();
			if ( start == -1 )
			{
				start = name_Transcript_map[trans].start;
				end = name_Transcript_map[trans].end;
			} else
			{
				if ( start > name_Transcript_map[trans].start )
					start = name_Transcript_map[trans].start;
				if ( end < name_Transcript_map[trans].end )
					end = name_Transcript_map[trans].end;
			}
		}
		gn.start = start;
		gn.end = end;
		name_Gene_map.insert( make_pair( genename, gn ) );
		chr_pos_GeneName_map[gn.chr][make_pair(start, end )] = genename;
	}
	cout<<"trans gene "<<name_Gene_map.begin()->first<<endl; 
}

void Genome::addgenomeseq(string infile)
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	getline( inf, line );
	if (line[0] != '>' )
	{
		cout<<"error wrong format fasta"<<endl; exit(1);
	}
	string header = line.substr(1);
	string seq = "";
	while (!inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		while( line[0] != '>' )
		{
			seq += line;
			getline(inf, line );
			if ( inf.eof() )
				break;
			if ( line.empty() )
				break;
		}
		if ( !seq.empty() )
		{
			genomeseq.insert( make_pair( header, seq ) );
		}
		if ( line.empty() )
			break;
		if ( line[0] == '>' )
		{
			header = line.substr(1);
			seq = "";
		} 
	}
	inf.close();
	
	cout<<"norm seq"<<endl;
	norm_seq();
}

void Genome::norm_seq()
{
	for ( map<string, string >::iterator ite = genomeseq.begin(); ite != genomeseq.end(); ++ite )
	{
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( ite->second[i] <= 'z' && ite->second[i] >= 'a' )
			{
				ite->second[i] -= 32;
			}
		}
	}
}

string Genome::getsubseq(string chr, int pos, int len )
{
	if ( genomeseq.find( chr ) == genomeseq.end() )
	{
		cout<<"error not find chr in genome: "<<chr<<endl; exit(1);
		
	} 
	if ( pos+len > (int)genomeseq[chr].size() )
	{
		cout<<"Warning coord exceeds chrom "<<chr<<": "<<genomeseq[chr].size()<<" "<<pos<<" "<<len<<endl; 
		exit(1);
	}
	return genomeseq[chr].substr(pos, len );
}

void Genome::get_kmer_distri( string kmer, map<string, vector<bool > > &plus, map<string, vector<bool > > &minus )
{
	int ks = (int)kmer.size();
	string rev_kmer = kmer;
	reverse(rev_kmer.begin(), rev_kmer.end() );
	for ( map<string, string >::iterator ite = genomeseq.begin(); ite != genomeseq.end(); ++ite )
	{
		string chr = ite->first;
		if ( chr.find("_") != std::string::npos || chr == "chrM" )
			continue;
		cout<<chr<<endl;
		int l = (int)ite->second.size();
		vector<bool > v;
		for ( int i = 0; i < l; ++i )
		{
			if ( i <= l - ks )
			{
				string subs = getsubseq( chr, i, ks );
				if ( subs == kmer )
				{
					v.push_back(true);
				} else
					v.push_back(false);
			} else
				v.push_back(false);
				
		}
		if ( v.size() != ite->second.size() )
		{
			cout<<"error v.size() != ite->second.size() in Genome::get_kmer_distri "<<v.size() << ite->second.size()<<endl;
			exit(1);
		}
		
		plus[chr] = v;
		
		vector<bool > vm;
		for ( int i = 0; i < ks-1; ++i )
		{
			vm.push_back( false );
		}
		for ( int i = 0; i <= l - ks; ++i )
		{
			string subs = getsubseq( chr, i, ks );
			if ( subs == rev_kmer )
				vm.push_back( true );
			else
				vm.push_back( false );
		}
		
		if ( vm.size() != ite->second.size() )
		{
			cout<<"error vm.size() != ite->second.size() in Genome::get_kmer_distri "<<vm.size() << ite->second.size()<<endl;
			exit(1);
		}
		
		minus[chr] = v;
	}
}

void Genome::readtranscriptfromucsc(string &infile )
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
		vector<string > parseditem = parse_string( line );
	/*	if ( parseditem.size() != 11 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		} */
		vector<string> exonStarts = parse_string( parseditem[8], ',');
		vector<string> exonEnds = parse_string( parseditem[9], ',');
		vector<int > t_exonStarts;
		vector<int > t_exonEnds;
		if ( exonStarts.size() != exonEnds.size() )
		{
			cout<<"error exonStarts size != exonEnds size: "<<line<<endl;
			exit(1);
		}
		for ( size_t i = 0; i < exonStarts.size(); ++i )
		{
			t_exonStarts.push_back( atoi(exonStarts[i].c_str() ) );
			t_exonEnds.push_back( atoi(exonEnds[i].c_str() ) );
		}
	//	Transcript tob();
		Transcript tob( parseditem[0], parseditem[1], parseditem[2][0], atoi(parseditem[3].c_str()), atoi(parseditem[4].c_str()),
			atoi( parseditem[5].c_str()), atoi( parseditem[6].c_str()), atoi( parseditem[7].c_str() ), t_exonStarts, t_exonEnds, parseditem[10] );
		string chr = parseditem[1];
		string name = parseditem[0];
		if (tob.strand != '+' && tob.strand != '-')
		{
			cout<<"error strand "<<line<<endl; exit(1);
		}
		int start = atoi( parseditem[3].c_str() );
		int end = atoi( parseditem[4].c_str() );
	/*	if ( start == 102513726 )
		{
			cout<<chr<<" "<<name<<" "<<start<<" "<<end<<endl;
		} */
		if ( name_Transcript_map.find( name ) != name_Transcript_map.end() )
			continue;
		name_Transcript_map.insert( make_pair(name, tob) );
		chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
		chr_pos_TranscriptNameS_map[chr][make_pair(start, end)].push_back( name );
	}
	inf.close();
}

void Genome::readchrlen(string &infile )
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
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int l = atoi(parseditem[1].c_str());
		chr_len_map.insert(make_pair(chr, l) );
		
	}
	inf.close();
}

void Genome::read3genegroups_ParsingIndex( string infile, vector<string > &lowPI_Ts, vector<string > &highPI_Ts )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	getline( inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		string gn = parseditem[0];
		double rk = atof(parseditem[2].c_str() );
		if ( rk <= 0.5 )
		{
			lowPI_Ts.push_back( gn );
		} else
		{
			highPI_Ts.push_back( gn );
		}
	}
	cout<<"low PI genes:"<<lowPI_Ts.size()<<" high PI genes:"<<highPI_Ts.size()<<endl;
}

void Genome::setgenicregion( vector<string > &active_Ts, vector<string > &silent_Ts,
		map<string, vector<pair<int, int> > > &active_TSS, 
		map<string, vector<pair<int, int > > > &active_GeneBody,
		map<string, vector<pair<int, int > > > &active_TES,
		map<string, vector<pair<int, int > > > &silent_TSS,
		map<string, vector<pair<int, int > > > &silent_GeneBody,
		map<string, vector<pair<int, int > > > &silent_TES,
		map<string, vector<pair<int, int > > > &intergenic )   // no direction of genes
{
	int upstream = 1000;
	int downstream = 1000;
	map<string, set<pair<int, int > > > chr_gene;
	for ( size_t i = 0; i < active_Ts.size(); ++i )
	{
		string gene = active_Ts[i];
		int tss = name_Transcript_map[gene].getTSS();
		char strand = name_Transcript_map[gene].strand;
		string chr = name_Transcript_map[gene].chr;
		int up = tss - upstream;
		int down = tss + downstream;
		active_TSS[chr].push_back( make_pair(up, down) );
		
		int tes = name_Transcript_map[gene].getTES();
		int tesup = tes - upstream;
		int tesdown = tes + downstream;
		active_TES[chr].push_back( make_pair(tesup, tesdown) );
		
		if ( strand == '+' )
		{
			if ( tesup - down > 1000 )
			{
				active_GeneBody[chr].push_back( make_pair(down, tesup) );
			}
			chr_gene[chr].insert( make_pair( up, tesdown ) );
		} else
		{
			if ( up - tesdown > 1000 )
			{
				active_GeneBody[chr].push_back( make_pair(tesdown, up ) );
			}
			chr_gene[chr].insert( make_pair( tesup, down ) );
		}
	}
	
	for ( size_t i = 0; i < silent_Ts.size(); ++i )
	{
		string gene = silent_Ts[i];
		int tss = name_Transcript_map[gene].getTSS();
		char strand = name_Transcript_map[gene].strand;
		string chr = name_Transcript_map[gene].chr;
		int up = tss - upstream;
		int down = tss + downstream;
		silent_TSS[chr].push_back( make_pair(up, down) );
		
		int tes = name_Transcript_map[gene].getTES();
		int tesup = tes - upstream;
		int tesdown = tes + downstream;
		silent_TES[chr].push_back( make_pair(tesup, tesdown) );
		
		if ( strand == '+' )
		{
			if ( tesup - down > 1000 )
			{
				silent_GeneBody[chr].push_back( make_pair(down, tesup) );
			}
			chr_gene[chr].insert( make_pair( up, tesdown ) );
		} else
		{
			if ( up - tesdown > 1000 )
			{
				silent_GeneBody[chr].push_back( make_pair(tesdown, up ) );
			}
			chr_gene[chr].insert( make_pair( tesup, down ) );
		}
	}
	
//	map<string, vector<pair<int, int > > > chr_merged_generegion;
	for ( map<string, set<pair<int, int > > >::iterator ite = chr_gene.begin(); ite != chr_gene.end(); ++ite )
	{
		string chr = ite->first;
		vector<pair<int, int > > merged_region = mergeragions( ite->second );
		for ( size_t i = 1; i < merged_region.size(); ++i )
		{
			if ( merged_region[i].first - merged_region[i-1].second > 1000 )
			{
				intergenic[chr].push_back( make_pair(merged_region[i-1].second, merged_region[i].first ) );
			}
		}
	}
	
}


void Genome::getnunredundanttss_Ts( vector<string >& Ts, vector<string > &nr_Ts )
{
	map<string, set<pair<int, char > > > app_tss;
	for ( size_t i = 0; i < Ts.size(); ++i )
	{
		string T = Ts[i];
		int tss = name_Transcript_map[T].getTSS();
		string chr = name_Transcript_map[T].chr;
		char strand = name_Transcript_map[T].strand;
		if ( app_tss[chr].find(make_pair(tss,strand) ) != app_tss[chr].end() )
			continue;
			
		app_tss[chr].insert( make_pair(tss, strand) );
		nr_Ts.push_back(T);
	}
}

void Genome::getnunredundanttss_Ts(  vector<string > &nr_Ts )
{
	map<string, set<pair<int, char > > > app_tss;
	for ( map<string, Transcript >::iterator ite = name_Transcript_map.begin(); ite != name_Transcript_map.end(); ++ite )
	{
		string T = ite->first;
		int tss = ite->second.getTSS();
		string chr = ite->second.chr;
		char strand = ite->second.strand;
		if ( app_tss[chr].find(make_pair(tss,strand) ) != app_tss[chr].end() )
			continue;
			
		app_tss[chr].insert( make_pair(tss, strand) );
		nr_Ts.push_back(T);
	}
}

void Genome::remove_promoterregion( map<string, vector<int > > &regions, int upstream, int downstream )
{
	map<string, set<pair<int, int > > > prom = get_chr_promoter( upstream, downstream );
	map<string, vector<int > > alt_regions;
	for ( map<string, vector<int > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			bool ov = false;
			for ( set<pair<int, int > >::iterator ci = prom[chr].begin(); ci != prom[chr].end(); ++ci )
			{
				if ( ci->second < *si )
					continue;
				if ( ci->first > *si )
					break;
				ov = true;
				break;
			}
			if ( !ov )
				alt_regions[chr].push_back( *si );
		}
	}
	regions = alt_regions;
}

pair<string, int> getnearestgene( string chr, int start, int end, Genome &gn)
{
	string res = "";
	int dis = -1;
	if ( gn.chr_pos_TranscriptName_map.find(chr) == gn.chr_pos_TranscriptName_map.end() )
	{
		cout<<"error chr pos transcript not find "<<chr<<endl; exit(1);
	}
	for ( map<pair<int, int >, string >::iterator ite = gn.chr_pos_TranscriptName_map[chr].begin(); ite != gn.chr_pos_TranscriptName_map[chr].end(); ++ite )
	{
		string name = ite->second;
		int tss = gn.name_Transcript_map[name].getTSS();
		int d = min( abs(tss- start), abs(tss-end) );
		if ( dis == -1 )
		{
			dis = d;
			res = gn.name_Transcript_map[name].geneName;
		} else if ( d < dis )
		{
			dis = d;
			res = gn.name_Transcript_map[name].geneName;
		}
		
		if ( tss > end )
			break; 
	}
	return make_pair( res, dis );
}

void GeneExpression_GW::read_rpkm_from_file(string infile )
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
		vector<string > parseditem = parse_string( line );
		vector<string > names = parse_string(parseditem[0], '+');
		double rpkm = atof( parseditem[1].c_str() );
		for ( size_t i = 0; i < names.size(); ++i )
		{
			name_exp_map.insert( make_pair(names[i], rpkm) );
		}
	}
	inf.close();
	
}

void TranscriptExpression_GW::read_rpkm_from_file(string infile )
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
		vector<string > parseditem = parse_string( line );
		vector<string > names = parse_string(parseditem[0], '+');
		double rpkm = atof( parseditem[1].c_str() );
		for ( size_t i = 0; i < names.size(); ++i )
		{
			name_exp_map.insert( make_pair(names[i], rpkm) );
		}
	}
	inf.close();
	
}

void TranscriptExpression_GW::filter_out_T(Genome &genome ) // with same tss or not in annotation
{
	map< string, set<int> > tss_pos;
	vector< string > fltnames;
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		string name = ite->first;
		if ( genome.name_Transcript_map.find(name) == genome.name_Transcript_map.end() )
		{
			fltnames.push_back( name );
			continue;
		}
		string chr = genome.name_Transcript_map[name].chr;
		int tss = genome.name_Transcript_map[name].getTSS();
		if ( tss_pos[chr].find(tss) != tss_pos[chr].end() )
		{
			fltnames.push_back( name );
		} else
		{
			tss_pos[chr].insert( tss );
		}
	}
	
	for ( size_t i = 0; i < fltnames.size(); ++i )
	{
		name_exp_map.erase( fltnames[i] );
	}
}

void TranscriptExpression_GW::separate_into_active_silent_gene( vector<string > &active_T, vector<string > &silent_T, double cutoff )
{
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		string gene = ite->first;
		double rpkm = ite->second;
		if ( rpkm >= cutoff )
			active_T.push_back( gene );
		else
			silent_T.push_back( gene ); 
	}
}

void TranscriptExpression_GW::separate_into_quantil_by_decending_rpkm( vector<vector<string > > &qual_Ts )
{
	qual_Ts.clear();
	
	multimap< double, string, greater<double > > sc_gene_map; 
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		sc_gene_map.insert( make_pair(ite->second, ite->first ) );
	}
	
	int total_g_size = (int)sc_gene_map.size();
	int q = total_g_size / 4;
	vector< string > tgv;
	int i = 0;
	vector<vector<string > > p_qual_ve;
	for ( multimap< double, string, greater<double > >::iterator ite = sc_gene_map.begin(); ite != sc_gene_map.end(); ++ite )
	{
		i += 1;
		double sc = ite->first;
		string gene = ite->second;
		tgv.push_back( gene );
		if ( i % q == 0 )
		{
			p_qual_ve.push_back( tgv );
			tgv.clear();
		}
		
	}
	if ( !tgv.empty() )
		p_qual_ve.push_back( tgv );
	if ( (int)p_qual_ve.size() == 4 )
		qual_Ts = p_qual_ve;
	else if ( (int)p_qual_ve.size() == 5 )
	{
		for ( size_t j = 0; j < 3; ++j )
		{
			qual_Ts.push_back( p_qual_ve[j] ); 
		}
		vector<string > cbgv = p_qual_ve[3];
		cbgv.insert( cbgv.end(), p_qual_ve[4].begin(), p_qual_ve[4].end() );
		qual_Ts.push_back( cbgv );
	} else
	{
		cout<<"error in TranscriptExpression_GW::separate_into_quantil_by_decending_rpkm"<<endl;
		cout<<"Unexpected p_qual_ve size() "<<p_qual_ve.size()<<endl;
		exit(1);
	}
	
}

void TranscriptExpression_GW::separate_into_4part_by_rpkm( vector<string > &silent_T, vector<vector<string > > &three_Ts, double cutoff )
{
	multimap< double, string > sc_gene_map;
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		string gene = ite->first;
		double rpkm = ite->second;
		if ( rpkm >= cutoff )
		{
			sc_gene_map.insert( make_pair(ite->second, ite->first ) );
		} else
			silent_T.push_back( gene ); 
	}
	
	
	int total_g_size = (int)sc_gene_map.size();
	int q = total_g_size / 3;
	vector< string > tgv;
	int i = 0;
	vector<vector<string > > p_qual_ve;
	for ( multimap< double, string >::iterator ite = sc_gene_map.begin(); ite != sc_gene_map.end(); ++ite )
	{
		i += 1;
		double sc = ite->first;
		string gene = ite->second;
		tgv.push_back( gene );
		if ( i % q == 0 )
		{
			p_qual_ve.push_back( tgv );
			tgv.clear();
		}
		
	}
	if ( !tgv.empty() )
		p_qual_ve.push_back( tgv );
	if ( (int)p_qual_ve.size() == 3 )
		three_Ts = p_qual_ve;
	else if ( (int)p_qual_ve.size() == 4 )
	{
		for ( size_t j = 0; j < 2; ++j )
		{
			three_Ts.push_back( p_qual_ve[j] ); 
		}
		vector<string > cbgv = p_qual_ve[2];
		cbgv.insert( cbgv.end(), p_qual_ve[3].begin(), p_qual_ve[3].end() );
		three_Ts.push_back( cbgv );
	} else
	{
		cout<<"error in TranscriptExpression_GW::separate_into_quantil_by_decending_rpkm"<<endl;
		cout<<"Unexpected p_qual_ve size() "<<p_qual_ve.size()<<endl;
		exit(1);
	}
	cout<<silent_T.size()<<" "<<three_Ts[0].size()<<" "<<three_Ts[1].size()<<" "<<three_Ts[2].size()<<endl;
	
	
}

void TranscriptExpression_GW::separate_into_4part_by_breaks( vector<vector<string > > &four_Ts, vector<double> &breaks )
{
	four_Ts.clear();
	for ( size_t i = 0; i < 4; ++i )
	{
		vector<string > nl;
		four_Ts.push_back( nl );
	}
	
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		if ( ite->second < breaks[0] )
			four_Ts[0].push_back( ite->first );
		else if ( ite->second < breaks[1] )
			four_Ts[1].push_back( ite->first );
		else if ( ite->second < breaks[2] )
			four_Ts[2].push_back( ite->first );
		else
			four_Ts[3].push_back( ite->first );
	}
	
}

void TranscriptExpression_GW::sortgene(vector<string > &sortedgene, double cutoff )
{
	multimap< double, string, greater<double> > sc_gene_map;
	for ( map<string, double >::iterator ite = name_exp_map.begin(); ite != name_exp_map.end(); ++ite )
	{
		string gene = ite->first;
		double rpkm = ite->second;
		if ( rpkm > cutoff )
		{
			sc_gene_map.insert( make_pair(ite->second, ite->first ) );
		} 
	}
	
	for ( multimap< double, string, greater<double> >::iterator ite = sc_gene_map.begin(); ite != sc_gene_map.end(); ++ite )
	{
		sortedgene.push_back( ite->second );
	}
	
	
}

void Sequence_Di::genecoding( )
{
	coding.insert( make_pair('a', 0));
	coding.insert( make_pair('A', 0));
	
	coding.insert( make_pair('T', 1));
	coding.insert( make_pair('t', 1));
	
	coding.insert( make_pair('c', 2));
	coding.insert( make_pair('C', 2));
	
	coding.insert( make_pair('g', 3));
	coding.insert( make_pair('G', 3));
	
	din_revcode.insert( make_pair(0, "AA") );
	din_revcode.insert( make_pair(1, "AT") );
	din_revcode.insert( make_pair(2, "AC") );
	din_revcode.insert( make_pair(3, "AG") );
	din_revcode.insert( make_pair(4, "TA") );
	din_revcode.insert( make_pair(5, "TT") );
	din_revcode.insert( make_pair(6, "TC") );
	din_revcode.insert( make_pair(7, "TG") );
	din_revcode.insert( make_pair(8, "CA") );
	din_revcode.insert( make_pair(9, "CT") );
	din_revcode.insert( make_pair(10, "CC") );
	din_revcode.insert( make_pair(11, "CG") );
	din_revcode.insert( make_pair(12, "GA") );
	din_revcode.insert( make_pair(13, "GT") );
	din_revcode.insert( make_pair(14, "GC") );
	din_revcode.insert( make_pair(15, "GG") );
}

int Sequence_Di::getdicode( char n1, char n2 )
{
	
	if ( n1 != 'a' && n1 != 'A' && n1 != 't' && n1 != 'T' && n1 != 'c' && n1 != 'C' && n1 != 'g' && n1 != 'G' )
	{
	//	cout<<"unexpect n "<<n1<<endl;
		return -1;
	}
	if ( n2 != 'a' && n2 != 'A' && n2 != 't' && n2 != 'T' && n2 != 'c' && n2 != 'C' && n2 != 'g' && n2 != 'G' )
	{	
	//	cout<<"unexpect n "<<n2<<endl;
		return -1;
	}
	
	int sc = coding[n1]*4+coding[n2];
	return sc;
}

vector<int > Sequence_Di::getdincode_fromseq( string seq )
{
	vector<int > c_ve;
	for ( size_t i = 0; i < seq.size() - 1; ++i )
	{
		char n1 = seq[i];
		char n2 = seq[i+1];
		int c = getdicode( n1, n2 );
		c_ve.push_back( c );
	}
	return c_ve;
}

vector<int > Sequence_Di::getdincode_aroundnucl( string chr, int pos, Genome &gm )
{
	int upstream = 200;
	int downstream = 200;
	string seq = gm.getsubseq( chr, pos-upstream, upstream+downstream );
	
	return getdincode_fromseq( seq );
}

double Sequence_Di::cal_AATTfreq_Nuclfrank( string chr, int pos, Genome &gm )
{

	int upstream = 200;
	int downstream = 200;
	int frank_ls = 20;
	int frank_le = 125;
	int frank_rs = 275;
	int frank_re = 381;
	
	string seq = gm.getsubseq( chr, pos-upstream, upstream+downstream );

	int aatt = 0;
	int total = 0;
	for ( size_t i = 0; i < seq.size() -1; ++i )
	{
		if ( ( i >= frank_ls && i <= frank_le ) || ( i >= frank_rs && i <= frank_re ) )
		{
			char n1 = seq[i];
			char n2 = seq[i+1];
			if ( coding[n1] == 0 && coding[n2] == 0 )
				aatt += 1;
			else if ( coding[n1] == 1 && coding[n2] == 1 )
				aatt += 1;
		
			total += 1;
		}
	}
	
	double fr = aatt*1.0/total;
	return fr;
	
}

double Sequence_Di::cal_AATTATTAfreq_Nuclfrank( string chr, int pos, Genome &gm )
{

	int upstream = 200;
	int downstream = 200;
	int frank_ls = 20;
	int frank_le = 125;
	int frank_rs = 275;
	int frank_re = 381;
	
	string seq = gm.getsubseq( chr, pos-upstream, upstream+downstream );

	int aatt = 0;
	int total = 0;
	for ( size_t i = 0; i < seq.size() -1; ++i )
	{
		if ( ( i >= frank_ls && i <= frank_le ) || ( i >= frank_rs && i <= frank_re ) )
		{
			char n1 = seq[i];
			char n2 = seq[i+1];
			if ( coding[n1] == 0 && coding[n2] == 0 )
				aatt += 1;
			else if ( coding[n1] == 1 && coding[n2] == 1 )
				aatt += 1;
			else if ( coding[n1] == 0 && coding[n2] == 1 )
				aatt += 1;
			else if ( coding[n1] == 1 && coding[n2] == 0 )
				aatt += 1;
				
			total += 1;
		}
	}
	
	double fr = aatt*1.0/total;
	return fr;
	
}

void Sequence_Di::cal_din_freq_from_dincode_ve( vector< vector<int > > &dincode_ve,
	map<int, vector<double > > &din_rate )
{
	map<int, vector<int > > din_freq;
	int col_n = (int)dincode_ve[0].size();
	for ( int i = 0; i < 16; ++i )
	{
		vector<int > nuls;
		for ( int j = 0; j < col_n; ++j )
			nuls.push_back(0);
		din_freq.insert( make_pair(i, nuls) );
	}
	
	for ( size_t i = 0; i < dincode_ve.size(); ++i )
	{
		for ( size_t j = 0; j < dincode_ve[i].size(); ++j )
		{
			int c = dincode_ve[i][j];
			if ( c == -1 )
				continue;
			if ( din_freq.find(c) == din_freq.end() )
			{
				cout<<"error din_freq cannot find c "<<c<<endl;
				exit(1);
			}
			if ( j >= col_n )
			{
				cout<<"error j >= col_n "<<j<<" "<<col_n<<endl;
				exit(1);
			}
			din_freq[c][j] += 1;
		}
	}
	
	int total_seq_n = (int)dincode_ve.size();
	double exp_freq = total_seq_n*1.0/16;
	
	for ( map<int, vector<int > >::iterator ite = din_freq.begin(); ite != din_freq.end(); ++ite )
	{
		vector<double > freq_ve;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			double fr = ite->second[i]*1.0 / exp_freq;
			freq_ve.push_back(fr);
		}
		din_rate[ite->first] = freq_ve;
	}
}

void Sequence_Di::cal_ave_freq( map<string, set<int > > &nucl_sets, 
	Genome &gm,
	map<int, vector<double > > &din_rate )
{
	vector< vector<int > > dincode_ve;
	for ( map<string, set<int > >::iterator ite = nucl_sets.begin(); ite != nucl_sets.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int pos = *si;

			vector<int > dinc = getdincode_aroundnucl(chr, pos, gm );
			dincode_ve.push_back( dinc );
		}
	}
	cout<<"cal freq from dincode "<<endl;
	cal_din_freq_from_dincode_ve( dincode_ve, din_rate );
	
}

void Sequence_Di::cal_ave_freq( map<string, map<int, char > > &nucl_dir_sets, 
	Genome &gm,
	map<int, vector<double > > &din_rate )
{
	vector< vector<int > > dincode_ve;
	for ( map<string, map<int, char > >::iterator ite = nucl_dir_sets.begin(); ite != nucl_dir_sets.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, char>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int pos = si->first;

			vector<int > dinc = getdincode_aroundnucl(chr, pos, gm );
			if ( si->second == '-' )
			{
				reverse( dinc.begin(), dinc.end() );
			}
			dincode_ve.push_back( dinc );
		}
	}
	cout<<"cal freq from dincode "<<endl;
	cal_din_freq_from_dincode_ve( dincode_ve, din_rate );
	
}

double Sequence_Di::cal_nucl_pref_score( string chr, int pos, Genome &gm )
{
	
	double b1 = 3;
	double p1 = -3;
	double b2 = 1;
	double p2 = -1;
	double sc = 0;
	
	int upstream = 200;
	int downstream = 200;
	int frank_ls = 100;
	int frank_le = 125;
	int frank_rs = 275;
	int frank_re = 300;
	
	string seq = gm.getsubseq( chr, pos-upstream, upstream+downstream );

	for ( size_t i = 0; i < seq.size() -1; ++i )
	{
		char n1 = seq[i];
		char n2 = seq[i+1];
		if ( ( i >= frank_ls && i <= frank_le ) || ( i >= frank_rs && i <= frank_re ) )
		{
			
			if ( coding[n1] == 0 && coding[n2] == 0 )
				sc += b1;
			else if ( coding[n1] == 1 && coding[n2] == 1 )
				sc += b1;
			else if ( coding[n1] == 0 && coding[n2] == 1 )
				sc += b1;
			else if ( coding[n1] == 1 && coding[n2] == 0 )
				sc += b1;
			else if ( coding[n1] == 2 && coding[n2] == 2 )
				sc += p1;
			else if ( coding[n1] == 3 && coding[n2] == 3 )
				sc += p1;
			else if ( coding[n1] == 2 && coding[n2] == 3 )
				sc += p1;
			else if ( coding[n1] == 3 && coding[n2] == 2 )
				sc += p1;
					
			
		}
		 else if ( i > frank_le && i < frank_rs )
		{
			if ( coding[n1] == 0 && coding[n2] == 0 )
				sc += p2;
			else if ( coding[n1] == 1 && coding[n2] == 1 )
				sc += p2;
			else if ( coding[n1] == 0 && coding[n2] == 1 )
				sc += p2;
			else if ( coding[n1] == 1 && coding[n2] == 0 )
				sc += p2;
			else if ( coding[n1] == 2 && coding[n2] == 2 )
				sc += b2;
			else if ( coding[n1] == 3 && coding[n2] == 3 )
				sc += b2;
			else if ( coding[n1] == 2 && coding[n2] == 3 )
				sc += b2;
			else if ( coding[n1] == 3 && coding[n2] == 2 )
				sc += b2;
		}  
	}
	
	return sc;
	
}

void Sequence_Di::cal_nucl_pref_score_region( Genome &gm, string chr, int start, int end, vector<double > &sc_ve )
{
	vector<double > sc_before;
	double agg = 0;
	double n = 0;
	for ( int pos = start; pos <= end; ++pos )
	{
		double sc = cal_nucl_pref_score( chr, pos, gm );
		sc_before.push_back( sc );
		agg += sc;
		n+=1;
	}
	double ave = agg / n;
	for ( size_t i = 0; i < sc_before.size(); ++i )
	{
		double sc = sc_before[i] - ave;
		sc_ve.push_back( sc );
	}
}

double Sequence_Di::cal_nucl_rotation_score( string chr, int pos, Genome &gm )
{
	double b1 = 1;
	int upstream = 73;
	int downstream = 73;
	string seq = gm.getsubseq( chr, pos-upstream, upstream+downstream );
	double sc = 0;
	for ( int i = 0; i < 14; ++i )
	{
		
		int j = (int)(i*10.3)+5;
		char n1 = seq[j];
		char n2 = seq[j+1];
		char n3 = seq[j+2];
		
		if ( ( coding[n1] == 0 && coding[n2] == 0 ) || (coding[n2] == 0 && coding[n3] == 0) )
			sc += b1;
		else if ( (coding[n1] == 1 && coding[n2] == 1) || (coding[n2] == 1 && coding[n3] == 1) )
			sc += b1;
		else if ( (coding[n1] == 0 && coding[n2] == 1) || (coding[n2] == 0 && coding[n3] == 1) )
			sc += b1;
		else if ( (coding[n1] == 1 && coding[n2] == 0) || (coding[n2] == 1 && coding[n3] == 0) )
			sc += b1;
		
	}
	
	return sc;
}

pair<int,double > Sequence_Di::realign_nucl_by_rotation_score( string chr, int pos, Genome &gm )
{
	map<double, int, greater<double > > score_pos;
	for ( int i = -5; i <= 5; ++i )
	{
		int rpos = pos + i;
		double sc = cal_nucl_rotation_score( chr, rpos, gm );
		score_pos.insert( make_pair(sc, rpos ) );
	/*	cout<<i<<" "<<sc<<endl;
		if ( rpos == pos )
		{
			cout<<"score at pos "<<pos<<" = "<<sc<<endl;
		}  */
	}
	
	int maxpos = score_pos.begin()->second;
//	cout<<"maxpos = "<<maxpos<<" maxsc = "<<score_pos.begin()->first<<endl;
	return make_pair( maxpos, score_pos.begin()->first);
}

void SNP_bank::readinsnp( string infile )
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
		vector<string > ps = parse_string( line );
		string chr = ps[1];
		int site = atoi(ps[2].c_str());
		string base = ps[9];
		if ( (int)base.size() < 3 )
		{
			cout<<"unexpected base in snp file "<<base<<endl;
			exit(1);
		}  
		vector<string> base_ve = parse_string( base, '/');
		vector<char> sve;
		for ( size_t i = 0; i < base_ve.size(); ++i )
		{
			char n = base_ve[i][0];
			sve.push_back( n );
		}
		
		chr_site_snp[chr][site] = sve; 
	}
}

void SNP_bank::readinsnp_2( string infile )
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
		vector<string > ps = parse_string( line );
		string chr = ps[1];
		int site = atoi(ps[2].c_str());
		if ( (int)ps[8].size() != 1 )
			continue;
		char refbase = ps[8][0];
		if ( ps[8][0] != ps[7][0] )
		{
			cout<<"unexpected ref base "<<ps[8][0]<<" "<<ps[7][0]<<endl;
			exit(1);
		}
		string base = ps[9];
		if ( (int)base.size() < 3 )
		{
			cout<<"unexpected base in snp file "<<base<<endl;
			exit(1);
		}  
		vector<string> base_ve = parse_string( base, '/');
		vector<char> sve;
		sve.push_back(refbase );
		for ( size_t i = 0; i < base_ve.size(); ++i )
		{
			char n = base_ve[i][0];
			if ( n != refbase ) 	
				sve.push_back( n );
		}
		
		chr_site_snp[chr][site] = sve; 
	}
}

void SNP_bank::readaltsnp( string infile )
{
	map<string, map<int, vector<char> > > alt_chr_site_snp; 
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
		vector<string > ps = parse_string( line );
		string chr = ps[0];
		int site = atoi(ps[1].c_str());
		if ( chr_site_snp[chr].find( site ) == chr_site_snp[chr].end( ) )
		{
			cout<<"error can not find snp "<<chr<<"\t"<<site<<endl; 
			exit(1);
		}
		alt_chr_site_snp[chr][site] = chr_site_snp[chr][site];
	}
	inf.close();
	chr_site_snp = alt_chr_site_snp;
}

void SNP_bank::generate_index()
{
	int L = 100000;
	for ( map<string, map<int, vector<char > > >::iterator ite = chr_site_snp.begin();
		ite != chr_site_snp.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, vector<char > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			int site = si->first;
			int index = site / L;
			chr_index_sites[chr][index].insert(site);
		}
	}
}

void SNP_bank::get_snp_withinregion( vector<int> &snps, string chr, int start, int end )
{
	int L = 100000;
	int c = start + (end-start)/2;
	int index = c / L;
	if ( chr_index_sites[chr].find(index) != chr_index_sites[chr].end() )
	{
		for ( set<int >::iterator ite = chr_index_sites[chr][index].begin(); 
			ite != chr_index_sites[chr][index].end(); ++ite )
		{
			if ( *ite < start )
				continue;
			if ( *ite > end )
				break;
			snps.push_back( *ite );
		}
	}
}

void SNP_bank::getsnp_atnucl( string chr, int c, vector<int > &pos )
{
	int index = c / 100000;
	int start = c - 200;
	int end = c + 200;
	if ( chr_index_sites[chr].find( index ) != chr_index_sites[chr].end() )
	{
		for ( set<int >::iterator ite = chr_index_sites[chr][index].begin();
			ite != chr_index_sites[chr][index].end(); ++ite )
		{
			if ( *ite < start )
				continue;
			if ( *ite > end )
				break;
			int rlp = *ite - c;
			pos.push_back( rlp );
		}
	}
}

void SNP_bank::getsnp_atnucl_set( map<string, set<int > > &nucl_sets, vector<double > &freq )
{
	map<int, int > pos_count;
	int nucl_num = 0;
	for ( map<string, set<int > >::iterator ite = nucl_sets.begin(); ite != nucl_sets.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			vector<int > pos;
			getsnp_atnucl( chr, *si, pos );
			nucl_num += 1;
			if ( !pos.empty() )
			{
				for ( size_t i = 0; i < pos.size(); ++i )
				{
					if ( pos_count.find(pos[i] ) == pos_count.end() )
					{
						pos_count[pos[i]] = 1;
					} else
					{
						pos_count[pos[i]] += 1;
					}
				}
			}
			if ( nucl_num >= 50000000 )
				break;
		}
		if ( nucl_num >= 50000000 )
			break;
	}
	
	for ( int i = -200; i <= 200; ++i )
	{
		int n = 0;
		if ( pos_count.find( i ) != pos_count.end() )
		{
			n = pos_count[i];
		}
		double fr = n*1000.0/nucl_num;
		freq.push_back( fr );
	}
	
}



void Conserv_bank::readinphastConsE( string infile )
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
		vector<string > ps = parse_string( line );
		string chr = ps[1];
		int start = atoi(ps[2].c_str() );
		int end = atoi(ps[3].c_str() );
		int sc = atoi(ps[5].c_str() );
		phastConsE[chr][make_pair(start, end)] = sc;
	}
	inf.close();
}

void Conserv_bank::readinphyloP( string infile )
{
	ifstream inf( infile.data() );
	if ( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read in "<<infile<<endl;
	string line;
	getline( inf, line );
	vector<string > parsed_items = parse_string(line);
	while(!inf.eof())
	{
		
		if ( line.empty() )
			break;
		vector<string> chrom_ve = parse_string( parsed_items[1], '=' );
		string chr = chrom_ve[1];
		vector<string> start_ve = parse_string( parsed_items[2], '=' );
		int start = atoi( start_ve[1].c_str() );
		
		vector<double > site_sc;
		getline( inf, line );
		parsed_items = parse_string(line);
		
		while ( (int)parsed_items.size() == 1 )
		{
			double sc = atof( parsed_items[0].c_str() );
			site_sc.push_back(sc);
			getline( inf, line );
			if ( line.empty() )
				break;
			if ( inf.eof() )
				break;
			parsed_items = parse_string(line);
			
		}
		
		
		if ( !site_sc.empty() )
			phyloP[chr][start]=site_sc;
		
		if ( line.empty() )
			break;
		if ( inf.eof() )
			break;
		if ( (int)parsed_items.size() != 4 )
		{
			cout<<"Unexpected line "<<line<<endl;
			exit(1);
		}
	}
	inf.close();
	
	
	
}

void Conserv_bank::readinphyloP_bunch( string infile )
{
	ifstream inf( infile.data() );
	if ( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	string line;
	
	while(!inf.eof())
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		
		readinphyloP( line );
		
	}
	
	inf.close();
		
}

void Conserv_bank::generate_index()
{
	int L = 100000;
	for ( map<string, map<int, vector<double > > >::iterator ite = phyloP.begin(); ite != phyloP.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, vector<double > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int start = si->first;
			int l = (int)si->second.size();
			int index1 = start / L;
			int index2 = (start+l) / L;
			for ( int index = index1; index <= index2; ++index )
			{
				chr_index_phyloPstart[chr][index].insert( start );
			}
		}
	}
}

bool Conserv_bank::getvalid_phyloPS( string chr, int c, vector<double > &sc_ve )
{
	bool fd = false;
	if ( phyloP.find( chr ) == phyloP.end() )
		return fd;
	
	int L = 100000;
	int index = c / L;
	if ( chr_index_phyloPstart[chr].find(index ) == chr_index_phyloPstart[chr].end() )
		return fd;
	
	for ( set<int>::iterator ite = chr_index_phyloPstart[chr][index].begin();
		ite != chr_index_phyloPstart[chr][index].end(); ++ite )
	{
		int start = *ite;
		int end = start + (int)phyloP[chr][start].size();
		if ( end <= c + 201 )
			continue;
		if ( start > c - 201 )
			break;
		
		int i = c - 200 - start;
		for ( ; i < c+200-start; ++i )
		{
			sc_ve.push_back( phyloP[chr][start][i] );
		}
		fd = true;
	}
	
	return fd;
}

bool Conserv_bank::get_core_flank_phyloPS( string chr, int c, double & sc_flank, double & sc_core )
{
	int upstream = 200;
	int downstream = 200;
	int frank_ls = 20;
	int frank_le = 125;
	int frank_rs = 275;
	int frank_re = 381;
	
	bool fd = false;
	vector< double > sc_ve;
	fd = getvalid_phyloPS( chr, c, sc_ve );
	
	if ( !fd )
		return fd;
	
	double flank_sum = 0;
	int flank_num = 0;
	for ( int i = frank_ls; i <= frank_le; ++i )
	{
		flank_sum += sc_ve[i-1];
		flank_num += 1;
	}
	for ( int i = frank_rs; i <= frank_re; ++i )
	{	
		flank_sum += sc_ve[i-1];
		flank_num += 1;
	}
	sc_flank = flank_sum / flank_num;
	double core_sum = 0;
	int core_num = 0;
	for ( int i = frank_le+1; i<= frank_rs-1; ++i )
	{
		core_sum += sc_ve[i-1];
		core_num += 1;
	}
	sc_core = core_sum / core_num;
	
	return fd;
}

void Conserv_bank::cal_phyloPS_nuclset( map<string, set<int > > &nucl_sets, vector<double > &ave_score )
{
	vector<vector<double > > mat;
	for ( map<string, set<int > >::iterator ite = nucl_sets.begin(); ite != nucl_sets.end(); ++ite )
	{
		string chr = ite->first;
		for ( set< int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			vector<double > sc_ve;
			bool fd = getvalid_phyloPS( chr, *si, sc_ve );
			if ( fd )
			{
				if ( mat.size() < 5000000 )
					mat.push_back( sc_ve );
					 
			}
				
		}
	}
	cout<<"valid nucl size "<<mat.size()<<endl;
	if ( !mat.empty() )
	{
		cal_column_average( mat, ave_score );
	}
}















