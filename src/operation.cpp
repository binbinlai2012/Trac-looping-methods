#include "operation.h"

void Position::addpos(int inst, int inend, string inchr)
{
	start = inst;
	end = inend;
	chr = inchr;
}

void Position::addpos(string str)
{
	vector<string > splstr1 = parse_string(str, ':');
	if ( splstr1.size() != 2 )
	{
		cout<<"error parse position: str "<< str<<endl;
		exit(1);
	}
	chr = splstr1[0];
	vector<string > splstr2 = parse_string(splstr1[1], '+');
	if ( splstr2.size() != 2 )
	{
		cout<<"error parse position: str "<< str<<endl;
		exit(1);
	}
	start = atoi( splstr2[0].c_str() );
	end = atoi( splstr2[1].c_str() );
}


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

double log_2(double r )
{
	double res = log(r)/log(2.0);
	return res;
}

void callendensity( vector<int > &len_ve, map<int, double > &len_den, int bin )
{
	int total_n = (int)len_ve.size();
	map<int, int > len_count;
	for ( size_t i = 0; i < len_ve.size(); ++i )
	{
		int tl = len_ve[i];
		int b = tl / bin;
		int r = b * bin;
		if ( len_count.find(r) == len_count.end() )
		{
			len_count[r] = 1;
		} else
		{
			len_count[r] += 1;
		}
	}
	
	for ( map<int, int >::iterator ite = len_count.begin(); ite != len_count.end(); ++ite )
	{
		double c = ite->second * 1.0 / total_n;
		len_den.insert(make_pair(ite->first, c) );
	}
}

bool overlaptest(Position &p1, Position &p2)
{
	if ( p1.chr != p2.chr )
		return false;

	if ( p1.start > p2.end || p1.end < p2.start )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, pair<int, int> p2)
{
	if ( p1.first > p2.second || p1.second < p2.first )
		return false;
	else
		return true;
}

bool overlaptest(Position &p1, Position &p2, int extension)
{
	if ( p1.chr != p2.chr )
		return false;

	if ( p1.start > (p2.end+extension) || p1.end < (p2.start-extension) )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, pair<int, int> p2, int extension)
{
	if ( p1.first > (p2.second+extension) || p1.second < (p2.first-extension) )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first )
			continue;
		else if ( ite->first.first > p1.second )
			break;
		else
			ovl = true;
			break;
		
	}
	return ovl;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first )
			continue;
		else if ( ite->first.first > p1.second )
			break;
		else
			ovl = true;
			idve.push_back( ite->second );
			
	}
	return ovl;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve, int extension )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first - extension )
			continue;
		else if ( ite->first.first > p1.second + extension )
			break;
		else
			ovl = true;
			idve.push_back( ite->second );
			
	}
	return ovl;
}

bool overlaptest( int site1, map<int, Region_id > &posmap, int extension, vector<int > &dis )
{
	bool ovl = false;
	for ( map<int, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first < site1 - extension )
			continue;
		else if ( ite->first > site1 + extension )
			break;
		else
		{
			ovl = true;
			dis.push_back( ite->first - site1);
		}
	}
	return ovl;
}


int matrixtrans( int stage, int elem, vector<int > &matrix )
{
	int n = 0;
	if ((int)matrix.size() != stage )
	{
		cout<<"error in marixtrans: "<<matrix.size()<<", "<<stage<<endl;
		exit(1);
	}
	for ( size_t i = 0; i < matrix.size(); ++i )
	{
		if ( matrix[i] >= elem )
		{
			cout<<"error in marixtrans: "<<matrix[i]<<", "<<elem<<endl;
			exit(1);
		}
		n += matrix[i] * (int)pow((double)elem, (double)i );
	}
	return n;
}

vector<int > revtransmatrix( int stage, int elem, int digit )
{
	vector<int> matrix;
	for ( int i = 0; i < stage; ++i )
		matrix.push_back(0);
	int maxn = (int)pow((double)elem, (double)stage );
	if ( digit >= maxn )
	{
		cout<<"error in revtransmatrix "<<stage<<","<<elem<<","<<digit<<endl;
		exit(1);
	}
	for (int i = 0; i < stage; ++i )
	{
		matrix[i] = digit % elem;
		digit -= matrix[i];
		digit /= elem;
	}

	return matrix;
}

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, string > > &region_map )
{
	region.clear();
	for ( map<string, map<pair<int, int>, string > >::iterator ite = region_map.begin(); ite != region_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int>, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			region[chr].push_back(seci->first);
		}
	}
}

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, Region_id > > &region_map )
{
	region.clear();
	for ( map<string, map<pair<int, int>, Region_id > >::iterator ite = region_map.begin(); ite != region_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int>, Region_id >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			region[chr].push_back(seci->first);
		}
	}
}

void filter_region_remain( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	filtered.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				continue;
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( overlap )
					filtered[chr].push_back(*seci);
			}
		}
	}
}

void filter_region_minus( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	filtered.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				filtered[chr].push_back( *seci );
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( !overlap )
					filtered[chr].push_back(*seci);
			}
		}
	}
}

void filter_region_minus_remain( map<string, vector< pair<int, int> > > &overlapped, 
				   map<string, vector< pair<int, int> > > &nonoverlapped,
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	overlapped.clear();
	nonoverlapped.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				nonoverlapped[chr].push_back( *seci );
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( !overlap )
					nonoverlapped[chr].push_back(*seci);
				else
					overlapped[chr].push_back(*seci);
			}
		}
	}
}


void region_merge_naive( vector< map<string, vector< pair<int, int> > > > &reg_ve,
	map<string, vector< pair<int, int> > > &merged )
{
	map<string, set<pair<int, int > > > chr_posset;
	for ( size_t i = 0; i < reg_ve.size(); ++i )
	{
		for ( map<string, vector<pair<int, int> > >::iterator ite = reg_ve[i].begin(); ite != reg_ve[i].end(); ++ite )
		{	
			string chr = ite->first;
			for ( size_t k = 0; k < ite->second.size(); ++k )
			{
				chr_posset[chr].insert(ite->second[k]);
			}
		}
	}
	
	for ( map<string, set<pair<int, int > > >::iterator ite = chr_posset.begin(); ite != chr_posset.end(); ++ite )
	{
		vector< pair<int, int > > m;
		m.push_back( *(ite->second.begin()) );
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( si->first <= m.back().second )
			{
				if ( si->second > m.back().second )
					m.back().second = si->second;
			} else
			{
				m.push_back( *si );
			}
		}
		merged.insert( make_pair(ite->first, m ) );
	}
}

void stitch_region( map<string, vector<vector<pair<int, int> > > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance )
{
	stitched.clear();
	for ( map<string, vector<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
	{
		string chr = ite->first;
		size_t i = 0;
		vector<pair<int, int > > clu;
		clu.push_back(ite->second[i]);
		i++;
		for ( ; i < ite->second.size(); ++i )
		{
			if ( ite->second[i].first - ite->second[i-1].second < distance )
			{
				clu.push_back( ite->second[i]);
			} else
			{
				stitched[chr].push_back( clu );
				clu.clear();
				clu.push_back( ite->second[i] );
			}
		}
		stitched[chr].push_back( clu );
	}

}

void transformregion( map<string, vector< pair<int, int> > > &newregion,
	map<string, set<pair<int, int > > > &oldregion )
{
	newregion.clear();
	for (map<string, set<pair<int, int > > >::iterator ite = oldregion.begin(); ite != oldregion.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			newregion[chr].push_back(*si);
		} 
	}
}

string inttostr(int i )
{
	string Res;
	ostringstream convert;
	convert << i;
	Res = convert.str();
	return Res;
}

void outputtable( ofstream &outf, vector<vector<int> > &mat, vector<string > &rowname, vector<string> &colname )
{
	
	for ( size_t i = 0; i < colname.size(); ++i )
	{
		outf<<"\t"<<colname[i];
	}
	outf<<endl;
	if ( rowname.size() != mat.size() )
	{
		cout<<"error output table rowname size != mat size "<<rowname.size()<<","<<mat.size()<<endl; exit(1);
	}
	for ( size_t i = 0; i < rowname.size(); ++i )
	{
		outf<<rowname[i];
		if ( mat[i].size() != colname.size() )
		{
			cout<<"error output table mat[i] size != colname size "<<colname.size()<<", "<<mat[i].size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j) 
		{
			outf<<"\t"<<mat[i][j];
		}
		outf<<endl;
	}
}

void outputtable( ofstream &outf, vector<vector<double> > &mat, vector<string > &rowname, vector<string> &colname )
{
	
	for ( size_t i = 0; i < colname.size(); ++i )
	{
		outf<<"\t"<<colname[i];
	}
	outf<<endl;
	if ( rowname.size() != mat.size() )
	{
		cout<<"error output table rowname size != mat size "<<rowname.size()<<","<<mat.size()<<endl; exit(1);
	}
	for ( size_t i = 0; i < rowname.size(); ++i )
	{
		outf<<rowname[i];
		if ( mat[i].size() != colname.size() )
		{
			cout<<"error output table mat[i] size != colname size "<<colname.size()<<", "<<mat[i].size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j) 
		{
			outf<<"\t"<<setprecision(3)<<mat[i][j];
		}
		outf<<endl;
	}
}

vector<vector<double> > transposmat( vector<vector<double> > &mat )
{
	if ( mat.empty() )
	{
		cout<<"error mat empty transposmat"<<endl; exit(1);
	}
	size_t il = mat.size();
	size_t jl = mat[0].size(); 
	if ( jl == 0 )
	{
		cout<<"error mat[0] empty transposmat"<<endl; exit(1);
	}
	for ( size_t i = 0; i < il; ++i )
	{
		if ( mat[i].size() != jl )
		{	
			cout<<"error mat i "<<i<<" != mat[0] "<<jl<<endl; exit(1);
		}
		
	}
	vector<vector<double > > res;
	for ( size_t j = 0; j < jl; ++j )
	{
		vector< double > subv;
		for ( size_t i = 0; i < il; ++i )
		{
			subv.push_back( mat[i][j] );
		}
		res.push_back( subv );
	}
	return res;
}

void normalizemat( vector<vector<double > > &mat, vector<double > &thr )
{
	
	for ( size_t i = 0; i < mat.size(); ++i )
	{
		if ( mat[i].size() != thr.size() )
		{
			cout<<"error in normalizemat: mat[i]size != thr size "<<mat[i].size()<<","<<thr.size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j )
		{
			if ( thr[j] > 0 )
			{
				mat[i][j] /= thr[j];
				if ( mat[i][j] > 1 )
					mat[i][j] = 1;
			}
		}
	}
}

int findnearTSSinchr( set<int > &tssset, pair<int, int > region )
{
	int distance = -1;
	int selectTSS = -1;
	for ( set<int>::iterator ite = tssset.begin(); ite != tssset.end(); ++ite )
	{
		int td = min(abs(region.first - *ite), abs(region.second - *ite) );
		if ( distance == -1 )
		{
			selectTSS = *ite;
			distance = td;
		} else
		{
			if ( td < distance )
			{
				distance = td;
				selectTSS = *ite;
			}
		}
		if ( *ite < region.first )
			continue;
		if ( *ite > region.second )
			break;
		
	}
	return selectTSS;
}

int findnearTSSinchr( set<int > &tssset, pair<int, int > region, vector< int > &containedTSS )
{
	int distance = -1;
	int selectTSS = -1;
	for ( set<int>::iterator ite = tssset.begin(); ite != tssset.end(); ++ite )
	{
		int td = min(abs(region.first - *ite), abs(region.second - *ite) );
		if ( distance == -1 )
		{
			selectTSS = *ite;
			distance = td;
		} else
		{
			if ( td < distance )
			{
				distance = td;
				selectTSS = *ite;
			}
		}
		if ( *ite < region.first )
			continue;
		if ( *ite > region.second )
			break;
		
		containedTSS.push_back( *ite );
	}
	return selectTSS;
}


set<Region_id > vetoset(vector<Region_id > &ve )
{
	set<Region_id > rs;
	for ( size_t i = 0; i < ve.size(); ++i )
		rs.insert( ve[i]);
	return rs;
}

vector<pair<int, int > > mergeragions( set<pair<int, int > > &r )
{
	vector<pair<int, int > > mr;
	if ( r.empty() )
		return mr;
	
	mr.push_back( *r.begin() );
	for ( set<pair<int, int > >::iterator ite = r.begin(); ite != r.end(); ++ite )
	{
		if ( mr.back().second >= ite->first )
		{
			mr.back().second = max( ite->second, mr.back().second );
		} else
		{
			mr.push_back( *ite );
		}
	}
	return mr;
}

int shift_tagpos( int start, int end, char strand, int seg_len )
{
	if ( seg_len == 0 )
	{
		return start + (end-start)/2;
	}
	int p = 0;
	if ( strand == '+' )
	{
		return start + seg_len/2;
	} else if ( strand == '-')
		return end - 1 - seg_len/2;
	else
	{
		cout<<"error strand "<<strand<<endl; 
		exit(1);
		return 0;
	}
}

void generate_wig( string outfile, map<string, map<int, int > > &chr_pos_sc, string name, int span )
{
	ofstream outf(outfile.data() );
	
	outf<<"track type=wiggle_0 name="+name<<endl;
	for ( map<string, map<int, int > >::iterator ite = chr_pos_sc.begin(); ite != chr_pos_sc.end(); ++ite )
	{
		string chr = ite->first;
		outf<<"variableStep chrom="+chr+" span="<<span<<endl;
		for ( map<int, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outf<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outf.close();
}

void generate_wig( string outfile, map<string, map<int, int > > &chr_pos_sc, string name, int span, string tchr )
{
	ofstream outf(outfile.data() );
	
	outf<<"track type=wiggle_0 name="+name<<endl;
	for ( map<string, map<int, int > >::iterator ite = chr_pos_sc.begin(); ite != chr_pos_sc.end(); ++ite )
	{
		string chr = ite->first;
		if (chr != tchr )
			continue;
		outf<<"variableStep chrom="+chr+" span="<<span<<endl;
		for ( map<int, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outf<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outf.close();
}

void generate_wig( string outfile, map<string, map<int, double > > &chr_pos_sc, string name, int span )
{
	ofstream outf(outfile.data() );
	
	outf<<"track type=wiggle_0 name="+name<<endl;
	for ( map<string, map<int, double > >::iterator ite = chr_pos_sc.begin(); ite != chr_pos_sc.end(); ++ite )
	{
		string chr = ite->first;
		outf<<"variableStep chrom="+chr+" span="<<span<<endl;
		for ( map<int, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outf<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outf.close();
}

void generate_wig( string outfile, map<string, map<int, double > > &chr_pos_sc, string name, int span, string tchr )
{
	ofstream outf(outfile.data() );
	
	outf<<"track type=wiggle_0 name="+name<<endl;
	for ( map<string, map<int, double > >::iterator ite = chr_pos_sc.begin(); ite != chr_pos_sc.end(); ++ite )
	{
		string chr = ite->first;
		if (chr != tchr )
			continue;
		outf<<"variableStep chrom="+chr+" span="<<span<<endl;
		for ( map<int, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outf<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outf.close();
}


void readwigfile( string infile, map<string, map<int, double > > &score_map )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line );
	vector<string > parsed_items;
	getline( inf, line );
	parsed_items = parse_string(line);
	while(!inf.eof())
	{
		
		if ( line.empty() )
			break;
		vector<string> chrom_ve = parse_string( parsed_items[1], '=' );
		string chr = chrom_ve[1];
		map<int, double > site_sc;
		getline( inf, line );
		parsed_items = parse_string(line);
		
		while ( parsed_items.size() == 2 )
		{
			int site = atoi( parsed_items[0].c_str() );
			double sc = atof( parsed_items[1].c_str() );
			site_sc.insert( make_pair(site, sc ) );
			
			getline( inf, line );
			if ( line.empty() )
				break;
			if ( inf.eof() )
				break;
			parsed_items = parse_string(line);
			
		}
		
		if ( !site_sc.empty() )
			score_map.insert( make_pair( chr, site_sc ) );
		
		if ( line.empty() )
			break;
		if ( inf.eof() )
			break;
		
	}
	
	inf.close();
	
	
}

void readwigandchangeres( string infile, map<string, map<int, double > > &score_map, int span )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
		
	string line;
	getline(inf, line );
	vector<string > parsed_items;
	getline( inf, line );
	parsed_items = parse_string(line);
	while(!inf.eof())
	{
		
		if ( line.empty() )
			break;
		vector<string> chrom_ve = parse_string( parsed_items[1], '=' );
		string chr = chrom_ve[1];
		map<int, double > site_sc;
		getline( inf, line );
		parsed_items = parse_string(line);
		
		while ( parsed_items.size() == 2 )
		{
			int site = atoi( parsed_items[0].c_str() );
			double sc = atof( parsed_items[1].c_str() );
			if ( sc > 0 )
				site_sc.insert( make_pair(site, sc ) );
			
			getline( inf, line );
			if ( line.empty() )
				break;
			if ( inf.eof() )
				break;
			parsed_items = parse_string(line);
			
		}
		
		if ( !site_sc.empty() )
		{
			map<int, double >::iterator ite = site_sc.begin();
			int site = ite->first;
			double sc = ite->second;
			int index = site/span;
			double aggsc = sc;
			++ite;
			for ( ; ite != site_sc.end(); ++ite )
			{
				site = ite->first;
				sc = ite->second;
				int tindex = (site) / span;
				if ( tindex == index )
				{
					aggsc += sc;
				} else
				{
					int bsite = index * span + 1;
					double bsc = aggsc / span;
					score_map[chr].insert( make_pair(bsite, bsc ) );
					index = tindex;
					aggsc = sc;
				}
			}
			int bsite = index * span + 1;
			double bsc = aggsc / span;
			score_map[chr].insert( make_pair(bsite, bsc) );
		}
		
		if ( line.empty() )
			break;
		if ( inf.eof() )
			break;
		
	}
	
	inf.close();
	
}

void readwigandchangeres_bunch( string infile, map<string, map<int, double > > &score_map, int span )
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
		
		readwigandchangeres( line, score_map, span );
		
	}
	inf.close();
}

void readwigfile_bunch( string infile, map<string, map<int, double > > &score_map )
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
		
		readwigfile( line, score_map );
		
	}
	inf.close();
}

void cal_column_average( vector<vector<double> > &mat, vector<double > &ave )
{
	vector< double > addve = mat[0];
	size_t rsize = addve.size();
	for ( size_t i = 1; i < mat.size(); ++i )
	{
		if ( mat[i].size() != rsize )
		{
			cout<<"error in cal_column_average: rows are not in the same size "<<rsize<<" "<<mat[i].size()<<endl;
			exit(1);
		}
		for ( size_t j = 0; j < rsize; ++j )
		{
			addve[j] += mat[i][j];
		}
	} 
	int totalr = (int)mat.size();
	ave.clear();
	for ( size_t j = 0; j < rsize; ++j )
	{
		double a = addve[j] / totalr;
		ave.push_back( a );
	}
	
}

void cal_column_sum( vector<vector<double> > &mat, vector<double > &addve )
{
	addve = mat[0];
	size_t rsize = addve.size();
	for ( size_t i = 1; i < mat.size(); ++i )
	{
		if ( mat[i].size() != rsize )
		{
			cout<<"error in cal_column_average: rows are not in the same size "<<rsize<<" "<<mat[i].size()<<endl;
			exit(1);
		}
		for ( size_t j = 0; j < rsize; ++j )
		{
			addve[j] += mat[i][j];
		}
	} 
	
	
}

int calsize( map<string, vector<pair<int, int > > >& regions )
{
	int n = 0;
	for ( map<string, vector<pair<int, int > > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		n += (int)ite->second.size();
	}
	return n;
}

void cal_summits_from_ve( vector<double > &ve, map<int, double > &summits )
{
	int tsite = 0;
	int tstate = 1;
	double tsc = ve[0];
	for ( int i = 1; i < ve.size(); ++i )
	{
		int site = i;
		double sc = ve[i];
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == -1 && tstate == 1 )
		{
			summits.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
	
}

void cal_summits_from_map( map<int, int> &m, map<int, int > &summits )
{
	map<int, int >::iterator ite = m.begin();
	int tsite = ite->first;
	int tstate = 1;
	int tsc = ite->second;
	++ite;
	for ( ; ite != m.end(); ++ite )
	{
		int site = ite->first;
		int sc = ite->second;
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == -1 && tstate == 1 )
		{
			summits.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
	
}

void cal_summits_from_map( map<int, double> &m, map<int, double > &summits )
{
	map<int, double >::iterator ite = m.begin();
	int tsite = ite->first;
	int tstate = 1;
	double tsc = ite->second;
	++ite;
	for ( ; ite != m.end(); ++ite )
	{
		int site = ite->first;
		double sc = ite->second;
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == -1 && tstate == 1 )
		{
			summits.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
	if ( summits.empty() )
	{
		summits.insert( make_pair( tsite, tsc ) );
	}
	
}

void cal_summits_from_map_win( map<int, int> &m, map<int, int > &summits, int win )
{
	map<int, vector<int> > index_score_ve;
	for ( map<int, int >::iterator ite = m.begin(); ite != m.end(); ++ite )
	{
		int site = ite->first;
		int index = site / win;
		int sc = ite->second;
		index_score_ve[index].push_back( sc );
		
	}
	map<int, double > index_score;
	for ( map<int, vector<int> >::iterator ite = index_score_ve.begin(); ite != index_score_ve.end(); ++ite )
	{
		double addv = 0;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			addv += ite->second[i];
		}
		double ave = addv/(int)ite->second.size();
		index_score.insert( make_pair( ite->first, ave ) );
	}
	map<int, double > summits_win;
	cal_summits_from_map( index_score, summits_win );
	for ( map<int, double >::iterator ite = summits_win.begin(); ite != summits_win.end(); ++ite )
	{
		int start = ite->first*win;
		int end = (ite->first+1)*win-1;
		int maxv = -1;
		int maxs = -1;
		for ( int k = start; k <= end; ++k )
		{
			if ( m.find(k) != m.end() )
			{
				if ( m[k] > maxv )
				{
					maxv = m[k];
					maxs = k;
				}
			}
		}
		if ( maxs == -1 )
		{
			cout<<"error in func cal_summits_from_map_win: cannot locat summit from max win "<<start<<" "<<end<<endl;
			exit(1);
		}
		summits.insert( make_pair(maxs, maxv ) );
	}
}

void cal_summits_from_map_win( map<int, double> &m, map<int, double > &summits, int win )
{
	map<int, vector<double> > index_score_ve;
	for ( map<int, double >::iterator ite = m.begin(); ite != m.end(); ++ite )
	{
		int site = ite->first;
		int index = site / win;
		double sc = ite->second;
		index_score_ve[index].push_back( sc );
		
	}
	map<int, double > index_score;
	for ( map<int, vector<double> >::iterator ite = index_score_ve.begin(); ite != index_score_ve.end(); ++ite )
	{
		double addv = 0;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			addv += ite->second[i];
		}
		double ave = addv/(int)ite->second.size();
		index_score.insert( make_pair( ite->first, ave ) );
	}
	map<int, double > summits_win;
	cal_summits_from_map( index_score, summits_win );
	for ( map<int, double >::iterator ite = summits_win.begin(); ite != summits_win.end(); ++ite )
	{
		int start = ite->first*win;
		int end = (ite->first+1)*win-1;
		double maxv = -1;
		int maxs = -1;
		for ( int k = start; k <= end; ++k )
		{
			if ( m.find(k) != m.end() )
			{
				if ( m[k] > maxv )
				{
					maxv = m[k];
					maxs = k;
				}
			}
		}
		if ( maxs == -1 )
		{
			cout<<"error in func cal_summits_from_map_win: cannot locat summit from max win "<<start<<" "<<end<<endl;
			exit(1);
		}
		summits.insert( make_pair(maxs, maxv ) );
	}
}

void cal_bottoms_from_ve( vector<double > &ve, map<int, double > &bottoms )
{
	int tsite = 0;
	int tstate = 1;
	double tsc = ve[0];
	for ( int i = 1; i < ve.size(); ++i )
	{
		int site = i;
		double sc = ve[i];
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == 1 && tstate == -1 )
		{
			bottoms.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
}

void cal_boundaries_for_peaks( map<int, double > &summits, map<int, double > &bottoms, 
	int leftsite, int rightsite,
	map<int, pair<int, int > > &peak_boundary )
{
	for ( map<int, double >::iterator ite = summits.begin(); ite != summits.end(); ++ite )
	{
		int site = ite->first;
		double sc = ite->second;
		
		int leftb = leftsite;
		for ( map<int, double >::iterator bi = bottoms.begin(); bi != bottoms.end(); ++bi )
		{
			int s = bi->first;
			if ( s < site )
			{
				leftb = s;
			} else
				break;
			
		}
		int rightb = rightsite;
		for ( map<int, double >::reverse_iterator rbi = bottoms.rbegin(); rbi != bottoms.rend(); ++rbi )
		{
			int s = rbi->first;
			if ( s > site )
				rightb = s;
			else
				break;
		}
		peak_boundary.insert( make_pair( site, make_pair(leftb, rightb ) ) );
			
	} 
}

void merge_tiny_peaks( map<int, pair<int, int > > &peak_boundary,
	map<int, double > &summit_score,
	map<int, double > &bottom_score,
	int win_str )  
{
	map<int, pair<bool, bool> > summit_type;
	for ( map<int, pair<int, int > >::iterator ite = peak_boundary.begin(); ite != peak_boundary.end(); ++ite )
	{
		int left_d = ite->first - ite->second.first;
		int right_d = ite->second.second - ite->first;
		bool left = false;
		if ( left_d < win_str )
			left = true;
		bool right = false;
		if ( right_d < win_str )
			right = true;
		summit_type.insert( make_pair(ite->first, make_pair(left, right) ) );
	}
	
	vector<pair<int, int > > cuts;
 	map<int, pair<bool, bool > >::iterator ite = summit_type.begin();
	int leftb = peak_boundary[ite->first].first;
	if ( ite->second.second )
	{
		++ite;
		while ( ite != summit_type.end() )
		{
			if ( !ite->second.second )
				break;
			++ite;
		}
	}
	if ( ite == summit_type.end() )
		--ite;
	int rightb = peak_boundary[ite->first].second;
	cuts.push_back( make_pair( leftb, rightb ) );
	 
	++ite;
	for ( ; ite != summit_type.end(); ++ite )
	{
		if ( ite->second.first )
		{
			++ite;
			while ( ite != summit_type.end() )
			{
				if ( !ite->second.first )
					break;
				++ite;
			}
		}
		if ( ite == summit_type.end() )
		{
			break;
		} 
		
		leftb = peak_boundary[ite->first].first;
		
		if ( ite->second.second )
		{
			++ite;
			while ( ite != summit_type.end() )
			{
				if ( !ite->second.second )
					break;
				++ite;
			}
		}
		if ( ite == summit_type.end() )
			--ite;
		rightb = peak_boundary[ite->first].second;
		cuts.push_back( make_pair( leftb, rightb ) );
	}
	
	if ( cuts.back().second != peak_boundary.rbegin()->second.second )
		cuts.back().second = peak_boundary.rbegin()->second.second;
	
	map<int, pair<int, int > > new_peak_boundary;
	for ( size_t i = 0; i < cuts.size(); ++i )
	{
		
		if ( i == 0 )
		{
			leftb = cuts[0].first;
		} else
		{
			int l = cuts[i-1].second;
			int r = cuts[i].first;
			if ( l == r )
				leftb = r;
			else
			{
				multimap< double, int > score_bottom;
				for ( map<int, double >::iterator ite = bottom_score.begin(); ite != bottom_score.end(); ++ite )
				{
					if ( ite->first >= l && ite->first <= r )
						score_bottom.insert( make_pair( ite->second, ite->first) );
				}
				if ( score_bottom.empty() )
				{
					cout<<"error in func merge_tiny_peaks: score_bottom empty "<<l<<" "<<r<<endl;
					exit(1);
				}
				leftb = score_bottom.begin()->second;
			}
		}
		
		multimap< double, int, greater<double> > score_summit;
		for ( map<int, double >::iterator ite = summit_score.begin(); ite != summit_score.end(); ++ite )
		{
			if ( ite->first >= cuts[i].first && ite->first <= cuts[i].second )
			{
				score_summit.insert( make_pair( ite->second, ite->first ) );
			}
		}
		if ( score_summit.empty() )
		{
			cout<<"error in func merge_tiny_peaks: score_summit empty "<<cuts[i].first<<" "<<cuts[i].second<<endl;
		}
		int summit = score_summit.begin()->second;
		
		if ( i == cuts.size()-1 )
		{
			rightb = cuts[i].second;
		} else
		{
			int l = cuts[i].second;
			int r = cuts[i+1].first;
			if ( l == r )
				rightb = r;
			else
			{
				multimap< double, int > score_bottom;
				for ( map<int, double >::iterator ite = bottom_score.begin(); ite != bottom_score.end(); ++ite )
				{
					if ( ite->first >= l && ite->first <= r )
						score_bottom.insert( make_pair( ite->second, ite->first) );
				}
				if ( score_bottom.empty() )
				{
					cout<<"error in func merge_tiny_peaks: score_bottom empty "<<l<<" "<<r<<endl;
					exit(1);
				}
				rightb = score_bottom.begin()->second;
			}
		}
		
		new_peak_boundary[summit] = make_pair( leftb, rightb ); 
	}
	
	peak_boundary = new_peak_boundary;
}

void get_boundary_summit_plus1( map<int, pair<int, int > > &peak_boundary,
	int oricenter,  // 1000
	int &upbound,
	int &downbound,
	int &summit )    // positions are relative to oricenter
{
	int tdis = 1000;
	upbound = 1000;
	downbound = 1000;
	summit = 1000;
	for ( map<int, pair<int, int > >::iterator ite = peak_boundary.begin(); ite != peak_boundary.end(); ++ite )
	{
		int left = ite->second.first;
		int d = abs( left - oricenter );
		if ( d < tdis )
		{
			upbound = ite->second.first - oricenter;
			downbound = ite->second.second - oricenter;
			summit = ite->first - oricenter;
			tdis = d;
		}
	}
	if ( summit == 1000 )
	{
		cout<<"error cannot get boundary summit plus1"<<endl;
		exit(1);
	}
}

void cal_density_and_win_smooth( map<int, int >& dis_count, map<int, double > &den, int halfw )
{
	int add = 0;
	for ( map<int, int >::iterator ite = dis_count.begin(); ite != dis_count.end(); ++ite )
	{	
		add += ite->second;
	}	
	map<int, double > p_d;
	for ( map<int, int >::iterator ite = dis_count.begin(); ite != dis_count.end(); ++ite )
	{	
		double d = ite->second*1.0/add;
		p_d.insert( make_pair(ite->first, d) );
	}	
	for ( map<int, double >::iterator ite = p_d.begin(); ite != p_d.end(); ++ite )
	{
		int s = ite->first;
		int start = s-halfw;
		int end = s + halfw;
		double agg = 0;
		int n = 0;
		for ( int k = start; k <= end; ++k )
		{
			if ( p_d.find(k) != p_d.end() )
			{
				agg += p_d[k];
				n += 1;
			}
		}
		double ave = agg/n;
		den.insert( make_pair(ite->first, ave ) );
	}
	
}

void smooth_ve( vector<double > &ve, int halfw )
{
	vector<double > newve;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		int startk = i-halfw;
		int endk = i+halfw;
		if ( startk < 0 )
			startk = 0;
		if ( endk > ve.size()-1 )
			endk = ve.size()-1;
		double agg = 0;
		for ( int k = startk; k <= endk; ++k )
		{
			agg += ve[k];
		}
		double ave = agg/(endk-startk+1);
		newve.push_back( ave );
	}
	ve = newve;
}

void cal_deviation_peak( map<int, pair<int, int > > &peak_boundary,
	vector< vector<double > > &score_mat,
	map<int, double > &peak_deviation,
	int length,
	int extend )  // extend=100. size of score_mat is 2201. valid size is 2000
{
	vector< double > addve_e;
	cal_column_sum( score_mat, addve_e );
	if ( addve_e.size() != length+2*extend )
	{
		cout<<"error addve_e size in func cal_deviation_peak "<<addve_e.size()<<endl;
		exit(1);
	}
	vector< double > addve;
	for ( int i = 0; i < length; ++i )
	{
		addve.push_back( addve_e[i+extend]);
	}
	
	map< int, vector<int > > peak_dis;
//	cout<<"size "<<addve.size()<<endl;
	for ( int i = 0; i < (int)addve.size(); ++i )
	{
	/*	if ( i % 100 == 0 )
		{
			cout<<i<<" "<<addve[i]<<endl;
		} */
		int n = (int)addve[i];
		for ( map<int, pair<int, int > >::iterator ite = peak_boundary.begin(); ite != peak_boundary.end(); ++ite )
		{
			int leftb = ite->second.first;
			int rightb = ite->second.second;
			if ( i > rightb )
				continue;
			if ( i < leftb )
				break;
			int d = abs( i - ite->first );
			for ( int m = 0; m < n; ++m )
			{
				peak_dis[ite->first].push_back( d );
			}
		}
	} 
	
	for ( map< int, vector<int > >::iterator ite = peak_dis.begin(); ite != peak_dis.end(); ++ite )
	{
		int s = ite->first;
		double d = 0;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			d += ( ite->second[i] * ite->second[i] );
		}
		double v = 0;
		if ( (int)ite->second.size() == 1 )
		{
			v = 1;
		} else if ( (int)ite->second.size() > 1 )
		{
			v = pow( d/((int)ite->second.size()), 0.5 );
		}
		peak_deviation.insert( make_pair(s, v) );
	} 
//	cout<<peak_dis.size()<<" "<<peak_deviation.size()<<endl;
}

map<int, double > find_max( vector<map<int, double > > close_map, int gap )
{
	map<int, double > final_map;
	for ( size_t i = 0; i < close_map.size(); ++i )
	{
		if ( (int)close_map[i].size() == 1 )
		{
			final_map.insert( close_map[i].begin(), close_map[i].end() );
			continue;
		}
		map<int, double >::iterator ite = close_map[i].begin();
		int maxp = ite->first;
		double maxv = ite->second;
		++ite;
		for ( ; ite != close_map[i].end(); ++ite )
		{
			if ( ite->second > maxv )
			{
				maxp = ite->first;
				maxv = ite->second;
			}
		}
		
		final_map.insert( make_pair(maxp, maxv ) );
		
		map<int, double > left_map;
		map<int, double > right_map;
		for ( ite = close_map[i].begin(); ite != close_map[i].end(); ++ite )
		{
			if ( ite->first < maxp - gap )
			{
				left_map.insert( make_pair(ite->first, ite->second ) );
			} else if ( ite->first > maxp + gap )
			{
				right_map.insert( make_pair(ite->first, ite->second ) );
			}
		}
		vector< map<int, double > > next_close_map;
		if ( !left_map.empty() )
		{
			next_close_map.push_back( left_map );
		} 
		if ( ! right_map.empty() )
		{
			next_close_map.push_back( right_map ); 
		}
		if ( !next_close_map.empty() )
		{
			map<int, double > next_final_map = find_max(next_close_map, gap );
			final_map.insert(  next_final_map.begin(), next_final_map.end() );
		}
	}
	
	return final_map;
}

void remove_nearby_summit( map<string, map<int, double > > &chr_summit_map,
	map<string, map<int, double > > &rm_chr_summit_map )
{
	int distance = 100;
	for ( map<string, map<int, double > >::iterator ite = chr_summit_map.begin(); 
		ite != chr_summit_map.end(); ++ite )
	{
		string chr = ite->first;
		
		vector< map<int, double > > separated_map;
		map<int, double > one_map;
		map<int, double >::iterator si = ite->second.begin();
		one_map.insert( make_pair( si->first, si->second ) );
		++si;
		while ( si != ite->second.end() )
		{
			if ( si->first < one_map.rbegin()->first + distance )
			{
				one_map.insert( make_pair( si->first, si->second ) );
				++si;
			} else
			{
				vector<map<int, double > > orimap;
				orimap.push_back( one_map );
				map<int, double > fmap = find_max( orimap, distance );
				rm_chr_summit_map[chr].insert( fmap.begin(), fmap.end() );
				one_map.clear();
				one_map.insert( make_pair( si->first, si->second ) );
				++si;
			}
			
		}
		if ( !one_map.empty() )
		{
			vector<map<int, double > > orimap;
			orimap.push_back( one_map );
			map<int, double > fmap = find_max( orimap, distance );
			rm_chr_summit_map[chr].insert( fmap.begin(), fmap.end() );
			
		}
	}
		
}

void remove_nearby_summit( map<int, double > &summit_map,
	map<int, double > &rm_summit_map, int distance )
{
	
	
		vector< map<int, double > > separated_map;
		map<int, double > one_map;
		map<int, double >::iterator si = summit_map.begin();
		one_map.insert( make_pair( si->first, si->second ) );
		++si;
		while ( si != summit_map.end() )
		{
			if ( si->first < one_map.rbegin()->first + distance )
			{
				one_map.insert( make_pair( si->first, si->second ) );
				++si;
			} else
			{
				vector<map<int, double > > orimap;
				orimap.push_back( one_map );
				map<int, double > fmap = find_max( orimap, distance );
				rm_summit_map.insert( fmap.begin(), fmap.end() );
				one_map.clear();
				one_map.insert( make_pair( si->first, si->second ) );
				++si;
			}
			
		}
		if ( !one_map.empty() )
		{
			vector<map<int, double > > orimap;
			orimap.push_back( one_map );
			map<int, double > fmap = find_max( orimap, distance );
			rm_summit_map.insert( fmap.begin(), fmap.end() );
			
		}
	
		
}

int getmean(vector<int> &ve )
{
	long a = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		a += ve[i];
	}
	int ave = a / (int)ve.size();
	return ave;
}

double getmean(vector<double > &ve )
{
	double a = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		a += ve[i];
	}
	double ave = a / (int)ve.size();
	return ave;
}

int getn50( vector<int> &ve )
{
	multiset<int, greater<int> > sortedve;
	int total = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		sortedve.insert(ve[i]);
		total += ve[i];
	}
	
	int hlf = total /2;
	int g = 0;
	int n = 0;
	for ( multiset<int, greater<int > >::iterator ite = sortedve.begin(); ite != sortedve.end(); ++ite )
	{
		g += *ite;
		if ( g >= hlf )
		{
			n = *ite;
			break;
		}
	}
	
	return n;

}

int getmedium( vector<int> &ve )
{
	multiset<int, greater<int> > sortedve;
	
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		sortedve.insert(ve[i]);
		
	}
	
	int s = (int)ve.size() / 2;
	int i = 0;
	int  m = 0;
	for ( multiset<int, greater<int > >::iterator ite = sortedve.begin(); ite != sortedve.end(); ++ite )
	{
		i+=1;
		if ( i == s )
		{
			m = *ite;
			break;
		}
	}
	return m; 
	
}

int getmax(vector<int > &ve )
{
	int m = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		if ( ve[i] > m )
			m = ve[i];
		
	}
	return m;
}

int getsum(vector<int > &ve )
{
	int s = 0; 
	for ( size_t i = 0; i < ve.size(); ++i )
		s += ve[i];
	return s;
}

pair<double, double > getmeanstd( vector<int > &ve )
{
	int sum = getsum(ve );
	double ave = sum*1.0 / (int)ve.size();
	double d = 0;
	for (size_t i = 0; i < ve.size(); ++i )
	{
		d += (ve[i]-ave) * (ve[i]-ave);
	}
	d /= (int)ve.size();
	double std = pow(d, 0.5);
	return make_pair(ave, std);
}

pair<double, double> getmeanstd(vector<double > &ve)
{
	double ave = getmean( ve );
	double agg = 0;
	for ( size_t i = 0; i < ve.size(); ++i  )
	{
		agg += (ve[i] - ave )*(ve[i] - ave );
	}
	double aggave = agg / (int)ve.size();
	double d = pow(aggave, 0.5);
	return make_pair(ave, d);
}

double getstd( vector<double > &ve )
{
	double add = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		add += ve[i];
	}
	int n = (int)ve.size();
	double ave = add/n;
	
	double agg = 0;
	for ( size_t i = 0; i < ve.size(); ++i  )
	{
		agg += (ve[i] - ave )*(ve[i] - ave );
		
	}
	double aggave = agg/n;
	double d = pow( aggave, 0.5 );
	return d;
}

double getCV( vector<double > &ve )
{
	
	double add = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		add += ve[i];
	}
	int n = (int)ve.size();
	double ave = add/n;
	if ( ave == 0 )
		return 0;
		
	double agg = 0;
	for ( size_t i = 0; i < ve.size(); ++i  )
	{
		agg += (ve[i] - ave )*(ve[i] - ave );
		
	}
	double aggave = agg/n;
	double d = pow( aggave, 0.5 );
	return d/ave;
}


double entropy( vector<double > &Pi )
{
	double s = 0;
	for ( size_t i = 0; i < Pi.size(); ++i )
	{
		s += Pi[i] * log(Pi[i]);
	}
	s *= (-1);
	return s;
}




