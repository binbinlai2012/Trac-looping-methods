#include "partition.h"
#include "pet.h"

void readinchr( string infile, map<string, int > &chr_len )
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
		int l = atoi(ps[1].c_str() );
		chr_len.insert( make_pair(chr, l ) );
	}
	inf.close();
}

void exit_with_help()
{
	cerr <<"call interaction backgroud model"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	cerr <<"-R		region file"<<endl;
	cerr <<"-L		chrlen file" <<endl;
	cerr <<"-b		bedpe file" <<endl;
	cerr <<"-w		bin size (default: 1K)"<<endl;
	cerr <<"-p		output prefix file input" <<endl;
	cerr <<"-l		lower bound pet length (default: 2K)"<<endl;
	cerr <<"-u		upper bound pet length (default: 2M)"<<endl;
	cerr <<"-c		minimal count in bin for consideration (default: 5)"<<endl;
	cerr <<"-d		minimal density(rpm) in bin for consideration (default: 0.1)"<<endl;
	exit(1);
}

void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}

int main( int argc, char* argv[] )
{
	
	string chrlenfile = "";
	string bedpefile = "";
	string outprefix = "";
	string regionfile = "";
	int win = 1000;
	int minl = 2000;
	int maxl = 2000000;
	int min_count = 5;
	double min_feaACC = 0.5;
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
		case 'L':
			chrlenfile = argv[i];
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
		case 'd':
			min_feaACC = atoi(argv[i]);
			break;
		default:
			exit_with_help();
		
		}
	}
	
	map<string, int > chr_len;
	readinchr( chrlenfile, chr_len );
	
	cout<<"Initiate partition from chr len file"<<endl;
	Partition_Bank p_bank;
	p_bank.ini_Partition_from_chrlen( chrlenfile, win );
	
	cout<<"Get par anchor from region file"<<endl;
	p_bank.add_partition_anchor_from_region( regionfile );
	
	cout<<"read in bedpe"<<endl;
	PET_bank pet_bank;
	pet_bank.readinPET( bedpefile, minl, maxl, chr_len );
	
	cout<<"assign bedpe to bin"<<endl;
	map<string, map<pair<size_t, size_t >, int > > par_pair_c;
	p_bank.assign_pet_to_par_pair2( pet_bank, minl, maxl, outprefix );
	p_bank.assign_pet_to_par_pair3( pet_bank, par_pair_c, minl, maxl );
	
	size_t minstep = (size_t)(minl / win);
	size_t maxstep = (size_t)(maxl / win);
	cout<<"generate valid bin pair space"<<endl;
	map<string, vector<pair<size_t, size_t > > > space;
	
	p_bank.generate_potential_space( space, minstep, maxstep, min_count, min_feaACC );
	
//	cout<<"cal space distance accessibility stat"<<endl;
//	p_bank.cal_space_dis_acc_stat( space, outprefix );

	cout<<"get space dist"<<endl;
	p_bank.space_dis_accm_distr( space, outprefix );
//	p_bank.space_dis_accm_distr( par_pair_c, minstep, maxstep, outprefix );
	cout<<"get PETs dist"<<endl;
	p_bank.pets_dis_accm_distr( space, par_pair_c, outprefix );
//	p_bank.pets_dis_accm_distr( par_pair_c, minstep, maxstep, outprefix );
	
	return 0;
	
}

