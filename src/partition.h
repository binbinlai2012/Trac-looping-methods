#ifndef PARTITION_H
#define PARTITION_H

#include "operation.h"
#include "gene.h"
#include "pet.h"

class Par_Bin
{
public:
	int start;
	double freq_GC;
	double mapa;
	double fea_ACC;
	int end_count;
	bool loop;
	bool peak_anchor;
};

class Partition
{
public:
	string chr;
	vector< Par_Bin > par_ve;
	int win;
	size_t get_bin_by_start( int start );

};

class Partition_Bank
{
public:
	map<string, Partition > chr_Par;
	string kmer;
	
	int total_valid_pets_in_anchor;
	
	void ini_Partition_from_genome( Genome &genome, int win );
	void add_Partition_freq_GC( Genome &genome, string inkmer );
	void add_Partition_mapa( string infile );
	
	void ini_Partition_from_chrlen( string infile, int win );
//	void ini_Partition_from_Region( string infile, int win );
	
	void write_partition( string outfile );
	
	void read_partition_with_genomic_fea( string infile );
	
	void add_Partition_fea_ACC_from_PET_bank( PET_bank & pet_bank );
	
	void write_partition_ext2( string outfile );
	
	void read_partition_from_par2( string infile );
	
	void add_partition_anchor_from_region( string infile );
	
	void assign_pet_to_par_pair( PET_bank & pet_bank, 
		map<string, vector<pair<size_t, size_t > > > &par_pair, int minl, int maxl, int &ck, bool anchor );
	void assign_pet_to_par_pair( PET_bank & pet_bank, 
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c, string prefix );
	void assign_pet_to_par_pair2( PET_bank & pet_bank,  
		int minl, int maxl, string prefix );
	void assign_pet_to_par_pair3( PET_bank & pet_bank, 
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		int minl, int maxl );
	
	
	void statistics_fea_ACC();
	void statistics_freq_GC();
	void statistics_mapa();
	
	void generate_simulate_space(map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > &space );
	// length, acc_m, pair<mapa_a, mapa_b>, pair<freq_a, freq_b>  pair_event_count
	void generate_simulate_space(map<int, map<int, map<int, map<int, int > > > > &space );
	// length, acc_m, mapa_m, freq_m,  pair_event_count
	void generate_simulate_space( map<int, vector< pair<int, int > > > &length_accm_map,
		map<int, vector<int > > &space, int mings);
	// length, acc_m_range;   length, acc_m_range_count;    minimal group size default 10000.
	
	void generate_observed_space( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > &space );
	// length, acc_m, pair<mapa_a, mapa_b>, pair<freq_a, freq_b>  pets_event_count
	void generate_observed_space( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		map<int, map<int, map<int, map<int, int > > > > &space );
	// length, acc_m, mapa_m, freq_m,  pair_event_count
	void generate_observed_space( map<int, vector< pair<int, int > > > &length_accm_map,
		map<int, vector<int > > &sm_space,
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		double ex_rate,
		map<int, vector<int > > &space );
	
	void simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > & par_pair_c,
		map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > & space,
		map<int, map<int, map<pair<int, int>, map<pair<int, int >, int > > > > & obsspace,
		int times,
		string outfile );
	void simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > & par_pair_c,
		map<int, map<int, map<int, map<int, int > > > > & space,
		map<int, map<int, map<int, map<int, int > > > > & obsspace,
		int times,
		string outfile );
	void simulation_and_cal_pv( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		map<int, vector< pair<int, int> > > &length_accm_map,
		map<int, vector<int > > &sm_space,
		map<int, vector<int > > &obsspace,
		int times,
		string outfile ); 	
		
	void generate_potential_space( map<string, vector<pair<size_t, size_t > > > &space,
		size_t minstep, size_t maxstep, int min_count, double min_feaACC );
	void cal_space_dis_acc_stat( map<string, vector<pair<size_t, size_t > > > &space,
		string prefix );
	
	void space_dis_accm_distr( map<string, vector<pair<size_t, size_t > > > &space,
		string prefix );
	void pets_dis_accm_distr( map<string, vector<pair<size_t, size_t > > > &space,
		map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		string prefix);
	void space_dis_accm_distr( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		size_t minstep, size_t maxstep, string prefix );
	void pets_dis_accm_distr( map<string, map<pair<size_t, size_t >, int > > &par_pair_c,
		size_t minstep, size_t maxstep, string prefix);
};

int transfer_acc_m( double acc_m );

void sim_func( int ve_size, int total_pets,  int times,
	vector< vector<int > > &mat );
	
pair<double, double > cal_pv( int c, vector< vector<int > > &mat );
	
double BC( int n, int k );
double cdf( int k, int n, double p );
// binomial distribution

#endif

