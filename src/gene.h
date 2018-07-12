#ifndef GENE_H
#define GENE_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <fstream>
#include "operation.h"
using namespace std;



class Transcript
{
public:
	int start;
	int end;
	string chr;
	char strand;
	string geneName;   // Symbol
	string name;       // Refseq id or Ensmbl name
	int cdsStart;
	int cdsEnd;
	int exonCount;
	vector<int > exonStarts;
	vector<int > exonEnds;
	map< Sample_id, double > expression_map;
	double epus_idx;
	bool active;
	
	Transcript()
	{

	}

	Transcript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName )
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
	
	Transcript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds )
	{
		start = instart;
		end = inend;
		chr = inchr;
		strand = instrand;
		name = inname;
		cdsStart = incdsStart;
		cdsEnd = incdsEnd;
		exonCount = inexonCount;
		exonStarts = inexonStarts;
		exonEnds = inexonEnds;
		
	}

	Transcript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName )
	{
		start = instart;
		end = inend;
		chr = inchr;
		strand = instrand;
		name = inname;
		geneName = ingeneName;
		
	} 
	
	void addTranscript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName);
	void addTranscript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName);

	int getTSS();
	int getTES();
	pair<int, int > getPromoter(int extention );
	pair<int, int > getPromoter();
	pair<int, int > getPromoter(int up, int down );
	pair<int, int > getgenebody(int up, int down );
	pair<int, int > getTESregion(int up, int down );
	string getgenetype();
};



class Gene
{
public:
	int start;
	int end;
	string chr;
	char strand;
	set<string> TranscriptName;
	map< Sample_id, double > expression_map;
	string type;
	double epus_idx;
	bool active;
	Gene()
	{
	}

	Gene(int instart, int inend, string inchr)
	{
		start = instart;
		end = inend;
		chr = inchr;
	}

	void addTranscript(string name);
	
};

class LincRNAs
{
public:
	int start;
	int end;
	string chr;
	int win;
	int id;
	
	LincRNAs()
	{
	}
	
	LincRNAs(int inid, string inchr, int instart, int inend, int inwin )
	{
		id = inid;
		chr = inchr;
		start = instart;
		end = inend;
		win = inwin;
	}
	
};

class Activepromoter
{
public:
	string chr;
	int TSS;
	double epus_idx;
	string name;
	map< Sample_id, double > k4me3_map;
};

class Genome
{
public:
	map<string, Transcript > name_Transcript_map;
	map<string, vector<pair<int, int > > > intergenicregion;
	map<string, int > chr_len_map;
	
	map<string, map<pair<int, int >, string > > chr_pos_TranscriptName_map;
	map<string, map<pair<int, int >, vector<string> > > chr_pos_TranscriptNameS_map;   // same pos has multiple distinct trancripts
	map<string, Gene > name_Gene_map;
	map<string, map<pair<int, int >, string > > chr_pos_GeneName_map;
	map<string, map<pair<int, int >, vector<string> > > chr_pos_GeneNameS_map;   
	map<int, LincRNAs > id_LincRNAs_map;
	map<string, map<pair<int, int >, int > > chr_pos_LincRNAsId_map;
	
	vector< Activepromoter > actpro_Vec;
	map<string, map< int, Region_id > > chr_pos_actpro_map;
	
	map<string, string > genomeseq;
	
	void norm_seq();  // transfer lower chara into upper chara
	void get_kmer_distri( string kmer, map<string, vector<bool > > &plus, map<string, vector<bool > > &minus );
	
	map<string, map<pair<int, int>, string > > get_chr_promoter_TranscriptName_map(int up, int down);
	map<string, map<pair<int, int>, string > > get_chr_TES_TranscriptName_map(int up, int down);
	map<string, set<pair<int, int > > > get_chr_genic( int tssdown, int tesup );
	map<string, set<pair<int, int > > > get_chr_promoter( int up, int down );
	map<string, set<pair<int, int > > > get_genebody( int up, int down );
	map<string, set<int > > get_chr_TSS();
	map<string, set<int > > get_chr_TES();
	vector<int > get_TSS_from_Gene( string gene);
	void transcripttogene();
	void addgenomeseq(string infile);
	string getsubseq(string chr, int pos, int len );
	
	void readtranscriptfromucsc(string &infile );
	void readchrlen(string &infile );
	
	void getintergenicregion( );
	
	void getnunredundanttss_Ts( vector<string >& Ts, vector<string > &nr_Ts );
	void getnunredundanttss_Ts( vector<string > &nr_Ts );
	
	void read3genegroups_ParsingIndex( string infile, vector<string > &lowPI_Ts, vector<string > &highPI_Ts );
	
	void setgenicregion( vector<string > &active_Ts, vector<string > &silent_Ts,
		map<string, vector<pair<int, int> > > &active_TSS, 
		map<string, vector<pair<int, int > > > &active_GeneBody,
		map<string, vector<pair<int, int > > > &active_TES,
		map<string, vector<pair<int, int > > > &silent_TSS,
		map<string, vector<pair<int, int > > > &silent_GeneBody,
		map<string, vector<pair<int, int > > > &silent_TES,
		map<string, vector<pair<int, int > > > &intergenic );
		
	void remove_promoterregion( map<string, vector<int > > &regions, int upstream, int downstream );
	
};


class GeneExpression_GW
{
public:
	map<string, double > name_exp_map;
	void read_rpkm_from_file( string infile );

};

class TranscriptExpression_GW
{
public:
	map<string, double > name_exp_map;
	void read_rpkm_from_file( string infile );
	void filter_out_T(Genome &genome);   // with same tss or not in annotation
	
	void separate_into_active_silent_gene( vector<string > &active_T, vector<string > &silent_T, double cutoff );
	void separate_into_quantil_by_decending_rpkm( vector<vector<string > > &qual_Ts );
	void separate_into_4part_by_rpkm( vector<string > &silent_T, vector<vector<string > > &three_Ts, double cutoff );
	void separate_into_4part_by_breaks( vector<vector<string > > &four_Ts, vector<double> &breaks );
	void sortgene( vector<string > &sortedgene, double cutoff );
};



class lincRNAexpression_GW
{
public:
	map<int, double> id_exp_map;
};

class Sequence_Di
{
public:
	map<char, int> coding;
	map<int, string > din_revcode;
	
	void genecoding( );
	int getdicode( char n1, char n2 );
	vector<int > getdincode_fromseq( string seq );
	vector<int > getdincode_aroundnucl( string chr, int pos, Genome &gm );
	double cal_AATTfreq_Nuclfrank( string chr, int pos, Genome &gm ); 
	double cal_AATTATTAfreq_Nuclfrank( string chr, int pos, Genome &gm );
	
	void cal_din_freq_from_dincode_ve( vector< vector<int > > &dincode_ve, map<int, vector<double > > &din_rate );
	void cal_ave_freq( map<string, set<int > > &nucl_sets, Genome &gm, map<int, vector<double > > &din_rate );
	void cal_ave_freq( map<string, map<int, char > > &nucl_dir_sets, Genome &gm, map<int, vector<double > > &din_rate );
	
	
	double cal_nucl_pref_score( string chr, int pos, Genome &gm );
	void cal_nucl_pref_score_region( Genome &gm, string chr, int start, int end, vector<double > &sc_ve );
	
	double cal_nucl_rotation_score( string chr, int pos, Genome &gm );
	pair<int, double > realign_nucl_by_rotation_score( string chr, int pos, Genome &gm );
};

class SNP_bank
{
public:
	map<string, map<int, vector<char> > > chr_site_snp;
	map<string, map<int, set<int> > > chr_index_sites;   // L == 100000
//	map<string, map<int, vector<char> > > alt_chr_site_snp;  
	
	void readinsnp( string infile );
	void readinsnp_2( string infile );   // refbase go to front
	void generate_index();
	void get_snp_withinregion( vector<int> &snps, string chr, int start, int end );
	
	
	void getsnp_atnucl( string chr, int c, vector<int > &pos );
	void getsnp_atnucl_set( map<string, set<int > > &nucl_sets, vector<double > &freq );
	
	void readaltsnp( string infile );
	
};

pair<string, int> getnearestgene( string chr, int start, int end, Genome &gn);

class Conserv_bank
{
public:
	map<string, map<pair<int, int >, int > > phastConsE;
	map<string, map<int, vector<double > > > phyloP;
	map<string, map< int, set<int > > > chr_index_phyloPstart;
	
	void readinphastConsE( string infile );
	void readinphyloP( string infile );
	void readinphyloP_bunch( string infile );
	void generate_index();
	
	bool getvalid_phyloPS( string chr, int c, vector<double > &sc_ve );
	bool get_core_flank_phyloPS( string chr, int c, double & sc_flank, double & sc_core );
	void cal_phyloPS_nuclset( map<string, set<int > > &nucl_sets, vector<double > &ave_score );
};


#endif

