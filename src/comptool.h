#include"burrows-wheeler-transform.h"
#include"induced-sorting.h"
#include"chaining.h"
#include<stdlib.h>
#include<stdint.h>
#include<string>
#include<sstream>
#include<fstream>

class CompTool{
    private:
    const int num_char_;
    int seq1_size_;
    int seq2_size_;
    int8_t* seq1_;
    int8_t* seq2_;

    void search(int argc, char** argv);
    int* create_SA(const string file, const int size);
    void search_forward_matches(const string seq1_name, const string seq2_name, int* SA, BWT& bwt, const int kmer_size, const int slide_letters, const int max_num_matches);
    void search_reverse_matches(const string seq1_name, const string seq2_name, int* SA, BWT& bwt, const int kmer_size, const int slide_letters, const int max_num_matches);
    void chain(int argc, char** argv);
    void run_chaining(string file, ofstream& ofs, const int near_dist, const bool forward);

    public:
    CompTool(): num_char_(5){};
    ~CompTool(){};
    void run_command(int argc, char** argv);
};
