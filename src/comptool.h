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
    int8_t* seq1_;
    int seq1_size_;
    int8_t* seq2_;
    int seq2_size_;

    int encode_char(const char c){
        switch(c){
        case '$': return 0;
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        case 'N': case 'n': return 1;  // Treat 'N' as 'A'
        default: return 1;             // Treat unknown character as 'A'
        }
    }

    int* read_fasta_and_create_int_array(const char* file, int& size);
    int8_t* read_fasta_and_create_int8_t_array(const string& file, int& size);
    void copy_int_array_to_int8_t_array(int* int_array, int8_t* int8_array, const int size);

    int* create_SA(const char* file);
    void search_alignment(int argc, char** argv, int* SA);

    void output_startpos(ofstream& ofs, const int sx, const int sy){
        ofs << sx << "\t" << sy << endl;
    }
    void output_alignment_forward(ofstream& ofs,const int sx, const int sy, const int k){
        ofs << sx << "\t" << sx+k << "\t" << sy << "\t" << sy+k << "\t" << 1 << endl;
    }
    void output_alignment_backward(ofstream& ofs,const int sx, const int sy, const int k){
        ofs << -sx << "\t" << -sx+k << "\t" << sy << "\t" << sy+k << "\t" << 1 << endl;
    }

    void chain_alignment(int argc, char** argv);
    void run_chaining(string file, ofstream& ofs, const int near_dist);

    public:
    CompTool(): num_char_(5){};
    ~CompTool(){};
    void run_command(int argc, char** argv);
};
