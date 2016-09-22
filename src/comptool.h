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

    int* create_SA(const string file, const int size);
    void search_alignment(int argc, char** argv);

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
