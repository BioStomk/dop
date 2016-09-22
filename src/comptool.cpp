#include"comptool.h"

int main(int argc, char** argv){
    CompTool comp;
    comp.run_command(argc, argv);
    return 0;
}


string basename(const string& path){
    return path.substr(path.find_last_of('/') + 1);
}


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


int get_seq_length(const string file){
    ifstream ifs;
    ifs.open(file);
    if(!ifs.is_open()) {cout << "Cannot find file:  " << file << endl; exit(1);}
    string header;
    getline(ifs, header);   // Skip header
    int size = 0;
    char buf;
    while(ifs >> buf) size++;
    ifs.close();
    return size;
}


int* read_fasta_and_create_int_array(const string file, const int size){
    ifstream ifs;
    ifs.open(file);
    if(!ifs.is_open()) {cout << "Cannot find file:  " << file << endl; exit(1);}
    string header;
    getline(ifs, header);
    int* s = new int[size];
    int* p = s;
    char buf;
    while(ifs >> buf) *p++ = encode_char(buf);
    *p++ = 0;   // Append '0'
    return s;
}


int8_t* read_fasta_and_create_int8_t_array(const string file, const int size){
    ifstream ifs;
    ifs.open(file);
    if(!ifs.is_open()) {cout << "Cannot find file:  " << file << endl; exit(1);}
    string header;
    getline(ifs, header);
    int8_t* s = new int8_t[size];
    int8_t* p = s;
    char buf;
    while(ifs >> buf) *p++ = (int8_t)encode_char(buf);
    *p++ = 0;   // Append '0'
    return s;
}


int8_t* convert_int_array_to_int8_t_array(int* int_array, const int size){
    int8_t* int8_t_array = new int8_t[size];
    for(int i = 0 ; i < size; i++)
        int8_t_array[i] = int_array[i];
    return int8_t_array;
}


void CompTool::run_command(int argc, char** argv){
    argc--; argv++;
    const string command = argv[0];
    if(command == "align")
        search_alignment(argc, argv);
    else if(command == "chain")
        chain_alignment(argc, argv);
    else
        cout << "Unknown command: \"" << command << "\"" << endl;
}


int* CompTool::create_SA(const string file, const int size){
    int* seq = read_fasta_and_create_int_array(file, size);
    int* SA = new int[size];
    IS is(seq, SA, size, num_char_);
    is.run();
    delete seq;
    return SA;
}


void CompTool::search_forward_matches(const string seq1_name, const string seq2_name, int* SA, BWT& bwt, const int kmer_size,
                                      const int slide_letters, const int max_num_matches, const bool outputs_start_pos){
    stringstream out_file;
    if(outputs_start_pos)
        out_file << "alignments-forward-startpos_" << seq1_name << "_" << seq2_name << ".tsv";
    else
        out_file << "alignments-forward-for-chaining_" << seq1_name << "_" << seq2_name << ".tsv";
    ofstream ofs(out_file.str().c_str());

    ofs << "#" << seq2_name << "\t" << seq1_name << endl;  // header
    for(int i = 0; i < seq2_size_ - kmer_size; i += slide_letters){
        int8_t* query = &seq2_[i];
        int lb, ub;  // lower- and upper-bound of matches in suffix array
        bwt.search(query, kmer_size, lb, ub);
        if(lb <= ub){
            for(int j = lb; j <= ub; j++){
                if(j == lb + max_num_matches) break;
                if(outputs_start_pos)
                    output_startpos(ofs, i, SA[j]);
                else
                    output_alignment_forward(ofs, i, SA[j], kmer_size);
            }
        }
    }
}


void CompTool::search_reverse_matches(const string seq1_name, const string seq2_name, int* SA, BWT& bwt, const int kmer_size,
                                      const int slide_letters, const int max_num_matches, const bool outputs_start_pos){
    stringstream out_file;
    if(outputs_start_pos)
        out_file << "alignments-backward-startpos_" << seq1_name << "_" << seq2_name << ".tsv";
    else
        out_file << "alignments-backward-for-chaining_" << seq1_name << "_" << seq2_name << ".tsv";
    ofstream ofs(out_file.str().c_str());

    ofs << "#" << seq2_name << "\t" << seq1_name << endl;  // header
    int8_t* query = new int8_t[kmer_size];
    for(int i = kmer_size-1; i < seq2_size_ - 1; i += slide_letters){
        // convert k-mers to reverse complements
        for(int j = 0; j < kmer_size; j++)
            query[j] = num_char_ - seq2_[i-j];
        int lb, ub;  // lower- and upper-bound of matches in suffix array
        bwt.search(query, kmer_size, lb, ub);
        if(lb <= ub){
            for(int j = lb; j <= ub; j++){
                if(j == lb + max_num_matches) break;
                if(outputs_start_pos)
                    output_startpos(ofs, i, SA[j]);
                else
                    output_alignment_backward(ofs, i, SA[j], kmer_size);
            }
        }
    }
    delete[] query;
}


void CompTool::search_alignment(int argc, char** argv){
    if(argc < 3) {cout << "File not enough" << endl; exit(1);}

    const string seq1_file = argv[1];
    const string seq2_file = argv[2];
    const string seq1_name = basename(seq1_file);
    const string seq2_name = basename(seq2_file);

    // Options
    int  kmer_size         = 15;
    int  slide_letters     = 1;
    int  bwt_interval      = 1;
    int  max_num_matches   = 1000000000;
    bool search_forward    = true;
    bool search_reverse    = true;
    bool outputs_start_pos = false;

    if(argc > 3){
        for(int i = 3; i < argc; i++){
            if     (argv[i][1] == 'k') kmer_size         = atoi(argv[++i]);
            else if(argv[i][1] == 'l') slide_letters     = atoi(argv[++i]);
            else if(argv[i][1] == 'i') bwt_interval      = atoi(argv[++i]);
            else if(argv[i][1] == 'm') max_num_matches   = atoi(argv[++i]);
            else if(argv[i][1] == 'f') search_reverse    = false;
            else if(argv[i][1] == 'r') search_forward    = false;
            else if(argv[i][1] == 's') outputs_start_pos = true;
        }
    }

    seq1_size_ = get_seq_length(seq1_file) + 1;
    seq2_size_ = get_seq_length(seq2_file) + 1;
    seq1_ = read_fasta_and_create_int8_t_array(seq1_file, seq1_size_);
    seq2_ = read_fasta_and_create_int8_t_array(seq2_file, seq2_size_);

    int* SA = create_SA(seq1_file, seq1_size_);
    BWT bwt(seq1_, SA, seq1_size_, num_char_, bwt_interval);

    if(search_forward) search_forward_matches(seq1_name, seq2_name, SA, bwt, kmer_size, slide_letters, max_num_matches, outputs_start_pos);
    if(search_reverse) search_reverse_matches(seq1_name, seq2_name, SA, bwt, kmer_size, slide_letters, max_num_matches, outputs_start_pos);

    delete SA;
}


// argv[1]: seq1, argv[2]: seq2
// (OPTION) [-n]: near_dist
// [-f]: runs forward only, [-b]: runs backward only
void CompTool::chain_alignment(int argc, char** argv){
    if(argc < 3) {cout << "Files not enough" << endl; exit(1);}
    // OPTION
    int near_dist = 50;
    bool runs_forward = true;
    bool runs_backward = true;
    if(argc > 3){
        for(int i = 3; i < argc; i++)
            if(argv[i][1] == 'n')  near_dist = atoi(argv[++i]);
            else if(argv[i][1] == 'f')  runs_backward = false;
            else if(argv[i][1] == 'b')  runs_forward = false;
    }

    const string seq1_file = argv[1];
    const string seq2_file = argv[2];

    const string seq1_name = basename(seq1_file);
    const string seq2_name = basename(seq2_file);

    stringstream out_file;
    out_file << "chains_" << seq1_name << "_" << seq2_name << ".tsv";
    ofstream ofs(out_file.str().c_str());
    ofs << "#" << seq2_name << "\t" << seq1_name << endl;

    if(runs_forward){
        stringstream ifs1_name;
        ifs1_name << "alignments-forward-for-chaining_" << seq1_name << "_" << seq2_name
                  << ".tsv";
        run_chaining(ifs1_name.str(), ofs, near_dist);
    }

    if(runs_backward){
        stringstream ifs2_name;
        ifs2_name << "alignments-backward-for-chaining_" << seq1_name << "_" << seq2_name
                  << ".tsv";
        run_chaining(ifs2_name.str(), ofs, near_dist);
    }
}


// Receive an alignment file, run chaining and output a chain file
void CompTool::run_chaining(string file, ofstream& ofs, const int near_dist){
    // Count alignments
    ifstream ifs;
    ifs.open(file.c_str());
    if(!ifs.is_open()) {cout << "File not found\n> " << file << endl; exit(1);}
    string line;
    int num_alignment = 0;
    getline(ifs, line);  // Skip header
    while(getline(ifs, line)) num_alignment++;
    ifs.close();

    // Load alignments
    ifs.open(file.c_str());
    getline(ifs, line);  // Skip header
    Alignment* alignments = new Alignment[num_alignment];
    int* buf = new int[5];
    for(int i = 0; i < num_alignment; i++){
        for(int j = 0; j < 5; j++){
            ifs >> buf[j];
        }
        alignments[i].set(buf);
    }
    delete[] buf;

    // Run chaining
    Chaining chaining(alignments, num_alignment, near_dist);
    chaining.run();
    chaining.output_major_chains(ofs);
}
