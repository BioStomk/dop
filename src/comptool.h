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
	// char('ACGT')をint(1234)に変換する
	int encode_char(const char c){
		switch(c){
		case '$': return 0;
		case 'A': case 'a': return 1;
		case 'C': case 'c': return 2;	  
		case 'G': case 'g': return 3;
		case 'T': case 't': return 4;
		case 'N': case 'n': return 1;	// 'N'は'A'として扱う
		default: return 1;
		}
	}
	// fastaファイルを読み込んで必要なメモリを確保し、'ACGTN'を'12340'に変換して格納する
	int* read_fasta_and_create_int_array(const char* file, int& size);
	int8_t* read_fasta_and_create_int8_t_array(const char* file, int& size);
	// intの配列の中身をint8_tの配列にコピーする
	void copy_int_array_to_int8_t_array(int* int_array, int8_t* int8_array, const int size);
	// 配列1をファイルから読み込んでsuffix-arrayを構築
	int* create_SA(const char* file);
	// 配列2からk-merをqueryとして取り出し、配列1のBWTを使ってマッチする位置を検索
	void search_alignment(int argc, char** argv, int* SA);
	// 座標のみを表示
	void output_startpos(ofstream& ofs, const int sx, const int sy){
		ofs << sx << "\t" << sy << endl;
	}
	// アラインメントの開始位置と終了位置の座標、およびスコア(1とする)を表示
	void output_alignment_forward(ofstream& ofs,const int sx, const int sy, const int k){
		ofs << sx << "\t" << sx+k << "\t" << sy << "\t" << sy+k << "\t" << 1 << endl;
	}
	void output_alignment_backward(ofstream& ofs,const int sx, const int sy, const int k){
		ofs << -sx << "\t" << -sx+k << "\t" << sy << "\t" << sy+k << "\t" << 1 << endl;
	}
	// アラインメントをChainingする
	void chain_alignment(int argc, char** argv);
	void run_chaining(string file, ofstream& ofs, const int near_dist);

	public:
	CompTool(): num_char_(5){};
	~CompTool(){};
	void run_command(int argc, char** argv);
};
