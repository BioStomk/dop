#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<stdint.h>
using namespace std;

//$:0, A:1, C:2, G:3, T:4

class BWT{
	private:
	const int* SA_;			//suffix_array
	const int length_;		//配列の長さ
	const int num_char_;	//文字の種類の数
	int8_t* bwt_;				//(SA_[i] - 1)番目の文字: (length_)
	int* cumm_;				//各文字の累積出現回数: (num_char_)
	int** occ_;				//occ_[x][i]: bwt_の[i*interval_]の位置までに文字xが出現する回数
							//(num_char_ * (length_/interval_)) ただし'$'の分のocc[0]は使用しない
	const int interval_;	//occ_内でデータを保持しておく間隔
	int8_t encode_char(const char c);	// char('ACGT')をint8_t(1234)に変換する
	public:
	BWT(const int8_t* sequence, const int* suffix_array, const int length, const int num_char, const int interval);
	~BWT();
	void interactive_search();	// コンソールを使って対話的に検索
	void search(const int8_t* query, const int query_length, int& lb, int& ub);	// queryにマッチする位置を探索
};
