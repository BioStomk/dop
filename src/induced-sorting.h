#include<iostream>
#include<fstream>
#include<stdlib.h>
using namespace std;

// $:0, A:1, C:2, G:3, T:4
// type: S:1 , L:0

/* InducedSorting
 * 配列sのsuffix-array SAを構築する */
class IS{
	private:
    const int* s_;		// 配列
	int* SA_;			// Suffix Array
	const int size_;	// 配列の長さ
	const int num_char_;// 使う文字の種類数 $ACGT:5
	int* type_;			// suffixのtype S:1 L:0 : (size_)
	int* bucket_;		// SA中で各文字を次に入れる位置 : (num_char_)
	bool isLMS(const int i){	// suffix[i]:S-type かつ suffix[i-1]:L-type ならLMS
		if(i > 0 && type_[i] && !type_[i-1]) return true;	
		else return false;
	}
	void get_buckets(bool computes_end);  // SA中で各文字が次に入る位置を計算
	void check_type();		// 各suffixを{S-type/L-type}に分類
	void sort_L_suffixes();
	void sort_S_suffixes();
	public:
	IS(const int* s, int* SA, const int size, const int num_char)
		:s_(s), SA_(SA), size_(size), num_char_(num_char){};
	~IS(){};
	void run();		// アルゴリズム本体を実行
	void disp_SA(){
		for(int i = 0; i < size_; i++)
			cout << SA_[i] << endl;
	}
	void disp_suffix(){
		for(int i = 0; i < size_ ; i++){
			cout << SA_[i] << "\t";
			for(int j = 0; j < 7; j++){
				if(SA_[i]+j >= size_) break;
				cout << s_[SA_[i]+j];
			}
			cout << endl;
		}
	}
};

