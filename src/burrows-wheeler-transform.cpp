#include"burrows-wheeler-transform.h"

/* コンストラクタ */
BWT::BWT(const int8_t* sequence, const int* suffix_array, const int length, const int num_char, const int interval)
	: SA_(suffix_array), length_(length), num_char_(num_char), interval_(interval)
{
	// bwt_[i]にsequence中の(SA_[i]-1)番目の文字をセット
	bwt_ = new int8_t[length_];
	for(int i = 0; i < length_; i++){
		if(SA_[i] == 0)
			bwt_[i] = 0;	//'$'
		else
			bwt_[i] = sequence[SA_[i]-1];
	}

	// cumm_[i]を計算
	cumm_ = new int[num_char_];
	for(int i = 0; i < num_char_; i++)
		cumm_[i] = 0;
	for(int i = 0; i < length_; i++)
		cumm_[bwt_[i]]++;
	cumm_[num_char_-1] = length_ - cumm_[num_char_-1];
	for(int i = num_char_-2; i >= 0; i--)
		cumm_[i] = cumm_[i+1] - cumm_[i];

	// occ_[i]を計算
	// ただし、'$'の分のocc[0]は使わない
	occ_ = new int*[num_char_];
	for(int i = 1; i < num_char_; i++){
		occ_[i] = new int[length_ / interval_];
		if(bwt_[0] == i)
			occ_[i][0] = 1;
		else
			occ_[i][0] = 0;
	}

	for(int j = 1; j < num_char_; j++){
		int tmp_count = 0;		//このinterval_行の間に文字が出現した回数
		for(int i = 1; i < length_; i++){
			if(bwt_[i] == j)
				tmp_count++;
			if(i % interval_ == 0){
				occ_[j][ i/interval_] = occ_[j][ i/interval_-1] + tmp_count;
				tmp_count = 0;
			}
		}
	}
}

/* デストラクタ */
BWT::~BWT(){
	delete bwt_;
	delete cumm_;
	for(int i = 1; i < num_char_; i++)
		delete occ_[i];
	delete occ_;
}


/* コンソールを使って対話的に検索する */
void BWT::interactive_search(){
	vector<int8_t> query;		//問い合わせ文字列
	string buf;				//コンソールから受け取る
	cout << "\n--To quit type 'q'--\n";
	
	while(true){
		cout << "query:  ";
		cin >> buf;
		if(buf == "q"){
			cout << endl;
			break;
		}
		bool received_valid_query = true;	//queryにACGT以外の文字が含まれていたらfalse
		for(unsigned int i = 0; i < buf.length(); i++){
			int8_t tmp = encode_char(buf[i]);
			if(tmp == -1){
				cout << "#Invalid query#" << endl;
				query.clear();
				received_valid_query = false;
				break;
			}else{
				query.push_back(tmp);
			}
		}
		if(received_valid_query){
			const int query_length = query.size();
			int8_t* p = new int8_t[query_length];
			for(int i = 0; i < query_length; i++){
				p[i] = query[i];
			}
			cout << "result: ";
			int lb, ub;
			search(p, query_length, lb, ub);
			for(int i = lb; i <= ub; i++){
				cout << SA_[i] << " ";
			}
			cout << endl;
			query.clear();
		}
	}
}

/* queryにマッチする位置を探索する */
void BWT::search(const int8_t* query, const int query_length, int& lb, int& ub){
	lb = 0;			//lower_bound in SA_
	ub = length_-1;	//upper_bound in SA_

	for(int i = query_length-1; i >= 0; i--){		//queryの後ろの文字から位置を絞り込んでいく
		//lbを計算
		if(lb == 0)									//ループの最初(末尾の文字)のとき
			lb = cumm_[query[i]];							
		else{
			int tmp_count = 0;						//lb-1以下で直近の足場からlb-1までの間にquery[i]が現れる回数
			for(int j = ((lb-1)/interval_)*interval_, k = 1; k <= (lb-1) % interval_; k++)
				if(bwt_[j+k] == query[i])
					tmp_count++;
			lb = cumm_[query[i]] + occ_[query[i]][(lb-1)/interval_] + tmp_count;
		}
		//ubを計算
		{
			int tmp_count  = 0;								
			for(int j = (ub/interval_)*interval_, k = 1; k <= ub % interval_; k++)
				if(bwt_[j + k] == query[i])
					tmp_count++;
			ub = cumm_[query[i]] + occ_[query[i]][ub/interval_] + tmp_count - 1;
		}
		if(lb > ub)		//queryはにマッチするものは存在せず
			break;
	}
}


/* char('ACGT')をint8_t(1234)に変換する */
int8_t BWT::encode_char(const char c){
	switch(c){
	case '$': return 0;
	case 'A': case 'a': return 1;
	case 'C': case 'c': return 2;	  
	case 'G': case 'g': return 3;
	case 'T': case 't': return 4;
	case 'N': case 'n': return 1;	// 'N'は'A'として扱う
	default: return -1;
	}
}
