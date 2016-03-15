#include "induced-sorting.h"

/* 各suffixを{S-type/L-type}に分類 */
void IS::check_type(){
	type_[size_-1] = 1;
	for(int i = size_-2; i >= 0; i--){
		if(s_[i] == s_[i+1])
			type_[i] = type_[i+1];
		else if(s_[i] < s_[i+1])
			type_[i] = 1;	// S
		else
			type_[i] = 0;	// L
	}
}

/* SA中で各文字が次に入る位置を計算
 * computes_endがtrueなら各バケットの末尾を
 * falseなら各バケットの先頭の位置を得る */
void IS::get_buckets(bool computes_end){
	for(int j = 0; j < num_char_; j++)
		bucket_[j] = 0;
	for(int i = 0; i < size_; i++)
		bucket_[s_[i]]++;
	int sum = 0;
	for(int j = 0; j < num_char_; j++){
		sum += bucket_[j];
		bucket_[j] = (computes_end ? sum - 1 : sum - bucket_[j]); 
	}
}

/* L-typeのsuffixをソート */
void IS::sort_L_suffixes(){
	get_buckets(false);	// 各文字のbucketの先頭の位置を得る
	for(int i = 0; i < size_; i++){
		int prev = SA_[i]-1;
		if(prev >= 0 && !type_[prev])  // type[prev] = 'L'
			SA_[bucket_[s_[prev]]++] = prev;
	}
}

/* S-typeのsuffixをソート */
void IS::sort_S_suffixes(){
	get_buckets(true);	// 各文字のbucketの末尾の位置を得る
	for(int i = size_-1 ;i >= 0; i--){
		int prev = SA_[i]-1;
		if(prev >= 0 && type_[prev])	 // type[prev] = 'S'
			SA_[bucket_[s_[prev]]--] = prev;
	}
}		


/* 本体 */
void IS::run(){
	/* (1) LMS-prefixesに半順序付け */
	type_ = new int[size_];
	check_type();
	bucket_ = new int[num_char_];
	get_buckets(true);
	for(int i = 0; i < size_; i++)
		SA_[i] = -1;
	for(int i = 0; i < size_; i++){
		if(isLMS(i))
			SA_[bucket_[s_[i]]--] = i;
	}
	sort_L_suffixes();
	sort_S_suffixes();

	/* (2) LMS-substringsのソートの準備 */
	int n1 = 0;		// LMSの数
	// LMSの数をカウントしながら、LMSをSAの前に詰めて入れる
	for(int i = 0; i < size_; i++){
		if(isLMS(SA_[i]))
			SA_[n1++] = SA_[i];
	}
	for(int i = n1; i < size_; i++)	// 残りは-1にする
		SA_[i] = -1;
	// LMS-substringsに順序を付ける
	int order = 0;	// LMSの順序
	int prev = -1;  // 1つ前の順序のLMS
	for(int i = 0; i < n1; i++){
		int pos = SA_[i];
		bool diff = false;	// posがprevより大きいか
		for(int d = 0; ; d++){
			if(prev == -1 || s_[pos+d] > s_[prev+d] || type_[pos+d] != type_[prev+d]){
				diff = true;
				break;
			}
			else if(d > 0 && isLMS(pos+d))	// 先に次のLMSに到着
				break;
		}
		if(diff){
			order++;
			prev = pos;
		}
		SA_[n1 + pos/2] = order-1;	// 
	}
	// LMS-substringsを順序に変換したものを後ろに寄せる
	for(int i = size_-1, j = size_-1; i >= n1; i--)
		if(SA_[i] != -1)
			SA_[j--] = SA_[i];
	
	delete type_;   // 再帰問題を呼ぶ前にメモリを解放して消費量を抑える
	delete bucket_;
	
	/* (3) "LMS-substringsを順序に変換したもの"を新たな配列と見なして再帰問題を解く */
	int* s1 = SA_+ size_- n1;	// 解くべき配列s1
	int* SA1 = SA_;		// SA1に新たなsuffix-arrayを作る(場所はSA_の先頭に上書きする)
	if(order < n1){
		IS reduced_prob(s1, SA1, n1, order);
		reduced_prob.run();
	}else{
		for(int i = 0; i < n1; i++)
			SA1[s1[i]] = i;
	}

	/* (4) LMSのsuffix-array SA1をもとに元の配列のSAを構築する */
	type_ = new int[size_];  // 元の配列のtypeを復元
	check_type();
	// s1上に元のLMSをリストアップ
	for(int i = 0, j = 0; i < size_; i++)
		if(isLMS(i)) s1[j++] = i;
	// SA1の値を元のLMSに書き換える(SA_の先頭にLMS-suffixが辞書順に並ぶ)
	for(int i = 0; i < n1; i++)
		SA_[i] = s1[SA1[i]];
	for(int i = n1; i < size_; i++)	// 残りは-1を入れる
		SA_[i] = -1;
	// LMSを各bucketに入れて残りのsuffixもソートする
	bucket_ = new int[num_char_];
	get_buckets(true);
	for(int i = n1-1; i >= 0; i--){
		int tmp = SA_[i];
		SA_[i] = -1;
		SA_[bucket_[s_[tmp]]--] = tmp;
	}
	sort_L_suffixes();
	sort_S_suffixes();
	
	delete bucket_;
	delete type_;
}

