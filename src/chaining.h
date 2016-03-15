#include<iostream>
#include<fstream>
#include<map>
using namespace std;

/* 2つの配列のアラインメント
 * 区間は[sx_,ex_)と[sy_,ey_) */
class Alignment{
	private:
	int sx_; int ex_;
	int sy_; int ey_;
	int score_;
	public:
	Alignment()
		:sx_(0), ex_(0), sy_(0), ey_(0), score_(0){};
	Alignment(int sx, int ex, int sy, int ey, int score)
		:sx_(sx), ex_(ex), sy_(sy), ey_(ey), score_(score) {};
	~Alignment(){};
	void set(int sx, int ex, int sy, int ey, int score){
		sx_ = sx; ex_ = ex;
		sy_ = sy; ey_ = ey;
		score_ = score;
	}
	void set(int* array){
		sx_ = array[0]; ex_ = array[1];
		sy_ = array[2]; ey_ = array[3];
		score_ = array[4];
	}
	void disp(){
		cout << "[" << sx_ << "," << ex_ << ")\t[" << sy_ << "," << ey_ << ")\t" << score_ << endl;
	}
	void output_startpos(ofstream& ofs) {ofs << ((sx_>0)? sx_: -sx_) << "\t" << sy_ << endl;}
	bool operator < (const Alignment& a){
		if(ex_ <= a.sx_ && ey_ <= a.sy_) return true;
		else return false;
	}
	bool operator > (const Alignment& a){
		if(ex_ >= a.sx_ && ey_ >= a.sy_) return true;
		else return false;
	}
	friend class Chain;
	friend class Chaining;
};


/* アラインメントをつなげた鎖 */
class Chain{
	private:
	Alignment* const end_al_;	//鎖の最後に来るアラインメント(固定)
	Chain* prev_;	//前につながる鎖
	int global_chain_score_;	//鎖に含まれるすべてのアラインメントのスコアの和
	int local_chain_score_;		//連続したひとまとまりのchainのスコア
	bool is_printed_;	// トラックバックの際にすでに出力したか
	public:
	Chain(Alignment& a)
		: end_al_(&a), prev_(NULL), global_chain_score_(a.score_),
		  local_chain_score_(a.score_), is_printed_(false){};
	~Chain(){};
	//受け取った鎖の後ろにアラインメントをつなげる
	void attach(Chain* prev_chain, const int near_dist){
		prev_ = prev_chain;
		global_chain_score_ += prev_chain->global_chain_score_;
		// 前の鎖が近かったらlocal_scoreも積算
		if(is_near(prev_chain, near_dist))
			local_chain_score_ += prev_chain->local_chain_score_;
	}
	// localなまとまりを判断するために使用
	bool is_near(Chain* prev_chain, const int near_dist){
		if(prev_chain->end_al_->ex_ > end_al_->sx_- near_dist
		   && prev_chain->end_al_->ey_ > end_al_->sy_- near_dist)
			return true;
		else return false;
	}
	friend class Chaining;
};


/* アラインメントの集合を受け取ってChainingを行う
 * 各アラインメントで終わるスコア最大の鎖を計算する */
class Chaining{
	private:
	Alignment* alignments_;     //アラインメントの集合: num_alignments_
	const int num_alignments_;	//アラインメントの本数
	Chain** chains_;	    	//Chainへのポインタの集合: num_alignments_
	multimap<int, Chain*> list_X_;	//各鎖Cの(C.sx,C)と(C.ex, C)を入れるリスト
	multimap<int, Chain*> list_Y_;	//(C.ey, C)を入れるリスト
	const int near_dist_;		// アラインメントが"近い"と判断する距離
	public:
	Chaining(Alignment* alignments, const int num_alignments, const int near_dist)	// コンストラクタ
		: alignments_(alignments), num_alignments_(num_alignments), near_dist_(near_dist)
    	{chains_ = new Chain*[num_alignments_];}
	~Chaining(){ delete chains_;}	// デストラクタ
	void run();  //各アラインメントで終わるスコア最大の鎖を計算
	//各鎖のスコアを表示
	void disp_global_scores(){
		cout << "global_chain_scores\n";
		for(int i = 0; i < num_alignments_; i++){
			cout << i << "\t" << chains_[i]->global_chain_score_ << endl;
		}
	}
	void disp_local_scores(){
		cout << "local_chain_scores\n";
		for(int i = 0; i < num_alignments_; i++){
			cout << i << "\t" << chains_[i]->local_chain_score_ << endl;
		}
	}
	// 連続したまとまりになっている鎖について、含まれるアラインメントの座標を出力
	void output_major_chains(ofstream& ofs) {
		for(int i = num_alignments_-1; i > 0; i--){
			Chain* c = chains_[i];
			if(c->local_chain_score_ > 5 && !(c->is_printed_)){
				do{
					c->end_al_->output_startpos(ofs);
					c->is_printed_ = true;
					c = c->prev_;
				}while(c->local_chain_score_ != 1);  // まとまりの先頭に来たらbreak
				c = chains_[i]; // 元の位置に戻す
			}
		}
	}
};
	
