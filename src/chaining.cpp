#include"chaining.h"

/* 各アラインメントで終わるスコア最大の鎖を計算 */
void Chaining::run(){
	//Chainを各アラインメントで初期化
	for(int i = 0; i < num_alignments_; i++){
		chains_[i] = new Chain(alignments_[i]);
	}

	//各鎖Cの(C.sx,C),(C.ex, C)をlist_X_に入れる
	for(int i = 0; i < num_alignments_; i++){
		list_X_.insert(pair<int, Chain*>(chains_[i]->end_al_->sx_, chains_[i]));
		list_X_.insert(pair<int, Chain*>(chains_[i]->end_al_->ex_, chains_[i]));
	}

	//list_X_から順番に取り出して以下を行う
	for(map<int, Chain*>::iterator itx = list_X_.begin(); itx!=list_X_.end(); itx++){
		Chain* chain_x = itx->second;
		if(itx->first == chain_x->end_al_->sx_){	//第一要素がsx_のとき
			for(map<int, Chain*>::reverse_iterator ity = list_Y_.rbegin();
				ity != list_Y_.rend(); ity++){
				Chain* chain_y = ity->second;
				if(*(chain_y->end_al_) < *(chain_x->end_al_)){
					chain_x->attach(chain_y, near_dist_);
					break;
				}
			}
		}else{  	//第一要素がex_のとき
			bool needs_to_insert = true;	//list_Y_に入れるかどうか
			//list_Y_の中に"自分より前でスコアが高いもの"がなければ入れる
			for(map<int, Chain*>::iterator ity = list_Y_.begin();
				ity != list_Y_.end(); ity++){
				Chain* chain_y = ity->second;
				if(chain_y->end_al_->ey_ <= chain_x->end_al_->ey_){
					if(chain_y->global_chain_score_ >= chain_x->global_chain_score_){
						needs_to_insert = false;
						break;
					}
				}
			}
			if(needs_to_insert){
				list_Y_.insert(pair<int, Chain*>(chain_x->end_al_->ey_, chain_x));
				//list_Y_の中に"自分より後でスコアが低いもの"があればすべて削除
				for(map<int, Chain*>::reverse_iterator ity = list_Y_.rbegin();
					ity != list_Y_.rend(); ity++){
					Chain* chain_y = ity->second;
					if(chain_y->end_al_->ey_ < chain_x->end_al_->ey_)
						break;
					if(chain_y->end_al_->ey_ >= chain_x->end_al_->ey_)
						if(chain_y->global_chain_score_ < chain_x->global_chain_score_){
							list_Y_.erase(--(ity.base()));
						}
				}
			}
		}
	}
}

