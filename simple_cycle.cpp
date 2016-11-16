// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <map>
#include <list>
#include <deque>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;




vector<vector<int> > simple_cycle( vector<vector<int> > & adjList );
bool  findCycles(int v, int s, vector<vector<int> >& adjList,vector<bool> & blocked,  deque<int>& stackLike,  vector<vector<int> >& B_Fruitless,vector<vector<int> > &  cycles);
void unblock(int  node, vector<bool> & blocked, vector<vector<int> > & B_Fruitless);

void displayVector( vector<int> vs){
  cout << "Vector: " ;
  for( vector<int>::iterator idx=vs.begin(); idx!=vs.end(); idx++ ){
    cout << *idx << ",";
  }
  cout << endl;
}

void displayVV( vector< vector<int> > vss){
  for( vector<vector<int> >::iterator id = vss.begin(); id !=vss.end(); id++){
    cout << "Vector: " ;
    vector<int> vs = *id;
    for( vector<int>::iterator idx=vs.begin(); idx!=vs.end(); idx++ ){
      cout << *idx << ",";
    }
    cout << endl;
  }
  cout<< "END " << endl;
}


// [[Rcpp::export]]
map<int, vector<int> > get_cycles(const mat & adjMat ){
  map<int, vector<int> > rlt;

  vector< vector<int> > adjList;
  for( unsigned int i = 0; i < adjMat.n_rows; i++)
  {
    vector<int> tmp;
    for( unsigned int j = 0; j < adjMat.n_rows; j++)
    {
      if( adjMat(i,j) > 0 )
      {
        tmp.push_back(j);
      }
    }
    adjList.push_back(tmp);
  }
  //displayVV(adjList);
  int count = 0;
  vector< vector<int> > cycles = simple_cycle(adjList);
  for( auto cy = cycles.begin(); cy != cycles.end(); cy++ ){
    rlt[count++] =  *cy;
  }
  return rlt;
}




bool  findCycles(int v, int s, vector<vector<int> >& adjList,vector<bool> & blocked,  deque<int>& stackLike,  vector<vector<int> >& B_Fruitless,vector<vector<int> > &  cycles) {

  bool f = false;
  stackLike.push_front(v);  // insert like a stack:  so at the front
  blocked[v] = true;

  /*cout << "BLOCKED list from visiting node" << v << ":::\t" ;
  for (unsigned int i =0; i<adjList.size() ;i++) {
    cout << blocked[i] << "\t";
  }
  cout << endl;

  */
  // explore all neighbours -- recursively
  for (unsigned int i = 0; i < adjList[v].size(); i++) {
    int w =  adjList[v][i];

    // found cycle through ANY of my neighbours w.
    if (w == s) {
      //cout << "----- FoUND cycle with s=" << s <<  "and v=" << v << " and Size=" << stackLike.size() << ": " ;
      vector<int>* cycle = new vector<int>;
      for (unsigned int j = 0; j < stackLike.size(); j++) {
        cycle->push_back( stackLike.at (stackLike.size() - j  - 1) );
        //cout << stackLike.at(stackLike.size() - j -1) << " ";
      }
      cycles.push_back(*cycle);
      //cout << endl;
      f = true;
    } else if (!blocked[w]) {
      if (findCycles(w, s, adjList, blocked, stackLike, B_Fruitless, cycles)) {
        f = true;
      }
    }
  } // for


  if (f) {
    //cout << "F is true! just found a cycle, and unblock starting with " << v << endl;
    unblock(v, blocked, B_Fruitless);
  } else {
    //cout << "--no cycles found for " << v << endl;

    // go through all my neighbors w.
    //  v is pushed on B_Fruitless[w].
    // looking at B_Fruitless[w] = [v1, v2, ..] , i know i can get from v1 to w, and v2 to w, etc.
    // later, whenever i block w, i recursively unblock v1, and v2, and v..

    for (unsigned int i = 0; i < adjList[v].size(); i++) {
      int w = adjList[v][i];

      // mark B_Fruitless[w] to point to v.  This says that going from v to w lead to an unfruitful search.
      // later when w is found to particiate in a cycle, i'd better get rid of this false assupmtion about
      // w not leading to fruitful cycles.

      vector<int>::iterator it;
      it = find(B_Fruitless[w].begin(), B_Fruitless[w].end(), v);
      if (it == B_Fruitless[w].end()) {
        //cout << "Pushing v=" << v << "on B_Fruitless list for w=" << w << "\n";
        B_Fruitless[w].push_back(v);
      }
    }
  }
  // find v and remove it from stack
  deque<int>::iterator it;

  it = find(stackLike.begin(), stackLike.end(), v);
  if (it != stackLike.end() ) {
    stackLike.erase(it);
  }
  return f;
} // bool  findCycles


vector<vector<int> > simple_cycle(vector<vector<int> > & adjList) {

  //cout << "Adj List:\n";
  //displayVV (adjList);


  vector<bool> blocked (adjList.size(), false);

  deque<int>  stackLike;

  // cycles
  vector<vector<int> > cycles;

  // B_Fruitless is the book keeping needed to avoid fruitless searches.  It is
  // referred to as B in Johnson's original algorithm


  // initialize B
  vector<vector<int> >  B_Fruitless ;

  for (unsigned int i=0; i< adjList.size()  ; i++ ) {
    vector<int>* k = new vector<int>;
    B_Fruitless.push_back(*k);
  }

  // loop to start new search from each node i
  for (unsigned int i=0; i< adjList.size()  ; i++ ) {
    // clear all book keeping
    for (unsigned int j =0; j< adjList.size(); j++) {
      blocked[j] = false;
      B_Fruitless[j].clear();
    }
    //cout << "START: i=" << i << endl;
    findCycles(i, i, adjList,  blocked, stackLike, B_Fruitless, cycles) ;
  }

  //cout << "Cycles:\n";
  //displayVV (cycles);

  return(cycles);

}

void unblock(int  node, vector<bool> & blocked, vector<vector<int> > & B_Fruitless)
{
  blocked[node] = false;
  while (B_Fruitless[node].size() > 0) {
    int w = B_Fruitless[node][0];
    B_Fruitless[node].erase(find (B_Fruitless[node].begin(), B_Fruitless[node].end(), w) ) ;
    if (blocked[w]) {
      unblock(w, blocked, B_Fruitless);
    }
  }
}
