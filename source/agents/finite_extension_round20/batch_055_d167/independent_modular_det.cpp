// Independent exact GMP Gaussian determinant on raw integer Phi rows.
// This file shares no reduction or FLINT determinant implementation.
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using Z=mpz_class;
using namespace std;
int intfield(const string&s,const string&k){auto p=s.find("\""+k+"\":");if(p==string::npos)throw runtime_error("missing integer field");return stoi(s.substr(p+k.size()+3));}
vector<Z> coeffs(const string&s){auto b=s.find("\"coefficients\":[");if(b==string::npos)throw runtime_error("no coefficients");b+=16;vector<Z>a;while(b<s.size()){while(b<s.size()&&(s[b]==' '||s[b]==','))b++;if(s[b]==']')break;if(s[b]!='\"')throw runtime_error("coefficient not quoted");auto e=s.find('"',b+1);if(e==string::npos)throw runtime_error("unterminated coefficient");a.emplace_back(s.substr(b+1,e-b-1));b=e+1;}return a;}
struct Answer{string status;Z determinant,factor;long swaps=0,pivots=0;};
Answer determinant_mod(const vector<vector<Z>>&src,const Z&p){
 const int n=src.size();vector<vector<Z>> a(n,vector<Z>(n));
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)mpz_mod(a[i][j].get_mpz_t(),src[i][j].get_mpz_t(),p.get_mpz_t());
 Answer out;out.status="nonzero_determinant";out.determinant=1;Z inv,factor,g;
 for(int k=0;k<n;k++){
  int r=k;while(r<n && a[r][k]==0)r++;
  if(r==n){out.status="zero_determinant";out.determinant=0;return out;}
  if(r!=k){swap(a[r],a[k]);out.determinant=-out.determinant;out.swaps++;}
  if(!mpz_invert(inv.get_mpz_t(),a[k][k].get_mpz_t(),p.get_mpz_t())){
   mpz_gcd(out.factor.get_mpz_t(),a[k][k].get_mpz_t(),p.get_mpz_t());
   if(out.factor<=1 || out.factor>=p)throw runtime_error("invalid nonunit gcd");
   out.status="composite_factor";return out;
  }
  out.determinant*=a[k][k];mpz_mod(out.determinant.get_mpz_t(),out.determinant.get_mpz_t(),p.get_mpz_t());out.pivots++;
  for(int i=k+1;i<n;i++)if(a[i][k]!=0){
   factor=a[i][k]*inv;mpz_mod(factor.get_mpz_t(),factor.get_mpz_t(),p.get_mpz_t());a[i][k]=0;
   for(int j=k+1;j<n;j++){
    mpz_submul(a[i][j].get_mpz_t(),factor.get_mpz_t(),a[k][j].get_mpz_t());
    mpz_mod(a[i][j].get_mpz_t(),a[i][j].get_mpz_t(),p.get_mpz_t());
   }
  }
 }
 return out;
}
int main(int argc,char**argv){
 if(argc!=5){cerr<<"Usage: independent_modular_det d raw.jsonl candidates.txt out.jsonl\n";return 2;}
 int d=stoi(argv[1]),N=4*d;map<pair<int,int>,vector<Z>> raw;ifstream f(argv[2]);string line;
 while(getline(f,line)){
  if(line.empty())continue;
  if(intfield(line,"d")!=d||intfield(line,"nd")!=N)throw runtime_error("mixed source identity");
  int n=intfield(line,"n"),m=intfield(line,"m");
  if(n<d||n>N||n%d||m<0||m>=n/d*(d-1)/2||!(n<N||m<=d))throw runtime_error("illegal source row");
  auto co=coeffs(line);if(co.size()!=size_t(n+1))throw runtime_error("coefficient length");
  if(co[0]!=0||co[1]!=12*co[2])throw runtime_error("low moments");
  if(!raw.emplace(make_pair(n,m),move(co)).second)throw runtime_error("duplicate row");
 }
 vector<vector<Z>> a;for(int b=1;b<=4;b++)for(int m=0;m<(b<4?b*(d-1)/2:d+1);m++){
  auto it=raw.find({b*d,m});if(it==raw.end())throw runtime_error("missing expected row");
  vector<Z> row(N-2);for(size_t k=3;k<it->second.size();k++)row[k-3]=it->second[k];a.push_back(move(row));
 }
 if(a.size()!=size_t(N-2)||raw.size()!=size_t(N-2))throw runtime_error("not square prefix");
 ifstream cf(argv[3]);ofstream out(argv[4]);if(!out)throw runtime_error("cannot write output");int count=0;
 while(getline(cf,line)){
  if(line.empty())continue;Z p(line);if(p<=2*N)throw runtime_error("characteristic gate");Z rem=p%(2*d);if(rem==1||rem==2*d-1)throw runtime_error("congruence gate");
  auto t=chrono::steady_clock::now();auto r=determinant_mod(a,p);double sec=chrono::duration<double>(chrono::steady_clock::now()-t).count();
  out<<"{\"candidate\":\""<<p<<"\",\"d\":"<<d<<",\"N\":"<<N<<",\"status\":\""<<r.status<<"\",\"determinant_mod_n\":\""<<r.determinant<<"\",\"factor\":\""<<r.factor<<"\",\"unit_pivots\":"<<r.pivots<<",\"row_swaps\":"<<r.swaps<<",\"seconds\":"<<sec<<"}\n";out.flush();if(!out)throw runtime_error("output write failed");
  cerr<<"{\"completed\":"<<++count<<",\"d\":"<<d<<",\"seconds\":"<<sec<<",\"status\":\""<<r.status<<"\"}"<<endl;
 }
}
