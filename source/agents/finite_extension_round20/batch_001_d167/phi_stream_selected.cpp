// Exact, low-memory transpose of the kappa=0 Phi reduction.
// Compile: g++ -O3 -std=c++17 -fopenmp phi_stream.cpp -lgmpxx -lgmp -o phi_stream
#include <gmpxx.h>
#include <omp.h>
#include <sys/resource.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
using namespace std;
using Z=mpz_class;
using Clock=chrono::steady_clock;
double secs(Clock::time_point t){return chrono::duration<double>(Clock::now()-t).count();}
long rss(){rusage r;getrusage(RUSAGE_SELF,&r);return r.ru_maxrss;}
vector<pair<int,int>> factor(int n){vector<pair<int,int>> f;for(int p=2;p*p<=n;p++)if(n%p==0){int e=0;while(n%p==0)n/=p,e++;f.emplace_back(p,e);}if(n>1)f.emplace_back(n,1);return f;}
vector<int> divisors(int n){vector<int>d{1};for(auto [p,e]:factor(n)){auto old=d;int pow=1;for(int i=1;i<=e;i++){pow*=p;for(int a:old)d.push_back(a*pow);}}return d;}
int phi(int n){int out=n;for(auto [p,e]:factor(n))out=out/p*(p-1);return out;}
int ipow(int a,int n){int r=1;while(n--)r*=a;return r;}
int valuation(int n,int p){int e=0;while(n%p==0)n/=p,e++;return e;}
vector<Z> gpoly(int r){
 // U_{r-1}(x/2), times x when needed to make it even, in T=x^2.
 vector<Z> a{1},b{0,1};int k=r-1;
 if(k==0)b=a;else for(int j=2;j<=k;j++){vector<Z> c(b.size()+1);for(size_t i=0;i<b.size();i++)c[i+1]+=b[i];for(size_t i=0;i<a.size();i++)c[i]-=a[i];a=move(b);b=move(c);}
 if(k%2)b.insert(b.begin(),Z(0));
 vector<Z> g;for(size_t i=0;i<b.size();i+=2)g.push_back(b[i]);return g;
}
struct Source {int n,mcount,gdegree;vector<Z> g;vector<vector<Z>> coeff;};
shared_ptr<Source> make_source(int d,int n,int q){
 auto src=make_shared<Source>();src->n=n;
 int raw=0;for(int i:divisors(2*n/(d*gcd(2,q))))raw+=phi(2*n/i);src->mcount=(raw+3)/4;
 int r=n*d/ipow(q,valuation(n,q)+1);src->g=gpoly(r);src->gdegree=src->g.size()-1;
 vector<Z> powers4(n+1);powers4[0]=1;for(int i=1;i<=n;i++)powers4[i]=-4*powers4[i-1];
 src->coeff.resize(n+1);
 for(int i=0;i<=n;i++){
  Z cni;mpz_bin_uiui(cni.get_mpz_t(),n+i,n-i);cni*=2*n;cni/=(n+i);if(i%2)cni=-cni;
  auto &row=src->coeff[i];row.resize(i+src->gdegree+1);Z bin;
  for(int j=0;j<=i;j++){
   mpz_bin_uiui(bin.get_mpz_t(),i,j);Z co=cni*bin*powers4[i-j];
   for(size_t k=0;k<src->g.size();k++)if(src->g[k]!=0)row[j+k]+=co*src->g[k];
  }
 }
 return src;
}
struct Prefix {
 vector<vector<int>> pre;vector<int> counts;
 Prefix(int maxs):pre(maxs+1),counts(maxs+1){for(int s=0;s<=maxs;s++){auto &p=pre[s];p.resize(s/2+1);int ct=0;for(int b=0;b<=s/2;b++){p[b]=ct;ct+=min(b,s-2*b)+1;}counts[s]=ct;}}
 inline int index(int s,int a,int b,int c)const{if(a<b)swap(a,b);if(b<c)swap(b,c);if(a<b)swap(a,b);return pre[s][b]+c;}
};
struct Workspace {vector<Z> a,b;Workspace(int n):a(n),b(n){}};
struct Result {vector<Z> coeff;uint64_t states=0,additions=0;int maxs;double seconds;};
Result reduce_row(const Source&src,int m,const Prefix&pr,Workspace&w){
 auto t0=Clock::now();int n=src.n,maxs=2*n+src.gdegree+m;Result res;res.coeff.resize(n+1);res.maxs=maxs;
 auto *cur=&w.a,*next=&w.b;
 for(int s=maxs;s>=0;s--){
  int j=s-n-m;
  if(j>=0)for(int i=0;i<=n;i++)if(j<(int)src.coeff[i].size()&&src.coeff[i][j]!=0){int a=s-i;if(a<0)throw runtime_error("negative input exponent");(*cur)[pr.index(s,a,i,0)]+=src.coeff[i][j];}
  for(int b=0;b<=s/2;b++)for(int c=0;c<=min(b,s-2*b);c++){
   int a=s-b-c;Z &co=(*cur)[pr.pre[s][b]+c];if(co==0)continue;res.states++;
   if(b==0){if(s>n)throw runtime_error("unexpected Phi degree above n");res.coeff[s]+=co;co=0;continue;}
   auto add=[&](int aa,int bb,int cc,unsigned long k){int ix=pr.index(s-1,aa,bb,cc);if(k==1)mpz_add((*next)[ix].get_mpz_t(),(*next)[ix].get_mpz_t(),co.get_mpz_t());else mpz_addmul_ui((*next)[ix].get_mpz_t(),co.get_mpz_t(),k);res.additions++;};
   if(c==0){add(a,b-1,0,2);add(a-1,b,0,2);add(a-1,b-1,1,2);}
   else{add(a+1,b-1,c-1,1);add(a-1,b+1,c-1,1);add(a-1,b-1,c+1,1);add(a,b,c-1,2);add(a,b-1,c,2);add(a-1,b,c,2);}
   co=0;
  }
  swap(cur,next);
 }
 if(res.coeff[0]!=0||res.coeff[1]!=12*res.coeff[2])throw runtime_error("low coefficient identities failed");
 Z check=0,bin;for(int i=0;i<=n;i++){mpz_bin_uiui(bin.get_mpz_t(),2*i,i);check+=bin*res.coeff[i];}if(check!=0)throw runtime_error("central binomial identity failed");
 res.seconds=secs(t0);return res;
}
int main(int argc,char **argv){
 if(argc<3){cerr<<"usage: phi_stream_selected d output.jsonl [threads=2] [nd=standard] [selection=prefix|all]\n";return 2;}
 int d=stoi(argv[1]),threads=argc>3?stoi(argv[3]):2;auto fac=factor(d);if(d<4||fac.size()!=1)throw runtime_error("d must be a prime power >=4");int q=fac[0].first;int nd=argc>4?stoi(argv[4]):d*(q==2?6:q==3?5:4);if(nd<3*d)throw runtime_error("nd<3d");
 string selection=argc>5?argv[5]:"prefix";
 if(selection!="prefix" && selection!="all")throw runtime_error("selection must be prefix or all");
 if(selection=="prefix" && (fac[0].second!=1 || q<5 || nd!=4*d))throw runtime_error("prefix requires odd prime d>=5 and nd=4d");
 auto selected=[&](int n,int m){return selection=="all" || n<4*d || (n==4*d && m<=d);};
 set<pair<int,int>> done;{ifstream old(argv[2]);string line;while(getline(old,line)){if(line.empty())continue;if(line.back()!='}')throw runtime_error("partial last row: repair JSONL before resuming");auto p=line.find("\"n\":");auto r=line.find("\"m\":");if(p==string::npos||r==string::npos)throw runtime_error("unrecognized checkpoint row");done.emplace(stoi(line.substr(p+4)),stoi(line.substr(r+4)));}}
 vector<shared_ptr<Source>> sources;vector<pair<shared_ptr<Source>,int>> jobs;int maxs=0;auto start=Clock::now();
 for(int n=d;n<=nd;n+=d){auto src=make_source(d,n,q);sources.push_back(src);for(int m=0;m<src->mcount;m++)if(selected(n,m) && !done.count({n,m})){jobs.emplace_back(src,m);maxs=max(maxs,2*n+src->gdegree+m);}}
 if(jobs.empty()){cerr<<"All rows already present.\n";return 0;}
 Prefix pref(maxs);ofstream out(argv[2],ios::app);if(!out)throw runtime_error("cannot open output");omp_set_num_threads(threads);
 cerr<<"{\"stage\":\"start\",\"d\":"<<d<<",\"nd\":"<<nd<<",\"jobs\":"<<jobs.size()<<",\"already_completed\":"<<done.size()<<",\"threads\":"<<threads<<",\"max_half_total_degree\":"<<maxs<<",\"max_layer_states\":"<<pref.counts[maxs]<<"}"<<endl;
 #pragma omp parallel
 {
  Workspace work(pref.counts[maxs]);
  #pragma omp for schedule(dynamic,1)
  for(int k=0;k<(int)jobs.size();k++){
   const auto &src=*jobs[k].first;int m=jobs[k].second;Result res=reduce_row(src,m,pref,work);
   ostringstream row;row<<"{\"d\":"<<d<<",\"nd\":"<<nd<<",\"n\":"<<src.n<<",\"m\":"<<m<<",\"coefficients\":[";for(size_t i=0;i<res.coeff.size();i++){if(i)row<<',';row<<'"'<<res.coeff[i]<<'"';}row<<"],\"max_half_total_degree\":"<<res.maxs<<",\"nonzero_states\":"<<res.states<<",\"additions\":"<<res.additions<<",\"seconds\":"<<res.seconds<<"}";
   #pragma omp critical
   {out<<row.str()<<'\n';out.flush();if(!out)throw runtime_error("output write failed");cerr<<"{\"stage\":\"row_complete\",\"d\":"<<d<<",\"n\":"<<src.n<<",\"m\":"<<m<<",\"seconds\":"<<res.seconds<<",\"states\":"<<res.states<<",\"maxrss_kib\":"<<rss()<<",\"total_elapsed\":"<<secs(start)<<"}"<<endl;}
  }
 }
 cerr<<"{\"stage\":\"complete\",\"d\":"<<d<<",\"rows_added\":"<<jobs.size()<<",\"seconds\":"<<secs(start)<<",\"maxrss_kib\":"<<rss()<<"}"<<endl;
}
