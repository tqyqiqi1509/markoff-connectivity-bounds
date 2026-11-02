// Independent-direction modular scalar Phi evaluation for odd prime d.
// Compile: g++ -O3 -std=c++17 phi_forward_scalar.cpp -o phi_forward_scalar
// No GMP, no polynomial cache and no backward row propagation are used.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sys/resource.h>
#include <utility>
#include <vector>
using namespace std;
using V=array<uint32_t,3>;
using Clock=chrono::steady_clock;
double elapsed(Clock::time_point t){return chrono::duration<double>(Clock::now()-t).count();}
uint64_t power(uint64_t a,uint64_t b,uint64_t p){uint64_t r=1;while(b){if(b&1)r=r*a%p;a=a*a%p;b>>=1;}return r;}
size_t triangle_index(int a,int b,int c){array<int,3> e{a,b,c};sort(e.begin(),e.end(),greater<int>());return size_t(e[1])*(e[1]+1)/2+e[2];}
bool isprime(uint64_t p){if(p<2)return false;for(uint64_t a=2;a*a<=p;a++)if(p%a==0)return false;return true;}
int main(int argc,char**argv){
 if(argc!=2){cerr<<"usage: phi_forward_scalar input.txt\n";return 2;}
 ifstream in(argv[1]);int d,count;in>>d>>count;if(!isprime(d)||d<5)throw runtime_error("odd prime parameter required");
 array<uint64_t,3> q,t;for(int h=0;h<3;h++){in>>q[h]>>t[h];if(!isprime(q[h])||q[h]<=8*d)throw runtime_error("invalid scalar modulus");}
 map<pair<int,int>,V> expected;for(int i=0;i<count;i++){int n,m;V v;in>>n>>m>>v[0]>>v[1]>>v[2];if(!in)throw runtime_error("input truncation");if(!expected.emplace(make_pair(n,m),v).second)throw runtime_error("duplicate row");}
 const int S=10*d-1,side=S+1;size_t width=size_t(S/2+2)*(S/2+3)/2;vector<V> before(width),current(width),pairs(size_t(side)*side);vector<V> tp(side);for(int h=0;h<3;h++)tp[0][h]=1;
 for(int s=1;s<=S;s++)for(int h=0;h<3;h++)tp[s][h]=uint64_t(tp[s-1][h])*t[h]%q[h];
 auto start=Clock::now();uint64_t states=0;
 for(int s=0;s<=S;s++){
  for(int b=0;b<=s/2;b++)for(int c=0;c<=min(b,s-2*b);c++){
   int a=s-b-c;size_t ix=size_t(b)*(b+1)/2+c;V value{};
   if(b==0)value=tp[s];
   else for(int h=0;h<3;h++){
    auto get=[&](int aa,int bb,int cc)->uint64_t{return before[triangle_index(aa,bb,cc)][h];};uint64_t v;
    if(c==0)v=2*(get(a,b-1,0)+get(a-1,b,0)+get(a-1,b-1,1));
    else v=get(a+1,b-1,c-1)+get(a-1,b+1,c-1)+get(a-1,b-1,c+1)+2*(get(a,b,c-1)+get(a,b-1,c)+get(a-1,b,c));
    value[h]=v%q[h];
   }
   current[ix]=value;if(c==0)pairs[size_t(a)*side+b]=value;states++;
  }
  before.swap(current);
  if(s%200==0)cerr<<"{\"stage\":\"forward_layer\",\"half_degree\":"<<s<<",\"states\":"<<states<<",\"seconds\":"<<elapsed(start)<<"}"<<endl;
 }
 double forward_time=elapsed(start);vector<vector<int>> filters{{1},{0,1},{-1,1},{0,-2,1}};
 uint64_t terms=0,checks=0;int failures=0;
 for(int k=1;k<=4;k++){
  int n=k*d,mc=k*(d-1)/2;const auto &G=filters[k-1];vector<vector<V>> coefficients(n+1);
  for(int i=0;i<=n;i++)coefficients[i].resize(i+G.size());
  for(int h=0;h<3;h++){
   uint64_t p=q[h];vector<uint64_t> fact(2*n+1),ifact(2*n+1),neg4(n+1);fact[0]=1;neg4[0]=1;
   for(int i=1;i<=2*n;i++)fact[i]=fact[i-1]*i%p;ifact[2*n]=power(fact[2*n],p-2,p);for(int i=2*n;i>=1;i--)ifact[i-1]=ifact[i]*i%p;for(int i=1;i<=n;i++)neg4[i]=neg4[i-1]*(p-4)%p;
   auto choose=[&](int a,int b){return fact[a]*ifact[b]%p*ifact[a-b]%p;};
   for(int i=0;i<=n;i++){
    uint64_t c=choose(n+i,n-i)*(2*n)%p*power(n+i,p-2,p)%p;if(i&1)c=(p-c)%p;
    for(int j=0;j<=i;j++){
     uint64_t co=c*choose(i,j)%p*neg4[i-j]%p;
     for(size_t a=0;a<G.size();a++)if(G[a]){uint64_t z=G[a]>0?co*G[a]%p:(p-co*uint64_t(-G[a])%p)%p;auto &out=coefficients[i][j+a][h];out=(out+z)%p;}
    }
   }
  }
  for(int m=0;m<mc;m++){
   array<uint64_t,3> answer{};
   for(int i=0;i<=n;i++)for(size_t j=0;j<coefficients[i].size();j++){
    int a=n-i+int(j)+m,b=i;if(a<b)swap(a,b);const V &v=pairs[size_t(a)*side+b];const V &co=coefficients[i][j];
    for(int h=0;h<3;h++)answer[h]=(answer[h]+uint64_t(co[h])*v[h])%q[h];terms++;
   }
   auto it=expected.find({n,m});if(it==expected.end())throw runtime_error("expected row missing");
   for(int h=0;h<3;h++){checks++;if(answer[h]!=it->second[h]){failures++;cerr<<"{\"mismatch_n\":"<<n<<",\"m\":"<<m<<",\"point\":"<<h<<",\"actual\":"<<answer[h]<<",\"expected\":"<<it->second[h]<<"}"<<endl;}}
  }
  cerr<<"{\"stage\":\"source_block_checked\",\"n\":"<<n<<",\"m_count\":"<<mc<<",\"checks\":"<<checks<<",\"failures\":"<<failures<<",\"seconds\":"<<elapsed(start)<<"}"<<endl;
 }
 rusage use;getrusage(RUSAGE_SELF,&use);
 cout<<"{\"d\":"<<d<<",\"rows\":"<<count<<",\"scalar_points\":3,\"checks\":"<<checks<<",\"failures\":"<<failures<<",\"max_half_degree\":"<<S<<",\"forward_states\":"<<states<<",\"source_terms_combined\":"<<terms<<",\"forward_seconds\":"<<forward_time<<",\"total_seconds\":"<<elapsed(start)<<",\"maxrss_kib\":"<<use.ru_maxrss<<",\"method\":\"independent forward scalar modular Phi plus direct GF_n reconstruction; diagnostic, not coefficient proof\"}"<<endl;
 return failures?1:0;
}
