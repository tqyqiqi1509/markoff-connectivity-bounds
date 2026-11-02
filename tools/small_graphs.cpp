#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
// Enumerate unordered triples, then all three Vieta edges in the S_3 quotient.
// The mathematical bridge from this quotient to pure Vieta is in the paper.
int main() {
  std::cout << "p,vertices_S3,visited_S3,connected\n";
  for(int p=5;p<=1304;++p) {
    bool prime=true; for(int d=2;d*d<=p;++d) if(p%d==0) prime=false;
    if(!prime) continue;
    std::vector<std::vector<int>> roots(p);
    for(int z=0;z<p;++z) roots[z*z%p].push_back(z);
    std::vector<std::array<int,3>> pts;
    std::unordered_map<uint64_t,int> ids;
    auto key=[p](std::array<int,3> a){std::sort(a.begin(),a.end());return (uint64_t(a[0])*p+a[1])*p+a[2];};
    int inv2=(p+1)/2;
    for(int x=0;x<p;++x) for(int y=x;y<p;++y) {
      int b=x*y%p, c=(x*x+y*y)%p;
      int disc=(int64_t(b)*b-4LL*c)%p; if(disc<0)disc+=p;
      for(int s:roots[disc]) {
        int z=int64_t(b+s)*inv2%p;
        if(z<y || (x==0&&y==0&&z==0))continue;
        std::array<int,3> a={x,y,z};
        if(ids.count(key(a))) return 2;
        ids[key(a)]=int(pts.size());pts.push_back(a);
      }
    }
    auto start=ids.find(key({3%p,3%p,3%p}));if(start==ids.end())return 3;
    std::vector<char> seen(pts.size());std::queue<int> todo;
    todo.push(start->second);seen[start->second]=1;size_t count=0;
    while(!todo.empty()) {
      int i=todo.front();todo.pop();++count;
      for(int j=0;j<3;++j) {
        auto a=pts[i];a[j]=(int64_t(a[(j+1)%3])*a[(j+2)%3]-a[j]+p)%p;
        auto it=ids.find(key(a));if(it==ids.end())return 4;
        if(!seen[it->second]){seen[it->second]=1;todo.push(it->second);}
      }
    }
    bool ok=count==pts.size();std::cout<<p<<','<<pts.size()<<','<<count<<','<<(ok?1:0)<<'\n';
    if(!ok)return 5;
  }
}
