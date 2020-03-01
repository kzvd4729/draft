#include<bits/stdc++.h>
#define long long long
using namespace std;
const int N=2e5;
//..............................................................................
/*
stirling number of seocond kind:
there are n different object(distinguishable) and k identical boxes
(indistinguishable). how many ways to put object into boxes such that each box
have atleast one object.
boxes are considered to be set. permutation of object doesnt matter.
if boxes are distinguishable multiply by k!
           1   k        
formula:  ---  ∑   (-1)^i kCi (k-i)^n
           k!  i=0
removing the constant and reformulating:
(-1)^i 1/(i!)  *  1/(k-i)! (k-i)^n
for fixed n we can calculte for all k using fft
*/
const int mod=880803841,g=26;vector<int>r;
void bitReverse(int n)
{
  r.resize(n);for(int i=0;i<n;i++)r[i]=0;int p=0;while((1<<p)<n)p++;
  for(int i=0;i<n;i++){for(int j=0;j<p;j++)if(i&(1<<j))r[i]|=(1<<(p-j-1));}
}
int big(int b,int p)
{
  int ret=1;
  while(p){if(p%2)ret=(1LL*ret*b)%mod;b=(1LL*b*b)%mod;p/=2;}return ret;
}
void dft(vector<int>&a,bool inv)
{
  int n=a.size();
  for(int i=0;i<n;i++)if(r[i]<i)swap(a[i],a[r[i]]);
  for(int ln=2;ln<=n;ln*=2)
  {
    int m=inv?big(g,mod-1-(mod-1)/ln):big(g,(mod-1)/ln);
    for(int i=0;i<n;i+=ln)
    {
      int r=1,u,v;
      for(int j=0;j<ln/2;j++)
      {
        u=a[i+j],v=(1LL*r*a[i+j+ln/2])%mod;a[i+j]=u+v<mod?u+v:u+v-mod;
        a[i+j+ln/2]=u-v>=0?u-v:u-v+mod;r=(1LL*r*m)%mod;
      }
    }
  }
  if(inv){int ni=big(n,mod-2);
  for(int i=0;i<a.size();i++)a[i]=(1LL*a[i]*ni)%mod;}
}
vector<int>multiply(vector<int>a,vector<int>b)
{
  int sz=a.size()+b.size();int n=1;while(n<sz)n<<=1;bitReverse(n);
  a.resize(n);b.resize(n);dft(a,false);dft(b,false);
  for(int i=0;i<n;i++)a[i]=(1LL*a[i]*b[i])%mod;dft(a,true);return a;
}
int fc[N+2];
vector<int>solve(int n,int k)
{
  fc[0]=1;
  for(int i=1;i<=k;i++)fc[i]=(1LL*fc[i-1]*i)%mod;
  vector<int>a(k+1),b(k+1);int sgn=1;
  for(int i=0;i<=k;i++)
  {
    a[i]=(1LL*sgn*big(fc[i],mod-2))%mod;
    b[i]=(1LL*big(fc[i],mod-2)*big(i,n))%mod;
    sgn*=-1;
  }
  vector<int>c=multiply(a,b);//c[i] stores stirling(n,i);
  return c;
}
int seg[5*N+2],lazy[5*N+2];
void build(int node,int lo,int hi)
{
  if(lo==hi)
  {
    seg[node]=1;return ;
  }
  int md=(lo+hi)/2;
  build(node*2,lo,md);build(node*2+1,md+1,hi);
  seg[node]=seg[node*2]+seg[node*2+1];
}
void tooLazy(int node,int lo,int hi)
{
  if(!lazy[node])return ;
  seg[node]=hi-lo+1-seg[node];
  if(lo!=hi)
  {
    lazy[node*2]^=1;lazy[node*2+1]^=1;
  }
  lazy[node]=0;
}
void upd(int node,int lo,int hi,int lt,int rt)
{
  tooLazy(node,lo,hi);
  if(lo>rt||hi<lt)return ;
  if(lo>=lt&&hi<=rt)
  {
    lazy[node]^=1;tooLazy(node,lo,hi);
    return ;
  }
  int md=(lo+hi)/2;
  upd(node*2,lo,md,lt,rt);upd(node*2+1,md+1,hi,lt,rt);
  seg[node]=seg[node*2]+seg[node*2+1];
}
int main()
{
  ios_base::sync_with_stdio(0);cin.tie(0);
  int t;cin>>t;
  while(t--)
  {
    int n,k,d;cin>>n>>k>>d;
    vector<int>c=solve(n,k);

    memset(lazy,0,sizeof(lazy));

    build(1,1,k);
    while(d--)
    {
      int l,r;cin>>l>>r;upd(1,1,k,l,r);
      cout<<c[seg[1]]<<"\n";
    }
  }
  return 0;
}
