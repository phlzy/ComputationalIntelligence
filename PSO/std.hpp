#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
using namespace std;
typedef long long ll;
const ll mod=9999973;
const double eps=1e-4;
int n;
double x[105],y[105],z[105];
double dis3(double a,double b,double c){
    double ans=0;
    for(int i=1;i<=n;i++){
        ans=max(ans,(x[i]-a)*(x[i]-a)+(y[i]-b)*(y[i]-b)+(z[i]-c)*(z[i]-c));
    }
    return ans;
}
double dis2(double a,double b){
    double l=-100000;
    double r=100000;
    double ans=0;
    while(r-l>=eps){
        double rmid=(r+l)/2;
        double lmid=(l+rmid)/2;
        if(dis3(a,b,lmid)<dis3(a,b,rmid)){
            r=rmid;
        }
        else l=lmid;
    }
    return dis3(a,b,l);
}
double dis(double a){
    double l=-100000;
    double r=100000;
    while(r-l>=eps){
        double rmid=(r+l)/2;
        double lmid=(l+rmid)/2;
        if(dis2(a,lmid)<dis2(a,rmid)){
            r=rmid;
        }
        else l=lmid;
    }
    return dis2(a,l);
}
int main(){
    scanf("%d",&n);
    for(int i=1;i<=n;i++)scanf("%lf%lf%lf",&x[i],&y[i],&z[i]);
    double l=-100000;
    double r=100000;
    while(r-l>=eps){
        double rmid=(r+l)/2;
        double lmid=(l+rmid)/2;
        if(dis(lmid)<dis(rmid)){
            r=rmid;
        }
        else l=lmid;
    }
    printf("%.8f\n",sqrt(dis(l)));
    return 0;
}
