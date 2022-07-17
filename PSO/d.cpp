#include <bits/stdc++.h>
using namespace std;

int main(){
    ifstream in("data1.txt");
    vector<double> p, s;
    vector<int> t1, t2;
    double b, c;
    int a, d, e;
    while(in >> a >> b >> c >> d >> e) {
        p.push_back(b);
        s.push_back(c);
        t1.push_back(d);
        t2.push_back(e);
    }
    ofstream out("datap.txt");
    out << "[";
    for (int i = 0; i < 28; ++i) {
        out << p[i] << ",]"[i == 27];
    }
    out << endl;
    out << "[";
    for (int i = 0; i < 28; ++i) {
        out << s[i] << ",]"[i == 27];
    }
    out << endl;
    out << "[";
    for (int i = 0; i < 28; ++i) {
        out << t1[i] << ",]"[i == 27];
    }
    out << endl;
    out << "[";
    for (int i = 0; i < 28; ++i) {
        out << t2[i] << ",]"[i == 27];
    }
    out << endl;
}