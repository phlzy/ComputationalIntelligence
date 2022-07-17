#include "sa.hpp"
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;
vector<int> p[305];
void work(int n, vector<double> init) {
    ComputationalIntelligence::SimulateAnneal sa;
    function<double(vector<double> &)> calc = [&](vector<double> &v) -> double {
        double ret = 0;
        int len = v.size();
        for (int i = 1; i <= n; ++i) {
            double sum = 0;
            for (int j = 0; j < len; ++j)
                sum += (v[j] - p[i][j]) * (v[j] - p[i][j]);
            ret = max(ret, sum);
        }
        return ret;
    };
    function<vector<double>(vector<double> &, double)> change =
        [](vector<double> &init, double T) -> vector<double> {
        mt19937 rd(time(nullptr));
        uniform_int_distribution<int> rnd(0, 200000);
        const double rmax = 1e5;
        vector<double> ret = init;
        for (auto &i : ret)
            i += (rnd(rd) - rmax) * T;
        return ret;
    };

    function<bool(double, double)> better = [](double x, double y) -> bool {
        return x < y;
    };

    pair<double, vector<double>> ans =
        sa.solve(init, calc, better, change, sa.defaultQuickCooling, 3000);
    cout << "radius: ";
    cout << fixed << setprecision(2) << sqrt(ans.first) << endl;
    cout << "centre of sphere: ";
    for (auto i : ans.second)
        cout << i << " ";
    cout << endl;
}
int main(int argc, char **argv) {
    cout << "input: dimension, point number" << endl;
    int n, m;
    cin >> n >> m;
    vector<double> init(m, 0);
    mt19937 rd(time(nullptr));
    for (int i = 1; i <= n; ++i) {
        p[i].resize(m);
        cout << "point " << i << endl;
        for (int j = 0; j < m; ++j) {
            // cin >> p[i][j];
            p[i][j] = rd() % 1000;
            cout << p[i][j] << " ";
            init[j] += p[i][j];
        }
        cout << endl;
    }
    for (auto &i : init)
        i /= n;
    work(n, init);
    return 0;
}