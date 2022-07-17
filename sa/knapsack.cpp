#include "sa.hpp"
#include <ctime>
#include <iostream>
#include <iomanip>
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
    cout << fixed << setprecision(2) << sqrt(ans.first) << endl;
    for (auto i : ans.second)
        cout << i << endl;
}
int main(int argc, char **argv) {
    int num = 100, capicity = 10000;
    
    return 0;
}