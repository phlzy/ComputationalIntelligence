#include "sa.hpp"
#include <bits/stdc++.h>
using namespace std;
ofstream out("data1.txt");
vector<int> point[305];
int point_num = 0;
int dimension = 3;
namespace PSO {
struct particle {
    vector<double> x;
    vector<double> v;
    double pBest;
    vector<double> pBest_pos;
    void init(int _size) {
        uniform_real_distribution<double> rnd(0, 1000);
        random_device rd;
        default_random_engine engine{rd()};
        x.resize(_size);
        for (auto &&i : x)
            i = rnd(engine);
        pBest_pos = x;
        v.resize(_size);
        for (auto &&i : v)
            i = rnd(engine);
    }
};
constexpr double c1 = 1.49445;
constexpr double c2 = 1.49445;
constexpr double w = 0.729;
constexpr double eps = 1e-3;
double gBest;
vector<double> gBest_pos;
vector<double> r1;
vector<double> r2;

vector<particle> particles;
double fitness(vector<double> &v) {
    double sum = 0;
    int sz = v.size();
    for (int i = 1; i <= point_num; ++i) {
        double tmp = 0;
        for (int j = 0; j < sz; ++j) {
            tmp += (point[i][j] - v[j]) * (point[i][j] - v[j]);
        }
        sum = max(sum, tmp);
    }
    return sqrt(sum);
}
void init(int particle_num) {
    gBest = 1e18;
    particles.resize(particle_num);
    r1.resize(dimension);
    r2.resize(dimension);
    uniform_real_distribution<double> rnd(0, 1);
    random_device rd;
    default_random_engine engine{rd()};
    for (auto &&i : r1)
        i = rnd(engine);
    for (auto &&i : r2)
        i = rnd(engine);
    for (auto &&i : particles) {
        i.init(dimension);
        i.pBest = fitness(i.x);
        if (i.pBest < gBest) {
            gBest = i.pBest;
            gBest_pos = i.x;
        }
    }
}
bool run() {
    double delta = 0, tmpgBest = gBest;
    for (auto &&p : particles) {
        for (int d = 0; d < dimension; ++d) {
            p.v[d] = w * p.v[d] + c1 * r1[d] * (p.pBest_pos[d] - p.x[d]) +
                     c2 * r2[d] * (gBest_pos[d] - p.x[d]);
            p.x[d] = p.x[d] + p.v[d];
        }
        double current_fitness = fitness(p.x);
        if (current_fitness < p.pBest) {
            p.pBest = current_fitness;
            p.pBest_pos = p.x;
        }
        if (current_fitness < gBest) {
            delta = max(delta, tmpgBest - current_fitness);
            gBest = current_fitness;
            gBest_pos = p.x;
        }
    }
    return delta < eps;
}
void work(int particle_num) {
    init(particle_num);
    constexpr int Gmax = 400;
    for (int i = 1; i <= Gmax; ++i) {
        bool flag = run();
        if (i > 100 && !flag) {
            break;
        }
    }
    cout << "PSO ans: ";
    cout << gBest << endl;
    out << gBest << " ";
    cout << "position: ";
    for (auto &&i : gBest_pos)
        cout << i << " ";
    cout << endl;
}
}; // namespace PSO
pair<double, vector<double>> work(int n, vector<double> init) {
    ComputationalIntelligence::SimulateAnneal sa;
    function<double(vector<double> &)> calc = [&](vector<double> &v) -> double {
        double ret = 0;
        int len = v.size();
        for (int i = 1; i <= n; ++i) {
            double sum = 0;
            for (int j = 0; j < len; ++j)
                sum += (v[j] - point[i][j]) * (v[j] - point[i][j]);
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
        sa.solve(init, calc, better, change, sa.defaultExponentCooling, 3000);
    cout << "sa ans: ";
    cout << fixed << setprecision(6) << sqrt(ans.first) << endl;
    return ans;
}
void checksa(int n, vector<double> &init) {
    for (auto &i : init)
        i /= n;
    vector<double> ans = init;
    double res = 1e18;
    for (int i = 1; i <= 2; ++i) {
        auto tmp = work(n, ans);
        if (tmp.first < res) {
            res = tmp.first;
            ans = tmp.second;
        }
    }
    cout << "sa ans: ";
    cout << fixed << setprecision(6) << sqrt(res) << endl;
    out << fixed << setprecision(2) << sqrt(res) << " ";
    cout << "position: ";
    for (auto i : ans)
        cout << i << " ";
    cout << endl;
}
class solver {
  public:
    int n;
    const double eps = 1e-6;
    double x[105], y[105], z[105];
    double dis3(double a, double b, double c) {
        double ans = 0;
        for (int i = 1; i <= n; i++) {
            ans = max(ans, (x[i] - a) * (x[i] - a) + (y[i] - b) * (y[i] - b) +
                               (z[i] - c) * (z[i] - c));
        }
        return ans;
    }
    double dis2(double a, double b) {
        double l = -100000;
        double r = 100000;
        while (r - l >= eps) {
            double rmid = (r + l) / 2;
            double lmid = (l + rmid) / 2;
            if (dis3(a, b, lmid) < dis3(a, b, rmid)) {
                r = rmid;
            } else
                l = lmid;
        }
        return dis3(a, b, l);
    }
    double dis(double a) {
        double l = -100000;
        double r = 100000;
        while (r - l >= eps) {
            double rmid = (r + l) / 2;
            double lmid = (l + rmid) / 2;
            if (dis2(a, lmid) < dis2(a, rmid)) {
                r = rmid;
            } else
                l = lmid;
        }
        return dis2(a, l);
    }
    void init(int _n = 100) {
        n = _n;
        for (int i = 1; i <= n; ++i)
            x[i] = point[i][0], y[i] = point[i][1], z[i] = point[i][2];
    }
    double solve() {
        init();
        double l = -100000;
        double r = 100000;
        while (r - l >= eps) {
            double rmid = (r + l) / 2;
            double lmid = (l + rmid) / 2;
            if (dis(lmid) < dis(rmid)) {
                r = rmid;
            } else
                l = lmid;
        }
        return sqrt(dis(l));
    }
};
int main(int argc, char **argv) {
    point_num = 100;
    for (dimension = 3; dimension <= 30; dimension++) {
        cout << dimension << endl;
        out << dimension << " ";
        int n = point_num;
        vector<double> init(dimension);
        mt19937 rd(time(nullptr));
        for (int i = 1; i <= n; ++i) {
            point[i].resize(dimension);
            for (int j = 0; j < dimension; ++j) {
                point[i][j] = rd() % 1000;
                init[j] += point[i][j];
            }
        }
        auto pso_start = chrono::system_clock::now();
        PSO::work(200);
        auto pso_end = chrono::system_clock::now();
        auto duration = (pso_end - pso_start);
        auto t1 = duration.count() / 1000000;
        cout << duration.count() / 1000000 << endl;
        auto sa_start = chrono::system_clock::now();
        checksa(point_num, init);
        auto sa_end = chrono::system_clock::now();
        duration = (sa_end - sa_start);
        auto t2 = duration.count() / 1000000;
        cout << duration.count() / 1000000 << endl;
        out << t1 << " " << t2 << endl;
    }
    // int n = point_num;
    // vector<double> init(dimension, 0);
    // mt19937 rd(time(nullptr));
    // for (int i = 1; i <= n; ++i) {
    //     point[i].resize(dimension);
    //     for (int j = 0; j < dimension; ++j) {
    //         point[i][j] = rd() % 1000;
    //         init[j] += point[i][j];
    //     }
    // }
    // PSO::work(100);
    // checksa(point_num, init);
    // solver s;
    // double ans = s.solve();
    // cout << "standard ans: ";
    // cout << ans << endl;
    return 0;
}