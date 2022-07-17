#include "sa.hpp"
#include <cassert>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;
using ll = long long;
int dis[105][105];
const int inf = 0x3f3f3f3f;
int work(int n, vector<ll> &init) {
    ComputationalIntelligence::SimulateAnneal sa;
    function<ll(vector<ll> &)> calc = [&](vector<ll> &v) -> ll {
        ll ret = 0, n = v.size();
        for (int i = 0; i < n; ++i) {
            ret += dis[v[i]][v[(i + 1) % n]];
        }
        return ret;
    };
    auto ans0 = calc(init);
    cout << "initial status " << ans0 << endl;
    function<vector<ll>(vector<ll> &, double)> change =
        [](vector<ll> &init, double T) -> vector<ll> {
        vector<ll> ret = init;
        int n = ret.size();
        mt19937 rd(time(nullptr));
        uniform_int_distribution<int> rnd(0, n - 1);
        int shuffle = round(T);
        shuffle = min(shuffle, n);
        for (int i = 0; i <= shuffle; ++i) {
            int a = rnd(rd), b = rnd(rd);
            while (a == b)
                b = rnd(rd);
            if (a > b)
                swap(a, b);
            random_shuffle(ret.begin() + a, ret.begin() + b);
            // swap(ret[a], ret[b]);
        }
        return ret;
    };
    function<bool(ll, ll)> better = [](ll x, ll y) -> bool { return x < y; };
    pair<ll, vector<ll>> ans = sa.solve(
        init, calc, better, change, sa.defaultExponentCooling, (double)3000);
    cout << "result " << ans.first << endl;
    // for (auto i : ans.second)
    //     cout << i << " ";
    // cout << endl;
    init = ans.second;
    return ans.first;
}
void randGraph(int n) {
    mt19937 rnd(time(nullptr));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dis[i][j] = rnd() % 1000 + n;
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dis[i][j] = min(dis[i][j], dis[j][i]);
        }
    }
    for (int i = 0; i < n; ++i)
        dis[i][i] = 0;
}
int dp[21][(1 << 21)];
int DP(int n) {
    int tot = (1 << n) - 1;
    memset(dp, 0x3f, sizeof(dp));
    dp[0][1] = 0;
    for (int st = 0; st <= tot; ++st) {
        for (int i = 0; i < n; ++i) {
            if ((st >> i & 1) == 0)
                continue;
            for (int j = 0; j < n; ++j) {
                if (st >> j & 1)
                    continue;
                int now = st | (1 << j);
                dp[j][now] = min(dp[j][now], dp[i][st] + dis[i][j]);
            }
        }
    }
    int ans = inf;
    for (int i = 1; i < n; ++i) {
        ans = min(ans, dp[i][tot] + dis[0][i]);
    }
    return ans;
}
void show(int n) {
    cout << n << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << dis[i][j] << " ";
        }
        cout << endl;
    }
}
int main(int argc, char **argv) {
    int n;
    // 10-25
    ofstream out("data3.txt");
    for (int cas = 5; cas <= 20; ++cas) {
        // cout << "input city number" << endl;
        // cin >> n;
        n = cas;
        cout << "city number " << cas << endl;
        randGraph(n);
        vector<ll> init(n);
        iota(init.begin(), init.end(), 0);
        random_shuffle(init.begin(), init.end());
        int sa_ans = inf;
        for (int i = 1; i <= 5; ++i) {
            cout << "sa time " << i << endl;
            sa_ans = work(n, init);
            out << sa_ans << " ";
        }
        out << endl;
        int dp_ans = DP(n);
        cout << "dp " << ' ' << dp_ans << endl;
        // show(n);
        out << cas << " " << sa_ans << " " << dp_ans << endl;
    }
    return 0;
}