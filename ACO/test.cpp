#include <Eigen/core>
#include <bits/stdc++.h>

// using namespace Eigen;
using namespace std;
using ll = long long;
const int N = 30;
const int inf = 0x3f3f3f3f;
int dis[N][N];

ofstream out("data1.txt");

namespace ACO {
double _t0;
const double _p = 0.9;
const double _alpha = 1;
const double _beta = 3;
struct ant {
    int loc;
    vector<int> vis;
    vector<int> path;
    // bool finished;
    int len;
    void init(int n) {
        loc = 0;
        vis.resize(n);
        fill(vis.begin(), vis.end(), 0);
        vis[0] = true;
        path = {loc};
        len = 0;
    }
};
struct city {
    int id;
    double prob;
    void init() {
        id = -1;
        prob = 0.0;
    }
};
vector<ant> ants;
vector<city> cities;
vector<double> cityProb;
double tot = 0;
double TAU[N][N];
double tmpTAU[N][N];
void init(int antNum, int cityNum) {
    ants.resize(antNum);
    for (auto &cur : ants) {
        cur.init(cityNum);
    }
    cities.resize(cityNum);
    for (auto &cur : cities) {
        cur.init();
    }
    cityProb.resize(cityNum);
    fill(cityProb.begin(), cityProb.end(), 0);
    for (int i = 0; i < cityNum; ++i) {
        for (int j = 0; j < cityNum; ++j) {
            TAU[i][j] = _t0;
            tmpTAU[i][j] = 0;
        }
    }
}

double calcProb(int antId, int cityId) {
    double res = 0;
    int cur = ants[antId].loc;
    if (dis[cur][cityId] < 0 || ants[antId].vis[cityId])
        return res;
    double eta = 1.0 / dis[cur][cityId];
    eta = pow(eta, _beta);
    double tau = pow(TAU[cur][cityId], _alpha);
    res = eta * tau / tot;
    return res;
}
void updateProb(int antId) {
    int cityNum = cities.size();
    fill(cityProb.begin(), cityProb.end(), 0);
    for (int i = 0; i < cityNum; ++i) {
        if (ants[antId].vis[i])
            continue;
        cityProb[i] = calcProb(antId, i);
    }
}
int select(int antId) {
    updateProb(antId);
    double sum = 0;
    int cityNum = cities.size();
    vector<pair<double, int>> tmp;
    for (int i = 0; i < cityNum; ++i) {
        if (ants[antId].vis[i])
            continue;
        sum += cityProb[i];
        tmp.push_back({cityProb[i], i});
    }
    if (tmp.empty())
        return -1;
    for (auto &[p, i] : tmp) {
        p /= sum;
    }
    random_shuffle(tmp.begin(), tmp.end());
    for (int i = 1; i < (int)tmp.size(); ++i) {
        tmp[i].first += tmp[i - 1].first;
    }
    std::uniform_real_distribution<double> P(0, 1.0);
    std::random_device rd;
    std::default_random_engine engine{rd()};
    double now = P(engine);
    for (auto &[p, i] : tmp) {
        if (p >= now)
            return i;
    }
    return tmp.back().second;
}
void updateAnt(int antId, int cityId) {
    ants[antId].len += dis[ants[antId].loc][cityId];
    // 对信息素的贡献
    tmpTAU[ants[antId].loc][cityId] += 1.0 / ants[antId].len;
    // 更新当前位置
    ants[antId].loc = cityId;
    ants[antId].vis[cityId] = 1;
    ants[antId].path.push_back(cityId);
}
void updateTAU() {
    int cityNum = cities.size();
    double sum = 0;
    for (int i = 0; i < cityNum; ++i) {
        for (int j = 0; j < cityNum; ++j) {
            if (i == j)
                continue;
            TAU[i][j] *= _p;
            TAU[i][j] += tmpTAU[i][j];
            tmpTAU[i][j] = 0;
            if (dis[i][j] <= 0)
                continue;
            double tmp = 1.0 / dis[i][j];
            sum += pow(TAU[i][j], _alpha) * pow(tmp, _beta);
        }
    }
    tot = sum;
}
void resetAnts() {
    int n = cities.size();
    std::uniform_int_distribution<int> rnd(0, n - 1);
    std::random_device rd;
    std::default_random_engine engine{rd()};
    for (auto &now : ants) {
        now.len = 0;
        now.loc = rnd(engine);
        now.path = {now.loc};
        fill(now.vis.begin(), now.vis.end(), 0);
        now.vis[now.loc] = 1;
    }
}
pair<int, vector<int>> work() {
    int antNum = ants.size();
    int cityNum = cities.size();
    for (int i = 0; i < cityNum; ++i) {
        for (int j = 0; j < antNum; ++j) {
            int now = select(j);
            if (now < 0)
                continue;
            updateAnt(j, now);
        }
    }
    int ans = inf;
    vector<int> path;
    for (auto &now : ants) {
        auto &u = now.path.front();
        auto &v = now.loc;
        if (now.len + dis[u][v] < ans) {
            ans = now.len + dis[u][v];
            path = now.path;
        }
    }
    updateTAU();
    resetAnts();
    return {ans, path};
}
void solve(int antNum, int cityNum, int epochNum) {
    init(antNum, cityNum);
    int ans = inf;
    vector<int> path;
    for (int i = 1; i <= epochNum; ++i) {
        out << "epoch " << i << " ";
        auto [res, resPath] = work();
        if (res < ans) {
            ans = res;
            path = resPath;
        }
        out << ans << endl;
    }
    out << ans << endl;
    for (auto i : path) {
        out << i << " ";
    }
    out << endl;
}
}; // namespace ACO
void genGraph(int n) {
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
int greedy(int pos, int vis, int num, int len) {
    if (__builtin_popcount(vis) == num) {
        return len + dis[pos][0];
    }
    int to, mn = inf;
    for (int i = 0; i < num; ++i) {
        if (vis >> i & 1)
            continue;
        if (dis[pos][i] < mn) {
            mn = dis[pos][i];
            to = i;
        }
    }
    return greedy(to, vis | (1 << to), num, len + mn);
}
void showGraph(int n) {
    cout << n << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << dis[i][j] << " ";
        }
        cout << endl;
    }
}
int main() {
    int n;
    for (int Case = 5; Case <= 20; ++Case) {
        n = Case;
        out << "city number " << n << endl;
        genGraph(n);
        int greedy_ans = greedy(0, 1, n, 0);
        ACO::_t0 = 1.0 * 40 / greedy_ans;
        int dp_ans = DP(n);
        ACO::solve(40, n, 10);
        out << "dp " << dp_ans << endl;
        out << "greedy " << greedy_ans << endl;
        cout << n << " done" << endl;
    }
    return 0;
}