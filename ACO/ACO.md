# 蚁群算法实验报告

## 实验名称

蚁群算法编程实现

## 算法思想

蚁群在寻找食物时，总能找到一条蚁穴到食物的最短路径，并且能随着环境的变化而更新最优路径。究其原因，是因为蚂蚁在运动过程中，会在所经过路径上留下一种称为信息素的物质，其他蚂蚁在运动中可以感知到这种物质，并以此来指导自己的运动方向。蚁群的这种行为表现出一种信息正反馈现象：某一路径上走过的蚂蚁越多，则后来者选择该路径的概率越大。

蚁群算法是一种基于群体智能的算法，用来在图中寻找优化路径的概率型。它由 Marco Dorigo 于 1992 年在他的博士论文中提出，其灵感来源就是上文中蚂蚁觅食的现象。蚁群算法被应用于解决大多数优化问题或者能够转化为优化求解的问题，已渗透到各个应用领域，展现出优异的性能和广阔的发展前景。


## 算法设计

蚂蚁在寻找食物的过程中往往是随机选择路径的，但它们能感知当前地面上的信息素浓度，并倾向于往信息素浓度高的方向行进。信息素由蚂蚁自身释放，是实现蚁群内间接通信的物质。由于较短路径上蚂蚁的往返时间比较短，单位时间内经过该路径的蚂蚁多，所以信息素的积累速度比较长路径快。因此，当后续蚂蚁在路口时，就能感知先前蚂蚁留下的信息，并倾向于选择一条较短的路径前行。这种正反馈机制使得越来越多的蚂蚁在巢穴与食物之间的最短路径上行进。由于其他路径上的信息素会随着时间蒸发，最终所有的蚂蚁都在最优路径上行进。

蚁群的觅食空间相当于问题的搜索空间，蚁群相当于一组有效解，蚁巢到食物的路径相当于一个有效解，找到的最短路径相当于最优解。

每只蚂蚁都随机选择一个城市作为其出发城市，并维护一个路径记忆向量，用来存放该蚂蚁依次经过的城市。蚂蚁在构建路径的每一步中，按照一个随机比例规则选择下一个要到达的城市。当所有蚂蚁构建完路径后，算法将会对所有的路径进行全局信息素的更新。若干轮迭代之后就可以找到最优解。

![Img](./FILES/ACO.md/img-20220529162858.png)


## 代码设计

这里以使用蚁群算法解决TSP问题为例：先用贪心算法求出一组解，计算初始情况下的信息素含量 $\tau_0$，然后将所有的边上的信息素都被初始化为 $\tau_0$。

接下来进行若干轮迭代，初始化每只蚂蚁的位置与禁忌表等信息，然后每只蚂蚁都根据受信息素浓度影响的概率来随机选择前进的方向，直到所有蚂蚁都构造出一个路径为止。在所有蚂蚁构造完路径之后，更新信息素浓度，所有路径上的信息素都会发生蒸发，我们为所有边上的信息素乘上一个小于 $1$ 的常数，然后蚂蚁根据自己构建的路径长度在它们本轮经过的边上释放信息素。蚂蚁构建的路径越短、释放的信息素就越多。一条边被蚂蚁爬过的次数越多，它所获得的信息素也越多。同时，每只蚂蚁的路径都对应一组解，一轮迭代结束后更新当前最优解。

迭代的过程结束后所获得的当前最优解就是蚁群算法找到的最优解。

## 源代码

蚁群算法求解TSP问题代码：

``` cpp
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
```

## 实验测试（包括测试数据、测试结果、结果分析等）

测试方式为随机生成点数 $n\in [5,20]$ 的图，对比蚁群算法的解和状态压缩动态规划求出的最优解之间的差距。

![Img](1.png)

下图为蚁群算法的收敛速度：

![Img](2.png)

可见蚁群算法找到的解与最优解十分接近，且收敛速度较快。

## 结论

蚁群算法可以高效的找出TSP问题的一组较为优秀的近似解，并在该问题上表现出比模拟退火算法更优秀的性能。