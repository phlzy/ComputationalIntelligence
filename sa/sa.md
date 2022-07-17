# 模拟退火算法实验报告

## 实验名称

编程实现模拟退火算法

## 算法思想

模拟退火算法（Simulated Annealing，SA）是一种模拟物理退火的过程而设计的优化算法。

模拟退火算法采用类似于物理退火的过程，先在一个高温状态下（相当于算法随机搜索），然后逐渐退火，在每个温度下（相当于算法的每一次状态转移）徐徐冷却（相当于算法局部搜索），最终达到物理基态（相当于算法找到最优解）。

模拟退火算法可以用于为 NP 复杂性问题提供有效近似算法，其优点在于可以克服优化过程陷入局部极小、克服初值依赖性。


## 算法设计

模拟退火算法从某一高温出发，在高温状态下计算初始解，然后以预设的邻域函数产生一个扰动量，从而得到新的状态，即模拟粒子的无序运动。

比较新旧状态下的能量，即目标函数的解。如果新状态的能量小于旧状态，则状态发生转化；如果新状态的能量大于旧状态，则以一定的概率准则发生转化（状态的接受和舍弃由 Metropolis 准则决定）。

当状态稳定后，便可以看作达到当前状态的最优解，便可以开始降温，在下一个温度继续迭代，最终达到低温的稳定状态，便得到了模拟退火算法产生的结果。


## 代码设计

为了便于在多种问题下使用模拟退火算法，我选择使用 C++ 的面向对象与泛型的机制实现这一算法，因此代码的可读性可能有一定程度的下降。

首先需要根据问题设置

- 初始温度 $T$
- 冷却调度函数 $\text{CoolingSchedule}$
- 初始状态 $\text{initial}$
- 能量计算函数 $\text{calc}$
- 能量比较函数 $\text{better}$
- 随机变化函数 $\text{change}$

然后就可以得到一个伪代码：

```
current := initial // 初始化状态
ans := calc(current) // 初始化能量
while T > lowerbound_of_T: // 温度不低于下界
    now := change(current) // 粒子随机运动，产生新状态
    value := calc(now) // 计算新状态的能量
    delta := value - ans
    if value is better than ans: // 能量变低了
        ans := value, current := now // 更新状态
    else if Metropolis(-delta / T): // 根据 Metropolis 准则接受状态
        current := now
    T := CoolingSchedule(T) // 降温
```

将其实现即可。

## 源代码

``` cpp
#include <algorithm>
#include <functional>
#include <random>
#include <utility>
#include <vector>
namespace ComputationalIntelligence {
class SimulateAnneal {
  private:
    double eps;
    double inf;
    int itercnt;

  public:
    SimulateAnneal() {
        eps = 1e-12;
        inf = 2e12;
        itercnt = 0;
    }
    double getEps() const { return eps; }
    double getInf() const { return inf; }
    int getItercnt() const { return itercnt; }
    void setEps(const double &_eps) { eps = _eps; }
    void setInf(const double &_inf) { inf = _inf; }

    std::function<double(double)> defaultExponentCooling =
        [](double temp) -> double { return temp * 0.998; };
    std::function<double(double)> defaultLinearCooling =
        [](double temp) -> double { return temp - 0.168; };

    std::function<double(double)> defaultInverseCooling =
        [](double temp) -> double { return temp / (0.0012 * temp + 1); };
    std::function<double(double)> defaultQuickCooling =
        [&](double temp) -> double { return temp / (getItercnt() + 1); };
    std::function<double(double)> defaultLogarithmicCooling =
        [&](double temp) -> double { return 3000.0 / (std::log1p(getItercnt()) + 1); };

    template <typename T>
    std::pair<T, std::vector<T>>
    solve(std::vector<T> &initial,
          const std::function<T(std::vector<T> &)> &calc,
          const std::function<bool(T, T)> &better,
          const std::function<std::vector<T>(std::vector<T> &, double)> &change,
          const std::function<double(double)> &CoolingSchedule,
          double temperature) {
        std::vector<T> current = initial;
        std::vector<T> ret = initial;
        T ans = calc(current);
        itercnt = 0;
        std::uniform_real_distribution<double> P(0, 1.0);
        std::random_device rd;
        std::default_random_engine engine{rd()};
        while (temperature > eps) {
            std::vector<T> now = change(current, temperature);
            T value = calc(now);
            double delta = -fabs(value - ans);
            if (better(value, ans)) {
                current = now;
                ret = now;
                ans = value;
            } else if (exp(delta / temperature) > P(engine)) {
                current = now;
            }
            temperature = CoolingSchedule(temperature);
            itercnt++;
        }
        return {ans, ret};
    }
};
}; // namespace ComputationalIntelligence

```

TSP 问题求解

``` cpp
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
    ofstream out("data.txt");
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
```

也可以用于求解其他的问题，如计算任意维空间内到所有点最大距离最小的点的位置

``` cpp
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
    cout << "position: ";
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
```

## 实验测试（包括测试数据、测试结果、结果分析等）

### 测试数据

本次实验以 TSP 问题为例，对模拟退火算法的性能进行测试。测试方法为先生成一定规模的随机图，然后用状态压缩动态规划算法求出最优解，与模拟退火算法得到的解进行比对。由于这里使用的状态压缩动态规划算法的时间复杂度为 $O(n^2 2^n)$（其中 $n$ 为城市数量），所以测试数据规模有限，$n$ 的范围是 $[5,20]$。

### 测试结果

测试结果以图表的形式呈现：

![](1.png)

上图比较了模拟退火算法得到的近似解与最优解的差距，可以看到在 $n\le 10$ 的时候得到的近似解就是最优解，而 $n$ 增加到一定程度后，几乎已经不能找到最优解了。

![](2.png)

上图列出了模拟退火算法得到的近似解与最优解的差值，由此可见随着 $n$ 的增加，近似解与最优解的差距可以认为是逐渐变大的。

![](3.png)

上图统计了模拟退火算法的收敛速度。可以发现绝大多数时候经过一轮退火（即温度从初温降低到最低温度一次）答案就会收敛，迭代最多的一次也只进行了 $4$ 轮就收敛了。这说明模拟退火算法的收敛速度是比较快的。

### 结果分析

本次实验采用了指数冷却的方式，即 $T_k=\alpha ^k T_0$，产生新状态的过程采用了根据温度高低随机改变区间内点的访问顺序的方法。

从结果上来看，虽然算法的收敛速度较快，但与最优解之间的误差还是很大的，说明很多参数都需要调整。

如果对冷却调度方式和产生新状态的过程进行优化，可能会得到更好的结果，但受时间所限我尚未尝试这些参数的调整。

## 结论

模拟退火算法可以高效的获得 NP-Hard 问题的比较优秀的可行解，但与最优解之间仍然存在一定的差距，如果需要找到更好的近似解，需要对迭代次数、冷却调度、随机状态生成方式等参数进行调整。