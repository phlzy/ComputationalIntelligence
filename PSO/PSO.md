# 粒子群算法实验报告

## 实验名称

粒子群算法编程实现

## 算法思想

粒子群优化算法是进化计算的一个分支，是一种模拟自然界的生物活动的随机搜索算法。PSO模拟了自然界鸟群捕食和鱼群捕食的过程。通过群体中的协作寻找到问题的全局最优解。


## 算法设计

先随机初始化若干的粒子，使它们位于随机的位置并具有一定初速度，然后进行若干轮迭代，每轮迭代根据粒子的历史最优解和全局最优解更新粒子的位置与速度，再更新历史最优解与全局最优解，直到答案收敛或迭代轮数达到上限为止。

## 代码设计

伪代码：

![Img](./FILES/PSO.md/img-20220605164025.png)

速度更新规则：

![Img](./FILES/PSO.md/img-20220605164101.png)

位置更新规则：

![Img](./FILES/PSO.md/img-20220605164121.png)

参数设置采用PPT中的推荐值：

``` cpp
constexpr double c1 = 1.49445;
constexpr double c2 = 1.49445;
constexpr double w = 0.729;
```

## 源代码

``` cpp
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
    return 0;
}
```

## 实验测试

本次实验没有像前两次一样选择TSP问题进行测试，而是选择了高维空间最小球覆盖问题进行测试。该问题的形式化描述如下：在 $k$ 维空间内存在 $n$ 个点 $p_1,p_2,\cdots,p_n$，定义函数 $dis(p)$ 表示点 $p$ 到这 $n$ 个点的距离的最大值，求 $dis(p)$ 的最小值。

与01背包问题类似，可以设计出与值域有关的多项式算法，设点的坐标的绝对值 $\in[0,R]$，存在复杂度为 $O(n\log^{k}{R})$ 的算法，即对每个维度进行三分搜索。然而在值域很大或维数很高的情况下，这种做法依然可能消耗非常多的时间，因此使用智能算法进行分析是有必要的。

### 测试数据

测试数据一共 $28$ 组，为 $3$ 到 $30$ 维空间内随机生成 $100$ 个点，对比粒子群算法与模拟退火算法的解的质量与运行效率。

### 测试结果

![](1.png)

![Img](2.png)

### 结果分析

粒子群算法得到的结果在绝大多数情况下都优于模拟退火算法的结果，而其时间开销在大多数情况下都比模拟退火算法少的多，这可能是因为粒子群算法的收敛速度更快。

## 结论

粒子群算法可以高效的找到质量较高的近似解。有些问题虽然存在多项式复杂度的算法可以求出精确解，但当数据规模较大时可能依然不适合在现实中使用，在这种情况下也可以考虑使用粒子群算法之类的随机搜索算法进行求解。