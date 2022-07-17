# 混合优化算法的实现与应用实验报告

## 实验名称

混合优化算法的实现与应用

## 算法思想

之前的实验表明，粒子群算法的收敛速度很快，但是也容易陷入局部最优解，经常在很少的迭代次数后就迅速收敛到某个局部最优解而不再发生变化了。因此，若在粒子群算法快速收敛后通过随机搜索获得一个更优秀的解，然后再采用粒子群算法在其邻域进行搜索，就可以使粒子群算法避免陷入之前的局部最优解。

然而，由于问题的解空间很大，普通的随机搜索算法难以在短时间内得到更好的解，因此采用模拟退火算法代替，以当前获得的局部最优解为初始状态进行搜索。模拟退火算法在温度较高时搜索范围广泛，在温度较低时则会在已经找到的最优解附近进行搜索，具有较高的搜索效率。


## 算法设计

最大割问题的数学描述如下：给定一个无向图 $G=(V,E)$，其中 $V$ 是图 $G$ 的点集，$E$ 是图 $G$ 的边集。设点集 $S$ 是点集 $V$ 的一个真子集，则其关于 $V$ 的补集 $T=V-S$。将边集 $E$ 中两个端点分别位于 $S,T$ 两个集合中的边的集合称为点集 $S$ 的割 $cut(S)$，即：

$$
cut(S)=\sum_{u\in S,v\in T,e_{u,v}\in E}e_{u,v}
$$

最大割问题是求一个点集 $S$ 使得其割集 $cut(S)$ 中的边的数量最大，即最大化 $|cut(S)|$。

运用粒子群算法求解该问题。粒子的位置对应向量中每个维度的值，每个粒子的位置对应一个可行解，速度表示从一个解转移到另一个解的倾向。

每个粒子的适应度直接使用其划分方案对应的割集的大小来表示。对于一个可行解，根据其每一位的值为 $0$ 还是 $1$ 可以将其下标分为两个集合，通过枚举两个集合中的顶点的所有可能组合统计割集的大小。显然，最后得到的适应度最高的粒子对应的划分方案就是求解出的最大割。

粒子群算法在随机初始化所有粒子的位置与速度后进行若干轮迭代，每轮迭代都会更新每个粒子的速度和位置，并更新每个粒子的历史最优解以及全局最优解，直到迭代次数达到上限或者答案收敛为止。在最大割问题中，每个粒子的速度更新公式如下：
$$
V_i=\omega V_i +c_1 r_1 * (pBest_i-X_i) +c_2 r_2 * (gBest_i-X_i)
$$

粒子的新位置 $X_i'=X_i+V_i$，由于速度是一个实数，此时得到的 $X_i'$ 并不能保证正好是 $0$ 或 $1$，因此根据 $X_i'$ 在数轴上距离 $0$ 和 $1$ 的距离，选择距离更近的数，将其修正为 $0$ 或 $1$。

粒子群算法的答案收敛后，使用模拟退火算法在当前答案的基础上进行新一轮的搜索。

在模拟退火算法中状态的编码方式与之前完全相同，在随机扰动时根据温度的高低调整状态改变的概率。首先对温度进行归一化，设初始温度为 $T$，当前温度为 $t$，则某一维状态发生变化的概率为 $\log_{T+1}(t+1)$。随着温度的降低，变化的概率也会逐渐降低。温度每发生一次改变，就对当前状态产生一次扰动，并计算新状态是否更优。若新状态更优，则更新当前状态为新状态，否则根据Metropolis准则接受新状态。

当粒子群算法的迭代次数达到上限，或得到的解经过若干轮模拟退火处理与粒子群算法的迭代流程依然不再发生改变时，说明该算法求得的解已经收敛，不再进行迭代。

## 代码设计

伪代码：

``` 
// 模拟退火算法
function simulate_annealing
    current_state = gBest
    while temprature > epsilon
        new_state = random_disturb(current_state)
        delta = fit(new_state) - fit(gBest)
        if delta > 0
            current_state = new_state
        else if random < exp(delta / temprature)
            current_state = new_state
        temprature = cooler(temprature)
    end while
    if fit(current_state) > fit(gBest)
        gBest = current_state
end function

// 粒子群算法
function PSO
    for each particle i
        initialize velocity Vi andposition Xi
        evaluate particle i and set pBesti= Xi
    end for
    gBest = max{pBest} 
    while not stop 
        for each particle i
            update Vi and Xi
            evaluate particle i
            if fit(Xi) > fit(pBesti)
                pBesti = Xi
            if fit(pBesti) > fit(gBest)
                gBest = pBesti
        end for
        if gBest convergence
            simulate_annealing
    end while
    result = gBest
end function
```

## 源代码

```cpp
#ifndef HYBRID_PSO_HPP
#define HYBRID_PSO_HPP 1
#include "sa.hpp"
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
namespace hybrid_pso {
class particle {
  public:
    std::vector<double> velo;
    std::vector<int> pos;
    int p_best;
    std::vector<int> p_best_pos;
    /**
     * @brief 初始化单个粒子
     *
     * @param size 图的顶点数量
     */
    void current_best(int size) {
        std::random_device rd;
        std::mt19937 rnd(rd());
        velo.resize(size);
        pos.resize(size);
        for (auto &&i : pos)
            i = (rnd() & 1);
        std::default_random_engine engine{rd()};
        std::uniform_real_distribution<double> drand(-1, 1);
        p_best_pos = pos;
        for (auto &&i : velo)
            i = drand(engine);
    }
};
constexpr double c1 = 1.49445;
constexpr double c2 = 1.49445;
constexpr double w = 0.729;
constexpr double eps = 1e-3;
int g_best;
std::vector<int> g_best_pos;
std::vector<double> r1;
std::vector<double> r2;
std::vector<std::vector<int>> graph;
std::vector<particle> particles;
/**
 * @brief 计算割集大小
 *
 * @return int 割集大小
 */
int fitness(std::vector<int> &pos) {
    // assert(graph.size() == pos.size());
    int n = pos.size();
    std::vector<int> left, right;
    int cnt = 0;
    for (int u = 0; u < n; ++u) {
        if (pos[u])
            left.push_back(u);
        else
            right.push_back(u);
    }
    if (left.empty() || right.empty())
        return -1;
    for (auto &&l : left) {
        for (auto &&r : right) {
            // assert(0 <= l && l < n);
            // assert(0 <= r && r < n);
            if (graph[l][r])
                cnt++;
        }
    }
    return cnt;
}
/**
 * @brief 初始化粒子群
 * 
 * @param particle_num 粒子数量
 * @param _graph 用于计算的图 (邻接矩阵表示)
 */
void init_pso(int particle_num, const std::vector<std::vector<int>> &_graph) {
    g_best = 0;
    particles.resize(particle_num);
    graph = _graph;
    int n = graph.size();
    assert(n > 1);
    // 随机初始化 r1, r2
    r1.resize(n), r2.resize(n);
    std::uniform_real_distribution<double> rnd(0, 1);
    std::random_device rd;
    std::default_random_engine engine{rd()};
    for (auto &&i : r1)
        i = rnd(engine);
    for (auto &&i : r2)
        i = rnd(engine);
    // 随机初始化粒子
    for (auto &&p : particles) {
        p.current_best(n);
        p.p_best = fitness(p.pos);
        if (g_best < p.p_best) {
            g_best = p.p_best;
            g_best_pos = p.pos;
        }
    }
}
/**
 * @brief 模拟退火算法
 * 
 */
void run_sa() {
    int n = graph.size();
    assert(n > 1);
    constexpr double EPS = 1e-12;
    constexpr double dT = 0.998;
    constexpr double MAX_TEMPERATURE = 3000;
    const double LOG_TEMPERATURE = log2(MAX_TEMPERATURE + 1);
    double T = MAX_TEMPERATURE;
    auto current_best = g_best_pos;
    auto res = g_best_pos;
    int max_cut = fitness(current_best);
    // 归一化函数
    std::function<double(double)> normalize = [&](const double &t) -> double {
        return log2(t + 1) / LOG_TEMPERATURE;
    };
    std::uniform_real_distribution<double> P(0, 1.0);
    std::random_device rd;
    std::default_random_engine engine{rd()};
    // 模拟退火开始迭代
    while (T > EPS) {
        auto tmp = current_best;
        double bound = normalize(T);
        // 随机扰动
        for (int i = 0; i < n; ++i) {
            if (P(engine) < bound)
                tmp[i] = !tmp[i];
        }
        int tmp_cut = fitness(tmp);
        double delta = tmp_cut - max_cut;
        // 接受更优解或根据 Metropolis 准则接受
        if (delta > 0) {
            current_best = tmp;
            max_cut = tmp_cut;
            res = tmp;
        } else if (exp(delta / T) > P(engine)) {
            current_best = tmp;
        }
        T *= dT;
    }
    // 模拟退火算法搜索到更优的解
    if (max_cut > g_best) {
        std::cout << "sa find a better solution!" << std::endl;
        g_best = max_cut;
        g_best_pos = res;
    }
}
/**
 * @brief 粒子群算法一次迭代
 * 
 */
void run_pso() {
    int n = graph.size();
    // 对每一个粒子进行更新
    for (auto &&p : particles) {
        // 对每个维度的速度进行更新
        for (int i = 0; i < n; ++i) {
            p.velo[i] = p.velo[i] * w +
                        c1 * r1[i] * (p.p_best_pos[i] - p.pos[i]) +
                        c2 * r2[i] * (g_best_pos[i] - p.pos[i]);
            double tmp_pos = p.pos[i] + p.velo[i];
            double dis0 = fabs(tmp_pos);
            double dis1 = fabs(tmp_pos - 1);
            p.pos[i] = (dis0 < dis1 ? 0 : 1);
        }
        // 计算适应度
        int current_fitness = fitness(p.pos);
        // 更新 pBest 与 gBest
        if (current_fitness > p.p_best) {
            p.p_best = current_fitness;
            p.p_best_pos = p.pos;
        }
        if (current_fitness > g_best) {
            g_best = current_fitness;
            g_best_pos = p.pos;
        }
    }
}
/**
 * @brief 混合算法驱动函数
 * 
 * @param particle_num 粒子数量
 * @param _graph 用于计算的图
 * @return int 达到收敛的迭代次数
 */
int work(int particle_num, const std::vector<std::vector<int>> &_graph) {
    // 初始化粒子群
    init_pso(particle_num, _graph);
    constexpr int Gmax = 400;
    int cnt = 0, tmp_g_best = -1, res = 0;
    // 开始迭代
    for (int i = 1; i <= Gmax; ++i) {
        run_pso();
        // 若粒子群算法收敛，开始模拟退火
        if (cnt >= 4)
            run_sa();
        if (g_best == tmp_g_best) {
            // 若模拟退火得到的答案也收敛
            cnt++;
        } else {
            // 模拟退火搜索到更优解，继续迭代
            tmp_g_best = g_best;
            cnt = 0;
            res = i;
        }
        // 答案收敛
        if (cnt >= 7)
            return res;
    }
    return Gmax;
}
/**
 * @brief 查看最后得到的解
 */
int print_result() {
    int n = graph.size();
    if (n < 2)
        return -1;
    std::vector<int> left, right;
    for (int u = 0; u < n; ++u) {
        if (g_best_pos[u])
            left.push_back(u);
        else
            right.push_back(u);
    }
    std::cout << "hybrid pso" << std::endl;
    // std::cout << "vertices set S:" << std::endl;
    // for (auto &v : left)
    //     std::cout << v << " ";
    // std::cout << std::endl;
    // std::cout << "vertices set T:" << std::endl;
    // for (auto &v : right)
    //     std::cout << v << " ";
    // std::cout << std::endl;
    std::cout << "max cut: " << g_best << std::endl;
    return g_best;
}
}; // namespace hybrid_pso

#endif // HYBRID_PSO_HPP
```

## 实验测试

本次实验生成了顶点数量 $5-100$ 的随机图与稀疏图用于测试，对比混合算法与粒子群算法在各方面的差异。

### 测试数据

测试数据一共 $192$ 组，为顶点数量 $5-100$ 的随机图与稀疏图，对比粒子群算法与混合要优化算法的解的质量与运行效率。

### 测试结果

#### 解的质量对比：

![Img](new/rand1.png)

![Img](new/rand2.png)

![Img](new/sp1.png)

#### 时间消耗对比：

![Img](new/rand3.png)

![Img](new/sp2.png)

#### 收敛轮数对比：

![Img](new/rand4.png)

![Img](new/sp3.png)

### 结果分析

混合粒子群模拟退火算法的收敛速度与解的质量朴素粒子群算法相比均有明显提升，求得的最大割的割集大小大约有 $10\%$ 的提高，但是时间消耗也是粒子群算法的 $45$ 倍左右。


## 结论

混合了模拟退火算法的粒子群算法能得到比朴素粒子群算法更优秀的解，收敛速度也比粒子群算法更快。然而，由于需要多次进行模拟退火，该算法需要消耗相较于粒子群算法数十倍的时间。如果能够将参数调整好，混合优化算法的效率有望进一步提升。