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