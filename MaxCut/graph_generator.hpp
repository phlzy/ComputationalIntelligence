#ifndef GRAPH_GENERATOR_HPP
#define GRAPH_GENERATOR_HPP 1

#include "dsu.hpp"
#include <functional>
#include <random>
namespace graph_generator {
/**
 * @brief 生成随机树 (邻接表形式存储)
 */
std::vector<std::vector<int>> random_tree(int _size) {
    assert(_size > 1);
    auto g = std::vector<std::vector<int>>(_size);
    std::vector<int> ord(_size);
    std::iota(ord.begin(), ord.end(), 0);
    std::random_shuffle(ord.begin(), ord.end());
    std::mt19937 rnd(std::random_device{}());
    for (int i = 1; i < _size; ++i) {
        int fa = rnd() % i;
        g[ord[i]].push_back(ord[fa]);
        g[ord[fa]].push_back(ord[i]);
    }
    return g;
}
/**
 * @brief 生成随机图 (邻接表形式存储)
 */
std::vector<std::vector<int>> random_graph(int _size, int extra_edge = 0) {
    assert(_size > 1);
    auto g = random_tree(_size);
    if (_size <= 4)
        return g;
    std::mt19937 rnd(std::random_device{}());
    if (extra_edge <= 0) {
        assert(_size <= 10000);
        int MAX_EDGE = _size * _size;
        extra_edge = rnd() % MAX_EDGE + 1;
    }
    basic_ds::dsu aux(_size);
    for (int i = 1; i <= extra_edge; ++i) {
        int u = rnd() % _size;
        int v = rnd() % _size;
        if (u == v)
            v = (v + 1) % _size;
        while (aux.same(u, v)) {
            bool lucky = rnd() % _size;
            if (lucky)
                break;
            u = rnd() % _size;
            v = rnd() % _size;
            if (u == v)
                v = (v + 1) % _size;
        }
        aux.merge(u, v);
        g[u].push_back(v), g[v].push_back(u);
    }
    for (int i = 0; i < _size; ++i) {
        std::sort(g[i].begin(), g[i].end());
        g[i].erase(std::unique(g[i].begin(), g[i].end()), g[i].end());
    }
    return g;
}
/**
 * @brief 生成完全图 (邻接表形式存储)
 */
std::vector<std::vector<int>> complete_graph(int _size) {
    assert(_size > 1);
    auto g = std::vector<std::vector<int>>(_size);
    for (int i = 0; i < _size; ++i) {
        for (int j = 0; j < _size; ++j) {
            if (i == j)
                continue;
            g[i].push_back(j);
        }
    }
    return g;
}
/**
 * @brief 生成随机稀疏图 (邻接表形式存储)
 */
std::vector<std::vector<int>> sparse_graph(int _size) {
    assert(_size > 1);
    if (_size <= 4)
        return random_tree(_size);
    const int SQRT_BOUND = _size * (int)sqrt(_size) - _size;
    const int LOG_BOUND = _size * (int)log2(_size) - _size;
    std::mt19937 rnd(std::random_device{}());
    int extra_edge = abs(SQRT_BOUND - LOG_BOUND) + 1;
    extra_edge = rnd() % extra_edge + LOG_BOUND + 1;
    return random_graph(_size, extra_edge);
}
/**
 * @brief 生成随机二分图 (邻接表形式存储)
 */
std::vector<std::vector<int>> binary_graph(int _size) {
    auto g = random_tree(_size);
    std::vector<int> left, right;
    std::function<void(int, int, int)> dfs = [&](int u, int fa,
                                                 int col) -> void {
        if (col)
            left.push_back(u);
        else
            right.push_back(u);
        for (auto v : g[u]) {
            if (v == fa)
                continue;
            dfs(v, u, !col);
        }
    };
    dfs(0, -1, 1);
    std::mt19937 rnd(std::random_device{}());
    const int LEFT_SIZE = left.size();
    const int RIGHT_SIZE = right.size();
    int extra_edge = LEFT_SIZE * RIGHT_SIZE;
    extra_edge = rnd() % extra_edge + 1;
    for (int i = 1; i <= extra_edge; ++i) {
        int u = rnd() % LEFT_SIZE;
        int v = rnd() % RIGHT_SIZE;
        g[left[u]].push_back(right[v]);
        g[right[v]].push_back(left[u]);
    }
    for (int i = 0; i < _size; ++i) {
        std::sort(g[i].begin(), g[i].end());
        g[i].erase(std::unique(g[i].begin(), g[i].end()), g[i].end());
    }
    return g;
}
/**
 * @brief 将邻接表转化为邻接矩阵
 */
std::vector<std::vector<int>> convert(std::vector<std::vector<int>> &graph) {
    int n = graph.size();
    auto res = std::vector<std::vector<int>>(n, std::vector<int>(n, 0));
    for (int u = 0; u < n; ++u) {
        for (auto &&v : graph[u])
            res[u][v] = 1;
    }
    return res;
}
}; // namespace graph_generator

#endif