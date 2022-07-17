#include "graph_generator.hpp"
#include "hybrid_pso.hpp"
#include "naive_pso.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
vector<int> ans_naive_pso, ans_hybrid_pso;
vector<int> time_naive_pso, time_hybrid_pso;
vector<int> epoch_naive_pso, epoch_hybrid_pso;
/**
 * @brief 运行粒子群算法，统计答案、运行时间、收敛轮数
 * 
 */
int run_naive_pso(vector<vector<int>> &graph) {
    auto t1 = chrono::system_clock::now();
    auto epoch = naive_pso::work(40, graph);
    epoch_naive_pso.push_back(epoch);
    auto t2 = chrono::system_clock::now();
    auto duration = t2 - t1;
    time_naive_pso.push_back(duration.count() / 1000000);
    return naive_pso::print_result();
}
/**
 * @brief 运行混合算法，统计答案、运行时间、收敛轮数
 * 
 */
int run_hybrid_pso(vector<vector<int>> &graph) {
    auto t1 = chrono::system_clock::now();
    auto epoch = hybrid_pso::work(40, graph);
    epoch_hybrid_pso.push_back(epoch);
    auto t2 = chrono::system_clock::now();
    auto duration = t2 - t1;
    time_hybrid_pso.push_back(duration.count() / 1000000);
    return hybrid_pso::print_result();
}
/**
 * @brief 对两种算法进行对比测试
 * @param size 图中顶点的数量
 */
void benchmark(int size) {
    // 生成随机稀疏图
    auto tmp = graph_generator::sparse_graph(size);
    auto graph = graph_generator::convert(tmp);
    ans_naive_pso.push_back(run_naive_pso(graph));
    ans_hybrid_pso.push_back(run_hybrid_pso(graph));
}
/**
 * @brief 测试驱动函数
 * 
 * @param n 图的最大规模
 */
void work(int n) {
    // 输出结果到文件
    ofstream out("data.txt");
    for (int i = 5; i <= n; ++i) {
        cout << "testing " << i << endl;
        benchmark(i);
    }
    cout << "done" << endl;
    out << "naive pso ans" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << ans_naive_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
    out << "hybrid pso ans" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << ans_hybrid_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
    out << "naive pso time" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << time_naive_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
    out << "hybrid pso time" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << time_hybrid_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
    out << "naive pso epoch" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << epoch_naive_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
    out << "hybrid pso epoch" << endl;
    out << "[";
    for (int i = 5; i <= n; ++i) {
        out << epoch_hybrid_pso[i - 5] << ",]"[i == n];
    }
    out << endl;
}
int main() {
    // 测试顶点数量 5-100 的图
    work(100);
    return 0;
}