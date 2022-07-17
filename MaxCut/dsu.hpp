#ifndef DSU_HPP
#define DSU_HPP 1

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

namespace basic_ds {
/**
 * @brief 并查集 (路径压缩+按秩合并)
 * @anchor 李喆元
 */
class dsu {
  protected:
    std::vector<int> fa;
    std::vector<int> rank;
    int component_cnt, size;

  public:
    void init() {
        assert(size > 0);
        component_cnt = size;
        fa.resize(size);
        std::iota(fa.begin(), fa.end(), 0);
        rank.resize(size);
        std::fill(rank.begin(), rank.end(), 1);
    }
    dsu(const int &_size) {
        assert(_size > 0);
        size = _size;
        init();
    }
    int get_size() const {
        assert(size > 0);
        return size;
    }
    bool is_connected() const {
        assert(component_cnt > 0);
        return component_cnt == 1;
    }
    int find(int u) {
        assert(0 <= u && u < size);
        if (u == fa[u])
            return u;
        return fa[u] = find(fa[u]);
    }
    /**
     * @brief 检查两个点是否在一个连通分量内
     */
    bool same(int u, int v) {
        assert(0 <= u && u < size);
        assert(0 <= v && v < size);
        return find(u) == find(v);
    }
    /**
     * @brief 合并两个连通分量
     * @return 合并后连通分量的根节点
     */
    int merge(int u, int v) {
        assert(0 <= u && u < size);
        assert(0 <= v && v < size);
        u = find(u), v = find(v);
        if (u == v)
            return u;
        if (rank[u] > rank[v])
            std::swap(u, v);
        rank[v] += rank[u];
        fa[u] = v;
        return v;
    }
};
}; // namespace basic_ds

#endif