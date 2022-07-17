#include <bits/stdc++.h>
using namespace std;

namespace t {
class foo {
  public:
    template <typename T>
    bool work(vector<T> vec, const function<bool(T, T)> &cmp,
              const function<vector<T>(vector<T> &, T)> &fuck) {
        double a = 114514, b = 1919810;
        return cmp(a, b);
    }
};
}; // namespace t

int main() {
    t::foo bar;
    function<bool(double, double)> cmp = [](double a, double b) -> bool {
        return a < b;
    };
    function<vector<double>(vector<double> &, double)> fuck =
        [](vector<double> &a, double b) -> vector<double> {
        auto v = a;
        for (auto &i : v)
            i += b;
        return v;
    };
    vector<double> v;
    auto ans = bar.work(v, cmp, fuck);
    cout << ans << endl;
}