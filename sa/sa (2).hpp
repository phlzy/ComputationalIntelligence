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
