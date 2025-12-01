#ifndef HEAT_EQUATION_INTEGRAL_H
#define HEAT_EQUATION_INTEGRAL_H

#include "simpson_integral.h"
#include <functional>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
class HeatEquationIntegral {
private:
    double x_0, y_0;
    
public:
    HeatEquationIntegral(double source_x = 0.5, double source_y = 0.5) 
        : x_0(source_x), y_0(source_y) {}
    
    // 二维热核
    double heat_kernel(double x, double y, double t) const {
        if (t <= 1e-12) return 0.0;
        double r_2 = (x - x_0) * (x - x_0) + (y - y_0) * (y - y_0);
        return std::exp(-r_2 / (4.0 * t)) / (4.0 * M_PI * t);
    }
    
    // 计算源项时间积分
    double compute_source_time_integral(double x, double y, double t,
                                       std::function<double(double)> source,
                                       int n_intervals = 10000) const {
        auto integrand = [&](double s) -> double {
            return source(s) * heat_kernel(x, y, t - s);
        };
        if(t < 1e-10) return 0.0;
        return SimpsonIntegral::integrate_adaptive(integrand, 0.0, t, 1e-12, 20);
    }
    
    // 计算初值空间卷积
    double compute_initial_convolution(double x, double y, double t,
                                      std::function<double(double, double)> initial_condition,
                                      double domain_scale = 5.0,
                                      int n_points = 1000) const {
        // 根据热核衰减确定积分区域
        double limit = domain_scale * std::sqrt(t);
        double ax = x - limit, bx = x + limit;
        double ay = y - limit, by = y + limit;
        
        auto integrand = [&](double xi, double eta) -> double {
            return heat_kernel(x - xi, y - eta, t) * initial_condition(xi, eta);
        };
        
        return SimpsonIntegral::integrate2d(integrand, ax, bx, ay, by, n_points, n_points);
    }
    
    // 完整解计算
    double compute_solution(double x, double y, double t,
                           std::function<double(double, double)> initial_condition,
                           std::function<double(double)> source) const {
        double initial_part = compute_initial_convolution(x, y, t, initial_condition);
        double source_part = compute_source_time_integral(x, y, t, source);
        return initial_part + source_part;
    }
};

#endif // HeatEquationIntegral_H