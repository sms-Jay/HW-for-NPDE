#ifndef SIMPSON_INTEGRAL_H
#define SIMPSON_INTEGRAL_H

#include <functional>
#include <cmath>
#include <stdexcept>

class SimpsonIntegral {
public:

    template<typename Func>
    static double integrate1d(Func f, double a, double b, int n = 1000) {
        if (n <= 0) throw std::invalid_argument("Number of intervals must be positive");
        if (a >= b) throw std::invalid_argument("Lower limit must be less than upper limit");
        

        if (n % 2 != 0) n++;
        
        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            if (i % 2 == 1) {
                sum += 4.0 * f(x);
            } else {
                sum += 2.0 * f(x);
            }
        }
        
        return sum * h / 3.0;
    }
    
    template<typename Func>
    static double integrate2d(Func f, double ax, double bx, 
                             double ay, double by, 
                             int nx = 100, int ny = 100) {
        if (nx <= 0 || ny <= 0) 
            throw std::invalid_argument("Number of intervals must be positive");
        if (ax >= bx || ay >= by) 
            throw std::invalid_argument("Lower limits must be less than upper limits");
        
        if (nx % 2 != 0) nx++;
        if (ny % 2 != 0) ny++;
        
        double hx = (bx - ax) / nx;
        double hy = (by - ay) / ny;
        double sum = 0.0;
        
        for (int i = 0; i <= nx; i++) {
            double x = ax + i * hx;
            double wx = get_simpson_weight(i, nx);
            
            for (int j = 0; j <= ny; j++) {
                double y = ay + j * hy;
                double wy = get_simpson_weight(j, ny);
                
                sum += wx * wy * f(x, y);
            }
        }
        
        return sum * hx * hy / 9.0;
    }
    
    template<typename Func>
    static double integrate_adaptive(Func f, double a, double b, 
                                   double tol = 1e-8, int max_depth = 20) {
        return adaptive_simpson_recursive(f, a, b, f(a), f(b), 
                                         f((a + b) / 2), tol, max_depth);
    }

private:
    static double get_simpson_weight(int index, int total_points) {
        if (index == 0 || index == total_points) {
            return 1.0;
        } else if (index % 2 == 1) {
            return 4.0;
        } else {
            return 2.0;
        }
    }
    
    template<typename Func>
    static double adaptive_simpson_recursive(Func f, double a, double b,
                                           double fa, double fb, double fc,
                                           double tol, int depth) {
        if (depth <= 0) {
            return (fa + 4 * fc + fb) * (b - a) / 6.0;
        }
        
        double c = (a + b) / 2;
        double d = (a + c) / 2;
        double e = (c + b) / 2;
        
        double fd = f(d);
        double fe = f(e);
        
        double left_simpson = (fa + 4 * fd + fc) * (c - a) / 6.0;
        double right_simpson = (fc + 4 * fe + fb) * (b - c) / 6.0;
        double whole_simpson = (fa + 4 * fc + fb) * (b - a) / 6.0;
        
        double error = std::abs(left_simpson + right_simpson - whole_simpson);
        
        if (error < 15 * tol) {
            return left_simpson + right_simpson + (left_simpson + right_simpson - whole_simpson) / 15.0;
        }
        
        return adaptive_simpson_recursive(f, a, c, fa, fc, fd, tol / 2, depth - 1) +
               adaptive_simpson_recursive(f, c, b, fc, fb, fe, tol / 2, depth - 1);
    }
};

#endif // SIMPSON_INTEGRAL_H