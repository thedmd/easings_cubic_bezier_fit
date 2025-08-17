#include <stdio.h>
#include <limits>
#include <utility>
#include <vector>
#include <span>
#include <algorithm>
#include <array>
#include <memory>

#include "FitCurves.hpp"

namespace math {

inline constexpr double PI         = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
inline constexpr double PI_DIV_TWO = PI / 2.0;

inline double simple_quad_in(double t)
{
    return t * t;
}

inline double simple_quad_out(double t)
{
    return -t * (t - 2.0);
}

inline double simple_quad_inout(double t)
{
    t *= 2.0;
    if (t < 1.0) return 1.0 / 2.0 * t * t;
    t -= 1.0;
    return -0.5 * (t * (t - 2.0) - 1.0);
}

inline double simple_quad_outin(double t)
{
    if (t < 0.5)
        return simple_quad_inout(0.5 + t) - 0.5;
    else
        return simple_quad_inout(t - 0.5) + 0.5;
}

inline double simple_sin_in(double t)
{
    return 1.0 - cos(t * PI_DIV_TWO);
}

inline double simple_sin_out(double t)
{
    return sin(t * PI_DIV_TWO);
}

inline double simple_sin_inout(double t)
{
    return -0.5 * (cos(PI * t) - 1.0);
}

inline double simple_sin_outin(double t)
{
    if (t < 0.5)
        return simple_sin_inout(0.5 + t) - 0.5;
    else
        return simple_sin_inout(t - 0.5) + 0.5;
}

inline double simple_exp_in(double t)
{
    return pow(2.0, 10.0 * (t - 1.0));
}

inline double simple_exp_out(double t)
{
    return -pow(2.0, -10.0 * t) + 1.0;
}

inline double simple_exp_inout(double t)
{
    t *= 2.f;
    if (t < 1.0) return 0.5 * pow(2.0, 10.0 * (t - 1.0));
    t -= 1.0;
    return 0.5 * (-pow(2.0, -10.0 * t) + 2.0);
}

inline double simple_exp_outin(double t)
{
    if (t < 0.5)
        return simple_exp_inout(0.5 + t) - 0.5;
    else
        return simple_exp_inout(t - 0.5) + 0.5;
}

inline double simple_circ_in(double t)
{
    return -(sqrt(1.0 - t * t) - 1.0);
}

inline double simple_circ_out(double t)
{
    t -= 1.0;
    return sqrt(1.0 - t * t);
}

inline double simple_circ_inout(double t)
{
    t *= 2.0;
    if (t < 1.0) return -0.5 * (sqrt(1.0 - t * t) - 1.0);
    t -= 2.0;
    return 0.5 * (sqrt(1.0 - t * t) + 1.0);
}

inline double simple_circ_outin(double t)
{
    if (t < 0.5)
        return simple_circ_inout(0.5 + t) - 0.5;
    else
        return simple_circ_inout(t - 0.5) + 0.5;
}

inline double simple_cube_in(double t)
{
    return t * t * t;
}

inline double simple_cube_out(double t)
{
    t -= 1.0;
    return t * t * t + 1.0;
}

inline double simple_cube_inout(double t)
{
    t *= 2.0;
    if (t < 1.0) return 0.5 * t * t * t;
    t -= 2.0;
    return 0.5 * (t * t * t + 2.0);
}

inline double simple_cube_outin(double t)
{
    if (t < 0.5)
        return simple_cube_inout(0.5 + t) - 0.5;
    else
        return simple_cube_inout(t - 0.5) + 0.5;
}

inline double simple_quartic_in(double t)
{
    return t * t * t * t;
}

inline double simple_quartic_out(double t)
{
    t -= 1.0;
    return 1.0 - t * t * t * t;
}

inline double simple_quartic_inout(double t)
{
    t *= 2.f;
    if (t < 1.0) return 0.5 * t * t * t * t;
    t -= 2.0;
    return -0.5 * (t * t * t * t - 2.0);
}

inline double simple_quartic_outin(double t)
{
    if (t < 0.5)
        return simple_quartic_inout(0.5 + t) - 0.5;
    else
        return simple_quartic_inout(t - 0.5) + 0.5;
}

inline double simple_quintic_in(double t)
{
    return t * t * t * t * t;
}

inline double simple_quintic_out(double t)
{
    t -= 1.0;
    return t * t * t * t * t + 1.0;
}

inline double simple_quintic_inout(double t)
{
    t *= 2.f;
    if (t < 1.0) return 0.5 * t * t * t * t * t;
    t -= 2.0;
    return 0.5 * (t * t * t * t * t + 2.0);
}

inline double simple_quintic_outin(double t)
{
    if (t < 0.5)
        return simple_quintic_inout(0.5 + t) - 0.5;
    else
        return simple_quintic_inout(t - 0.5) + 0.5;
}

inline double simple_back_in(double t)
{
    const double c1 = 1.70158;
    const double c3 = c1 + 1.0;

    return c3 * t * t * t - c1 * t * t;
}

inline double simple_back_out(double t)
{
    const double c1 = 1.70158;
    const double c3 = c1 + 1.0;

    return 1.0 + c3 * pow(t - 1.0, 3.0) + c1 * pow(t - 1.0, 2.0);
}

inline double simple_back_inout(double t)
{
    const double c1 = 1.70158;
    const double c2 = c1 * 1.525;

    if (t < 0.5)
        return (pow(2.0 * t, 2.0) * ((c2 + 1.0) * 2.0 * t - c2)) / 2.0;
    else
        return (pow(2.0 * t - 2.0, 2.0) * ((c2 + 1.0) * (t * 2.0 - 2.0) + c2) + 2.0) / 2.0;
}

inline double simple_back_outin(double t)
{
    if (t < 0.5)
        return simple_back_inout(0.5 + t) - 0.5;
    else
        return simple_back_inout(t - 0.5) + 0.5;
}

inline double simple_elastic_in(double t)
{
    if (t <= 0.0)
        return 0.0;
    else if (t >= 1.0)
        return 1.0;

    const double c4 = (2.0 * PI) / 3.0;

    return -pow(2.0, 10.0 * t - 10.0) * sin((t * 10.0 - 10.75) * c4);
}

inline double simple_elastic_out(double t)
{
    if (t <= 0.0)
        return 0.0;
    else if (t >= 1.0)
        return 1.0;

    const double c4 = (2.0 * PI) / 3.0;

    return pow(2.0, -10.0 * t) * sin((t * 10.0 - 0.75) * c4) + 1.0;
}

inline double simple_elastic_inout(double t)
{
    if (t <= 0.0)
        return 0.0;
    else if (t >= 1.0)
        return 1.0;

    const double c5 = (2.0 * PI) / 4.5;

    if (t < 0.5)
        return -(pow(2.0, 20.0 * t - 10.0) * sin((20.0 * t - 11.125) * c5)) / 2.0;
    else
        return (pow(2.0, -20.0 * t + 10.0) * sin((20.0 * t - 11.125) * c5)) / 2.0 + 1.0;
}

inline double simple_elastic_outin(double t)
{
    if (t < 0.5)
        return simple_elastic_inout(0.5 + t) - 0.5;
    else
        return simple_elastic_inout(t - 0.5) + 0.5;
}

inline double simple_bounce_out(double t)
{
    const double n1 = 7.5625;
    const double d1 = 2.75;

    if (t < 1.0 / d1)
    {
        return n1 * t * t;
    }
    else if (t < 2.0 / d1)
    {
        t -= 1.5 / d1;
        return n1 * t * t + 0.75;
    }
    else if (t < 2.5 / d1)
    {
        t -= 2.25 / d1;
        return n1 * t * t + 0.9375;
    }
    else
    {
        t -= 2.625 / d1;
        return n1 * t * t + 0.984375;
    }
}

inline double simple_bounce_in(double t)
{
    return 1.0 - simple_bounce_out(1.0 - t);
}

inline double simple_bounce_inout(double t)
{
    if (t < 0.5)
        return (1.0 - simple_bounce_out(1.0 - 2.0 * t)) / 2.0;
    else
        return (1.0 + simple_bounce_out(2.0 * t - 1.0)) / 2.0;
}

inline double simple_bounce_outin(double t)
{
    if (t < 0.5)
        return simple_bounce_inout(0.5 + t) - 0.5;
    else
        return simple_bounce_inout(t - 0.5) + 0.5;
}

} // namespace math

using source_func_t = double(*)(double);

struct cubic_bezier_t { Point2 CP1, CP2; };

template <size_t N>
using sampled_curve_t = std::array<double, N>;

constexpr auto linear_sample(std::span<const Point2> points, double t) -> double
{
    if (points.empty())
        return 0.0;

    if (t <= points.front().x)
        return points.front().y;
    else if (t >= points.back().x)
        return points.back().y;

    auto it = std::lower_bound(points.begin(), points.end(), t, [](const auto& point, auto t) { return point.x < t; });
    if (it == points.end())
        return points.back().y;
    else if (it == points.begin())
        return points.front().x;

    const auto [t0, v0] = *(it - 1);
    const auto [t1, v1] = *it;

    const auto segment_duration = t1 - t0;

    auto blendTime = segment_duration ? (t - t0) / segment_duration : t0;

    return v0 + (v1 - v0) * blendTime;
}

template <typename T, T p1 = T{0}, T p2 = T{1}>
constexpr auto bezier(T cp1, T cp2, T t) -> T
{
    const auto a = T{1} - t;
    const auto b = a * a * a;
    const auto c = t * t * t;

    return b * p1 + T{3} * t * a * a * cp1 + T{3} * t * t * a * cp2 + c * p2;
}

template <typename T, T p1 = T{0}, T p2 = T{1}>
constexpr auto bezierDt(T cp1, T cp2, T t) -> T
{
    const auto a = T{1} - t;
    const auto b = a * a;
    const auto c = t * t;
    const auto d = T{2} * t * a;

    return T{-3} * p1 * b + T{3} * cp1 * (b - d) + T{3} * cp2 * (d - c) + T{3} * c * p2;
}

template <double p1x = 0.0, double p1y = 0.0, double p2x = 1.0, double p2y = 1.0>
constexpr auto sample_bezier_at_x(const cubic_bezier_t& curve, double linear_t) -> double
{
#if 1
    constexpr auto newton_raphson_iteration_limit = 10;
    constexpr auto newton_raphson_max_error       = 1e-6f;

    auto curve_t = static_cast<double>(linear_t);
    for (int i = 0; i < newton_raphson_iteration_limit; ++i)
    {
        const auto x  = bezier  <double, p1x, p2x>(curve.CP1.x, curve.CP2.x, curve_t);
        const auto dx = bezierDt<double, p1x, p2x>(curve.CP1.x, curve.CP2.x, curve_t);
        if (dx == 0.0f) [[unlikely]] // avoid division by zero
            break;
        const auto diff = x - linear_t;
        if (fabs(diff) < newton_raphson_max_error)
            break;
        curve_t -= diff / dx;

        if (curve_t < 0.0f) [[unlikely]]
        {
            curve_t = 0.0f;
            break;
        }
        else if (curve_t > 1.0f) [[unlikely]]
        {
            curve_t = 1.0f;
            break;
        }
    }
#else
    Point2 mapping[1001];
    for (int i = 0; i < 1001; ++i)
    {
        const auto t = static_cast<double>(i) / 1000.0;
        mapping[i] = { t, bezier<double, 0.0, 1.0>(curve.CP1.x, curve.CP2.x, t) };
    }

    auto curve_t = linear_sample(mapping, linear_t);

#endif

    return bezier<double, p1y, p2y>(curve.CP1.y, curve.CP2.y, curve_t);
}

template <size_t N, typename F>
    requires std::is_invocable_r_v<double, F, double>
constexpr void sample_curve_into(sampled_curve_t<N>& result, F&& f)
{
    double step = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i)
        result[i] = f(i * step);
}

template <size_t N, double p1x = 0.0, double p1y = 0.0, double p2x = 1.0, double p2y = 1.0>
constexpr void sample_curve_into(sampled_curve_t<N>& result, const cubic_bezier_t& curve)
{
    double step = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i)
        result[i] = sample_bezier_at_x<p1x, p1y, p2x, p2y>(curve, i * step);
}

template <size_t N>
constexpr auto fit_error(sampled_curve_t<N>& reference, sampled_curve_t<N>& curve) -> double
{
    double error = 0.0;
    double step = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i)
    {
        const auto t = i * step;
        const auto curveValue     = curve[i];
        const auto referenceValue = reference[i];
        const auto diff = static_cast<double>(curveValue - referenceValue);
        error += diff * diff;
    }
    return error;
}

template <size_t N, double p1x = 0.0, double p1y = 0.0, double p2x = 1.0, double p2y = 1.0>
constexpr auto fit_error(sampled_curve_t<N>& reference, const cubic_bezier_t& curve) -> double
{
    sampled_curve_t<N> curveSamples;
    sample_curve_into<N, p1x, p1y, p2x, p2y>(curveSamples, curve);
    return fit_error(reference, curveSamples);
}

template <size_t N>
constexpr auto fit_error(source_func_t referenceFunction, const cubic_bezier_t& curve) -> double
{
    sampled_curve_t<N> reference;
    sample_curve_into(reference, referenceFunction);

    sampled_curve_t<N> curveSamples;
    sample_curve_into(curveSamples, curve);
    return fit_error(reference, curveSamples);
}

struct function_t
{
    const char*     enumName;
    const char*     name;
    source_func_t   func;
    cubic_bezier_t  curve;
    double          error = fit_error<1000>(func, curve);
};

template <size_t N = MAXPOINTS, typename F>
    requires std::is_invocable_r_v<double, F, double>
auto graphics_gems_fit_curve(F&& f, double max_error) -> cubic_bezier_t
{
    Point2 data[N];
    double step = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i)
    {
        data[i].x = i * step;
        data[i].y = f(data[i].x);
    }

    Point2 curve[4];
    FitCurve(
        data,
        static_cast<int>(N),
        max_error,
        [](void* ctx, int n, const BezierCurve curve)
        {
            auto* result = static_cast<Point2*>(ctx);
            result[0] = curve[0];
            result[1] = curve[1];
            result[2] = curve[2];
            result[3] = curve[3];
        },
        curve
    );

    cubic_bezier_t result;
    result.CP1 = curve[1];
    result.CP2 = curve[2];
    return result;
}

template <size_t N = 100>
void graphics_gems_fit_curve_many(std::span<function_t> functions, double max_error)
{
    for (auto& function : functions)
    {
        auto curve = graphics_gems_fit_curve<N>(function.func, max_error);
        auto error = fit_error<N>(function.func, function.curve);

        if (error < function.error)
        {
            function.curve = curve;
            function.error = error;
        }
    }
}

template <size_t N = 100>
void brute_force_fit_curve_many(std::span<function_t> functions, int NX, int NY)
{
    std::vector<sampled_curve_t<N>> references(functions.size());

#pragma omp parallel for
    for (int i = 0; i < std::ssize(functions); ++i)
        sample_curve_into(references[i], functions[i].func);

    const auto x_min   = 0.0;
    const auto x_max   = 1.0;
    const auto x_units = NX;
    const auto x_step  = (x_max - x_min) / x_units;

    const auto y_min   = -1.0;
    const auto y_max   =  2.0;
    const auto y_units = NY;
    const auto y_step  = (y_max - y_min) / y_units;

#pragma omp parallel for
    for (int x1 = 0; x1 < x_units; ++x1)
    {
#pragma omp parallel for
        for (int y1 = 0; y1 < y_units; ++y1)
        {
//#pragma omp parallel for
            for (int x2 = 0; x2 < x_units; ++x2)
            {
//#pragma omp parallel for
                for (int y2 = 0; y2 < y_units; ++y2)
                {
                    cubic_bezier_t curve;
                    curve.CP1.x = x_min + x1 * x_step;
                    curve.CP1.y = y_min + y1 * y_step;
                    curve.CP2.x = x_min + x2 * x_step;
                    curve.CP2.y = y_min + y2 * y_step;

                    sampled_curve_t<N> curveValues;
                    sample_curve_into(curveValues, curve);

                    for (auto& function : functions)
                    {
                        const auto index = &function - functions.data();
                        const auto error = fit_error(references[index], curveValues);
                        if (error < function.error)
                        {
                            function.curve = curve;
                            function.error = error;
                        }
                    }
                }
            }
        }
    }
}

template <size_t N = 100, double p1x = 0.0, double p1y = 0.0, double p2x = 1.0, double p2y = 1.0>
auto brute_force_refine_curve(sampled_curve_t<N> reference, int NX, int NY, cubic_bezier_t base, double rangeX, double rangeY, int iterations_left) -> cubic_bezier_t
{
    if (iterations_left <= 0)
        return base;

    const auto x_units = NX + (NX % 2 ? 0 : 1); // always use odd number of units to ensure original point is included
    const auto y_units = NY + (NX % 2 ? 0 : 1); // always use odd number of units to ensure original point is included

    const auto x1_min   = std::max(0.0, base.CP1.x - rangeX * 0.5);
    const auto x1_max   = std::min(1.0, base.CP1.x + rangeX * 0.5);
    const auto x1_step  = (x1_max - x1_min) / x_units;

    const auto x2_min   = std::max(0.0, base.CP2.x - rangeX * 0.5);
    const auto x2_max   = std::min(1.0, base.CP2.x + rangeX * 0.5);
    const auto x2_step  = (x1_max - x1_min) / x_units;

    const auto y1_min   = std::max(-1.0, base.CP1.y - rangeY * 0.5);
    const auto y1_max   = std::min( 2.0, base.CP1.y + rangeY * 0.5);
    const auto y1_step  = (y1_max - y1_min) / y_units;

    const auto y2_min   = std::max(-1.0, base.CP2.y - rangeY * 0.5);
    const auto y2_max   = std::min( 2.0, base.CP2.y + rangeY * 0.5);
    const auto y2_step  = (y1_max - y1_min) / y_units;

    cubic_bezier_t best_fit_curve = base;
    double best_error = fit_error<N, p1x, p1y, p2x, p2y>(reference, base);
    for (int x1 = 0; x1 < x_units; ++x1)
    {
        for (int y1 = 0; y1 < y_units; ++y1)
        {
            for (int x2 = 0; x2 < x_units; ++x2)
            {
                for (int y2 = 0; y2 < y_units; ++y2)
                {
                    cubic_bezier_t curve;
                    curve.CP1.x = x1_min + x1 * x1_step;
                    curve.CP1.y = y1_min + y1 * y1_step;
                    curve.CP2.x = x2_min + x2 * x2_step;
                    curve.CP2.y = y2_min + y2 * y2_step;

                    const auto error = fit_error<N, p1x, p1y, p2x, p2y>(reference, curve);
                    if (error < best_error)
                    {
                        best_fit_curve = curve;
                        best_error = error;
                    }
                }
            }
        }
    }

    return brute_force_refine_curve<N, p1x, p1y, p2x, p2y>(reference, NX, NY, best_fit_curve, rangeX * 0.25, rangeY * 0.25, iterations_left - 1);
}

template <size_t N = 100>
void brute_force_refine_curve_many(std::span<function_t> functions, int NX, int NY, double x_step, double y_step, int recurstions_left)
{
#pragma omp parallel for
    for (int i = 0; i < std::ssize(functions); ++i)
    {
        auto& function = functions[i];

        sampled_curve_t<N> reference;
        sample_curve_into(reference, function.func);

        function.curve = brute_force_refine_curve<N>(reference, NX, NY, function.curve, x_step, y_step, recurstions_left);
        function.error = fit_error<N>(reference, function.curve);
    }
}

template <size_t N = 100>
void brute_force_refine_curve_many_until(std::span<function_t> functions, int NX, int NY, double x_step, double y_step, int recurstions_left, double max_error, int iteration_limit)
{
    for (int iteration = 0; iteration < iteration_limit; ++iteration)
    {
        brute_force_refine_curve_many<N>(functions, NX, NY, x_step, y_step, recurstions_left);
        bool all_below_max_error = true;
        for (const auto& function : functions)
        {
            if (function.error > max_error)
            {
                all_below_max_error = false;
                break;
            }
        }
        if (all_below_max_error)
            break;
    }
}


// initial values from https://gist.github.com/zz85/2a0e4a0b944ec89aa5eb
//
// QuadIn:     { { 0.26,  0.00 }, { 0.60,  0.20 } } }, // 0.0000014195846674133613
// QuadOut:    { { 0.40,  0.80 }, { 0.74,  1.00 } } }, // 0.0000014195846674136283
// QuadInOut:  { { 0.48,  0.04 }, { 0.52,  0.96 } } }, // 0.00021443512293854952
// SineIn:     { { 0.32,  0.00 }, { 0.60,  0.36 } } }, // 0.00003239261842406147
// SineOut:    { { 0.40,  0.64 }, { 0.68,  1.00 } } }, // 0.00003239261842406522
// SineInOut:  { { 0.36,  0.00 }, { 0.64,  1.00 } } }, // 0.000042771492108870344
// ExpoIn:     { { 0.62,  0.02 }, { 0.84, -0.08 } } }, // 0.00019519295753793093
// ExpoOut:    { { 0.16,  1.08 }, { 0.38,  0.98 } } }, // 0.00019519295753793044
// ExpoInOut:  { { 0.84, -0.12 }, { 0.16,  1.12 } } }, // 0.002911151068100008
// CircIn:     { { 0.54,  0.00 }, { 1.00,  0.44 } } }, // 0.000051489431433657843
// CircOut:    { { 0.00,  0.56 }, { 0.46,  1.00 } } }, // 0.00005148943143365869
// CircInOut:  { { 0.88,  0.14 }, { 0.12,  0.86 } } }, // 0.003209319580927181
// CubicIn:    { { 0.32,  0.00 }, { 0.66, -0.02 } } }, // 0.0000837264957522941
// CubicOut:   { { 0.34,  1.02 }, { 0.68,  1.00 } } }, // 0.00008372649575229691
// CubicInOut: { { 0.62, -0.04 }, { 0.38,  1.04 } } }, // 0.00020290989758759337
// QuartIn:    { { 0.46,  0.00 }, { 0.74, -0.04 } } }, // 0.00004920119056999463
// QuartOut:   { { 0.26,  1.04 }, { 0.54,  1.00 } } }, // 0.00004920119056999559
// QuartInOut: { { 0.70, -0.10 }, { 0.30,  1.10 } } }, // 0.0007318300363503209
// QuintIn:    { { 0.52,  0.00 }, { 0.78, -0.10 } } }, // 0.00015727214523005402
// QuintOut:   { { 0.22,  1.10 }, { 0.48,  1.00 } } }, // 0.0001572721452300539
// QuintInOut: { { 0.76, -0.14 }, { 0.24,  1.14 } } }, // 0.002051715614688728


// ///////////////////////////////////////////
// // Animation Cubic Bezier easings        //
// // see - http://matthewlein.com/ceaser   //
// ///////////////////////////////////////////
//
// https://gist.github.com/ugonnanwosu/5176c8892f32c55dd5ac
//
// defaults
// $ease-linear:          "cubic-bezier(0.250, 0.250, 0.750, 0.750)";
// $ease-default:         "cubic-bezier(0.250, 0.100, 0.250, 1.000)";
// $ease-in:              "cubic-bezier(0.420, 0.000, 1.000, 1.000)";
// $ease-out:             "cubic-bezier(0.000, 0.000, 0.580, 1.000)";
// $ease-in-out:          "cubic-bezier(0.420, 0.000, 0.580, 1.000)";
//
// // penner equations (approximated)
// $ease-in-quad:         "cubic-bezier(0.550, 0.085, 0.680, 0.530)";
// $ease-in-cubic:        "cubic-bezier(0.550, 0.055, 0.675, 0.190)";
// $ease-in-quart:        "cubic-bezier(0.895, 0.030, 0.685, 0.220)";
// $ease-in-quint:        "cubic-bezier(0.755, 0.050, 0.855, 0.060)";
// $ease-in-sine:         "cubic-bezier(0.470, 0.000, 0.745, 0.715)";
// $ease-in-expo:         "cubic-bezier(0.950, 0.050, 0.795, 0.035)";
// $ease-in-circ:         "cubic-bezier(0.600, 0.040, 0.980, 0.335)";
// $ease-in-back:         "cubic-bezier(0.600, -0.280, 0.735, 0.045)";
//
// $ease-out-quad:        "cubic-bezier(0.250, 0.460, 0.450, 0.940)";
// $ease-out-cubic:       "cubic-bezier(0.215, 0.610, 0.355, 1.000)";
// $ease-out-quart:       "cubic-bezier(0.165, 0.840, 0.440, 1.000)";
// $ease-out-quint:       "cubic-bezier(0.230, 1.000, 0.320, 1.000)";
// $ease-out-sine:        "cubic-bezier(0.390, 0.575, 0.565, 1.000)";
// $ease-out-expo:        "cubic-bezier(0.190, 1.000, 0.220, 1.000)";
// $ease-out-circ:        "cubic-bezier(0.075, 0.820, 0.165, 1.000)";
// $ease-out-back:        "cubic-bezier(0.175, 0.885, 0.320, 1.275)";
//
// $ease-in-out-quad:     "cubic-bezier(0.455, 0.030, 0.515, 0.955)";
// $ease-in-out-cubic:    "cubic-bezier(0.645, 0.045, 0.355, 1.000)";
// $ease-in-out-quart:    "cubic-bezier(0.770, 0.000, 0.175, 1.000)";
// $ease-in-out-quint:    "cubic-bezier(0.860, 0.000, 0.070, 1.000)";
// $ease-in-out-sine:     "cubic-bezier(0.445, 0.050, 0.550, 0.950)";
// $ease-in-out-expo:     "cubic-bezier(1.000, 0.000, 0.000, 1.000)";
// $ease-in-out-circ:     "cubic-bezier(0.785, 0.135, 0.150, 0.860)";
// $ease-in-out-back:     "cubic-bezier(0.680, -0.550, 0.265, 1.550)";

function_t functions[] =
{
    { "Linear",           "linear",       [](auto t) { return t; }    , { { 1.0/3.0,   1.0/3.0  }, { 2.0/3.0,  2.0/3.0  } } }, // 0.0000000000000000
    { "EaseInQuad",       "quadIn",       &math::simple_quad_in       , { { 0.26,      0.00     }, { 0.60,     0.20     } } }, // 0.0000014195846674133613
    { "EaseOutQuad",      "quadOut",      &math::simple_quad_out      , { { 0.40,      0.80     }, { 0.74,     1.00     } } }, // 0.0000014195846674136283
    { "EaseInOutQuad",    "quadInOut",    &math::simple_quad_inout    , { { 0.48,      0.04     }, { 0.52,     0.96     } } }, // 0.00021443512293854952
    { "EaseOutInQuad",    "quadOutIn",    &math::simple_quad_outin    , { { 0.829607,  0.677567 }, { 0.922548, 0.845253 } } },
    { "EaseInSine",       "sinIn",        &math::simple_sin_in        , { { 0.32,      0.00     }, { 0.60,     0.36     } } }, // 0.00003239261842406147
    { "EaseOutSine",      "sinOut",       &math::simple_sin_out       , { { 0.40,      0.64     }, { 0.68,     1.00     } } }, // 0.00003239261842406522
    { "EaseInOutSine",    "sinInOut",     &math::simple_sin_inout     , { { 0.36,      0.00     }, { 0.64,     1.00     } } }, // 0.000042771492108870344
    { "EaseOutInSine",    "sinOutIn",     &math::simple_sin_outin     , { { 0.794529,  0.644224 }, { 0.899130, 0.841548 } } },
    { "EaseInExpo",       "expIn",        &math::simple_exp_in        , { { 0.62,      0.02     }, { 0.84,    -0.08     } } }, // 0.00019519295753793093
    { "EaseOutExpo",      "expOut",       &math::simple_exp_out       , { { 0.16,      1.08     }, { 0.38,     0.98     } } }, // 0.00019519295753793044
    { "EaseInOutExpo",    "expInOut",     &math::simple_exp_inout     , { { 0.84,     -0.12     }, { 0.16,     1.12     } } }, // 0.002911151068100008
    { "EaseOutInExpo",    "expOutIn",     &math::simple_exp_outin     , { { 0.333333,  0.500000 }, { 0.666667, 0.500000 } } },
    { "EaseInCirc",       "circIn",       &math::simple_circ_in       , { { 0.54,      0.00     }, { 1.00,     0.44     } } }, // 0.000051489431433657843
    { "EaseOutCirc",      "circOut",      &math::simple_circ_out      , { { 0.00,      0.56     }, { 0.46,     1.00     } } }, // 0.00005148943143365869
    { "EaseInOutCirc",    "circInOut",    &math::simple_circ_inout    , { { 0.88,      0.14     }, { 0.12,     0.86     } } }, // 0.003209319580927181
    { "EaseOutInCirc",    "circOutIn",    &math::simple_circ_outin    , { { 0.988885,  0.879079 }, { 0.998044, 0.938129 } } },
    { "EaseInCubic",      "cubeIn",       &math::simple_cube_in       , { { 0.32,      0.00     }, { 0.66,    -0.02     } } }, // 0.0000837264957522941
    { "EaseOutCubic",     "cubeOut",      &math::simple_cube_out      , { { 0.34,      1.02     }, { 0.68,     1.00     } } }, // 0.00008372649575229691
    { "EaseInOutCubic",   "cubeInOut",    &math::simple_cube_inout    , { { 0.62,     -0.04     }, { 0.38,     1.04     } } }, // 0.00020290989758759337
    { "EaseOutInCubic",   "cubeOutIn",    &math::simple_cube_outin    , { { 0.856523,  0.576503 }, { 0.942109, 0.826672 } } },
    { "EaseInQuart",      "quarticIn",    &math::simple_quartic_in    , { { 0.798788,  0.184878 }, { 0.923863, 0.695910 } } },
    { "EaseOutQuart",     "quarticOut",   &math::simple_quartic_out   , { { 0.46,      0.00     }, { 0.74,    -0.04     } } }, // 0.00004920119056999463
    { "EaseInOutQuart",   "quarticInOut", &math::simple_quartic_inout , { { 0.26,      1.04     }, { 0.54,     1.00     } } }, // 0.00004920119056999559
    { "EaseOutInQuart",   "quarticOutIn", &math::simple_quartic_outin , { { 0.70,     -0.10     }, { 0.30,     1.10     } } }, // 0.0007318300363503209
    { "EaseInQuint",      "quinticIn",    &math::simple_quintic_in    , { { 0.52,      0.00     }, { 0.78,    -0.10     } } }, // 0.00015727214523005402
    { "EaseOutQuint",     "quinticOut",   &math::simple_quintic_out   , { { 0.22,      1.10     }, { 0.48,     1.00     } } }, // 0.0001572721452300539
    { "EaseInOutQuint",   "quinticInOut", &math::simple_quintic_inout , { { 0.76,     -0.14     }, { 0.24,     1.14     } } }, // 0.002051715614688728
    { "EaseOutInQuint",   "quinticOutIn", &math::simple_quintic_outin , { { 0.939403,  0.644893 }, { 0.971176, 0.856456 } } },
    { "EaseInBack",       "backIn",       &math::simple_back_in       , { { 0.600,    -0.280    }, { 0.735,    0.045    } } }, // 0.034971136882655643
    { "EaseOutBack",      "backOut",      &math::simple_back_out      , { { 0.175,     0.885    }, { 0.320,    1.275    } } },
    { "EaseInOutBack",    "backInOut",    &math::simple_back_inout    , { { 0.680,    -0.550    }, { 0.265,    1.550    } } }, // 0.48314638943436461
    { "EaseOutInBack",    "backOutIn",    &math::simple_back_outin    , { { 0.916716,  0.557008 }, { 0.969700, 0.857933 } } },

};

double parabola(double t)
{
    return (1.0 - (t - 0.5) * (t - 0.5) * 4.0);
}

cubic_bezier_t fit_parabola_curve()
{
    cubic_bezier_t curve;
    curve.CP1.x = 1.0 / 3.0;
    curve.CP1.y = 1.0 / 3.0 + 1.0;
    curve.CP2.x = 2.0 / 3.0;
    curve.CP2.y = 1.0 / 3.0 + 1.0;

    return curve;
}

template <size_t N = 100>
void fit_elastic_in()
{
    // math::simple_elastic_in

    constexpr const double zero_points[] = { 0.0, 0.025, 0.175, 0.325, 0.475, 0.625, 0.775, 0.925, 1.0 };

    auto fit_segment = []<double p2y>(double x1, double x2) -> std::pair<cubic_bezier_t, double>
    {
        const auto offset = x1;
        const auto scale  = (x2 - offset);

        auto function = [&](double t) -> double
        {
            return math::simple_elastic_in(t * scale + offset);
        };

        auto initial_guess = graphics_gems_fit_curve<N>(function, 0.0000001);

        sampled_curve_t<N> reference;
        sample_curve_into(reference, function);

        auto curve = initial_guess;
        int iteration_limit = 10;
        while (fit_error<N, 0, 0, 1, p2y>(reference, curve) > 0.0000001)
        {
            if (iteration_limit-- <= 0)
                break;

            curve = brute_force_refine_curve<N, 0, 0, 1, p2y>(reference, 5, 10, curve, 0.5, 0.5, 10);
        }

        auto error = fit_error<N, 0, 0, 1, p2y>(reference, curve);

        curve.CP1.x = curve.CP1.x * scale + offset;
        curve.CP2.x = curve.CP2.x * scale + offset;

        return { curve, error };
    };

    printf("// Fitting elastic-in curve segments:\n");
    printf("{\n");
    printf("    { %9.6f, %9.6f },\n", 0.0, 0.0);

    constexpr auto segment_count = std::size(zero_points) - 1;
    constexpr auto last_segment_index = segment_count - 1;
    for (size_t i = 0; i < segment_count; ++i)
    {
        const auto offset = zero_points[i];
        const auto scale  = (zero_points[i + 1] - offset);

        const auto p2y = i == last_segment_index ? 1.0 : 0.0;

        cubic_bezier_t curve;
        float error = 0.0;
        if (p2y < 1)
        {
            auto result = fit_segment.operator()<0>(offset, zero_points[i + 1]);
            curve = result.first;
            error = result.second;
        }
        else
        {
            auto result = fit_segment.operator()<1>(offset, zero_points[i + 1]);
            curve = result.first;
            error = result.second;
        }

        printf("    { %9.6f, %9.6f }, // error: %.20lf\n", curve.CP1.x, curve.CP1.y, error);
        printf("    { %9.6f, %9.6f },\n", curve.CP2.x, curve.CP2.y);
        printf("    { %9.6f, %9.6f },\n", zero_points[i + 1], p2y);
    }
    printf("}\n");
}

template <size_t N = 100>
void fit_elastic_in_out_first_half()
{
    // first half of math::simple_elastic_in

    constexpr const double roots[] = { 0.0, 0.10625 * 2.0, 0.21875 * 2.0, 0.33125 * 2.0, 0.44375 * 2.0, 1.0 };

    auto fit_segment = []<double p1y, double p2y>(double x1, double x2) -> std::pair<cubic_bezier_t, double>
    {
        const auto offset = x1;
        const auto scale  = (x2 - offset);

        auto function = [&](double t) -> double
        {
            return math::simple_elastic_inout((t * scale + offset) * 0.5) * 2.0;
        };

        auto initial_guess = graphics_gems_fit_curve<N>(function, 0.0000001);

        sampled_curve_t<N> reference;
        sample_curve_into(reference, function);

        auto curve = initial_guess;
        int iteration_limit = 10;
        while (fit_error<N, 0, p1y, 1, p2y>(reference, curve) > 0.0000001)
        {
            if (iteration_limit-- <= 0)
                break;

            curve = brute_force_refine_curve<N, 0, p1y, 1, p2y>(reference, 5, 10, curve, 0.5, 0.5, 10);
        }

        auto error = fit_error<N, 0, p1y, 1, p2y>(reference, curve);

        curve.CP1.x = curve.CP1.x * scale + offset;
        curve.CP2.x = curve.CP2.x * scale + offset;

        return { curve, error };
    };

    printf("// Fitting elastic-in-out (first half) curve segments:\n");
    printf("{\n");
    printf("    { %9.6f, %9.6f },\n", 0.0, 0.0);

    constexpr auto segment_count = std::size(roots) - 1;
    constexpr auto last_segment_index = segment_count - 1;
    for (size_t i = 0; i < segment_count; ++i)
    {
        const auto offset = roots[i];
        const auto scale  = (roots[i + 1] - offset);

        const auto p1y = 0.0;
        const auto p2y = i == last_segment_index ? 1.0 : 0.0;

        std::pair<cubic_bezier_t, double> result;

        auto& [curve, error] = result;

        if (p1y < 1 && p2y < 1)
            result = fit_segment.operator()<0, 0>(offset, roots[i + 1]);
        else if (p1y >= 1 && p2y < 1)
            result = fit_segment.operator()<1, 0>(offset, roots[i + 1]);
        else if (p1y < 1 && p2y >= 1)
            result = fit_segment.operator()<0, 1>(offset, roots[i + 1]);
        else
            result = fit_segment.operator()<1, 1>(offset, roots[i + 1]);

        printf("    { %9.6f, %9.6f }, // error: %.20lf\n", curve.CP1.x, curve.CP1.y, error);
        printf("    { %9.6f, %9.6f },\n", curve.CP2.x, curve.CP2.y);
        printf("    { %9.6f, %9.6f },\n", roots[i + 1], p2y);
    }
    printf("}\n");
}

int main()
{
    constexpr size_t SAMPLES = 100;

    // Convert elastic-in to set of cubic bezier curves
    fit_elastic_in();

    // Convert elastic-in-out (first half) to set of cubic bezier curves
    fit_elastic_in_out_first_half();

    // Fit the curves using Philip J. Schneider piecewise cubic fitting from "Graphics Gems"
    graphics_gems_fit_curve_many<SAMPLES>(functions, 0.0000001);

    // Fit the curves using brute force method
    // brute_force_fit_curve_many(functions, 50, 100);

    // Refine the curves stepping down narrowing fit window until the error is below a certain threshold
    brute_force_refine_curve_many_until<SAMPLES>(functions, 5, 10, 0.5f, 0.5f, 10, 0.0000001, 100);

    // Print the results
    printf("// Penner's equations approximations:\n");
    for (const auto& function : functions)
    {
        printf("%-20s (%9.6f, %9.6f, %9.6f, %9.6f) // max error: %.20lf\n",
            function.name, function.curve.CP1.x, function.curve.CP1.y, function.curve.CP2.x, function.curve.CP2.y, function.error);
    }

    printf("\n\n// Table:\n");
    for (const auto& function : functions)
    {
        char enumName[64];
        snprintf(enumName, sizeof(enumName), "%s:", function.enumName);

        printf("    case %-17s return make_pair(vec2(%9.6ff, %9.6ff), vec2(%9.6ff, %9.6ff)); // max error: %.20lf\n",
            enumName, function.curve.CP1.x, function.curve.CP1.y, function.curve.CP2.x, function.curve.CP2.y, function.error);
    }

    system("pause");

    return 0;
}