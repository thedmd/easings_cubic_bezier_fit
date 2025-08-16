#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Point2Struct {	/* 2d point */
    double x, y;
} Point2;
typedef Point2 Vector2;

#define MAXPOINTS	1000		/* The most points you can have */

typedef Point2 *BezierCurve;

/* Forward declarations */
void FitCurve(const Point2 *d, int nPts, double error, void(*result)(void* ctx, int n, const BezierCurve curve), void* ctx);

#ifdef __cplusplus
} // extern "C"
#endif