#ifndef FUNDAMENTALS_H
#define FUNDAMENTALS_H

#include <bits/stdc++.h>

using namespace std;

constexpr long double pi();

class angle
{
    private:
    long double value, x, y;
    int quad;
    public:
    angle(long double);
    angle(long double, long double);
    void set_angle(long double);
    void set_angle(long double, long double);
    long double get_value(bool) const;
    pair< long double, long double > get_ratio() const;
    int get_quadrant() const;
};

bool operator<(const angle&, const angle&);
bool operator>(const angle&, const angle&);
bool operator==(const angle&, const angle&);

class point
{
    private:
    long double x, y, r_sq, ox, oy;
    angle theta;
    public:
    point(long double, long double, long double, long double);
    point(long double, angle);
    void set_point(long double, long double, long double, long double);
    void set_point(long double, angle);
    pair< long double, long double>  get_point_cartesian() const;
    pair< long double, angle>  get_point_polar() const;
    pair< long double, long double>  get_point_unnormalized() const;
    long double get_radius_squared() const;
};

bool less_than_cartesian(const point&, const point&);
bool less_than_polar(const point&, const point&);
template <class RandomAccessIterator>
void sort_points(RandomAccessIterator, RandomAccessIterator, bool polar = false);

long double twice_triangle_signed_area(point a, point b, point c);
#endif
