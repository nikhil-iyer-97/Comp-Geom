#include "fundamentals.h"

constexpr long double pi()
{
    return atan((long double)(1.0)) * (long double)(4.0);
}

/*************** Class: angle ***************/

angle::angle(long double rad = 0.0)
{
    set_angle(rad);
}

angle::angle(long double y, long double x)
{
    set_angle(y, x);
}

void angle::set_angle(long double rad)
{
    if(rad >= 2*pi() or rad < 0.0)
    {
        rad = fmod((180.0*(rad/pi())), 360);
        if(rad < 0.0)
            rad += 360;
        rad = (pi()*rad)/180.0;
    }
    value = rad;
    if(((value*2.0)/pi())==1.0 or ((value*2.0)/pi())==3.0)
    {
        y = -1.0*(((value*2.0)/pi())-2.0);
        x = 0;
    }
    else{
        y = tan(value);
        x = 1.0;
    }
    quad = (int)fmod(((2.0*rad)/pi()), 4)+1;
}

void angle::set_angle(long double opp, long double adj)
{
    y = opp;
    x = adj;
    value = atan2(y, x);
    if(value < 0.0) value += 2.0*pi();
    if(x>0 and y>0) quad = 1;
    if(x<=0 and y>0) quad = 2;
    if(x<=0 and y<=0) quad = 3;
    if(x>0 and y<=0) quad = 4;
}

long double angle::get_value(bool type = true) const
{
    if(type)
        return value;
    else
        return ((value*180.0)/pi());
}

pair< long double, long double > angle::get_ratio() const
{
    return make_pair(y, x);
}

int angle::get_quadrant() const
{
    return quad;
}

/*********************************************/

/************ Operators : angle **************/

bool operator < (const angle& a, const angle& b)
{
    if(a.get_quadrant() != b.get_quadrant())
        return a.get_quadrant() < b.get_quadrant();
    return a.get_value() < b.get_value();
}

bool operator > (const angle& a, const angle& b)
{
    return b < a;
}

bool operator == (const angle& a, const angle& b)
{
    return ((!(b < a)) and (!(b > a)));
}

/*********************************************/

/*************** Class: point ***************/

point::point(long double x=0.0, long double y=0.0, long double ox = 0.0, long double oy = 0.0)
{
    set_point(x, y, ox, oy);
}

point::point(long double r, angle theta)
{
    set_point(r, theta);
}

void point::set_point(long double x, long double y, long double ox = 0.0, long double oy = 0.0)
{
    this->x = x - ox;
    this->y = y - oy;
    this->ox = ox;
    this->oy = oy;
    r_sq = (this->x)*(this->x) + (this->y)*(this->y);
    theta = angle((this->y), (this->x));
}

void point::set_point(long double r, angle theta)
{
    this->r_sq = r*r;
    this->theta = angle(theta.get_value());
    this->ox = 0.0;
    this->oy = 0.0;
    x = r*cos(theta.get_value());
    y = r*sin(theta.get_value());
}

pair< long double, long double>  point::get_point_cartesian() const
{
    return make_pair(x, y);
}

pair< long double, angle>  point::get_point_polar() const
{
    return make_pair(sqrt(r_sq), theta);
}

pair< long double, long double>  point::get_point_unnormalized() const
{
    return make_pair(x+ox, y+oy);
}

long double point::get_radius_squared() const
{
    return r_sq;
}

/*********************************************/

/**** Sort Function and Related : point ******/

bool less_than_cartesian(const point& a, const point& b)
{
    return a.get_point_cartesian() < b.get_point_cartesian();
}

bool less_than_polar(const point& a, const point& b)
{
    if(!(a.get_point_polar().second == b.get_point_polar().second))
        return a.get_point_polar().second < b.get_point_polar().second;
    return a.get_radius_squared() < b.get_radius_squared();
}

template <class RandomAccessIterator>
void sort_points(RandomAccessIterator first, RandomAccessIterator last, bool polar)
{
    if(!is_same< typename iterator_traits<RandomAccessIterator>::value_type, point >::value)
        return;
    if(polar)
        sort(first, last, less_than_polar);
    else
        sort(first, last, less_than_cartesian);
}

/*********************************************/

/*************** Miscellaneous ***************/

long double twice_triangle_signed_area(point a, point b, point c)
{
    long double x0 = a.get_point_cartesian().first;
    long double y0 = a.get_point_cartesian().second;
    long double x1 = b.get_point_cartesian().first;
    long double y1 = b.get_point_cartesian().second;
    long double x2 = c.get_point_cartesian().first;
    long double y2 = c.get_point_cartesian().second;
    return (x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1));
}

/*********************************************/
