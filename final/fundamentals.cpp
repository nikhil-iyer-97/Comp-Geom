#include "fundamentals.h"

constexpr long double pi()
{
    /*
    This function evaluates the value of pi at compile time.

    Return Value:
    The (approximated) decimal value of pi.*/
    return atan((long double)(1.0)) * (long double)(4.0);
}

/*************** Class: angle ***************/

angle::angle(long double rad = 0.0)
{
    /*
    Constructor to initialize an angle instance.

    Parameters:
    rad - the value of the angle in radians. */
    set_angle(rad);
}

angle::angle(long double y, long double x)
{
    /*
    Constructor to initialize an angle instance,
    where the angle is represented by an opposite and adjecent ratio.

    Parameters:
    y - length of opposite edge of the angle.
    x - length of adjecent edge of the angle. */

    set_angle(y, x);
}

void angle::set_angle(long double rad)
{
    /*
    Sets the value of an angle instance.

    Parameters:
    rad - the value of the angle in radians. */

    if(rad > 2*pi() or rad < 0.0 or fabs(rad-2*pi()) < eps)
    {
        rad = fmod((180.0*(rad/pi())), 360);
        if(rad < 0.0)
            rad += 360;
        rad = (pi()*rad)/180.0;
    }
    value = rad;
    if((fabs(((value*2.0)/pi())-1.0) < eps) or
      ((fabs(((value*2.0)/pi())-3.0) < eps)))
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
    /*
    Sets the value of an angle instance,
    where the angle is represented by an opposite and adjecent ratio.

    Parameters:
    y - length of opposite edge of the angle.
    x - length of adjecent edge of the angle. */

    y = opp;
    x = adj;
    value = atan2(y, x);
    if(value < 0.0) value += 2.0*pi();
    quad = (int)fmod(((2.0*value)/pi()), 4)+1;
}

long double angle::get_value(bool type = true) const
{
    /*
    Returns the value of the angle in radians or degrees.

    Parameters:
    type - the angle is in degrees if this parameter is false, else in radian.

    Return Value:
    The value in the desired units.
    */

    if(type)
        return value;
    else
        return ((value*180.0)/pi());
}

pair< long double, long double > angle::get_ratio() const
{
    /*
    Returns the value of the angle as an opposite and ajecent ratio.

    Return Value:
    A pair representing the opposite and adjecent lengths respectively.
    */
    return make_pair(y, x);
}

int angle::get_quadrant() const
{
    /*
    Returns the quadrant of the angle.

    Return Value:
    1, 2, 3 or 4, based on the quadrant.
    */

    return quad;
}

long double angle::get_tan() const
{
    /*
    Returns the tangent of the angle.

    Return Value:
    The tangent of the angle.
    */

    return y/x;
}

/*********************************************/

/************ Operators : angle **************/

bool operator < (const angle& a, const angle& b)
{
    /*
    Less than comparison of angles.

    Parameters:
    a - 1st angle.
    b - second angle.

    Return Value:
    boolean true is a < b, else boolean false.
    */
    if(a.get_quadrant() != b.get_quadrant())
        return a.get_quadrant() < b.get_quadrant();
    return a.get_value() < b.get_value();
}

bool operator > (const angle& a, const angle& b)
{
    /*
    Greater than comparison of angles.

    Parameters:
    a - 1st angle.
    b - second angle.

    Return Value:
    boolean true is a > b, else boolean false.
    */
    return b < a;
}

bool operator == (const angle& a, const angle& b)
{
    /*
    Equal to comparison of angles.

    Parameters:
    a - 1st angle.
    b - second angle.

    Return Value:
    boolean true is a == b, else boolean false.
    */
    return ((!(b < a)) and (!(b > a)));
}

bool operator != (const angle& a, const angle& b)
{
    /*
    Not Equal to comparison of angles.

    Parameters:
    a - 1st angle.
    b - second angle.

    Return Value:
    boolean true is a != b, else boolean false.
    */
    return !((!(b < a)) and (!(b > a)));
}

/*********************************************/

/*************** Class: point ***************/

point::point(long double x = 0.0, long double y = 0.0,
             long double ox = 0.0, long double oy = 0.0)
{
    /*
    Constructor of an instance of point.

    Parameters:
    x - x coordinate of point
    y - y coordinate of point
    ox - x coordinate of origin
    oy - x coordinate of origin
    */
    set_point(x, y, ox, oy);
}

point::point(long double r, angle theta)
{
    /*
    Constructor of an instance of point in polar form.

    Parameters:
    r - radius of point
    theta - angle of point
    */
    set_point(r, theta);
}

void point::set_point(long double xx, long double yy,
                      long double oox = 0.0, long double ooy = 0.0)
{
    /*
    Sets the value of a point.

    Parameters:
    x - x coordinate of point
    y - y coordinate of point
    ox - x coordinate of origin
    oy - x coordinate of origin
    */
    this->x = xx - oox;
    this->y = yy - ooy;
    r_sq = (this->x)*(this->x) + (this->y)*(this->y);
    theta.set_angle((this->y), (this->x));
}

void point::set_point(long double r, angle theta)
{
    /*
    Sets the value of a point using polar coordinates.

    Parameters:
    r - radius of point
    theta - angle of point
    */
    this->r_sq = r*r;
    this->theta.set_angle(theta.get_value());
    x = r*cos(theta.get_value());
    y = r*sin(theta.get_value());
}

pair< long double, long double>  point::get_point_cartesian() const
{
    /*
    Returns the Cartesian coordinates of the point.

    Return Value:
    (x coordinate, y coordinate)
    */
    return make_pair(x, y);
}

pair< long double, angle>  point::get_point_polar() const
{
    /*
    Returns the polar coordinates of the point.

    Return Value:
    (radius of point, angle of point)
    */
    return make_pair(sqrt(r_sq), theta);
}

pair< long double, long double>  point::get_point_unnormalized(point p) const
{
    /*
    Returns the original cartesian coordinates of the point,
    whose origin is shifted.

    Parameters:
    p - the origin

    Return Value:
    (x coordinate, y coordinate)
    */
    auto xx = p.get_point_cartesian();
    long double a = x + xx.first, b = y + xx.second;
    return make_pair(a, b);
}

long double point::get_radius_squared() const
{
    /*
    Returns squared radius of the point.

    Return Value:
    Squared radius of the point.
    */
    return r_sq;
}

/*********************************************/

/**** Sort Function and Related : point ******/

bool less_than_cartesian(const point& a, const point& b)
{
    /*
    less than comparison of points based on cartesian coordinates.

    Parameters:
    a - 1st point.
    b - second point.

    Return Value:
    boolean true is a < b, else boolean false.
    */
    return a.get_point_cartesian() < b.get_point_cartesian();
}

bool less_than_polar(const point& a, const point& b)
{
    /*
    less than comparison of points based on polar coordinates.

    Parameters:
    a - 1st point.
    b - second point.

    Return Value:
    boolean true is a < b, else boolean false.
    */
    if(a.get_point_polar().second != b.get_point_polar().second)
        return a.get_point_polar().second < b.get_point_polar().second;
    return a.get_radius_squared() < b.get_radius_squared();
}

bool operator==(const point& p, const point& q)
{
    /*
    Equal to comparison of points.

    Parameters:
    a - 1st point.
    b - second point.

    Return Value:
    boolean true is a == b, else boolean false.
    */
    return ((!less_than_cartesian(p, q)) && (!less_than_cartesian(q, p)));
}

bool operator!=(const point& p, const point& q)
{
    /*
    Not Equal to comparison of points.

    Parameters:
    a - 1st point.
    b - second point.

    Return Value:
    boolean true is a != b, else boolean false.
    */
    return !(p == q);
}

void swap(point& p, point& q)
{
    /*
    Swap two points.

    Parameters:
    a - 1st point.
    b - second point.
    */
    point t(p.get_point_cartesian().first, p.get_point_cartesian().second);
    p.set_point(q.get_point_cartesian().first, q.get_point_cartesian().second);
    q.set_point(t.get_point_cartesian().first, t.get_point_cartesian().second);
}

template <class RandomAccessIterator>
void sort_points(RandomAccessIterator f, RandomAccessIterator l, bool pol)
{
    /*
    Sorts a container of points.

    Parameters:
    f - iterator to 1st point.
    l - iterator to last point.
    pol - performs polar sort if true else cartesian sort..
    */
    if(!is_same< typename iterator_traits<RandomAccessIterator>::value_type,
        point >::value)
        return;
    if(pol)
        sort(f, l, less_than_polar);
    else
        sort(f, l, less_than_cartesian);
}

/*********************************************/

/*************** Miscellaneous ***************/

long double twice_triangle_signed_area(point a, point b, point c)
{
    /*
    Returns 2*(signed area of triangle).

    Parameters:
    a - 1st vertex
    b - 2nd vertex
    c - 3rd vertex

    Return Value:
    2*(signed area of triangle)
    */
    long double x0 = a.get_point_cartesian().first;
    long double y0 = a.get_point_cartesian().second;
    long double x1 = b.get_point_cartesian().first;
    long double y1 = b.get_point_cartesian().second;
    long double x2 = c.get_point_cartesian().first;
    long double y2 = c.get_point_cartesian().second;
    return (x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1));
}

/*********************************************/
