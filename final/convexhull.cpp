#include "convexhull.h"

vector<point> graham_scan(const vector<point>& pts)
{
    /*
    This function calculates the convex hull of a given set of points,
    based on Graham Scan Algorithm.

    Parameters:
    points - a vector of points whose convex hull is to be found.

    Return Value:
    A vector of points representing the convex hull in anti-clockwise order.*/

	if(pts.size() <= 2)
        return pts;
	vector<point> answer;

	/* Finds leftmost, bottommost point. */

    point lefbtm(pts[0].get_point_cartesian().first,
                 pts[0].get_point_cartesian().second);
    for (int i = 0; i < pts.size(); ++i)
    {
    	if(lefbtm.get_point_cartesian().second >
           pts[i].get_point_cartesian().second)
    	      lefbtm.set_point(pts[i].get_point_cartesian().first,
                               pts[i].get_point_cartesian().second);
        else if(fabs(lefbtm.get_point_cartesian().second -
                     pts[i].get_point_cartesian().second) < eps and
                     lefbtm.get_point_cartesian().first >
                     pts[i].get_point_cartesian().first)
            lefbtm.set_point(pts[i].get_point_cartesian().first,
                             pts[i].get_point_cartesian().second);
    }
    vector<point> shifted_pts;
    point tmp(lefbtm.get_point_cartesian().first,
              lefbtm.get_point_cartesian().second,
              lefbtm.get_point_cartesian().first,
              lefbtm.get_point_cartesian().second);

    answer.push_back(tmp);	    // The first point is obviously the origin.

    for (int i = 0; i < pts.size(); ++i)
    {
     	if(pts[i] != lefbtm)
        {
     		tmp.set_point(pts[i].get_point_cartesian().first,
                           pts[i].get_point_cartesian().second,
     		               lefbtm.get_point_cartesian().first,
                           lefbtm.get_point_cartesian().second);
     		shifted_pts.push_back(tmp);
     	}
    }

    /* Sorts points according to polar coordinates. */

    sort_points(shifted_pts.begin(), shifted_pts.end(), true);
    /* Traverses boundary counterclockwise, in sets of three points.
    Finds convex hull points based on sign of triangular area. */

    int j = 0, k = 0;
    while(j<shifted_pts.size())
    {
        if((j+1) < shifted_pts.size() and
           shifted_pts[j].get_point_polar().second ==
           shifted_pts[(j+1)].get_point_polar().second)
        {
            j++;
            continue;
        }
        if(k < 2)
        {
            answer.push_back(shifted_pts[j]);
            j++,k++;
            continue;
        }
        if(twice_triangle_signed_area(answer[k-1],answer[k],shifted_pts[j]) > 0)
        {
            answer.push_back(shifted_pts[j]);
            j++,k++;
            continue;
        }
        else
        {
            k--;
            answer.pop_back();
        }
    }
    vector<point> ret;
    for (int i = 0; i < answer.size() ; i++)
    {
        auto xx = answer[i].get_point_unnormalized(lefbtm);
        point aa(xx.first, xx.second);
    	ret.push_back(aa);
    }
    return ret;
}

vector<point> jarvis_march(const vector<point>& points)
{
    /*
    This function calculates the convex hull of a given set of points,
    based on Jarvis March Algorithm.

    Parameters:
    points - a vector of points whose convex hull is to be found.

    Return Value:
    A vector of points representing the convex hull in anti-clockwise order.*/

    vector<point> jarvis;
    // Vector to store the points that are part of the convex hull
    int flag = 0;
    point first_point;
    first_point.set_point(0.0, 0.0);
    // stores the coordinates for the first point
    pair< long double, long double> mp_cord = first_point.get_point_cartesian();

    for(int i = 0; i < points.size(); i++){
        if(flag == 0){  // if detected for the first time
            first_point = points[i];
            flag=1;
            mp_cord = points[i].get_point_cartesian();
            continue;
        }
        else{
            auto temp = points[i].get_point_cartesian();
            if(temp.second < mp_cord.second){
                first_point = points[i];
                mp_cord = temp;
            }
            else if(fabs(temp.second - mp_cord.second) < eps &&
                    temp.first < mp_cord.first){
                first_point = points[i];
                mp_cord = temp;
            }
        }
    }
    // push the bottomost leftmost points into the vector
    jarvis.push_back(first_point);

    long double firstpt_x = mp_cord.first, firstpt_y = mp_cord.second;
    flag = 0 ;
    point cur_point, origin;
    cur_point.set_point(0.0 ,0.0);
    angle min_angle;
    min_angle.set_angle(0.0);
    origin.set_point(firstpt_x, firstpt_y);
    for(int i = 0; i < points.size(); i++){
        if(points[i] != first_point){

            long double x = points[i].get_point_cartesian().first;
            long double y = points[i].get_point_cartesian().second;

            point temp = points[i];
            // to check for 2nd point of the convex hull
            temp.set_point(x, y, firstpt_x, firstpt_y);

            if(flag == 0){      // if encountered for the first time
                min_angle.set_angle(temp.get_point_cartesian().second,
                                    temp.get_point_cartesian().first);
                cur_point = temp;
                flag = 1;
            }

            else{
                angle temp_angle;
                temp_angle.set_angle(temp.get_point_cartesian().second,
                                    temp.get_point_cartesian().first);

                if(temp_angle < min_angle){
                    min_angle = temp_angle;
                    cur_point = temp;
                }
                else if (temp_angle == min_angle &&
                temp.get_radius_squared() > cur_point.get_radius_squared()){
                    // collinear points
                    min_angle = temp_angle;
                    cur_point = temp;
                }
            }
        }
    }
    // change coordinates wrt first_point
    cur_point.set_point(cur_point.get_point_unnormalized(origin).first,
    cur_point.get_point_unnormalized(origin).second);
    jarvis.push_back(cur_point);    // push second point of the hull

    while(true) // keep looping till all points of the hull are detected
    {
        int check = 0,flag = 0;
        point nxt_point;  // to check for next point of the hull
        nxt_point.set_point(0.0, 0.0);
        angle min_angle;
        min_angle.set_angle(0.0);
        origin.set_point(cur_point.get_point_cartesian().first,
                        cur_point.get_point_cartesian().second);
        for(int i=0; i<points.size();i++){
            if(points[i]!= first_point && points[i] != cur_point){

                if(twice_triangle_signed_area(first_point,cur_point,points[i])
                    > 0){   // to check if they are in anti-clockwise order
                    point temp;
                    long double x = points[i].get_point_cartesian().first;
                    long double y = points[i].get_point_cartesian().second;
                    temp.set_point(x, y, cur_point.get_point_cartesian().first,
                                    cur_point.get_point_cartesian().second);
                    check = 1;

                    if(flag == 0)       // if encountered for the first time
                    {
                        min_angle.set_angle(temp.get_point_cartesian().second,
                                            temp.get_point_cartesian().first);
                        nxt_point = temp;
                        flag = 1;
                    }

                    else{
                        angle temp_angle;
                        temp_angle.set_angle(temp.get_point_cartesian().second,
                                            temp.get_point_cartesian().first);

                        if(temp_angle < min_angle){
                            min_angle = temp_angle;
                            nxt_point = temp;
                        }
                        else if (temp_angle == min_angle &&
                                temp.get_radius_squared() >
                                nxt_point.get_radius_squared()){
                            // if collinear points
                            min_angle = temp_angle;
                            nxt_point = temp;
                        }
                    }
                }
            }
        }
        //no point left to be detected
        //and all the detected points satisfy the signed triangle rule
        if(check == 0)   break;
        // nxt_point to become cur_point for next iteration
        cur_point.set_point(nxt_point.get_point_unnormalized(origin).first,
                        nxt_point.get_point_unnormalized(origin).second);
        jarvis.push_back(cur_point);
    }
    return jarvis;  // return the points that are part of the convex hull
}

pair<point, point> upper_bridge(const vector<point>& points, long double a)
{
    /*
    This function calculates the upper bridge of a given set of points,
    across a vertical line, based on Kirkpatrick Seidel Algorithm.
	Since this function is for internal use, the user should use with care.

    Parameters:
    points - a vector of points whose upper bridge is to be found.
    a - The vertical line across which bridge is to be found is x = a.

    Return Value:
    A pair of points representing the upper bridge.*/

    vector<point> pts;
    for(auto p:points) pts.push_back(p);
    if(pts.size() == 2)
    {
        sort_points(pts.begin(), pts.end());
        return make_pair(pts[0], pts[1]);
    }
    vector<point> cand;
    auto it = pts.begin();
    vector< pair<point, point> > pairs;
    //Pairing up points
    while(true)
    {
        auto p = *(it++), q = *it;
        it++;
        if(less_than_cartesian(q, p))
            swap(p, q);
        pairs.push_back(make_pair(p, q));
        if(it==pts.end())
            break;
        it++;
        if(it==pts.end())
        {
            it--;
            cand.push_back(*it);
            break;
        }
        it--;
    }
    vector< pair< angle, pair<point, point> > > slopes;
    //Calculating slopes of line segment formed using each pair of points.
    for(auto p:pairs)
    {
        if(fabs(p.first.get_point_cartesian().first -
                p.second.get_point_cartesian().first) < eps)
        {
            if(p.first.get_point_cartesian().second >
               p.second.get_point_cartesian().second)
                cand.push_back(p.first);
            else
                cand.push_back(p.second);
        }
        else
        {
            long double x1 = p.first.get_point_cartesian().first;
            long double y1 = p.first.get_point_cartesian().second;
            long double x2 = p.second.get_point_cartesian().first;
            long double y2 = p.second.get_point_cartesian().second;
            angle ang;
            if(fabs(y2-y1) < eps)
                ang.set_angle(0.0, abs(x2-x1));
            else if(y1 - y2 < 0)
                ang.set_angle(y2-y1, x2-x1);
            else
                ang.set_angle(y1-y2, x1-x2);
            slopes.push_back(make_pair(ang, p));
        }
    }
    if(slopes.empty())
        return upper_bridge(cand, a);
    //Finding median slope.
    nth_element(slopes.begin(), slopes.begin()+slopes.size()/2, slopes.end(),
    [](pair< angle, pair<point, point> > a,
       pair< angle, pair<point, point> > b) -> bool
    {
        return a.first.get_tan() < b.first.get_tan();
    });
    angle median(slopes[slopes.size()/2].first.get_ratio().first,
                 slopes[slopes.size()/2].first.get_ratio().second);
    //Partitioning pairs based on their slopes.
    vector< pair< angle, pair<point, point> > > small, equal, large;
    for(auto p:slopes)
    {
        if(p.first.get_tan() < median.get_tan())
            small.push_back(p);
        else if(fabs(p.first.get_tan()-median.get_tan()) < eps)
            equal.push_back(p);
        else
            large.push_back(p);
    }
    //Finding supporting line.
    long double mx = pts[0].get_point_cartesian().second -
                     median.get_tan()*pts[0].get_point_cartesian().first;
    for(auto p: pts)
    {
        long double tmp = p.get_point_cartesian().second -
                          median.get_tan()*p.get_point_cartesian().first;
        mx = max(mx, tmp);
    }
    vector< point > maxpts;
    for(auto p:pts)
    {
        long double tmp = p.get_point_cartesian().second -
                          median.get_tan()*p.get_point_cartesian().first;
        if(fabs(tmp - mx) < eps)
            maxpts.push_back(p);
    }
    point pk = maxpts[0], pm = maxpts[0];
    for(auto p:maxpts)
    {
        if(less_than_cartesian(p, pk))
            pk = p;
        if(!less_than_cartesian(p, pm))
            pm = p;
    }
    if((!(pk.get_point_cartesian().first > a)) and
        pm.get_point_cartesian().first > a)
        return make_pair(pk, pm);
    //Eliminating invalid points.
    if(pk.get_point_cartesian().first > a)
    {
        for(auto p: small)
            cand.push_back(p.second.first);
        for(auto p: equal)
            cand.push_back(p.second.first);
        for(auto p: large)
        {
            cand.push_back(p.second.first);
            cand.push_back(p.second.second);
        }
    }
    else if(!(pm.get_point_cartesian().first > a))
    {
        for(auto p: large)
            cand.push_back(p.second.second);
        for(auto p: equal)
            cand.push_back(p.second.second);
        for(auto p: small)
        {
            cand.push_back(p.second.first);
            cand.push_back(p.second.second);
        }
    }
    //Recursively finding upper bridge from remaining candidates.
    return upper_bridge(cand, a);
}

vector<point> upper_hull(const vector<point>& points, point pmin, point pmax)
{
    /*
    This function calculates the upper hull of a given set of points,
    based on Kirkpatrick Seidel Algorithm.
	Since this function is for internal use, the user should use with care.

    Parameters:
    points - a vector of points whose upper hull is to be found.
    pmin - left extreme point.
    pmax - right extreme point.

    Return Value:
    A vector of points representing the upper hull from left to right.*/

    if(pmin == pmax)
    {
        vector<point> tmp;
        tmp.push_back(pmin);
        return tmp;
    }
    vector<point> pts;
    for(auto p:points)
        pts.push_back(p);
    if(pts.size() <= 2)
    {
        sort_points(pts.begin(), pts.end());
        return pts;
    }
    //Finding middle point.
    nth_element(pts.begin(), pts.begin() +
               (pts.size())/2, pts.end(), less_than_cartesian);
    long double a = pts[pts.size()/2].get_point_cartesian().first;
    //Finding the bridge.
    auto br = upper_bridge(pts, a);
    vector<point> tleft, tright;
    tleft.push_back(br.first);
    if(pmin != br.first)
        tleft.push_back(pmin);
    tright.push_back(br.second);
    if(pmax != br.second)
        tright.push_back(pmax);
    //Removing points not on hull.
    for(auto p:pts)
    {
        long double aa, bb, cc;
        aa = twice_triangle_signed_area(pmin, br.first, p);
        bb = twice_triangle_signed_area(br.first, br.second, p);
        cc = twice_triangle_signed_area(br.second, pmax, p);;
        if(!(aa > 0 or bb > 0 or cc > 0)) continue;
        if(p.get_point_cartesian().first > a)
            tright.push_back(p);
        else
            tleft.push_back(p);
    }
    //Recursively finding the upper hull and joining the two parts.
    auto h1 = upper_hull(tleft, pmin, br.first);
    auto h2 = upper_hull(tright, br.second, pmax);
    for(auto p:h2)
        h1.push_back(p);
    return h1;
}

vector<point> kirkpatrick_seidel(const vector<point>& points)
{
    /*
    This function calculates the convex hull of a given set of points,
    based on Kirkpatrick Seidel Algorithm.

    Parameters:
    points - a vector of points whose convex hull is to be found.

    Return Value:
    A vector of points representing the convex hull in anti-clockwise order.*/

    vector<point> pts;
    //Finding extreme points for upper hull.
    for(auto p: points)
        pts.push_back(p);
    point pmin(pts[0].get_point_cartesian().first,
               pts[0].get_point_cartesian().second);
    point pmax(pts[0].get_point_cartesian().first,
               pts[0].get_point_cartesian().second);
    for(auto p:pts)
    {
        auto tmp1 = pmin.get_point_cartesian();
        auto tmp2 = p.get_point_cartesian();
        if(tmp1.first > tmp2.first)
            pmin.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
        else if(fabs(tmp1.first - tmp2.first) < eps and
                tmp1.second < tmp2.second)
            pmin.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
        if(less_than_cartesian(pmax, p))
            pmax.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
    }
    //Finding the upper hull.
    auto uh = upper_hull(pts, pmin, pmax);
    pts.clear();
    //Finding extreme points for lower hull.
    for(auto p: points)
        pts.push_back(point(p.get_point_cartesian().first,
                      -1.0 * p.get_point_cartesian().second));
    pmin.set_point(pts[0].get_point_cartesian().first,
               pts[0].get_point_cartesian().second);
    pmax.set_point(pts[0].get_point_cartesian().first,
               pts[0].get_point_cartesian().second);
    for(auto p:pts)
    {
        auto tmp1 = pmin.get_point_cartesian();
        auto tmp2 = p.get_point_cartesian();
        if(tmp1.first > tmp2.first)
            pmin.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
        else if(fabs(tmp1.first - tmp2.first) < eps and
                tmp1.second < tmp2.second)
            pmin.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
        if(less_than_cartesian(pmax, p))
            pmax.set_point(p.get_point_cartesian().first,
                           p.get_point_cartesian().second);
    }
    //Finding the lower hull. This is performed by negating the points'
    //y coordinates and finding the upper hull.
    auto lh = upper_hull(pts, pmin, pmax);
    for(int i = 0; i < lh.size(); i++)
    {
        lh[i].set_point(lh[i].get_point_cartesian().first,
                      -1.0 * lh[i].get_point_cartesian().second);
    }
    reverse(uh.begin(), uh.end());
    //Joining upper and lower hulls.
    for(auto p: uh)
    {
        if((p != (*(--lh.end()))) && (p != (*lh.begin())))
            lh.push_back(p);
    }
    return lh;
}
