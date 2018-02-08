#include <bits/stdc++.h>
#include "fundamentals.h"

using namespace std;

pair<point, point> upper_bridge(const vector<point>& points, long double a)
{
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
    nth_element(slopes.begin(), slopes.begin()+slopes.size()/2, slopes.end(),
    [](pair< angle, pair<point, point> > a,
       pair< angle, pair<point, point> > b) -> bool
    {
        return a.first.get_tan() < b.first.get_tan();
    });
    angle median(slopes[slopes.size()/2].first.get_ratio().first,
                 slopes[slopes.size()/2].first.get_ratio().second);
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
    return upper_bridge(cand, a);
}

vector<point> upper_hull(const vector<point>& points, point pmin, point pmax)
{
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
    nth_element(pts.begin(), pts.begin() +
               (pts.size())/2, pts.end(), less_than_cartesian);
    long double a = pts[pts.size()/2].get_point_cartesian().first;
    auto br = upper_bridge(pts, a);
    vector<point> tleft, tright;
    tleft.push_back(br.first);
    if(pmin != br.first)
        tleft.push_back(pmin);
    tright.push_back(br.second);
    if(pmax != br.second)
        tright.push_back(pmax);
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
    auto h1 = upper_hull(tleft, pmin, br.first);
    auto h2 = upper_hull(tright, br.second, pmax);
    for(auto p:h2)
        h1.push_back(p);
    return h1;
}

vector<point> kirkpatrick_seidel(const vector<point>& points)
{
    vector<point> pts;
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
    auto uh = upper_hull(pts, pmin, pmax);
    pts.clear();
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
    auto lh = upper_hull(pts, pmin, pmax);
    for(int i = 0; i < lh.size(); i++)
    {
        lh[i].set_point(lh[i].get_point_cartesian().first,
                      -1.0 * lh[i].get_point_cartesian().second);
    }
    reverse(uh.begin(), uh.end());
    for(auto p: uh)
    {
        if((p != (*(--lh.end()))) && (p != (*lh.begin())))
            lh.push_back(p);
    }
    return lh;
}

vector<point> graham_scan(const vector<point>& pts)
{
	if(pts.size()<=2) return pts;

	vector<point> answer;			// Final vector for answer. Contains set of points in the Hull, arranged CCW.

	/* Finds leftmost, bottommost point. */

    point lefbtm(pts[0].get_point_cartesian().first, pts[0].get_point_cartesian().second);
    for (int i = 0; i < pts.size(); ++i)
    {
    	if(lefbtm.get_point_cartesian().second > pts[i].get_point_cartesian().second){
    		lefbtm.set_point(pts[i].get_point_cartesian().first, pts[i].get_point_cartesian().second);
    	}
        else if(fabs(lefbtm.get_point_cartesian().second-pts[i].get_point_cartesian().second) < eps and
                lefbtm.get_point_cartesian().first > pts[i].get_point_cartesian().first)
            lefbtm.set_point(pts[i].get_point_cartesian().first, pts[i].get_point_cartesian().second);
    }

    print_point(lefbtm);
    /* Shifts origin to this point, finds all other points with respect to this. */

    vector<point> shifted_pts;
    point temp(lefbtm.get_point_cartesian().first, lefbtm.get_point_cartesian().second,
     lefbtm.get_point_cartesian().first, lefbtm.get_point_cartesian().second);
    shifted_pts.push_back(temp); // The first point is obviously the origin.
    print_point(temp);
    for (int i = 0; i < pts.size(); ++i)
    {
     	if(pts[i] != lefbtm ){
     		temp.set_point(pts[i].get_point_cartesian().first, pts[i].get_point_cartesian().second,
     		lefbtm.get_point_cartesian().first, lefbtm.get_point_cartesian().second);
     		shifted_pts.push_back(temp);
            //cout<<"***********\n";
            //print_point(pts[i]);
            //print_point(temp);
     	}
    }


    /* Sorts points according to polar coordinates. */

    sort_points(++shifted_pts.begin(), shifted_pts.end(), true);


    /* Traverses boundary counterclockwise, in sets of three points. */
    /* Finds and adds convex hull points to the answer, based on sign of triangular area. */

    answer.push_back(shifted_pts[0]);							// Initializes answer vector.

    int j = 1;			// Assumes that n >= 3.

    int k;
    while( j < shifted_pts.size() ){
        //cout<<j<<" "<<answer.size()<<endl;
        //print_point(shifted_pts[j]);
    	k = (int)answer.size() - 1;
        if(j != shifted_pts.size()-1)
        {
            if(shifted_pts[j].get_point_polar().second == shifted_pts[j+1].get_point_polar().second)
            {
                cout<<"same angle\n";
                j++;
                continue;
            }
        }
        if(k < 1)
        {
            cout<<"compulsary push\n";
            answer.push_back(shifted_pts[j++]);
            continue;
        }
    	if(twice_triangle_signed_area(answer[k-1], answer[k], shifted_pts[j]) > 0){			// Accepts if CCW turn.
            cout<<"valid pt\n";
            answer.push_back(shifted_pts[j++]);
    	}

    	else
        {
            cout<<"pop\n";
            answer.pop_back(); // Pops till it finds CCW turn.

        }

    }

    for (int i = 0; i < answer.size() ; ++i)
    {
    	answer[i].set_point(answer[i].get_point_unnormalized().first, answer[i].get_point_unnormalized().second);
    }

    return answer;
}


//vector<point> jarvis_march(const vector<point>&);
