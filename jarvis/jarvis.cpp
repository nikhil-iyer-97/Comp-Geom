vector<point> jarvis_march(const vector<point>& points)
{
    vector<point> jarvis;
    int flag = 0;
    point first_point;
    first_point.set_point(0.0, 0.0);
    pair< long double, long double> mp_cord = first_point.get_point_cartesian();

    for(int i = 0; i < points.size(); i++){
        if(flag == 0){
            first_point = points[i];
            flag=1;
            mp_cord = points[i].get_point_cartesian();
            continue;
        }
        else{
            pair< long double, long double> temp = points[i].get_point_cartesian();
            if(temp.second < mp_cord.second){
                first_point = points[i];
                mp_cord = temp;
            }
            else if(temp.second == mp_cord.second && temp.first < mp_cord.first){
                first_point = points[i];
                mp_cord = temp;
            }
        }
    }
    jarvis.push_back(first_point);

    long double firstpt_x = mp_cord.first, firstpt_y = mp_cord.second;
    flag = 0 ;
    point cur_point;
    cur_point.set_point(0.0 ,0.0);
    angle min_angle;
    min_angle.set_angle(0.0);

    for(int i = 0; i < points.size(); i++){
        if(points[i].get_point_cartesian() != first_point.get_point_cartesian()){

            long double x = points[i].get_point_cartesian().first;
            long double y = points[i].get_point_cartesian().second;

            point temp = points[i];
            temp.set_point(x, y, firstpt_x, firstpt_y);

            if(flag == 0){
                min_angle.set_angle(temp.get_point_cartesian().second, temp.get_point_cartesian().first);
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
                    min_angle = temp_angle;
                    cur_point = temp;
                }
            }
        }
    }
    cur_point.set_point(cur_point.get_point_unnormalized().first,
    cur_point.get_point_unnormalized().second);
    jarvis.push_back(cur_point);

    while(true)
    {
        int check = 0,flag = 0;
        point nxt_point;
        nxt_point.set_point(0.0, 0.0);
        angle min_angle;
        min_angle.set_angle(0.0);

        for(int i=0; i<points.size();i++){
            if(points[i].get_point_cartesian() != first_point.get_point_cartesian()
            && points[i].get_point_cartesian() != cur_point.get_point_cartesian()){

                if(twice_triangle_signed_area(first_point,cur_point,points[i]) > 0){
                    point temp;
                    long double x = points[i].get_point_cartesian().first;
                    long double y = points[i].get_point_cartesian().second;
                    temp.set_point(x, y, cur_point.get_point_cartesian().first,
                                    cur_point.get_point_cartesian().second);
                    check = 1;

                    if(flag == 0)
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
                                temp.get_radius_squared() > nxt_point.get_radius_squared()){
                            min_angle = temp_angle;
                            nxt_point = temp;
                        }
                    }
                }
            }
        }
        if(check == 0)   break;
        cur_point.set_point(nxt_point.get_point_unnormalized().first,
                        nxt_point.get_point_unnormalized().second);
        jarvis.push_back(cur_point);
    }

    return jarvis;
}
