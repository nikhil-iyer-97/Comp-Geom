#include <bits/stdc++.h>
#include "fundamentals.h"
#include "convexhull.h"
#include <GL/glut.h>
using namespace std;

vector<point> points,lines;

int rint()
{
    return rand()%100000-50000;
}

/*int main()
{
    srand(time(NULL));
    vector<angle> v;
    for(int i=0;i<10;i++) v.push_back(angle(rint(),rint()));
    for(auto a:v) cout<<a.get_value(false)<<" ";
    cout<<endl;
    //sort(v.begin(),v.end(),'a');
    for(auto a:v) cout<<a.get_value(false)<<" ";
    cout<<endl;
    vector<point> v;
    for(int i=0;i<3;i++) v.push_back(point(3,5));
    for(int i=0;i<3;i++) v.push_back(point(3,4));
    for(auto p:v) print_point(p);
    cout<<endl;
    sort_points(v.begin(),v.end(),true); //Polar sort
    cout<<"Sorted\n";
    for(auto p:v) print_point(p);
    cout<<endl;
    cout<<twice_triangle_signed_area(point(0,3),point(),point(4,0))<<endl;
    vector<point> v;
    //srand(time(NULL));
    //for(int i=0;i<5;i++) v.push_back(point(rint(), rint()));
    for(int i=0;i<10000;i++) v.push_back(point(rint()*1.1,rint()*1.1));
    for(point p:v) cout<<p.get_point_cartesian().first<<" "<<p.get_point_cartesian().second<<endl;
    cout<<"Hull\n";
    //for(auto p:v) print_point(p);
    auto x = kirkpatrick_seidel(v);
    //cerr<< x.size() <<endl;
    for(point p:x) cout<<p.get_point_cartesian().first<<" "<<p.get_point_cartesian().second<<endl;
    cout<<"Complete\n";
}
*/
void init2D(float r, float g, float b)
{
	glClearColor(r,g,b,0.0);
	glMatrixMode (GL_PROJECTION);
	gluOrtho2D (-25000.0, 25000.0, -25000.0, 25000.0);
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 1.0, 0.0);
    GLfloat thickness = 3.0;

    glPointSize(thickness);
    glBegin(GL_POINTS);
	for(int i = 0; i < points.size(); i++)
		glVertex2f(points[i].get_point_cartesian().first, points[i].get_point_cartesian().second);
	glEnd();


    glLineWidth(thickness);
	glBegin(GL_LINES);
	for(int i=0; i<lines.size(); i++)
	{
		glVertex2f(lines[i%lines.size()].get_point_cartesian().first, lines[i%lines.size()].get_point_cartesian().second);
		glVertex2f(lines[(i+1)%lines.size()].get_point_cartesian().first, lines[(i+1)%lines.size()].get_point_cartesian().second);
	}
	glEnd();
	glFlush();
}
int main(int argc, char** argv)
{
    //srand(time(NULL));
    for(int i=0; i<20; i++)
        points.push_back(point(rint()/2.61,rint()/4.22));
    /*points.push_back(point(1,1));
    points.push_back(point(2,2));
    points.push_back(point(3,3));
    points.push_back(point(2,4));
    points.push_back(point(-1,5));*/
    lines = graham_scan(points);
    //while(lines.size()>3) lines.pop_back();
     glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (1500, 1500);
	glutInitWindowPosition (1000, 1000);
	glutCreateWindow ("Convex Hull");
	init2D(0.0,0.0,0.0);
	glutDisplayFunc(display);
	glutMainLoop();
}
