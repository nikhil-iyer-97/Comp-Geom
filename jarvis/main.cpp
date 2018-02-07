#include <bits/stdc++.h>
#include "fundamentals.cpp"
#include "jarvis.cpp"
#include <GL/glut.h>

using namespace std;

vector<point> points,lines;

int rint()
{
    return rand()%100000-50000;
}

void print_angle(angle a)
{
    cout<<"Degrees: "<<a.get_value(false)<<endl;
    cout<<"Radians: "<<a.get_value()<<endl;
    cout<<"Ratio (y x): "<<a.get_ratio().first<<" "<<a.get_ratio().second<<endl;
}

void print_point(point p)
{
    cout<<"Cartesian: "<<p.get_point_cartesian().first<<" "<<p.get_point_cartesian().second<<endl;
    cout<<"R squared: "<<p.get_radius_squared()<<endl;
    cout<<"Angle: \n";
    print_angle(p.get_point_polar().second);
    cout<<"Unnormalized: "<<p.get_point_unnormalized().first<<" "<<p.get_point_unnormalized().second<<endl;
}

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
    for(int i=0; i<1000; i++)
        points.push_back(point(rint()/2.61,rint()/4.22));
    lines = jarvis_march(points);
    glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (1500, 1500);
	glutInitWindowPosition (1000, 1000);
	glutCreateWindow ("Convex Hull");
	init2D(0.0,0.0,0.0);
	glutDisplayFunc(display);
	glutMainLoop();
}
