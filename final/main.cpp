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
		glVertex2f(points[i].get_point_cartesian().first,
                   points[i].get_point_cartesian().second);
	glEnd();


    glLineWidth(thickness);
    glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for(int i=0; i<lines.size(); i++)
	{
		glVertex2f(lines[i%lines.size()].get_point_cartesian().first,
                   lines[i%lines.size()].get_point_cartesian().second);
		glVertex2f(lines[(i+1)%lines.size()].get_point_cartesian().first,
                   lines[(i+1)%lines.size()].get_point_cartesian().second);
	}
	glEnd();
	glFlush();
}

int main(int argc, char** argv)
{
    //Compile as g++ std=c++1y -lpthread -lGL -lGLU -lglut
    srand(time(NULL));
    for(int i=0; i<1000; i++)
        points.push_back(point(rint()/2.61,rint()/4.22));
    lines = kirkpatrick_seidel(points);
    glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (1500, 1500);
	glutInitWindowPosition (1000, 1000);
	glutCreateWindow ("Convex Hull");
	init2D(0.0,0.0,0.0);
	glutDisplayFunc(display);
	glutMainLoop();
}
