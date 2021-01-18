/* $Id: graphics.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef GRAPHICS_H
#define GRAPHICS_H

using namespace std; //df

class Point;
class Line;
class Lattice;
class Hyperplane;

void message(const char* s);
void message(double d);
void add_orbit(const Point&);
void set_line(const Line&);
void show_lattice();
void set_lattice(Lattice*);
void show_hyperplane(const Hyperplane& H);
void check_event();
void wait_if_pause();

#endif
