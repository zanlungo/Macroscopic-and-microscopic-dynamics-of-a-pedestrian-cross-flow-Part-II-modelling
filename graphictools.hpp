#ifndef __GRAPHICTOOLS_H__
#define __GRAPHICTOOLS_H__

#include <SDL2/SDL.h>

using namespace std;

extern int PIXEL_X, PIXEL_Y;


void SDL_Inizio(); //graphics
void SDL_Update();

void DrawBallnp(int r,int g,int b,int rg,int x1,int y1);
void DrawEllnp(int r,int g,int b,double A,double B,double theta,double x0[2],double XMAX,double YMAX);
void DrawLine(int r,int g,int b,double x0,double x1,double v0,double v1,double XMAX,double YMAX);
void DrawRectangle(int r,int g,int b,int x, int y,int w, int h);

#endif
