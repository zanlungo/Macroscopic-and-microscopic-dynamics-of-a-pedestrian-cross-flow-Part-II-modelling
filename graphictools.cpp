#include "graphictools.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <iostream> 
#include <math.h>  

using namespace std;

SDL_Window *window;
SDL_Renderer *renderer;
int PIXEL_X=1000;
int PIXEL_Y=1000;





void SDL_Inizio() //inizializza con dim X,Y
{
  if ( SDL_Init(SDL_INIT_AUDIO|SDL_INIT_VIDEO) < 0 )
    {
      printf("Unable to init SDL: %s\n", SDL_GetError());
      exit(1);
    }
  atexit(SDL_Quit);

  window = SDL_CreateWindow("Realistic Simulator", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, PIXEL_X, PIXEL_Y, SDL_WINDOW_OPENGL);
  renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
  SDL_RenderClear(renderer);
  SDL_RenderPresent(renderer);
  SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");
  SDL_RenderSetLogicalSize(renderer, PIXEL_X, PIXEL_Y);

  /*if ( TTF_Init() < 0 ) {
    printf("TTFcould not initialize! TTF_Error: %s\n", TTF_GetError());
    exit(1);
  }
  font = TTF_OpenFont("VeraMono.ttf", 36);*/
}





void DrawBallnp(int r,int g,int b,int rg,int x1,int y1)
{
  SDL_SetRenderDrawColor(renderer,r,g,b,255);
  for (int x=x1-rg;x<x1+rg;x++)
  {
    for (int y=y1-rg;y<y1+rg;y++) 
    {
      if ((((x-x1)*(x-x1)+(y-y1)*(y-y1))<(rg*rg))&&(((x>=0)&&(x<PIXEL_X)&&(y>=0)&&(y<PIXEL_Y)))) SDL_RenderDrawPoint(renderer,(x+PIXEL_X)%PIXEL_X,(y+PIXEL_Y)%PIXEL_Y);
    }   
  }     
}


void DrawEllnp(int r,int g,int b,double A,double B,double theta,double x0[2],double XMAX,double YMAX)
{
  double c=A/B;
  double ct=cos(theta);
  double st=sin(theta);
  int x1=int(x0[0]*PIXEL_X/(XMAX));
  int y1=int(x0[1]*PIXEL_Y/(YMAX));
  double ratio=B/XMAX;
  int rg=int(ratio*PIXEL_X)+1;
  if(c>1) {rg*=c;rg++;}
  SDL_SetRenderDrawColor(renderer,r,g,b,255);
  for (int x=x1-rg;x<x1+rg;x++)
    {
      for (int y=y1-rg;y<y1+rg;y++) 
	{
	  if((x>=0)&&(x<PIXEL_X)&&(y>=0)&&(y<PIXEL_Y)) 
	    {
	      double dx=x-x1;
	      double dy=y-y1;
	      double cx=dx*ct-dy*st;
	      double cy=dx*st+dy*ct;
	      cx/=c;
	      double l=sqrt(cx*cx+cy*cy);
	      l/=PIXEL_X;
	      if(l<ratio)
		{
		  SDL_RenderDrawPoint(renderer,x,y);
		}
	    }	        
	}
    }    
}

void DrawLine(int r,int g,int b,double x0,double x1,double v0,double v1,double XMAX,double YMAX)
{
  int X1=int(x0*PIXEL_X/XMAX);
  int Y1=int(x1*PIXEL_Y/YMAX);
  int X2=X1+int(v0*PIXEL_X/XMAX);
  int Y2=Y1+int(v1*PIXEL_Y/YMAX);
  SDL_SetRenderDrawColor(renderer,r,g,b, 255);
  SDL_RenderDrawLine(renderer, X1, Y1, X2, Y2);
}


void DrawRectangle(int r,int g,int b,int x, int y,int w, int h)
{
  SDL_SetRenderDrawColor(renderer,r,g,b, 255);
  SDL_Rect rect = {x, y, w, h};
  SDL_RenderFillRect(renderer, &rect);
}

void SDL_Update()
{
  SDL_RenderPresent(renderer);
}
