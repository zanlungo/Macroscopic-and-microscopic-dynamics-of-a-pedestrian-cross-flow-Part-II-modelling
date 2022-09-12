#include "group.hpp"
#include "ellipse.hpp"
#include "graphictools.hpp"
#include <iostream>  
#include <fstream> 
#include<math.h>                                                    
#include<cstdlib>
#include<SDL2/SDL.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <chrono>
#include <thread>


using namespace std;

struct stat st = {0};

//to compile (SDL2 needed)
/*g++ -O3 -o simulator_cir simulator_cir.cpp graphictools.cpp `sdl2-config --cflags --libs`*/

double EA=sqrt(0.0225); //ellipse major axis
double EB=EA; //ellipse minor
double ROB=0.02;   //size of an obstacle scan in graphics
double RHO=1.5;  //the following parameters are used to define the environment, but many of them are just for graphics purpose
double L;
double W;
double L_2;
double W_2;
double LM;
double LP;
double LS;
double x_max;  
double y_max;
double x_min;
double y_min;
double XMAX;
double YMAX;
double XMAX_2;
double YMAX_2;
double x_max_g;   
double y_max_g;
double x_min_g;
double y_min_g;
double XMAX_g;
double YMAX_g;
double XMAX_2_g;
double YMAX_2_g;
double DMC_x,DMC_y;
int ripet;  //number of experiment repetitions
double T_int=1.;

class Colore //ball colour for graphics
{
public:
  int r,g,b;//red green blue
  Colore()
  {
    r=0;
    g=0;
    b=0;
  }
  Colore(int r_,int g_,int b_)//initializes
  {
    r=r_;
    g=g_;
    b=b_;
  }
  void Set(int r_,int g_,int b_)//sets the values
  {
    r=r_;
    g=g_;
    b=b_;
  }
  void Copy(Colore C)//copies
  {
    r=C.r;
    g=C.g;
    b=C.b;
  }
}; 

void SDL_Inizio();     //graphics 

Colore *colori;   //other functions related to graphics
void DrawBallnp(int r,int g,int b,int rg,int x1,int y1);     //graphics with nonperiodic conditions (just a ball)
void DrawEllnp(int r,int g,int b,double A,double B,double theta,double x0[2],double XMAX,double YMAX);
void DrawLine(int r,int g,int b,double x0,double x1,double v0,double v1,double XMAX,double YMAX);
void DrawRectangle(int r,int g,int b,int x, int y,int w, int h);
void SDL_Update();

void PAINT()//this is for graphics
{
  double W_2_w=W_2+0.1; 
  DrawRectangle(0, 0, 0, 0, 0, PIXEL_X, PIXEL_Y);
  DrawRectangle(255, 255, 255, 0, (YMAX_2_g-W_2_w)*PIXEL_Y/YMAX_g, PIXEL_X, 2*W_2_w*PIXEL_Y/YMAX_g);
  DrawRectangle(255, 255, 255, (XMAX_2_g-W_2_w)*PIXEL_X/XMAX_g, 0, 2*W_2_w*PIXEL_X/XMAX_g, PIXEL_Y);
}

void Drawcorpo(Ellisse E,int g)  //this is for graphics  (draws robot body)
{
  double px[2];
  px[0]=E.x[0]-x_min_g;
  px[1]=E.x[1]-y_min_g;   
  DrawEllnp(colori[g].r,colori[g].g,colori[g].b,E.A,E.B,E.theta,px,XMAX_g,YMAX_g); 
}


void Drawcorpo(Vector2D S,double raggio)  //this is for graphics (draws obstacles) 
{
  int RA=int(raggio*PIXEL_Y/YMAX_g);//body dimension in pixel 
  int xi=int((S.x-x_min_g)*PIXEL_X/XMAX_g);
  int yi=int((YMAX_g-S.y+y_min_g)*PIXEL_Y/YMAX_g); 
  DrawBallnp(255,255,0,RA,xi,yi);
}  

void DrawScene(int n,Ellisse *E,int *ids,int sizeo,Vector2D *O)     //this is for graphics  (full scene)                
{
  PAINT();
  for(int k=0;k<n;k++) Drawcorpo(E[k],ids[k]);  //pedestrains
  for(int k=0;k<sizeo;k++) Drawcorpo(O[k],ROB);  //obstacles      
}                   



class Ped_min //classes with minimal information for pedestrians
{
public:
  Vector2D r,v;
  bool direction;
  bool gyro;
  double theta;
};

class Ped_min_s //used to handle initial conditions
{
public:
  Ped_min Pm;
  double t;
  int ts;
};

Ped_min_s **Starting;


double Check_lt(double tt)
{
  if(tt>=M_PI) tt-=2*M_PI;
  else if(tt<-M_PI) tt+=2*M_PI;
  return tt;
}


int n_peds;   //number of pedestrians outside of group
int nobs;  //hard-coded number of scan points in the map
double DT; //time step,time interval 
double TI=40;  //time interval 
Vector2D *Obstacles;   //position of obstacle scans



//functions for random noise

void start();

double random_0_1()  //real random between 0 and 1
{
  double a;
  a=(double)rand()/RAND_MAX;
  return a;
}

int random(int i);//random integer up to i

double Try_Gauss(double sigma)//gaussian distributed number with square deviation sigma and mean 0
{
  float x1,x2,w,y1,y2;
  do 
    {
      x1=2.0*random_0_1()-1.0;
      x2=2.0*random_0_1()-1.0;
      w=x1*x1+x2*x2;
    } 
  while (w>=1.0);
  w=sqrt((-2.0*log(w))/w);
  y1=x1*w;
  y2=x2*w;
  return y1*sigma;
}

Wall P[8];//walls in the environment

void build_map()//builds the map with walls
{
  double delta=0.05;
  double delta_2=delta/2;
  double npx=L/delta;
  ofstream map("map.dat");
  nobs=0;  
  for(int i=0;i<npx;i++)
    for(int j=0;j<4;j++)
      {
	double px=i*delta+x_min+delta_2;
	if((px<LM)||(px>LP)){
	  double py=j*delta+delta_2+W_2;
	  map << px << " " << py << endl;
	  nobs++;
	  py=-j*delta+delta_2-W_2;
	  map << px << " " << py << endl;
	  nobs++;
	}
      }
  for(int j=0;j<npx;j++)
    for(int i=0;i<4;i++)
      {
	double py=j*delta+y_min+delta_2;
	if((py<(LM-L_2))||(py>(LP-L_2))){
	  double px=i*delta+delta_2+W_2+L_2;
	  map << px << " " << py << endl;
	  nobs++;
	  px=-i*delta+delta_2-W_2+L_2;
	  map << px << " " << py << endl;
	  nobs++;
	}
      }
  Obstacles=new Vector2D[nobs];
}

void Init_walls()//builds the walls in the environment
{
  LS=n_peds/(2*RHO*3);
  L=3+2*(LS);
  L_2=L/2.;
  W_2=W/2.;
  LP=L_2+W_2;
  LM=L_2-W_2;
  x_max=L;
  x_max_g=L_2+1.5*W;
  y_max_g=1.5*W;
  y_max=L_2;
  x_min=0;
  x_min_g=L_2-1.5*W;
  y_min=-L_2;
  y_min_g=-1.5*W;  
  XMAX=x_max-x_min;
  YMAX=y_max-y_min;
  XMAX_g=x_max_g-x_min_g;
  YMAX_g=y_max_g-y_min_g;
  XMAX_2=XMAX/2;
  YMAX_2=YMAX/2;
  XMAX_2_g=XMAX_g/2;
  YMAX_2_g=YMAX_g/2;
  DMC_x=L_2-1.5;
  DMC_y=-1.5;  
  build_map();
  ifstream map("map.dat");//reads position of obstacle scans from the map
  for(int i=0;i<nobs;i++) 
    {
      double x,y;
      map >> x;
      map >> y;
      Obstacles[i].Init(x,y);
    }
  double W_w=W+0.2;
  double W_2_w=W_2+0.1;
  double LM_w=LM-0.1;
  double LP_w=LP+0.1;
  colori=new Colore[n_peds];   
  P[0].Init(0,0,-W_2_w,0,LM_w);
  P[1].Init(1,0,W_2_w,0,LM_w);
  P[2].Init(0,0,-W_2_w,LP_w,L);
  P[3].Init(1,0,W_2_w,LP_w,L);
  P[4].Init(0,1,LM_w,-L_2,-W_2_w);
  P[5].Init(1,1,LP_w,-L_2,-W_2_w);
  P[6].Init(0,1,LM_w,W_2_w,L_2);
  P[7].Init(1,1,LP_w,W_2_w,L_2);  
}


void Run(Genome &genome,int np)//simulation code
{
  bool stop=false;
  Ellisse *Ellipses;
  Ellipses=new Ellisse[np];//pedestrian bodies
  Pedestrian *singleped;
  singleped=new Pedestrian[np]; //pedestrians  
  Vector2D *sgoals;//pedestrian goals
  State2D *pass2s;//to pass information on other pedestrians
  pass2s=new State2D[np-1];
  int *pass_id;//to handle pedestrian ids (flow they belong)
  pass_id=new int[np];
  int *sdirections;
  sdirections=new int[np];
  sgoals=new Vector2D[np];
  chrono::system_clock::time_point starttime;//for playing in real time
  for(int rip=0;rip<ripet;rip++)//repetitions in the experiment
    {     
      double phys_coll=0;    //amount of collision  
      char name[200];
      if(stop) {cout << "stop!" << endl;break;}	  
      Init_walls();	  
      for(int i=0;i<np;i++)//initialises pedestrian (velocity, flow they belong to)
	{
	  singleped[i].active=false;
	  singleped[i].Init(genome);
	  singleped[i].genome.Fix_Vel(singleped[i].genome.vp+Try_Gauss(singleped[i].genome.GS));
	  sdirections[i]=!(Starting[rip][i].Pm.direction);
	  if(sdirections[i]) colori[i].Set(0,255,0);
	  else colori[i].Set(255,0,255);
	}
      double t;
      for(t=0;t<TI;)//for the integration time
	{
	  int active_ell=0;
	  for(int i=0;i<np;i++)  //for all  pedestrians: if actvive (in the simulation area) go on simulating
	    {
	      if(singleped[i].active)
		{	      
		  Ellipses[active_ell].Copy_id(singleped[i].E);
		  pass_id[active_ell]=i;
		  active_ell++;
		}
	    }	  
	  for(int i=0;i<np;i++)  //for all  pedestrians: if not active but entering in the next step: take initial conditions (checks for possible overlapping)
	    {
	      if(!singleped[i].active)
		{
		  if(fabs(Starting[rip][i].t-t)<(2*DT/3))
		    {
		      singleped[i].active=true;
		      double x=Starting[rip][i].Pm.r.x+DMC_x;
		      double y=Starting[rip][i].Pm.r.y+DMC_y;
		      double cx=0;
		      double cy=0;
		      Ellisse E;
		      while(1)
			{
			  E.Initialise(x,y,Starting[rip][i].Pm.v.x,Starting[rip][i].Pm.v.y,Starting[rip][i].Pm.theta,0,EA,EB,i);
			  bool wall=false;
			  double r,rp,rm;
			  double l[2];
			  for(int k=0;k<8;k++) {wall=wall||E.Overlap_parete(P[k],r,rp,rm);}
			  bool again=wall;
			  for(int k=0;k<active_ell;k++) again=again||Overlap(l,E,Ellipses[k]);
			  if(!again) break;
			  else
			    {
			      if(sdirections[i])
				{
				  if(!wall) {cy=-0.1;cx=Try_Gauss(0.1);}
				  if(wall)
				    {
				      if(x>L_2) {cy=-0.1;cx=-0.1;}
				      else {cy=-0.1;cx=0.1;}
				    }
				}
			      else
				{
				  if(!wall) {cx-=0.1;cy=Try_Gauss(0.1);}
				  if(wall)
				    {
				      if(y>0) {cy=-0.1;cx=-0.1;}
				      else {cy=0.1;cx=-0.1;}
				    }
				}
			      x+=cx;
			      y+=cy;
			    }
			}			  
		      singleped[i].Init(x,y,Starting[rip][i].Pm.v.x,Starting[rip][i].Pm.v.y,EA,EB,singleped[i].genome.EBV,Starting[rip][i].Pm.theta,0,i);
		      Ellipses[active_ell].Copy_id(singleped[i].E);
		      active_ell++;
		    }
		}
	    }		  
	  active_ell=0;
	  for(int i=0;i<np;i++)  //for all independent pedestrians: builds copies of bodies for physical simulation
	    {
	      if(singleped[i].active)
		{	      
		  Ellipses[active_ell].Copy_id(singleped[i].E);
		  pass_id[active_ell]=i;
		  active_ell++;
		}
	    }
	  for(int i=0;i<np;i++)
	    {
	      if(singleped[i].active)
		{
		  if(sdirections[i]) //updates goals
		    {	  
		      if((singleped[i].Self.r.x<(LM+0.5))&&(singleped[i].Self.r.y<W_2)&&(singleped[i].Self.r.y>-W_2))
			{
			  double dx=LM+0.5-singleped[i].Self.r.x;
			  double dy=W_2-singleped[i].Self.r.y;
			  sgoals[i].Init(2*(dx+Try_Gauss(singleped[i].genome.GD2)),dy+Try_Gauss(singleped[i].genome.GD2));
			}
		      else if((singleped[i].Self.r.x>(LP-0.5))&&(singleped[i].Self.r.y<W_2)&&(singleped[i].Self.r.y>-W_2))
			{
			  double dx=LP-0.5-singleped[i].Self.r.x;
			  double dy=W_2-singleped[i].Self.r.y;
			  sgoals[i].Init(2*(dx+Try_Gauss(singleped[i].genome.GD2)),dy+Try_Gauss(singleped[i].genome.GD2));
			}
		      else sgoals[i].Init(0+Try_Gauss(singleped[i].genome.GD2),1+Try_Gauss(singleped[i].genome.GD2));
		    }
		  else
		    {
		      if((singleped[i].Self.r.y<(-W_2+0.5))&&(singleped[i].Self.r.x<LP)&&(singleped[i].Self.r.x>LM))
			{
			  double dy=-W_2+0.5-singleped[i].Self.r.y;
			  double dx=LP-singleped[i].Self.r.x;
			  sgoals[i].Init(dx+Try_Gauss(singleped[i].genome.GD2),2*(dy+Try_Gauss(singleped[i].genome.GD2)));
			}
		      else if((singleped[i].Self.r.y>(W_2-0.5))&&(singleped[i].Self.r.x<LP)&&(singleped[i].Self.r.x>LM))
			{
			  double dy=W_2-0.5-singleped[i].Self.r.y;
			  double dx=LP-singleped[i].Self.r.x;
			  sgoals[i].Init(dx+Try_Gauss(singleped[i].genome.GD2),2*(dy+Try_Gauss(singleped[i].genome.GD2)));
			}	
		      else sgoals[i].Init(1+Try_Gauss(singleped[i].genome.GD2),0+Try_Gauss(singleped[i].genome.GD2));
		    }	  
		  int pass=0;//a counter to pass correctly information (positions of others)
		  for(int j=0;j<np;j++)
		    {
		      if(singleped[j].active)
			{
			  if(j!=i)  {pass2s[pass]=singleped[j].Self;pass++;}
			}
		    } 
		  singleped[i].NextV(nobs,Obstacles,pass,pass2s,sgoals[i],DT,3,5,t,Ellipses,active_ell,P,8);//computes new velocity
		}
	    }	  
	  active_ell=0;
	  for(int i=0;i<np;i++)
	    {
	      if(singleped[i].active)
		{ 
		  singleped[i].IR(DT,2*M_PI,2*M_PI);//updates angular velocity in pedestrians
		  singleped[i].Update_onlyV(DT);    //updates velocity in pedestrians
		  Ellipses[active_ell].Copy(singleped[i].E);  //copies in physical bodies
		  pass_id[active_ell]=i;
		  active_ell++;
		}
	    }	  
	  bool done=false;
	  SDL_Event event;   //these next lines are for graphics
	  while ( SDL_PollEvent(&event) )        
	    {
	      if ( event.type == SDL_QUIT )
		{ 
		  done = true; 
		} 
	    }		
	  if(done) break;
	  DrawScene(active_ell,Ellipses,pass_id,nobs,Obstacles); //draws the scene  
	  if(t == 0)
	    starttime = chrono::system_clock::now();
	  else
	    while(true)
	      {
		const auto now = chrono::system_clock::now();
		const int dt = t*1000-chrono::duration_cast<chrono::milliseconds>(now-starttime).count();
		if(dt > 0)
		  this_thread::sleep_for(chrono::milliseconds(dt));
		else
		  break;
	      }
	  SDL_Update();
	  bool stop;
	  double phys_coll_memo=phys_coll;
	  stop=Physical_evolution(active_ell,Ellipses,phys_coll,8,P,DT,t);//simulation of the physical system
	  if(stop) {cout << "stop!" << endl;break;}	  
	  int step=int((t-DT)/T_int);	  
	  for(int i=0;i<active_ell;i++)
	    {
	      singleped[pass_id[i]].Init(Ellipses[i]);//copies back info on position and velocity after physical simulation
	    }
	  for(int i=0;i<np;i++)
	    {
	      if(singleped[i].active)
		{
		  if((singleped[i].Self.r.x>x_max_g)||(singleped[i].Self.r.y>y_max_g)) singleped[i].active=false;
		}		  
	    }
	  if((!active_ell)&&(t>10)) break;// stops simulation when everybody has passed the crossing area
	}
      delete [] Obstacles;	  
    }
  delete [] singleped;
  delete [] pass2s;
  delete [] sgoals;
  delete [] sdirections;
  delete [] Ellipses; 
  delete [] pass_id;
}


void Initialise(Genome &genome)
{
  string token;
  W=3.4;//environment size
  n_peds=56;//number of subjects
  DT=0.05;//integration step
  ripet=6;//number of experiment repetitions
  char nome[100];
  ifstream read;
  sprintf(nome,"best_cir.dat");
  read.open(nome);  //gets model parameters (from a GA optimisation) 
  read >> token;  //we keep the original form of the "genome" that includes parameters not used by this version of the model (e.g., for group behaviour)
  read >> genome.r_0;//not used
  read >> token;
  read >>genome.Cr;//not used
  read >> token;
  read >>genome.Ct;//not used
  read >> token;
  read >>genome.eta;//not used
  read >> token;
  read >>genome.k;
  read >> token;
  read >>genome.vp;
  read >> token;
  read >>genome.S0;//not used
  read >> token;
  read >>genome.CC;//not used
  read >> token;
  read >>genome.coll_rad;
  read >> token;
  read >>genome.velover;//not used
  read >> token;
  read >>genome.rmax;//not used
  read >> token;
  read >>genome.body_radius;
  read >> token;
  read >> genome.ext_disk;
  read >> token;
  read >>genome.body_diameter;
  read >> token;
  read >> genome.ext_disk_2;
  read >> token;
  read >>genome.extern_radius;
  read >> token;
  read >>genome.extern_diameter;
  read >> token;
  read >>genome.radius_view;
  read >> token;
  read >> genome.KT;
  read >> token;
  read >> genome.KTT;  
  read >> token;
  read >>genome.CBO;
  read >> token;
  read >>genome.MU;  
  read >> token; 
  read >>genome.TMAX;
  read >> token;
  read >>genome.BETA;
  read >> token;
  read >> genome.KVT;    
  read >> token;
  read >> genome.KSTEP;     
  read >> token;
  read >>genome.GD2;
  read >> token;
  read >> genome.GS;    
  read >> token;
  read >> genome.EBV;      
  read >> token;  
  read >> genome.SIG_DRIFT;//not used
  read.close();
  ifstream map("map.dat");//reads position of obstacle scans from the map
  for(int i=0;i<nobs;i++) 
    {
      double x,y;
      map >> x;
      map >> y;
      Obstacles[i].Init(x,y);
    }
  Starting=new Ped_min_s *[6];  
  for(int j=0;j<6;j++)
    {
      Starting[j]=new Ped_min_s[n_peds];      
    }  
  char name[200];
  ifstream leggi;
  double ll;
  for(int i=0;i<6;i++)
    {  
      sprintf(name,"tracking_data/start_6_%d.dat",i);
      leggi.open(name);
      for(int j=0;j<56;j++)
	{
	  leggi >> ll;
	  leggi >> Starting[i][j].ts;
	  leggi >> Starting[i][j].t;
	  Starting[i][j].t=int(Starting[i][j].t/DT)*DT;	  
	  leggi >> Starting[i][j].Pm.gyro;
	  leggi >> Starting[i][j].Pm.direction;
	  Starting[i][j].Pm.direction=Starting[i][j].Pm.direction;	  
	  leggi >> Starting[i][j].Pm.r.x;
	  leggi >> Starting[i][j].Pm.r.y; 
	  leggi >> Starting[i][j].Pm.v.x;
	  leggi >> Starting[i][j].Pm.v.y;
	  leggi >> Starting[i][j].Pm.theta;
	  leggi >> ll;	  
	}
      leggi.close();      
    }      
}


int main()
{
  Genome genome;  
  start();
  SDL_Inizio();//starts graphics  
  Initialise(genome);//reads solution 
  Run(genome,n_peds);//simulates
  return 0;
}

