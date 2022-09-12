#ifndef __GROUP_H__    
#define __GROUP_H__ 
#include <iostream>  
#include <fstream> 
#include<math.h>                                                    
#include<cstdlib>
#include "ellipse.hpp"

//pedestrian and robot behaviour class (group and robot classes not used in the SAFETY 2022 paper)

using namespace std;

double NEAR=0.3;//not used
double NEAR2=0.6;
double NEAROV=0.75;

//random number calls

int roundmio(double r)     //rounds a double
{
  int i=int(r);
  if((r-i)>0.5) i++;
  return i;
}

double roundmio2(double r,int dig)     //rounds a double up to dig digits (used in order to limit the solution space of GA)
{
  dig--;
  dig=pow(10,dig);
  double l10=log10(r);
  double mult;
  if(l10<0)  mult=pow(10,-int(l10)+1);
  else mult=pow(10,-int(l10));
  r=r*mult;
  r=double(roundmio(r*dig))/(mult*dig);
  return r;
}


double random_01()  //real random between 0 and 1
{
  double a;
  a=(double)rand()/RAND_MAX;
  return a;
}

double randommio(double g)  //real random between 0 and 1
{
  double a;
  a=(double)rand()/RAND_MAX;
  return g*a;
}

int randomi(int i) /*random integer between 0 and i-1*/
{
  int a;
  double r;
  r=RAND_MAX/double(i);
  a=rand();
  if(a!=RAND_MAX) a=int(a/r);
  else a=randomi(i);
  return a;
} 


double GaussN(double sigma)   //gaussian distributed number with square deviation sigma and mean 0
{
  float x1,x2,w,y1,y2;
  do 
  {
    x1=2.0*random_01()-1.0;
    x2=2.0*random_01()-1.0;
    w=x1*x1+x2*x2;
  } 
  while (w>=1.0);
  w=sqrt((-2.0*log(w))/w);
  y1=x1*w;
  y2=x2*w;
  return y1*sigma;
}

double Check_th(double theta)   //keeps theta in a given range
{
  double th=theta;
  if(th>M_PI) th=-2*M_PI+th;
  else if(th<-M_PI) th=2*M_PI+th;
  return th;
}

class Vector2D //for operations on 2d vectors
{
public:
  double x,y,m,th;//cartesian and polar
  Vector2D()
  {
    x=0;
    y=0;
    m=0;
    th=0;
  }
  Vector2D(double xp,double yp)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=0;
    if(x<0) th=-th;
  }
  Vector2D(double xp,double yp,double theta)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=Check_th(theta);
    if(x<0) th=-th;
  }
  void BCX(double L)//boundary conditions
  {
    if(x>=L) x-=L;
    else if(x<0) x+=L; 
  }
  void BCX(Vector2D V,double L)
  {
    double L_2=L/2;
    if((x-V.x)>=L_2) x-=L;
    else if((x-V.x)<-L_2) x+=L; 
  }  
  bool Uguale(Vector2D U)//checks if equal
  {
    bool u=true;
    u=u&&(x==U.x);
    u=u&&(y==U.y);
    u=u&&(m==U.m);
    u=u&&(th==U.th);
    return u;
  }
  void Init(double xp,double yp)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=0;
    if(x<0) th=-th;
  }
  void Init(double xp,double yp,double theta)
  {
    x=xp;
    y=yp;
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    else th=Check_th(theta);
    if(x<0) th=-th;
  }
  void Check_theta()
  {
    th=Check_th(th);
  }
  void Set_zero(double theta)
  {
    x=0;
    y=0;
    m=0;
    th=theta;
    Check_theta();
  }
  void Set_zero()
  {
    x=0;
    y=0;
    m=0;
    th=0;
  }
  void Magnitude()
  {
    m=sqrt(x*x+y*y);
    if(m) th=acos(y/m);
    if(x<0) th=-th;
  }
  void Scale(double s)
  {
    x*=s;
    y*=s;
    Magnitude();
  }
  void Add(Vector2D A)
  {
    x+=A.x;
    y+=A.y;
    Magnitude();
  }
  void Noise(double n)
  {
    x+=GaussN(n);
    y+=GaussN(n);
    Magnitude();
  }
  void Print()
  {
    cout << "x=" << x << " y=" << y << " m=" << m << " th=" << th << endl;
  }
  void Rotate_Clock(double alpha)
  {
    double xp=x;
    double yp=y;
    x=xp*cos(alpha)+yp*sin(alpha);
    y=-xp*sin(alpha)+yp*cos(alpha);
    Magnitude();
  }
  void Rotate_CClock(double alpha)
  {
    double xp=x;
    double yp=y;
    x=xp*cos(alpha)-yp*sin(alpha);
    y=xp*sin(alpha)+yp*cos(alpha);
    Magnitude();
  }
};

double Scalar(Vector2D v1,Vector2D v2)//scalar product
{
  return v1.x*v2.x+v1.y*v2.y;
}

Vector2D Average(Vector2D v1,Vector2D v2)  //computes  average
{
  Vector2D Av;
  Av.x=(v1.x+v2.x)/2;
  Av.y=(v1.y+v2.y)/2;
  Av.Magnitude();
  return Av;
}

Vector2D Average(Vector2D v1,Vector2D v2,double w) //weighted average
{
  Vector2D Av;
  double wm=1-w;
  Av.x=v1.x*w+v2.x*wm;
  Av.y=v1.y*w+v2.y*wm;
  Av.Magnitude();
  return Av;
}

Vector2D Average(Vector2D v1,Vector2D v2,Vector2D v3)//average of 3
{
  Vector2D Av;
  Av.x=(v1.x+v2.x+v3.x)/3;
  Av.y=(v1.y+v2.y+v3.y)/3;
  Av.Magnitude();
  return Av;
}

Vector2D Average(Vector2D v1,Vector2D v2,Vector2D v3,double w)//weighted average
{
  Vector2D Av;
  double wm=(1-w)/2;
  Av.x=v1.x*w+v2.x*wm+v3.x*wm;
  Av.y=v1.y*w+v2.y*wm+v3.y*wm;
  Av.Magnitude();
  return Av;
}

Vector2D Sum(Vector2D v1,Vector2D v2)//sums
{
  Vector2D S;
  S.x=v1.x+v2.x;
  S.y=v1.y+v2.y;
  S.Magnitude();
  return S;
}

Vector2D Versor(Vector2D v)//normalises
{
  Vector2D R;
  R=v;
  if(R.m) R.Scale(1./R.m);
  return R;
}

Vector2D Ort_vers(Vector2D v)//normalises and rotates of pi/2
{
  Vector2D R(v.y,-v.x);
  if(R.m) R.Scale(1./R.m);
  return R;
}

double Projection(Vector2D v1,Vector2D v2)//of v1 along v2
{
  if(v2.m) return Scalar(v1,v2)/v2.m;
  else return 0;
}

Vector2D Clock_Rot(Vector2D v)//rotates clockwise of pi/2
{
  Vector2D R(v.y,-v.x);
  return R;
}

double Clock_Ort(Vector2D v1,Vector2D v2)//projection of v1 along v2 rotated
{
  if(v2.m) return Scalar(v1,Clock_Rot(v2))/v2.m;
  else return 0;
}

double Dist(Vector2D v1,Vector2D v2)
{
  double dx=v1.x-v2.x;
  double dy=v1.y-v2.y;
  return sqrt(dx*dx+dy*dy);
}

class State2D//position and velocity
{
public:
  Vector2D r,v;
  State2D(){}
  State2D(double rx,double ry,double vx,double vy)
  {
    r.Init(rx,ry);
    v.Init(vx,vy);
  }
  State2D(Vector2D rp,Vector2D vp)
  {
    r=rp;
    v=vp;
  }
  bool Uguale(State2D U)
  {
    bool u=true;
    u=u&&r.Uguale(U.r);
    u=u&&v.Uguale(U.v);
    return u;
  }
  void Init(double rx,double ry,double vx,double vy)
  {
    r.Init(rx,ry);
    v.Init(vx,vy);
  }
  void Init(double rx,double ry,double vx,double vy,double theta)
  {
    r.Init(rx,ry);
    v.Init(vx,vy,theta);
  }
  void Init(Vector2D rp,Vector2D vp)
  {
    r=rp;
    v=vp;
  }
  void Noise(double np,double nv)
  {
    r.Noise(np);
    v.Noise(nv);
  }
  void Print()
  {
    cout << "r:" << endl;
    r.Print();
    cout << "v:" << endl;
    v.Print();
  }
};

State2D CM(State2D S1,State2D S2)//centre of mass
{
  State2D C;
  C.r=Average(S1.r,S2.r);
  C.v=Average(S1.v,S2.v);
  return C;
}

State2D CM(State2D S1,State2D S2,double w)//weighted centre of mass
{
  State2D C;
  C.r=Average(S1.r,S2.r,w);
  C.v=Average(S1.v,S2.v,w);
  return C;
}

State2D CM(State2D S1,State2D S2,State2D S3)//centre of mass of 3 
{
  State2D C;
  C.r=Average(S1.r,S2.r,S3.r);
  C.v=Average(S1.v,S2.v,S3.v);
  return C;
}

State2D CM(State2D S1,State2D S2,State2D S3,double w)//weighted centre of mass of 3
{
  State2D C;
  C.r=Average(S1.r,S2.r,S3.r,w);
  C.v=Average(S1.v,S2.v,S3.v,w);
  return C;
}

State2D CM(State2D S,int no,State2D *So)//centre of mass of #no entities
{
  State2D C;
  C=S;
  for(int i=0;i<no;i++)
    {
      C.r.Add(So[i].r);
      C.v.Add(So[i].v);
    }
  C.r.Scale(1./(no+1));
  C.v.Scale(1./(no+1));
  return C;
}

State2D CM(State2D S,int no,State2D *So,double w)//weighted centre of mass of #no entities
{
  State2D C=S;
  C.r.Scale(w);
  C.v.Scale(w);
  double wm=(1-w)/no;
  State2D P;
  for(int i=0;i<no;i++)
    {
      P=So[i];
      P.r.Scale(wm);
      P.v.Scale(wm);
      C.r.Add(P.r);
      C.v.Add(P.v);
    }
  return C;
}

Vector2D Rel(Vector2D V1,Vector2D V2)//relative position of 1 with respect to 2
{
  Vector2D r;
  r.x=V1.x-V2.x;
  r.y=V1.y-V2.y;
  r.Magnitude();
  return r;
}

Vector2D Rel_CM(State2D S1,State2D S2)//relative position of 1 with respect to 2 in the CM frame (y along velocity)
{
  Vector2D r;
  State2D C=CM(S1,S2);
  r.x=S1.r.x-S2.r.x;
  r.y=S1.r.y-S2.r.y;
  r.Magnitude();
  double tx,ty;
  ty=Projection(r,C.v);
  tx=Clock_Ort(r,C.v);
  r.Init(tx,ty);
  r.Magnitude();
  return r;
}

Vector2D Rel_goal(State2D S1,State2D S2,Vector2D goal)//relative position of 1 with respect to 2 in the goal frame (y along goal)
{
  Vector2D r;
  r.x=S1.r.x-S2.r.x;
  r.y=S1.r.y-S2.r.y;
  r.Magnitude();
  double tx,ty;
  ty=Projection(r,goal);
  tx=Clock_Ort(r,goal);
  r.Init(tx,ty);
  r.Magnitude();
  return r;
}


void Order_CM(State2D S0,int no,State2D *S,int *order)//decides who is first neighbour based on the group velocity (not used in the corrent version in which when the goal is unknown I pass as a goal the velocity of the others) THIS IS FOR GROUP BEHAVIOUR!
{
  State2D C=CM(S0,no,S);
  double *x;
  x=new double[no+1];
  x[0]=Clock_Ort(S0.r,C.v);
  for(int i=1;i<no+1;i++) x[i]=Clock_Ort(S[i-1].r,C.v);
  for(int i=0;i<no+1;i++)
    {
      order[i]=0;
      for(int j=0;j<no+1;j++)
	{
	  if(j!=i) if(x[i]>x[j]) order[i]++;
	}
    }
  delete [] x;
}

void Order_goal(State2D S0,int no,State2D *S,Vector2D goal,int *order)//decides who is first neighbour used on the goal (that may be a weighted velocity, although this possibility is not used in this version) THIS IS FOR GROUP BEHAVIOUR!
{
  double *x;
  x=new double[no+1];
  x[0]=Clock_Ort(S0.r,goal);
  for(int i=1;i<no+1;i++) x[i]=Clock_Ort(S[i-1].r,goal);
  for(int i=0;i<no+1;i++)
    {
      order[i]=0;
      for(int j=0;j<no+1;j++)
	{
	  if(j!=i) if(x[i]>x[j]) order[i]++;
	}
    }
  delete [] x;
}


Vector2D Rotate(Vector2D f,Vector2D g)//goes back to the original frame
{
  Vector2D nf;
  Vector2D yv=g;
  if(g.m) yv.Scale(1./g.m);
  Vector2D xv=Clock_Rot(yv);
  nf.x=f.y*yv.x+f.x*xv.x;
  nf.y=f.y*yv.y+f.x*xv.y;
  nf.Magnitude();
  return nf;
}

class Gene    //used for ga optimisation
{
public:
  bool real;
  double value;
  bool bvalue;
  bool change;   //if one it is changed by the GA algorithm
  double max;          //max and min that value may assume
  double min;
  double delta;
  double range;          //of mutation
  double pmutate;   //probability of mutation
  Gene()
  {
    real=true;
    bvalue=false;
    value=0;
    max=0;
    change=false;
    min=0;
    delta=0;
    pmutate=0;
    range=0;
  }
  Gene(bool isreal,bool c,double mx,double mn,double p,double r,bool v)    //initialises real gene
  {
    real=isreal;
    change=c;
    if(!isreal)
      {
	if(!c)
	  {
	    bvalue=v;
	    pmutate=0;
	  }
	else
	  {
	    pmutate=p;
	    bvalue=randomi(2);
	  }
      }
    else
      {
	max=mx;
	if(c) 
	  {
	    min=mn;
	    delta=max-min;
	    pmutate=p;
	    range=delta*r;
	    value=min+random_01()*delta;
	    if(value!=0) value=roundmio2(value,4);
	  }
	else //fixed value
	  {
	    min=max;
	    delta=0;
	    pmutate=0;
	    range=0;
	    value=max;
	  }
      }
  }
  Gene(bool c,bool v,double p)    //initialises real gene
  {
    real=false;
    change=c;
    if(!c)
      {
	bvalue=v;
	pmutate=0;
      }
    else
      {
	pmutate=p;
	bvalue=randomi(2);
      }
  }
  Gene(bool c,double mx,double mn,double p,double r)    //initialises real gene
  {
    real=true;
    max=mx;
    change=c;
    if(c) 
      {
	min=mn;
	delta=max-min;
	pmutate=p;
	range=delta*r;
	value=min+random_01()*delta;
	if(value!=0) value=roundmio2(value,4);
      }
    else //fixed value
      {
	min=max;
	delta=0;
	pmutate=0;
	range=0;
	value=max;
      }
  }
  void Init(bool isreal,bool c,double mx,double mn,double p,double r,bool v)    //initialises real gene
  {
    real=isreal;
    change=c;
    if(!isreal)
      {
	if(!c)
	  {
	    bvalue=v;
	    pmutate=0;
	  }
	else
	  {
	    pmutate=p;
	    bvalue=randomi(2);
	  }
      }
    else
      {
	max=mx;
	if(c) 
	  {
	    min=mn;
	    delta=max-min;
	    pmutate=p;
	    range=delta*r;
	    value=min+random_01()*delta;
	    if(value!=0) value=roundmio2(value,4);
	  }
	else //fixed value
	  {
	    min=max;
	    delta=0;
	    pmutate=0;
	    range=0;
	    value=max;
	  }
      }
  }
  void Init(bool c,bool v,double p)    //initialises real gene
  {
    real=false;
    change=c;
    if(!c)
      {
	bvalue=v;
	pmutate=0;
      }
    else
      {
	pmutate=p;
	bvalue=randomi(2);
      }
  }
  void Init(bool c,double mx,double mn,double p,double r)    //initialises real gene
  {
    real=true;
    max=mx;
    change=c;
    if(c) 
      {
	min=mn;
	delta=max-min;
	pmutate=p;
	range=delta*r;
	value=min+random_01()*delta;
	if(value!=0) value=roundmio2(value,4);
      }
    else //fixed value
      {
	min=max;
	delta=0;
	pmutate=0;
	range=0;
	value=max;
      }
  }
  void Init()   //random start between max and min
  {
    if(change) 
      {
	if(real)
	  {
	    value=min+random_01()*delta;
	    if(value!=0) value=roundmio2(value,4);    //I only use 4 digit precision (for technical reasons)
	  }
	else bvalue=randomi(2);
      }
  }
  void Mutate()
  {
    if(change)
      {	
	double t=random_01();
	if(t<pmutate)
	  {	
	    if(real)
	      {
		value+=GaussN(range);    //mutates with sigma=range
		if(value>max) value=max;
		else if(value<min) value=min;
		if(value!=0) value=roundmio2(value,4);    //I only use 4 digit precision (for technical reasons)
	      }
	    else if(bvalue) bvalue=false;
	    else bvalue=true;
	  }
      }
  }
};

class Generic_genome          //a collection of genes
{
public:
  int size;                 //size genes
  Gene *genes;
  double fitness;
  Generic_genome()
  {
    size=0;
    fitness=0;
  }
  Generic_genome(int s)
  {
    size=s;
    genes=new Gene[size];    //initialises the genes array
    fitness=0;
  }
  void Init(int s)
  {
    size=s;
    genes=new Gene[size];
    fitness=0;
  }
  ~Generic_genome()
  {
    if(size) delete [] genes;
    size=0;
  }
  void Free()
  {
    if(size) delete [] genes;
    size=0;
  }
  void Copy(Generic_genome &P)
  {
    if(size!=P.size) {cout << "error in passing genome" << endl;return;}
    for(int i=0;i<size;i++) genes[i]=P.genes[i];
    fitness=P.fitness;
  }
  void Mutate()     //mutation for ga
  {
    for(int i=0;i<size;i++)
      {
	genes[i].Mutate();
      }
    fitness=0;
  }
};

void Cross(Generic_genome &C,Generic_genome &G1,Generic_genome &G2)   //cross operator for ga
{
  if((G1.size!=G2.size)||(G1.size!=C.size)) {cout << "error in crossing genome" << endl;return;}
  int cp=randomi(C.size);
  for(int i=0;i<cp;i++) C.genes[i]=G1.genes[i];
  for(int i=cp;i<C.size;i++) C.genes[i]=G2.genes[i];
}

void Select(Generic_genome &S,int genomes,Generic_genome *Genomes,int tour)   //selection operator for ga (uses tournament selection of size tour)
{
  int selected;
  double champion=-1e12;
  int tr;
  while(1)
    {
      tr=randomi(genomes);
      selected=tr;
      champion=Genomes[tr].fitness;
      if(Genomes[tr].fitness==Genomes[tr].fitness) break;
    }
  for(int i=1;i<tour;i++)
    {
      tr=randomi(genomes);
      if(Genomes[tr].fitness>champion)
	{
	  selected=tr;
	  champion=Genomes[tr].fitness;
	}
    }
  S.Copy(Genomes[selected]);
}



class Genome    //this is the pedestrian genome, i.e. the model parameters
{                //NB the original model is simple, but describes behaviour in "standard" situations. Many of the following parameters are used to handle more "unusual" behaviours while leaving the "usual" behaviour unchanged
public:
  double r_0;      //the distance between pedestrians (2014 PRE paper) GROUP BEHAVIOUR! 
  double Cr;       //the weight of the radial potential (2014 PRE paper) GROUP BEHAVIOUR! 
  double S0;       //the intensity of second neighbour force (2015 EPL paper) GROUP BEHAVIOUR! 
  double Ct;       //the weight of the angular potential (2014 PRE paper) GROUP BEHAVIOUR! 
  double eta;      //non-newtonianity  (2014 PRE paper) GROUP BEHAVIOUR! 
  double k;        //k_{v_p} in 2022 SAFETY paper
  double vp;       //preferred velocity determined by mu_v and sigma_v in 2022 SAFETY paper
  double CC;   //collision force intensity C in 2022 SAFETY paper
  double coll_rad;   //collision force maximum radius R_{max} in 2022 SAFETY paper
  double velover;     //NOT USED
  double rmax;      //when the distance goes over this value a stronger force is applied  GROUP BEHAVIOUR!  
  double body_radius;    
  double ext_disk;
  double ext_disk_2;   
  double body_diameter;//d_{int} in 2022 SAFETY paper
  double extern_radius;
  double extern_diameter;//d_{ext} in 2022 SAFETY paper
  double KVT;  //k^v_{\theta} in 2022 SAFETY paper
  double KT;//k_{\omega} in 2022 SAFETY paper
  double KTT;//k^{\omega}_{\theta} in 2022 SAFETY paper
  double KSTEP;//k_s in 2022 SAFETY paper
  double CBO;//gamma in 2022 SAFETY paper
  double TMAX;//tau_2 in 2022 SAFETY paper
  double BETA;//tau_1 in 2022 SAFETY paper
  double MU;
  double GS;//sigma_v in 2022 SAFETY paper
  double EBV;//B_l in 2022 SAFETY paper
  double GD2;//sigma_g in 2022 SAFETY paper
  bool SIG_DRIFT;
  double radius_view;    //no interaction if distance goes beyond this
  double fitness;        //overall fitness for GA, not used in simulator
  double fitness_pot;    //fitness concerning the ability to keep formation
  double fitness_coll; //fitness concerning the ability to avoid collisions
  double fitness_walk; //fitness concerning the ability to avoid collisions
  double fitness_dir;
  double fitness_near;
  int count_pot;       //counter used in the computation of fitness_pot
  int count_walk;
  int count_dir;
  int count_near;
  Genome()   //initialises to zero
  {
    r_0=0;
    Cr=0;
    S0=0; 
    Ct=0;
    eta=0;
    k=0;
    vp=0;
    CC=0;
    KT=0;
    KTT=0;
    KVT=0;
    KSTEP=0;
    CBO=0;
    TMAX=0;
    BETA=0;
    MU=0;
    GD2=0;
    GS=0;
    EBV=0;
    coll_rad=0;
    velover=0;
    rmax=0.;
    body_radius=0;
    ext_disk=0;
    ext_disk_2=0;
    body_diameter=0;
    radius_view=0;
    extern_radius=0;
    extern_diameter=0;
    SIG_DRIFT=false;
  }
  Genome(double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd) //initialises to passed values at creation
  {
    r_0=rp;
    Cr=Crp;
    S0=Sp; 
    Ct=Ctp;
    eta=ep;
    k=kp;
    vp=v;
    CC=cc;
    KT=kt;
    KTT=ktt;
    KVT=kvt;
    KSTEP=kstep;
    CBO=cbo;
    MU=mu;
    BETA=beta;
    GD2=gd2;
    GS=gs;
    EBV=ebv;    
    TMAX=MU+BETA;
    coll_rad=cr;
    velover=vo;
    rmax=rm;
    body_radius=br;
    body_diameter=2*br;
    ext_disk=ed;
    ext_disk_2=ext_disk*2;
    radius_view=rv;
    extern_radius=body_radius+ext_disk;
    extern_diameter=body_diameter+ext_disk*2;
    SIG_DRIFT=sd;    
  }
  bool Uguale(Genome U)
  {
    bool u=true;
    u=u&&(r_0==U.r_0);    
    u=u&&(Cr==U.Cr);      
    u=u&&(S0==U.S0);    
    u=u&&(Ct==U.Ct);    
    u=u&&(eta==U.eta);  
    u=u&&(k==U.k);     
    u=u&&(vp==U.vp);    
    u=u&&(CC==U.CC); 
    u=u&&(coll_rad==U.coll_rad);  
    u=u&&(velover==U.velover);    
    u=u&&(rmax==U.rmax);     
    u=u&&(body_radius==U.body_radius);  
    u=u&&(ext_disk==U.ext_disk);
    u=u&&(ext_disk_2==U.ext_disk_2);   
    u=u&&(body_diameter==U.body_diameter);
    u=u&&(extern_radius==U.extern_radius);
    u=u&&(extern_diameter==U.extern_diameter);
    u=u&&(KVT==U.KVT);
    u=u&&(KSTEP==U.KSTEP);    
    u=u&&(KT==U.KT);
    u=u&&(KTT==U.KTT);    
    u=u&&(CBO==U.CBO);
    u=u&&(TMAX==U.TMAX);
    u=u&&(BETA==U.BETA);
    u=u&&(MU==U.MU);
    u=u&&(GD2==U.GD2);
    u=u&&(GS==U.GS);
    u=u&&(EBV==U.EBV);   
    u=u&&(SIG_DRIFT==U.SIG_DRIFT);
    u=u&&(radius_view==U. radius_view);
    u=u&&(fitness==U.fitness);   
    u=u&&(fitness_pot==U.fitness_pot);   
    u=u&&(fitness_coll==U.fitness_coll);
    u=u&&(fitness_walk==U.fitness_walk);
    u=u&&(fitness_dir==U.fitness_dir);
    u=u&&(fitness_near==U.fitness_near);
    u=u&&(count_pot==U.count_pot); 
    u=u&&(count_walk==U.count_walk);
    u=u&&(count_dir==U.count_dir);
    u=u&&(count_near==U.count_near);
    return u;
  }
  void Init(double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd) //initialises to passed values at given time
  {
    r_0=rp;
    Cr=Crp;
    S0=Sp; 
    Ct=Ctp;
    eta=ep;
    k=kp;
    vp=v;
    CC=cc;
    KT=kt;
    KTT=ktt;
    KVT=kvt;
    KSTEP=kstep;
    CBO=cbo;
    MU=mu;    
    BETA=beta;
    GD2=gd2;
    GS=gs;
    EBV=ebv;    
    TMAX=MU+BETA;
    coll_rad=cr;
    velover=vo;
    rmax=rm;
    body_radius=br;
    body_diameter=2*br;
    ext_disk=ed;
    ext_disk_2=ext_disk*2;
    radius_view=rv;
    extern_radius=body_radius+ext_disk;
    extern_diameter=body_diameter+ext_disk*2;
    SIG_DRIFT=sd;    
  }
  void Init(Generic_genome &G)                  //initialises from genome
  {
    if(G.size!=25) cout << "error in converting genome" << endl;
    r_0=G.genes[0].value;
    Cr=G.genes[1].value;
    Ct=G.genes[2].value;
    eta=G.genes[3].value;
    k=G.genes[4].value;
    vp=G.genes[5].value;
    S0=G.genes[6].value;
    CC=G.genes[7].value;
    coll_rad=G.genes[8].value;
    velover=G.genes[9].value;
    rmax=G.genes[10].value;
    body_radius=G.genes[11].value;
    body_diameter=2*body_radius;
    ext_disk=G.genes[12].value;
    ext_disk_2=ext_disk*2;
    radius_view=G.genes[13].value;
    extern_radius=body_radius+ext_disk;
    extern_diameter=body_diameter+ext_disk*2;
    KT=G.genes[14].value;
    KTT=G.genes[15].value;
    CBO=G.genes[16].value;
    MU=G.genes[17].value;
    BETA=G.genes[18].value;
    TMAX=BETA+MU;        
    KVT=G.genes[19].value;
    KSTEP=G.genes[20].value;
    GD2=G.genes[21].value;
    GS=G.genes[22].value;
    EBV=G.genes[23].value;    
    SIG_DRIFT=G.genes[24].bvalue;    
  }
  void Copy_In(Generic_genome &G)                  //initialises from genome
  {
    if(G.size!=25) cout << "error in converting genome" << endl;
    G.genes[0].value=r_0;
    G.genes[1].value=Cr;
    G.genes[2].value=Ct;
    G.genes[3].value=eta;
    G.genes[4].value=k;
    G.genes[5].value=vp;
    G.genes[6].value=S0;
    G.genes[7].value=CC;
    G.genes[8].value=coll_rad;
    G.genes[9].value=velover;
    G.genes[10].value=rmax;
    G.genes[11].value=body_radius;
    G.genes[12].value=ext_disk;
    G.genes[13].value=radius_view;
    G.genes[14].value=KT;
    G.genes[15].value=KTT;
    G.genes[16].value=CBO;
    G.genes[17].value=MU;
    G.genes[18].value=BETA;        
    G.genes[19].value=KVT;
    G.genes[20].value=KSTEP;
    G.genes[21].value=GD2;
    G.genes[22].value=GS;
    G.genes[23].value=EBV;  
    G.genes[24].bvalue=SIG_DRIFT;    
  }
  void Fix_Vel(double v)    //if the preferred velocity is different than pedestrian other values are scaled too
  {
    //double ratio=vp/1.336;
    vp=v;
    //Cr=0.62*ratio;
    //Ct=0.08*ratio;
  }
  void Init_Fit()   //initialises fitness
  {
    fitness=0;
    fitness_pot=0;
    fitness_coll=0;
    fitness_walk=0;
    fitness_dir=0;
    fitness_near=0;
    count_pot=0;
    count_walk=0;
    count_dir=0;
    count_near=0;
  }
  void Compute_Fit()   //computes fitness
  {
    if(count_pot) fitness_pot/=count_pot;
    if(count_walk) fitness_walk/=count_walk;
    if(count_dir) fitness_dir/=count_dir;
    if(count_near) fitness_near/=count_near;
    //fitness=fitness_walk-(fitness_pot+fitness_coll+fitness_dir);
  }
};

class Pedestrian      //a pedestrian (robot) class
{
public:
  Genome genome;      //parameters
  State2D Self;            //position and velocity
  Vector2D nextv;      //next velocity
  Vector2D acc;        //current total acceleratio 
  Vector2D drift;     //acceleration due to local goal
  Vector2D drift2;     //acceleration due to body orientation
  Vector2D drift3;  //acceleration due to overlapping
  Vector2D centriped;//to simulate behaviour in curved corridors
  double omega_centr;
  Vector2D group;     //acceleration due to group coesive behaviour
  Vector2D group_coll; //acceleration due to collision avoidance inside group
  Vector2D collision; //acceleration due collision avoidance with obstacles 
  Vector2D ped_coll; //acceleration due collision avoidance with pedestrians outside group
  Vector2D par_step;  //used for the simulation of the robot motion (v omega from vectors)
  Vector2D ort_step;  //used for the simulation of the robot motion (v omega from vectors)
  double domega;//angular acceleration
  Vector2D f_ell;//acceleration term due to ellipse collision prediction
  double vlin;  //linear velocity
  double omega;  //angular velocity
  double radius;  //radius of curvature, used for the simulation of the robot motion
  int rotating;  //a counter saying how much the robot should rotate before assuming the correct position and start moving
  double arr_time;//used to compute fitness is some environment/setting
  Vector2D arr_place;//used to compute fitness is some environment/setting
  Ellisse E;//body
  Ellisse EV; //for overlapping  
  double v_actual;//used to compute fitness is some environment/setting
  double v_ratio;//used to compute fitness is some environment/setting
  double sig;//to compute the "beta" function in the paper
  double msig;//to compute the "beta" function in the paper
  bool active;
  Pedestrian()
  {
    domega=0;
    vlin=0;  //linear velocity
    omega=0;  //angular velocity
    radius=0;  //radius of curvature, used for the simulation of the robot motion
    rotating=0;  //a counter saying how much the robot should rotate before assuming the correct position and start moving
    arr_time=0;
    v_actual=0;
    sig=0;
    msig=0;
    active=false;
  }
  Pedestrian(State2D S)   //initialises passing position and velocity together
  {
    Self=S;
    rotating=0;  //rotating should be initialised to 0
    arr_time=0;
    arr_place=S.r;
    active=true;
  }
  Pedestrian(Vector2D r,Vector2D v)  //initialises passing position and velocity explicitely
  { 
    Self.r=r;
    Self.v=v;
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;   
  }  
  Pedestrian(Genome g)  //initialises passing genome
  {
    genome=g;
    rotating=0;
    arr_time=0;
    active=false;
  }
  Pedestrian(double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)   //initialises passing parameters
  {
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=v;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;   
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;    
    genome.SIG_DRIFT=sd;      
    rotating=0;
    arr_time=0;
    active=false;
  }
  Pedestrian(State2D S,Genome g)//initialises passing genome and position and velocity
  {
    Self=S;
    genome=g;
    rotating=0;
    arr_time=0;
    arr_place=S.r;
    active=true;
  }
  Pedestrian(State2D S,double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)   //initialises passing parameters and position and velocity
  {
    Self=S;
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=v;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;    
    genome.CBO=cbo;
    genome.MU=mu;
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.SIG_DRIFT=sd;
    genome.KVT=kvt;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;     
    genome.KSTEP=kstep;    
    rotating=0;
    arr_time=0;
    arr_place=S.r;
    active=true;
  }
  Pedestrian(Vector2D r,Vector2D v,Genome g)//initialises passing genome and position and velocity (as separate vectors)
  {
    Self.r=r;
    Self.v=v;
    genome=g;
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;
  }
  Pedestrian(Vector2D r,Vector2D v,double rp,double Crp,double Sp,double Ctp,double ep,double kp,double vps,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)   //initialises passing parameters and position and velocity (as separate vectors)
  {
    Self.r=r;
    Self.v=v;
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=vps;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;       
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;    
    genome.SIG_DRIFT=sd;    
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;
  }
  void MK()//not used
  {
    genome.KT*=2;
    genome.KTT*=2;
    genome.KVT*=2;
  }
  bool Uguale(Pedestrian U)
  {
    bool u=true;
    u=u&&genome.Uguale(U.genome);
    u=u&&Self.Uguale(U.Self);
    u=u&&nextv.Uguale(U.nextv);
    u=u&&acc.Uguale(U.acc);
    u=u&&drift.Uguale(U.drift);
    u=u&&drift2.Uguale(U.drift2);
    u=u&&drift3.Uguale(U.drift3);
    u=u&&centriped.Uguale(U.centriped);
    u=u&&(omega_centr==U.omega_centr);
    u=u&&group.Uguale(U.group);    
    u=u&&group_coll.Uguale(U.group_coll); 
    u=u&&collision.Uguale(U.collision); 
    u=u&&ped_coll.Uguale(U.ped_coll);
    u=u&&par_step.Uguale(U.par_step);
    u=u&&ort_step.Uguale(U.ort_step);
    u=u&&(domega==U.domega);   
    u=u&&f_ell.Uguale(U.f_ell);
    u=u&&(vlin==U.vlin);   
    u=u&&(omega==U.omega);  
    u=u&&(radius==U.radius); 
    u=u&&(rotating==U.rotating);
    u=u&&(arr_time==U.arr_time);
    u=u&&arr_place.Uguale(U.arr_place);
    u=u&&E.Uguale(U.E);
    u=u&&EV.Uguale(U.EV);
    u=u&&(v_actual==U.v_actual);
    u=u&&(sig==U.sig);  
    u=u&&(msig==U.msig);
    u=u&&(active==U.active);
    return u;
  }
  void Init(State2D S)  //Init functions are similar to above but can be called after creation of the object
  {
    Self=S;
    rotating=0;
    arr_time=0;
    arr_place=S.r;
    active=true;
  }
  void Init(Vector2D r,Vector2D v)
  {
    Self.r=r;
    Self.v=v;
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;
  }  
  void Init(Genome g)
  {
    genome=g;
    rotating=0;
    arr_time=0;
    active=false;
  }
  void Init(double x,double y,double vx,double vy,double A,double B,double BV,double theta,double omega,int id)
  {
    Self.r.Init(x,y);
    Self.v.Init(vx,vy,theta);
    E.Initialise(x,y,vx,vy,theta,omega,A,B,id);
    EV.Initialise(x,y,0,0,theta,0,A,BV,id);    
    arr_time=0;
    arr_place.Init(x,y);
    active=true;
  }
  void Init(Ellisse Ep)
  {
    Self.r.Init(Ep.x[0],Ep.x[1]);
    Self.v.Init(Ep.v[0],Ep.v[1]);
    E.Copy(Ep);
    EV.Update(Self.r.x,Self.r.y,Self.v.x,Self.v.y,Self.v.th,0);
    active=true;
  }  
  void Init(double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)
  {
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=v;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;  
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;    
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;    
    genome.SIG_DRIFT=sd;     
    rotating=0;
    arr_time=0;
    active=false;
  }
  void Init(State2D S,Genome g)
  {
    Self=S;
    genome=g;
    rotating=0;
    arr_time=0;
    arr_place=S.r;
    active=true;
  }
  void Init(State2D S,double rp,double Crp,double Sp,double Ctp,double ep,double kp,double v,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)
  {
    Self=S;
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=v;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;   
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;    
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;     
    genome.SIG_DRIFT=sd;     
    rotating=0;
    arr_time=0;
    arr_place=S.r;
    active=true;    
  }
  void Init(Vector2D r,Vector2D v,Genome g)
  {
    Self.r=r;
    Self.v=v;
    genome=g;
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;    
  }
  void Init(Vector2D r,Vector2D v,double rp,double Crp,double Sp,double Ctp,double ep,double kp,double vps,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)
  {
    Self.r=r;
    Self.v=v;
    genome.r_0=rp;
    genome.Cr=Crp;
    genome.S0=Sp; 
    genome.Ct=Ctp;
    genome.eta=ep;
    genome.k=kp;
    genome.vp=vps;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt;    
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;    
    genome.SIG_DRIFT=sd;      
    rotating=0;
    arr_time=0;
    arr_place=r;
    active=true;    
  }
  void Standard_parameters(double vp,bool s,double r0,double cc,double cr,double vm,double vo,double vs,double rm,double br,double ed,double rv,double kt,double ktt,double kvt,double kstep,double cbo,double mu,double beta,double gd2,double gs,double ebv,bool sd)   //part of these parameters are not passed because "standard"
  {
    genome.r_0=0.75;
    genome.Cr=0.62;
    genome.Ct=0.08;
    genome.eta=-0.43;
    genome.k=1.52;
    genome.vp=1.336;
    genome.S0=0.5;
    if(vp>0)
      {
	double ratio=vp/genome.vp;
	genome.Cr*=ratio;
	genome.Ct*=ratio;
	genome.vp=vp;	
      }
    else
      {
	genome.vp=0;
	genome.Cr=0;
	genome.Ct=0;
      }
    if(!s) genome.S0=0;
    genome.r_0=r0;
    genome.CC=cc;
    genome.coll_rad=cr;
    genome.velover=vo;
    genome.rmax=rm;
    genome.body_radius=br;
    genome.body_diameter=2*br;
    genome.ext_disk=ed;
    genome.ext_disk_2=genome.ext_disk*2;
    genome.radius_view=rv;
    genome.extern_radius=genome.body_radius+genome.ext_disk;
    genome.extern_diameter=genome.body_diameter+genome.ext_disk_2;
    genome.KT=kt;
    genome.KTT=ktt; 
    genome.KVT=kvt;
    genome.KSTEP=kstep;
    genome.CBO=cbo;
    genome.MU=mu;  
    genome.BETA=beta;
    genome.TMAX=mu+beta;
    genome.GD2=gd2;
    genome.GS=gs;
    genome.EBV=ebv;      
    genome.SIG_DRIFT=sd;   
    rotating=0;
    arr_time=0;
    active=false;
  }
  void Active()
  {
    active=true;
  }
  void De_Active()
  {
    active=false;
  }  
  void v2vp(int gs)//not one used for groups, not in this simulator
  {
    if(genome.k>0)
      {
	if(gs==1)
	  {
	    v_actual=genome.vp;
	  }
	else
	  {
	    double phi=genome.eta*M_PI/2;
	    double acc;
	    if(gs==3) acc=16/(3*genome.r_0)*genome.Ct*phi;
	    else if(gs==2) acc=4/genome.r_0*genome.Ct*phi; 
	    double delta_v=acc/genome.k;
	    v_actual=genome.vp+delta_v;
	  }
      }
  }     
  void Sigmoid(double t)//originally was a sigmoid function but it has been simplified to linear
  {
    sig=0;
    if(t<genome.BETA) sig=1;
    else if((t<genome.TMAX)&&(genome.MU>0)) sig=(t-genome.BETA)/genome.MU;
    msig=1-sig;
  }
  void Ext_Fitness_Coll(double ext)//fitness functions not used in this simulator
  {
    {
      genome.fitness_coll+=ext;
    }
  }
  void Ext_Fitness_Walk(double time)
  {
    double deltat=time-arr_time;
    if(deltat>0)
      {
	double disp=Dist(Self.r,arr_place);
	double vel=disp/deltat;
	double fitq=vel/(v_actual);
	if(v_actual>0) genome.fitness_walk+=fitq;
	//if(fitq!=fitq) cout << "ATTENZIONE! " << Self.r.x << " " << Self.r.y << " " << arr_place.x << " " << arr_place.y << " " << disp << " " << deltat << " " << v_actual << endl;
	genome.count_walk++;
	arr_time=time;
	arr_place=Self.r;
      }
  }
  void Ext_Fitness_Walk(Ellisse E,Vector2D g,double dt)
  {
    Vector2D disp(E.x[0]-Self.r.x,E.x[1]-Self.r.y);
    double vel=Projection(disp,g)/dt;
    if(v_actual)
      {
	double fitq=vel/(v_actual*v_ratio);
	genome.fitness_walk+=fitq;
	genome.count_walk++;
      }
  }
  double EXT_Fitness_Walk(Ellisse E,Vector2D g,double dt)
  {
    Vector2D disp(E.x[0]-Self.r.x,E.x[1]-Self.r.y);
    double vel=Projection(disp,g)/dt;
    double fitq=vel/(v_actual*v_ratio);
    genome.fitness_walk+=fitq;
    genome.count_walk++;
    return fitq;
  }
  bool Ext_Fitness_Out(double iw,double ew)
  {
    double fx=fabs(Self.r.x);
    double fy=fabs(Self.r.y);
    bool out=((fx>ew)||(fy>ew));
    out=out||((fx<iw)&&(fy<iw));
    if(out) return true;
    else return false;
  }
  Vector2D Drift(Vector2D preferred)   //computes drift to local goal (eg 2014 PRE) preferred is the direction vector given by path planner (it has been normalised to preferred velocity)
  {
    Vector2D d;
    d.x=genome.k*(preferred.x-Self.v.x);
    d.y=genome.k*(preferred.y-Self.v.y);
    d.Magnitude();
    return d;
  }
  Vector2D Drift2()//drift force due to body orientation terms
  {
    Vector2D d;
    Vector2D vv(Self.v.m*sin(E.theta),Self.v.m*cos(E.theta));
    double dx,dy;
    dx=vv.x-Self.v.x;
    dy=vv.y-Self.v.y;
    d.x=genome.KVT*dx;
    d.y=genome.KVT*dy;
    d.Magnitude();
    genome.fitness_dir+=sqrt(dx*dx+dy*dy);
    genome.count_dir++;
    return d;
  }
  Vector2D Drift3(int nped,State2D *Peds,double DT)//drift force due to overlapping
  {
    double intensity=0;
    for(int i=0;i<nped;i++)
      {
	Vector2D rel_pos=Rel(Self.r,Peds[i].r);
	if((rel_pos.m<2*EV.B)&&(Scalar(Self.v,rel_pos)<0))
	  {
	    Ellisse Echeck(Peds[i].r.x,Peds[i].r.y,0,0,Peds[i].v.th,0,EV.A,EV.B,i);
	    double l[2];
	    if(Overlap(l,EV,Echeck))
	      {
		intensity+=1-Contact(l,EV,Echeck);
	      }
	  }
      }
    Vector2D d;
    if(intensity>0)
      {
	d.x=-genome.KSTEP*intensity*Self.v.x;
	d.y=-genome.KSTEP*intensity*Self.v.y;
	d.Magnitude();
	double bstep=d.m*DT;
	if(bstep>Self.v.m)
	  {
	    double scale=Self.v.m/bstep;
	    d.Scale(scale);
	  }
      }
    return d;
  }
  Vector2D CPW(Vector2D r,double v,double t)    //interaction with walls (obstacles)
  {
    Vector2D f;
    if(r.m==0)
      {
	double theta=randommio(2*M_PI);
	r.Init(genome.body_radius/2*cos(theta),genome.body_radius/2*sin(theta));
      }
    if(r.m<genome.body_radius)   //then maxim force opposed to the obstacle position (normalised to 1, then scaled below)
      {
	f.x=r.x/r.m;
	f.y=r.y/r.m;
      }
    else if(r.m<genome.extern_radius)  //decreases linearly to be 0 at genome.extern_radius
      {
	double er=1-(r.m-genome.body_radius)/genome.ext_disk;
	f.x=r.x/r.m*er;
	f.y=r.y/r.m*er;
      }
    f.Magnitude();
    f.Scale(genome.CC*v/t);//scales force based also on velocity, using velmin if too small
    return f;
  }
  Vector2D CPM(Vector2D r,double v,double t)//collosion avoidance inside group, as above CPW but using diameters and not using minimum velocity
  {
    Vector2D f;
    if(r.m==0)
      {
	double theta=randommio(2*M_PI);
	r.Init(genome.body_radius/2*cos(theta),genome.body_radius/2*sin(theta));
      }
    if(r.m<genome.body_diameter)
      {
	f.x=r.x/r.m;
	f.y=r.y/r.m;
      }
    else if(r.m<genome.extern_diameter)
      {
	double er=1-(r.m-genome.body_diameter)/genome.ext_disk_2;
	f.x=r.x/r.m*er;
	f.y=r.y/r.m*er;
      }
    f.Magnitude();
    f.Scale(genome.CC*v/t);  //intensity is different for obstacles, group member,etc
    return f;
  }
  Vector2D CPP(Vector2D r,double v,double t)//collision avoidance term
  {
    Vector2D f;
    if(r.m==0)
      {
	double theta=randommio(2*M_PI);
	r.Init(genome.body_radius/2*cos(theta),genome.body_radius/2*sin(theta));
      }    
    if(r.m<genome.body_diameter)
      {
	f.x=r.x/r.m;
	f.y=r.y/r.m;
      }
    else if(r.m<genome.extern_diameter)
      {
	double er;
	if(genome.ext_disk_2>0) er=1-(r.m-genome.body_diameter)/genome.ext_disk_2;
	f.x=r.x/r.m*er;
	f.y=r.y/r.m*er;
      }
    f.Magnitude();
    f.Scale(genome.CC*v/t);
    return f;
  }
  Vector2D CollisionForce_av(int nobs,Vector2D *Obs,double DT)   //returns average force if vector of obstacles is passed. It also computes fitness!
  {
    Vector2D fi; //details of the computation are found in EPL 2011, computes the "first collision time" tmin, minor modifications are present with respect to original paper
    Vector2D f;  //as explained in ICSR, we compute an average force to deal with obstacles
    Vector2D *r;
    if(Self.v.m)
      {
	double tmin=1e10;
	r=new Vector2D [nobs];
	int count=0;
	bool findnear=false;
	for(int i=0;i<nobs;i++)
	  {
	    r[i]=Rel(Obs[i],Self.r);
	    if(r[i].m<NEAR) findnear=true;
	    if(r[i].m<genome.coll_rad)
	      {
		double distalvel=Projection(r[i],Self.v);
		if(r[i].m>0)
		  {
		    if((distalvel/r[i].m)>cos(M_PI/4))
		      {
			double ortdist=fabs(Clock_Ort(r[i],Self.v));
			if(ortdist<genome.extern_radius)
			  {
			    if(Self.v.m>0)
			      {
				double lt=distalvel/Self.v.m;
				if(lt<tmin) tmin=lt;
			      }
			  }
		      }
		  }
	      }
	  }
	if(findnear) genome.fitness_near+=1;
	if(tmin<1e10)
	  {
	    Vector2D newself=Self.r;
	    Vector2D step=Self.v;
	    step.Scale(tmin);
	    newself.Add(step);
	    for(int i=0;i<nobs;i++)
	      {
		if((r[i].m<genome.coll_rad)&&(Scalar(r[i],Self.v)>0))
		  {
		    Vector2D newrel=Rel(newself,Obs[i]);
		    double timinpass=tmin;
		    if(timinpass<DT) timinpass=DT;
		    fi=CPW(newrel,Self.v.m,timinpass);
		    if(fi.m) {f.Add(fi);count++;}
		  }
	      }
	    if(count) f.Scale(1./count);
	  }
	delete [] r;
      }
    return f;
  }
 
  Vector2D PedCollision(int nped,State2D *Peds,double DT)  //here collision avoidance outside group (EPL 2011)
  {
    Vector2D fi; 
    Vector2D f;
    double tmin=1e10;
    Vector2D *rel_vel;
    Vector2D *rel_pos;
    rel_vel= new Vector2D [nped];
    rel_pos= new Vector2D [nped];
    for(int i=0;i<nped;i++)
      {
	rel_vel[i]=Rel(Peds[i].v,Self.v);
	rel_pos[i]=Rel(Self.r,Peds[i].r);
	if(rel_pos[i].m<NEAR2) genome.fitness_near+=1;  //including collision
	if(rel_pos[i].m<genome.radius_view)
	  {
	    double distalvel=Projection(rel_pos[i],rel_vel[i]);
	    if((distalvel>0)&&(Scalar(rel_pos[i],Self.v)<0)&&(rel_vel[i].m>0))
	      {
		double ortdist=fabs(Clock_Ort(rel_pos[i],rel_vel[i]));
		if(ortdist<genome.extern_diameter)
		  {
		    if(rel_vel[i].m>0)
		      {
			double lt=distalvel/rel_vel[i].m;
			if(lt<tmin) tmin=lt;
		      }
		  }
	      }
	  }
      }
    if(tmin<1e10)
      {
	Vector2D newself=Self.r;
	Vector2D step=Self.v;
	step.Scale(tmin);
	newself.Add(step);
	for(int i=0;i<nped;i++)
	  {
	    if((rel_pos[i].m<genome.radius_view)&&(Scalar(rel_pos[i],Self.v)<0))
	      {
		Vector2D newoth=Peds[i].r;
		step=Peds[i].v;
		step.Scale(tmin);
		newoth.Add(step);		
		Vector2D newrel=Rel(newself,newoth);
		double timinpass=tmin;
		if(timinpass<DT) timinpass=DT;
		fi=CPP(newrel,rel_vel[i].m,timinpass);
		f.Add(fi);
	      }
	  }
      }
    delete [] rel_vel;
    delete [] rel_pos;
    return f;
  }
  Vector2D EllCollision(double &tcoll,double t,double dt,Ellisse *El,int nes,double ER,double IR)//term due to ellipse overlapping
  {
    domega=0;
    Vector2D f;
    double dv[2];
    dv[0]=0;
    dv[1]=0;
    int fc=E.Collision(tcoll,dv,domega,t,genome.TMAX,dt,El,nes,ER,IR);
    if(fc)
      {
	double tdiv=tcoll;
	if(tdiv<dt) tdiv=dt;
	dv[0]*=genome.CBO/tdiv;
	dv[1]*=genome.CBO/tdiv;
	domega*=genome.CBO/tdiv;
      }
    f.Init(dv[0],dv[1]);
    return f;
  }  
  Vector2D EllCollision(double &tcoll,double t,double dt,Ellisse *El,int nes, Wall *P,int nw)//term due to ellipse overlapping (walls)
  {
    domega=0;    
    Vector2D f;
    double dv[2];
    dv[0]=0;
    dv[1]=0;
    int fc=E.Collision(tcoll,dv,domega,t,genome.TMAX,dt,El,nes, P,nw);
    if(fc)
      {
	double tdiv=tcoll;
	if(tdiv<dt) tdiv=dt;
	dv[0]*=genome.CBO/tdiv;
	dv[1]*=genome.CBO/tdiv;
	domega*=genome.CBO/tdiv;
      }
    f.Init(dv[0],dv[1]);
    return f;
  }
  void NextV(int nobs,Vector2D *Obstacles,int nped,State2D *Peds,Vector2D preferred,double DT,double vmax,double accmax,double t,Ellisse *El,int nes, Wall *P,int nw)
  {//computes next velocity
    msig=1;//if smaller than 1 scales drifts
    sig=0;
    domega=0;
    omega_centr=0;    
    genome.count_near++;
    if((preferred.m!=genome.vp)&&(preferred.m>0)) preferred.Scale(genome.vp/preferred.m);      //scales goal to preferred speed
    acc.Init(0,0);
    if(nes)
      {
	double tcoll;
	if(genome.TMAX>0.)
	  {
	    f_ell=EllCollision(tcoll,t,DT,El,nes,P,nw);
	    if(tcoll<genome.TMAX)
	      {
		Sigmoid(tcoll);
		f_ell.Scale(sig);	    
		acc.Add(f_ell);
		domega*=sig;
	      }
	  }
      }
    drift=Drift(preferred);     //computes drift (local goal)
    drift2=Drift2();
    drift3=Drift3(nped,Peds,DT);    
    if(msig<1)
      {
	if(genome.SIG_DRIFT)
	  {
	    drift.Scale(msig);
	    drift2.Scale(msig);
	  }
	drift3.Scale(msig);
      }
    acc.Add(drift);   //adds to acceleration
    acc.Add(drift2);
    acc.Add(drift3);   
    if(nped)
      {
	ped_coll=PedCollision(nped,Peds,DT);     //computes outside-group collision avoidance
	if(msig<1)
	  {	    
	    ped_coll.Scale(msig);	    
	  }	
	acc.Add(ped_coll);   //adds to acceleration
      } 
    if(nobs)
      {
	collision=CollisionForce_av(nobs,Obstacles,DT);    //computes obstacle collision avoidance
	if(msig<1)
	  {    
	    collision.Scale(msig);
	  }		
	acc.Add(collision);  //adds to acceleration
      }
    if(acc.m>accmax) acc.Scale(accmax/acc.m);//scales acceleration if too much
    nextv=acc;
    nextv.Scale(DT);//computes variation in velocity
    nextv.Add(Self.v);//adds to previous velocity
    if(nextv.m>vmax) nextv.Scale(vmax/nextv.m);   //scales if too much
  }
  void Update_onlyV(double DT)//UPDATES VELOCITY AND POSITION USING VECTORS (SIMULATOR)
  {
    Self.v=nextv;
    E.Update(Self.r.x,Self.r.y,Self.v.x,Self.v.y);
    EV.Update(Self.r.x,Self.r.y,Self.v.x,Self.v.y,Self.v.th,0);   
  }  
  void IR(double DT,double omega_max,double a_omega_max)//computes angular velocity
  {
    double thn=E.theta;
    ModmPpP(thn);
    double delo=Self.v.th-thn; 
    ModmPpP(delo);
    double a_omega;
    if(genome.SIG_DRIFT) a_omega=((-E.omega)*genome.KT+delo*genome.KTT)*msig+domega;
    else
      {
	a_omega=(-E.omega)*genome.KT+delo*genome.KTT+domega;
      }
    double faomega=fabs(a_omega);
    if(faomega>a_omega_max) a_omega*=a_omega_max/faomega;
    E.omega+=a_omega*DT;
    double fomega=fabs(E.omega);
    if(fomega>omega_max) E.omega*=omega_max/fomega;
  }
};



#endif

	    



