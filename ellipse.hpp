#ifndef __ELLIPSE_H__    
#define __ELLIPSE_H__ 
#include <iostream>  
#include <fstream> 
#include<math.h>                                                    
#include<cstdlib>
 
 
using namespace std;
  

double TMP=2*M_PI;

class Wall//defines a wall
{
public:
  bool up;
  bool vert;
  double level;
  double min;
  double max;
  Wall(){}
  void Init(bool u,bool v,double l,double mn,double mx)
  {
    up=u;
    vert=v;
    level=l;
    min=mn;
    max=mx;
  }
};


double randommio()  //real random number between 0 and 1
{
  double a;
  a=(double)rand()/RAND_MAX;
  return a;
}

double Try_Gauss()//normally distributed random number
{
  float x1,x2,w,y1,y2;
  do 
  {
    x1=2.0*randommio()-1.0;
    x2=2.0*randommio()-1.0;
    w=x1*x1+x2*x2;
  } 
  while (w>=1.0);
  w=sqrt((-2.0*log(w))/w);
  y1=x1*w;
  y2=x2*w;
  return y1;
}

void start()    //ititialises random number
{
  int stime;
  long ltime;
  ltime=time(NULL);
  stime=(unsigned) ltime/2;
  srand(stime);
  int i=rand();
}

int random(int i) //random integer between 0 and i-1
{
  int a;
  double r;
  r=RAND_MAX/double(i);
  a=rand();
  if(a!=RAND_MAX) a=int(a/r);
  else a=random(i);
  return a;
} 


void Mod(double &theta,double l)//modulo function
{
  if(theta>=l) {int np=int(theta/l);theta-=np*l;}
  else if(theta<0) {int np=int(theta/l)-1;theta-=np*l;}
}

void ModmPpP(double &theta)//modulo function
{
  if(theta>=M_PI) theta-=TMP;
  else if(theta<-M_PI) theta+=TMP;
}


double Scalar(double x1[2],double x2[2])//scalar product
{
  return x1[0]*x2[0]+x1[1]*x2[1];
}

double Vprod(double x1[2],double x2[2])//vector product
{
  return x1[0]*x2[1]-x2[0]*x1[1];
}

double Norm2(double x[2])//square of norm
{
  return Scalar(x,x);
}

double Norm(double x[2])//norm
{
  return sqrt(Norm2(x));
}

void Normalize(double x[2])//normalises
{
  double nx=Norm(x);
  if(nx>0)
    {
      double n=1./nx;
      x[0]*=n;
      x[1]*=n;
    }
}

double Normalized(double n[2],double x[2])//normalises
{
  double nor=Norm(x);
  if(nor>0)
    {
      double nm=1./nor;
      n[0]=x[0]*nm;
      n[1]=x[1]*nm;
    }
  return nor;
}

void Eigen(double D[2][2],double v[2],double &theta,double M[2][2])//computes an eigenvalue/vector
{
  double A=M[0][0];
  double B=M[1][1];
  double C=M[1][0];
  if(C!=0)
    {
      double sq=A*A+B*B-2*A*B+4*C*C;
      if(sq<0) sq=0;
      sq=sqrt(sq);
      double AB=A+B;
      D[0][0]=(AB-sq)/2;
      D[1][1]=(AB+sq)/2;
    }
  else
    {
      D[0][0]=M[0][0];
      D[1][1]=M[1][1];
    }
  if(D[0][0]>D[1][1])
    {
      double temp=D[1][1];
      D[1][1]=D[0][0];
      D[0][0]=temp;
    } 
  if(C!=0) v[1]=(D[0][0]-A)/C;
  else v[1]=0;
  v[0]=1./sqrt(1+v[1]*v[1]);
  v[1]*=v[0];
  theta=asin(v[1]);
}


void Matprod(double P[2][2],double M1[2][2],double M2[2][2])//product of matrices
{
  for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
      {
	P[i][j]=0;
      for(int k=0;k<2;k++)
	{
	  P[i][j]+=M1[i][k]*M2[k][j];
	}
      }
}

void R(double x[2],double t,double x0[2])//rotates
{
  double xt=cos(t)*x[0]+sin(t)*x[1];
  x[1]=cos(t)*x[1]-sin(t)*x[0];
  x[0]=xt;
  x[0]+=x0[0];
  x[1]+=x0[1];
}

void TR(double x[2],double &t,double rt,double x0[2])//rotates
{
  x[0]+=x0[0];
  x[1]+=x0[1];
  t+=rt;
  Mod(t,TMP);
  double xt=cos(rt)*x[0]+sin(rt)*x[1];
  x[1]=cos(rt)*x[1]-sin(rt)*x[0];
  x[0]=xt;
}

void TRM(double x[2],double &t,double rt,double x0[2])//rotates
{
  x[0]-=x0[0];
  x[1]-=x0[1];
  t+=rt;
  Mod(t,TMP);
  double xt=cos(rt)*x[0]+sin(rt)*x[1];
  x[1]=cos(rt)*x[1]-sin(rt)*x[0];
  x[0]=xt;
}

void TR(double x[2],double &t,double rt)//rotates
{
  t+=rt;
  Mod(t,TMP);  
  double xt=cos(rt)*x[0]+sin(rt)*x[1];
  x[1]=cos(rt)*x[1]-sin(rt)*x[0];
  x[0]=xt;
}

void TR(double x[2],double rt)//rotates
{
  double xt=cos(rt)*x[0]+sin(rt)*x[1];
  x[1]=cos(rt)*x[1]-sin(rt)*x[0];
  x[0]=xt;
}


class Ellisse//class for ellipse dynamics
{
public:
  double x[2];//position 
  double v[2];//velocity
  double theta;//orientation
  double omega;//angular velocity
  double A,B,c;//axes
  double A2,B2;//squared
  double iA2,iB2;//inverse
  double M,I;//area
  int id;
  Ellisse()
  {
     for(int i=0;i<2;i++)
      {
	x[i]=0;
	v[i]=0;
      }
    theta=0;
    omega=0;
    A=0;
    B=0;
    A2=0;
    B2=0;
    iA2=0;
    iB2=0;
    c=0;
    M=0;
    I=0;
    id=0;   
  }
  Ellisse(double xi[2],double vi[2],double t,double o,double a,double b,int id_)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=xi[i];
	v[i]=vi[i];
      }
    theta=t;
    Mod(theta,TMP);
    omega=o;
    A=a;
    B=b;
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);
    id=id_;
  }
  Ellisse(double xi,double yi,double vxi,double vyi,double t,double o,double a,double b,int id_)
  {
    x[0]=xi;
    x[1]=yi;    
    v[0]=vxi;
    v[1]=vxi;    
    theta=t;
    Mod(theta,TMP);
    omega=o;
    A=a;
    B=b;
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);
    id=id_;
  }
  bool Uguale(Ellisse E)//=
  {
    bool u=true;
    for(int i=0;i<2;i++)
      {
	u=u&&(x[i]==E.x[i]);	
	u=u&&(v[i]==E.v[i]);
      }
    u=u&&(theta==E.theta);
    u=u&&(omega==E.omega);
    u=u&&(A==E.A);
    u=u&&(B==E.B);
    u=u&&(A2==E.A2);
    u=u&&(B2==E.B2);
    u=u&&(iA2==E.iA2);
    u=u&&(iB2==E.iB2);
    u=u&&(c==E.c);
    u=u&&(M==E.M);
    u=u&&(I==E.I);
    u=u&&(id==E.id);
    return u;
  }	
  void Initialise(double xi[2],double vi[2],double t,double o,double a,double b,int id_)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=xi[i];
	v[i]=vi[i];
      }
    theta=t;
    Mod(theta,TMP);
    omega=o;
    A=a;
    B=b;
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);
    id=id_;
  }
  void Initialise(double xi[2],double t,double a,double b,int id_)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=xi[i];
      }
    theta=t;
    Mod(t,TMP);
    A=a;
    B=b;
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);
    id=id_;    
  }
  void Initialise(double xi,double yi,double vxi,double vyi,double t,double o,double a,double b,int id_)
  {
    x[0]=xi;
    x[1]=yi;    
    v[0]=vxi;
    v[1]=vxi;    
    theta=t;
    Mod(theta,TMP);
    omega=o;
    A=a;
    B=b;
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);
    id=id_;
  }
    void Update(double xi,double yi,double vxi,double vyi)
  {
    x[0]=xi;
    x[1]=yi;    
    v[0]=vxi;
    v[1]=vyi;    
  }  
  void Update(double xi,double yi,double vxi,double vyi,double t,double o)
  {
    x[0]=xi;
    x[1]=yi;    
    v[0]=vxi;
    v[1]=vyi;    
    theta=t;
    Mod(theta,TMP);
    omega=o;
  }  
  double Vprod(double x1[2],double x2[2])
  {
    return x1[0]*x2[1]-x2[0]*x1[1];
  }
  void AB(double a,double b)//passes axes
  {
    A=a;
    B=b;    
    A2=a*a;
    B2=b*b;
    iA2=1./A2;
    iB2=1./B2;
    c=A/B;
    M=M_PI*A*B;
    I=M/4*(A2+B2);     
  }
  double NormE2(double v[2])
  {
    double dx=v[0];
    double dy=v[1];
    return dx*dx*iA2+dy*dy*iB2;
  }
  double NormE(double v[2])
  {
    return sqrt(NormE2(v));
  }
  void NormalizeE(double v[2])
  {
    double ne=NormE(v);
    if(ne>0)
      {
	double n=1./ne;
	v[0]*=n;
	v[1]*=n;
      }
  }
  void Versor(double e[2],double phi)
  {
    e[0]=sin(phi)*iA2;
    e[1]=cos(phi)*iB2; 
    Normalize(e);
  }
  void Versor(double e[2],double l[2])
  {
    double dx=l[0]-x[0];
    double dy=l[1]-x[1];
    double phi=atan2(dx,dy);
    phi-=theta;
    e[0]=sin(phi)*iA2;
    e[1]=cos(phi)*iB2; 
    Normalize(e);
    TR(e,theta);
  }
  void Puntellisse(double l[2],double phi)
  {
    l[0]=sin(phi);
    l[1]=cos(phi);
    NormalizeE(l);
    l[0]+=x[0];
    l[1]+=x[1];
  }
  double Vprod(double l[2],double &d,double phi)
  {
    Puntellisse(l,phi);
    double er[2];
    d=Normalized(er,l);
    double v[2];
    Versor(v,phi);
    return ::Vprod(v,er);
  }
  double Distellisse_max(double l[2])
  {
    double pu,pm,pd,d;
    if(x[1]>=0)
      {
	if(x[0]>=0) pd=0;
	else pd=1.5*M_PI;
      }
    else
      {
	if(x[0]>=0) pd=0.5*M_PI;
	else pd=M_PI;
      }
    pu=pd+0.5*M_PI;
    double diff=0.5*M_PI;
    double vu=Vprod(l,d,pu);
    while(diff>1e-12)
      {
	pm=(pd+pu)/2.;
	double vm=Vprod(l,d,pm);
	if(vu*vm<=0) pd=pm;
	else
	  {
	    pu=pm;
	    vu=vm;
	  }
	diff=pu-pd;
      }
    return d;
  }  
  double Distellisse(double l[2])
  {
    double pu,pm,pd,d;
    if(x[1]>=0)
      {
	if(x[0]>=0) pd=M_PI;
	else pd=0.5*M_PI;
      }
    else
      {
	if(x[0]>=0) pd=1.5*M_PI;
	else pd=0;
      }
    pu=pd+0.5*M_PI;
    double diff=0.5*M_PI;
    double vu=Vprod(l,d,pu);
    while(diff>1e-12)
      {
	pm=(pd+pu)/2.;
	double vm=Vprod(l,d,pm);
	if(vu*vm<=0) pd=pm;
	else
	  {
	    pu=pm;
	    vu=vm;
	  }
	diff=pu-pd;
      }
    return d;
  }
  void BuildM(double M[2][2],double gamma)
  {
    double c2=cos(theta);
    double s2=sin(theta);
    double cs=c2*s2;
    double gamma2=gamma*gamma;
    c2*=c2;
    s2*=s2;
    M[0][0]=gamma2*(iA2*c2+iB2*s2);
    M[1][1]=(iA2*s2+iB2*c2);
    M[0][1]=M[1][0]=gamma*cs*(iA2-iB2);
  }
  void Copy_id(Ellisse C)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=C.x[i];
	v[i]=C.v[i];
      }
    theta=C.theta;
    omega=C.omega;
    A=C.A;
    B=C.B;
    A2=C.A2;
    B2=C.B2;
    iA2=C.iA2;
    iB2=C.iB2;
    c=C.c;
    M=C.M;
    I=C.I;
    id=C.id;
  }
  void Copy(Ellisse C)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=C.x[i];
	v[i]=C.v[i];
      }
    theta=C.theta;
    omega=C.omega;
    A=C.A;
    B=C.B;
    A2=C.A2;
    B2=C.B2;
    iA2=C.iA2;
    iB2=C.iB2;
    c=C.c;
    M=C.M;
    I=C.I;
  }
  void Scale(Ellisse C,double s)
  {
    for(int i=0;i<2;i++)
      {
	x[i]=C.x[i];
	v[i]=C.v[i];
      }
    theta=C.theta;
    omega=C.omega;
    A=C.A*s;
    B=C.B*s;
    double s2=s*s;
    A2=C.A2*s2;
    B2=C.B2*s2;
    iA2=C.iA2/s2;
    iB2=C.iB2/s2;
    c=C.c;
    M=C.M*s2;
    I=C.I*s2;
    }
  void CopyAB(Ellisse C)
  {
    A=C.A;
    B=C.B;
    A2=C.A2;
    B2=C.B2;
    iA2=C.iA2;
    iB2=C.iB2;
    c=C.c;
    M=C.M;
    I=C.I;    
  }
  void Updatev(double a[2],double n,double delta)
  {
    for(int i=0;i<2;i++) v[i]+=a[i]*delta;
    omega+=n*delta;
  }
  void Updatev(double a[2],double n)
  {
    for(int i=0;i<2;i++) v[i]+=a[i];
    omega+=n;
  }
  void Update(double dt)
  {
    for(int i=0;i<2;i++) x[i]+=v[i]*dt;
    theta+=omega*dt;
    Mod(theta,TMP);
  }
  void Update(double ex[2],double ev[2],double etheta,double eomega,double dt)
  {
    for(int i=0;i<2;i++) x[i]=ex[i]+ev[i]*dt;
    theta=etheta+eomega*dt;
    Mod(theta,TMP); 
  }
    void Update(Ellisse Ex,double dt)
  {
    for(int i=0;i<2;i++) x[i]=Ex.x[i]+Ex.v[i]*dt;
    theta=Ex.theta+Ex.omega*dt;
    Mod(theta,TMP);     
  }
  bool Overlap_parete(Wall P,double &root,double &rootp,double &rootm) //overlapping with wall 
  {
    double level,max,min;
    if(P.vert)
      {
	level=P.level-x[0];
	max=P.max-x[1];
	min=P.min-x[1];
      }
    else
      {
	level=P.level-x[1];
	max=P.max-x[0];
	min=P.min-x[0];
      }
    bool over=false;
    double KA,KB,KC;
    double s=sin(theta);
    double c=cos(theta);
    double c2=c*c;
    double s2=s*s;
    double cs=c*s;
    if(P.vert)
      {
	KA=s2*iA2+c2*iB2;
	KB=level*cs*(iB2-iA2);
	KC=level*level*(c2*iA2+s2*iB2)-1;
      }
    else
      {
	KA=c2*iA2+s2*iB2;
	KB=level*cs*(iB2-iA2);
	KC=level*level*(s2*iA2+c2*iB2)-1;
      }
    double D=KB*KB-KA*KC;
    if(D>=0) 
      {
	if(KA!=0)
	  {
	    rootp=(-KB+sqrt(D))/KA;
	    rootm=(-KB-sqrt(D))/KA;
	    if(((rootp<=max)&&(rootp>=min))||((rootm<=max)&&(rootm>=min)))
	      {
		over=true;
		if(P.vert)
		  {
		    rootp+=x[1];
		    rootm+=x[1];
		  }
		else
		  {
		    rootp+=x[0];
		    rootm+=x[0];
		  }
		root=(rootp+rootm)/2;
	      }
	  }
      }
    return over;
  }  
  bool Collision(double &dt,double l[2],Wall P)//collision with wall
  {
    double root,rootp,rootm;
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    E1fin.Update(x,v,theta,omega,dt);
    bool over=E1fin.Overlap_parete(P,root,rootp,rootm); 
    if(over)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double in=0;
	double fin=dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);
	    if(E1m.Overlap_parete(P,root,rootp,rootm))
	      {
		fin=middle;
		if(P.vert)
		  {
		    l[0]=P.level;
		    l[1]=root;
		  }
		else
		  {
		    l[1]=P.level;
		    l[0]=root;
		  }	      
	      }
	    else in=middle;
	    delta=fin-in;
	  }
	dt=middle;
	if(fabs(rootp-rootm)>1e-6)
	  {
	    if((fabs(rootp-P.min)<1e-8)||(fabs(rootp-P.max)<1e-8))
	      {
		if(P.vert) {l[0]=P.level;l[1]=rootp;}
		else {l[1]=P.level;l[0]=rootp;}
	      }
	    else if((fabs(rootm-P.min)<1e-8)||(fabs(rootm-P.max)<1e-8))
	      {
		if(P.vert) {l[0]=P.level;l[1]=rootm;}
		else {l[1]=P.level;l[0]=rootm;}
	      }
	  }
      }
    return over;
  }
  bool Collision_INT(double &dt,double l[2],double IR)
  {
    double root,rootp,rootm;
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    E1fin.Update(x,v,theta,omega,dt);
    bool over=Overlap(l,IR,E1fin);
    if(over)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double in=0;
	double fin=dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);
	    if(Overlap(l,IR,E1m))
	      {
		fin=middle;	      
	      }
	    else in=middle;
	    delta=fin-in;
	  }
	dt=middle;
	}
    return over;
  }
  bool Collision_EXT(double &dt,double l[2],double ER)
  {
    double root,rootp,rootm;
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    E1fin.Update(x,v,theta,omega,dt);
    bool inter=Internal(l,ER,E1fin);
    if(!inter)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double in=0;
	double fin=dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);
	    if(!Internal(l,ER,E1m))
	      {
		fin=middle;	      
	      }
	    else in=middle;
	    delta=fin-in;
	  }
	dt=middle;
      }
    return inter;
  }  
  bool Collision(double &t,double tmax,double dt,double l[2],Wall P)
  {
    double root,rootp,rootm;    
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    int steps=1;
    if(dt!=0) steps=int(tmax/dt)+1;
    bool over=false;
    int step;
    for(step=1;step<=steps;step++)
      {
	E1fin.Update(x,v,theta,omega,dt*step);
	over=E1fin.Overlap_parete(P,root,rootp,rootm);
	if(over) break;
      }
    if(over)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double fin=dt*step;  
	double in=fin-dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);
	    if(E1m.Overlap_parete(P,root,rootp,rootm))
	      {
		fin=middle;
		if(P.vert)
		  {
		    l[0]=P.level;
		    l[1]=root;
		  }
		else
		  {
		    l[1]=P.level;
		    l[0]=root;
		  }
	      }
	    else in=middle;
	    delta=fin-in;
	  }
	t=middle;
	if(fabs(rootp-rootm)>1e-6)
	  {
	    if((fabs(rootp-P.min)<1e-8)||(fabs(rootp-P.max)<1e-8))
	      {
		if(P.vert) {l[0]=P.level;l[1]=rootp;}
		else {l[1]=P.level;l[0]=rootp;}
	      }
	    else if((fabs(rootm-P.min)<1e-8)||(fabs(rootm-P.max)<1e-8))
	      {
		if(P.vert) {l[0]=P.level;l[1]=rootm;}
		else {l[1]=P.level;l[0]=rootm;}
	      }		  
	  }
      }
    return over;
  }
  bool Internal(double l[2],Ellisse E1,Ellisse E2) //assumes E1 much larger
  {
    bool inter;
    Ellisse E1r;
    Ellisse E2r;
    E1r.Copy(E1);
    E2r.Copy(E2);
    TRM(E1r.x,E1r.theta,-E1.theta,E1.x);
    TRM(E2r.x,E2r.theta,-E1.theta,E1.x); 
    double gamma=E1.A/E1.B;
    double M[2][2];
    E2r.BuildM(M,gamma);
    double D[2][2];
    double v[2];
    Eigen(D,v,E2r.theta,M);
    E1r.A=E1.B;
    E2r.AB(1./sqrt(D[0][0]),1./sqrt(D[1][1]));
    E2r.x[0]/=gamma;
    double mt2=E2r.theta;
    TR(E2r.x,E2r.theta,-mt2);
    double ln=Norm(E2r.x);
    double max=E2r.B;
    if(E2r.A>E2r.B) max=E2r.A;
    if(ln<E1r.A-max) {cout << "eccomi" << endl;inter=true;return inter;}    
    else
      {
	double distel=E2r.Distellisse_max(l);
	inter=(distel<E1r.A);
      }
    TR(l,mt2);
    l[0]*=gamma;
    R(l,E1.theta,E1.x);
    return inter;
  }
  bool Internal(double l[2],double R,Ellisse E2) //assumes E1 much larger
  {
    bool inter;
    Ellisse E2r;
    E2r.Copy(E2);
    double mt2=E2r.theta;
    TR(E2r.x,E2r.theta,-mt2);
    double ln=Norm(E2r.x);
    double max=E2r.B;
    if(E2r.A>E2r.B) max=E2r.A;
    if(ln<R-max) {inter=true;return inter;}  
    else
      {
	double distel=E2r.Distellisse_max(l);
	inter=(distel<R);
      }
    TR(l,mt2);
    return inter;
  }   
  bool Overlap(double l[2],Ellisse E1,Ellisse E2)//overlap between ellipses
  {
    bool over;
    Ellisse E1r;
    Ellisse E2r;
    E1r.Copy(E1);
    E2r.Copy(E2);
    TRM(E1r.x,E1r.theta,-E1.theta,E1.x);
    TRM(E2r.x,E2r.theta,-E1.theta,E1.x); 
    double gamma=E1.A/E1.B;
    double M[2][2];
    E2r.BuildM(M,gamma);
    double D[2][2];
    double v[2];
    Eigen(D,v,E2r.theta,M);
    E1r.A=E1.B;
    E2r.AB(1./sqrt(D[0][0]),1./sqrt(D[1][1]));
    E2r.x[0]/=gamma;
    double mt2=E2r.theta;
    TR(E2r.x,E2r.theta,-mt2);
    double ln=Norm(E2r.x);
    double max=E2r.B;
    if(E2r.A>E2r.B) max=E2r.A;
    if(ln>max+E1r.A) {over=false;return over;}  
    else
      {
	double ln=E2r.NormE(E2r.x);
	if(ln<1)
	  {
	    l[0]=E2r.x[0];
	    l[1]=E2r.x[1];
	    over=true;
	  }
	else
	  {
	    double distel=E2r.Distellisse(l);
	    over=(distel<E1r.A);
	  }
      }
    TR(l,mt2);
    l[0]*=gamma;
    R(l,E1.theta,E1.x);
    return over;
  }
  bool Overlap(double l[2],double R,Ellisse E2)//overlap with other ellipse
  {
    bool over;
    Ellisse E2r;
    E2r.Copy(E2);
    double mt2=E2r.theta;
    TR(E2r.x,E2r.theta,-mt2);
    double ln=Norm(E2r.x);
    double max=E2r.B;
    if(E2r.A>E2r.B) max=E2r.A;
    if(ln>max+R) {over=false;return over;}  
    else
      {
	double ln=E2r.NormE(E2r.x);
	if(ln<1)
	  {
	    l[0]=E2r.x[0];
	    l[1]=E2r.x[1];
	    over=true;
	  }
	else
	  {
	    double distel=E2r.Distellisse(l);
	    over=(distel<R);
	  }
      }
    TR(l,mt2);
    return over;
  }  
  bool Collision(double &t,double tmax,double dt,double l[2],Ellisse E2)//collision with ellipse E2
  {
    Ellisse E1fin;
    Ellisse E2fin;
    Ellisse E1m;
    Ellisse E2m;
    E1fin.AB(A,B);
    E2fin.CopyAB(E2);
    int steps=1;
    if(dt!=0) steps=int(tmax/dt)+1;
    bool over=false;
    int step;
    for(step=1;step<=steps;step++)
      {
	E1fin.Update(x,v,theta,omega,dt*step);
	E2fin.Update(E2,dt*step);
	over=Overlap(l,E1fin,E2fin);
	if(over) break;
      }
    if(over)
      {
	E1m.AB(A,B);
	E2m.CopyAB(E2);
	double delta=dt;
	double middle;
	double fin=dt*step;  
	double in=fin-dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);
	    E2m.Update(E2,middle);      
	    if(Overlap(l,E1m,E2m)) fin=middle;
	    else in=middle;
	    delta=fin-in;
	    //cout << "middle " << middle << endl;
	  }
	t=middle;
      }
    return over;
  }
  bool Collision_INT(double &t,double tmax,double dt,double l[2],double R)
  {
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    int steps=1;
    if(dt!=0) steps=int(tmax/dt)+1;
    bool over=false;
    int step;
    for(step=1;step<=steps;step++)
      {
	E1fin.Update(x,v,theta,omega,dt*step);
	over=Overlap(l,R,E1fin);
	if(over) break;
      }
    if(over)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double fin=dt*step;  
	double in=fin-dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);    
	    if(Overlap(l,R,E1m)) fin=middle;
	    else in=middle;
	    delta=fin-in;
	  }
	t=middle;
      }
    return over;
  }
  bool Collision_EXT(double &t,double tmax,double dt,double l[2],double R)
  {
    Ellisse E1fin;
    Ellisse E1m;
    E1fin.AB(A,B);
    int steps=1;
    if(dt!=0) steps=int(tmax/dt)+1;
    bool inter=true;
    int step;
    for(step=1;step<=steps;step++)
      {
	E1fin.Update(x,v,theta,omega,dt*step);
	inter=Internal(l,R,E1fin);
	if(!inter) break;
      }
    if(!inter)
      {
	E1m.AB(A,B);
	double delta=dt;
	double middle;
	double fin=dt*step;  
	double in=fin-dt;
	while(delta>1e-12)
	  {
	    middle=in+delta/2;
	    E1m.Update(x,v,theta,omega,middle);    
	    if(!Internal(l,R,E1m)) fin=middle;
	    else in=middle;
	    delta=fin-in;
	  }
	t=middle;
      }
    return (!inter);
  }    
  void Collision(double &tcoll,int &fc,double dv[2],double &domega,double t,double tmax,double dt,Ellisse *E,int nes, Wall *P,int nw)//collisions with environment
  {
    fc=0;
    double l[2];
    double ll[2];
    tcoll=tmax;
    int fcj;
    double tc;
    for(int j=0;j<nw;j++)
      {    
	if(Collision(tc,tcoll,dt,ll,P[j]))
	  {
	    if((tc<tcoll)&&(tc>0))
	      {
		tcoll=tc;
		fc=1;
		fcj=j;
		for(int i=0;i<2;i++) l[i]=ll[i];
	      }
	  }
      }
    for(int j=0;j<nes;j++)
      {
	if(E[j].id!=id)
	  {
	    if(Collision(tc,tcoll,dt,ll,E[j]))
	      {
		if((tc<tcoll)&&(tc>0))
		  {
		    tcoll=tc;
		    fcj=j;
		    fc=2;
		    for(int i=0;i<2;i++) l[i]=ll[i];		  
		  }
	      }
	  }
      }
    cout << "I am " << id << endl;
    if(fc==2) cout << "Now " << t << " future collision with " << fcj << " at t=" << t+tcoll << endl;
    else if(fc==1) cout << "Now " << t << " future collision with wall at t=" << t+tcoll << endl;
    else cout << "No future collision between " << t << " and " << t+tcoll << endl;
    if(fc)
      {
	Ellisse E1dummy;
	E1dummy.AB(A,B);
	E1dummy.Update(x,v,theta,omega,tcoll);
	double versor[2];
	double dp;
	double ps1;
	double pv1;
	double dr1[2];
	if(fc==2)
	  {
	    double dp;
	    double ps2;
	    double pv2;
	    double dr2[2];	    
	    Ellisse E2dummy;
	    E2dummy.CopyAB(E[fcj]);
	    E2dummy.Update(E[fcj],tcoll);	    
	    E2dummy.Versor(versor,l);	  
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
		dr2[i]=l[i]-E2dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    ps2=Scalar(E[fcj].v,versor);
	    pv1=Vprod(dr1,versor);
	    pv2=Vprod(dr2,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+1./E2dummy.M+pv1*pv1/E1dummy.I+pv2*pv2/E2dummy.I;
	    C2=-2*(ps1-ps2-omega*pv1+E[fcj].omega*pv2);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;    
	      } 
	    domega=-(dp*pv1)/E1dummy.I;	    
	  }
	else
	  {
	    E1dummy.Versor(versor,l);
	    versor[0]=-versor[0];
	    versor[1]=-versor[1];
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    pv1=Vprod(dr1,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+pv1*pv1/E1dummy.I;
	    C2=-2*(ps1-omega*pv1);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;	      
	      } 
	    domega=-(dp*pv1)/E1dummy.I;
	  }
	cout << "Velocity variation at collision dvx=" << dv[0] << " dvy=" << dv[1] << endl;
	cout << "Angular velocity variation at collision domega=" << domega << endl;	
      }
  }
  int Collision(double &tcoll,double dv[2],double &domega,double t,double tmax,double dt,Ellisse *E,int nes,double ER,double IR)//collisions with environment
  {
    int fc=0;
    double l[2];
    double ll[2];
    tcoll=tmax;
    int fcj;
    double tc;
    if(Collision_EXT(tc,tcoll,dt,ll,ER))
      {
	if((tc<tcoll)&&(tc>0))
	  {
	    tcoll=tc;
	    fc=1;
	    fcj=0;
	    for(int i=0;i<2;i++) l[i]=ll[i];
	  }
      }
    if(Collision_INT(tc,tcoll,dt,ll,IR))
      {
	if((tc<tcoll)&&(tc>0))
	  {
	    tcoll=tc;
	    fc=1;
	    fcj=1;
	    for(int i=0;i<2;i++) l[i]=ll[i];
	  }
      }
    for(int j=0;j<nes;j++)
      {
	if(E[j].id!=id)
	  {
	    if(Collision(tc,tcoll,dt,ll,E[j]))
	      {
		if((tc<tcoll)&&(tc>0))
		  {
		    tcoll=tc;
		    fcj=j;
		    fc=2;
		    for(int i=0;i<2;i++) l[i]=ll[i];		  
		  }
	      }
	  }
      }
    if(fc)
      {
	Ellisse E1dummy;
	E1dummy.AB(A,B);
	E1dummy.Update(x,v,theta,omega,tcoll);
	double versor[2];
	double dp;
	double ps1;
	double pv1;
	double dr1[2];
	if(fc==2)
	  {
	    double dp;
	    double ps2;
	    double pv2;
	    double dr2[2];	    
	    Ellisse E2dummy;
	    E2dummy.CopyAB(E[fcj]);
	    E2dummy.Update(E[fcj],tcoll);	    
	    E2dummy.Versor(versor,l);	  
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
		dr2[i]=l[i]-E2dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    ps2=Scalar(E[fcj].v,versor);
	    pv1=Vprod(dr1,versor);
	    pv2=Vprod(dr2,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+1./E2dummy.M+pv1*pv1/E1dummy.I+pv2*pv2/E2dummy.I;
	    C2=-2*(ps1-ps2-omega*pv1+E[fcj].omega*pv2);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;    
	      } 
	    domega=-(dp*pv1)/E1dummy.I;	    
	  }
	else
	  {
	    E1dummy.Versor(versor,l);
	    versor[0]=-versor[0];
	    versor[1]=-versor[1];
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    pv1=Vprod(dr1,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+pv1*pv1/E1dummy.I;
	    C2=-2*(ps1-omega*pv1);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;	      
	      } 
	    domega=-(dp*pv1)/E1dummy.I;
	  }
      }
    return fc;
  }
  int Collision(double &tcoll,double dv[2],double &domega,double t,double tmax,double dt,Ellisse *E,int nes, Wall *P,int nw)//collisions with environment
  {
    int fc=0;
    double l[2];
    double ll[2];
    tcoll=tmax;
    int fcj;
    double tc;
    for(int j=0;j<nw;j++)
      {    
	if(Collision(tc,tcoll,dt,ll,P[j]))
	  {
	    if((tc<tcoll)&&(tc>0))
	      {
		tcoll=tc;
		fc=1;
		fcj=j;
		for(int i=0;i<2;i++) l[i]=ll[i];
	      }
	  }
      }
    for(int j=0;j<nes;j++)
      {
	if(E[j].id!=id)
	  {
	    if(Collision(tc,tcoll,dt,ll,E[j]))
	      {
		if((tc<tcoll)&&(tc>0))
		  {
		    tcoll=tc;
		    fcj=j;
		    fc=2;
		    for(int i=0;i<2;i++) l[i]=ll[i];		  
		  }
	      }
	  }
      }
    if(fc)
      {
	Ellisse E1dummy;
	E1dummy.AB(A,B);
	E1dummy.Update(x,v,theta,omega,tcoll);
	double versor[2];
	double dp;
	double ps1;
	double pv1;
	double dr1[2];
	if(fc==2)
	  {
	    double dp;
	    double ps2;
	    double pv2;
	    double dr2[2];	    
	    Ellisse E2dummy;
	    E2dummy.CopyAB(E[fcj]);
	    E2dummy.Update(E[fcj],tcoll);	    
	    E2dummy.Versor(versor,l);	  
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
		dr2[i]=l[i]-E2dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    ps2=Scalar(E[fcj].v,versor);
	    pv1=Vprod(dr1,versor);
	    pv2=Vprod(dr2,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+1./E2dummy.M+pv1*pv1/E1dummy.I+pv2*pv2/E2dummy.I;
	    C2=-2*(ps1-ps2-omega*pv1+E[fcj].omega*pv2);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;    
	      } 
	    domega=-(dp*pv1)/E1dummy.I;	    
	  }
	else
	  {
	    E1dummy.Versor(versor,l);
	    versor[0]=-versor[0];
	    versor[1]=-versor[1];
	    for(int i=0;i<2;i++)
	      {
		dr1[i]=l[i]-E1dummy.x[i];
	      }
	    ps1=Scalar(v,versor);
	    pv1=Vprod(dr1,versor);
	    double C1,C2;
	    C1=1./E1dummy.M+pv1*pv1/E1dummy.I;
	    C2=-2*(ps1-omega*pv1);
	    dp=C2/C1;
	    for(int i=0;i<2;i++)
	      {
		dv[i]=(versor[i]*dp)/E1dummy.M;	      
	      } 
	    domega=-(dp*pv1)/E1dummy.I;
	  }
      }
    return fc;
    } 
  double E()
  {
    return (Norm2(v)*M+I*omega*omega)/2.;
    }
 };
  
bool Overlap(double l[2],Ellisse E1,Ellisse E2)//overlap between two ellipses
{
  bool over;
  Ellisse E1r;
  Ellisse E2r;
  E1r.Copy(E1);
  E2r.Copy(E2);
  TRM(E1r.x,E1r.theta,-E1.theta,E1.x);
  TRM(E2r.x,E2r.theta,-E1.theta,E1.x); 
  double gamma=E1.A/E1.B;
  double M[2][2];
  E2r.BuildM(M,gamma);
  double D[2][2];
  double v[2];
  Eigen(D,v,E2r.theta,M);
  E1r.A=E1.B;
  E2r.AB(1./sqrt(D[0][0]),1./sqrt(D[1][1]));
  E2r.x[0]/=gamma;
  double mt2=E2r.theta;
  TR(E2r.x,E2r.theta,-mt2);
  double ln=Norm(E2r.x);
  double max=E2r.B;
  if(E2r.A>E2r.B) max=E2r.A;
  if(ln>max+E1r.A) {over=false;return over;}  
  else
    {
      double ln=E2r.NormE(E2r.x);
      if(ln<1)
	{
	  l[0]=E2r.x[0];
	  l[1]=E2r.x[1];
	  over=true;
	}
      else
	{
	  double distel=E2r.Distellisse(l);
	  over=(distel<E1r.A);
	}
    }
  TR(l,mt2);
  l[0]*=gamma;
  R(l,E1.theta,E1.x);
  return over;
}

double Contact(double l[2],Ellisse E1,Ellisse E2)//contact between two ellipses
{
  double over=1;
  double nover=0.;
  double delta=1.;
  double middle;
  Ellisse Es1;
  Ellisse Es2;
  while(delta>1e-6)
    {
      middle=nover+delta/2;
      Es1.Scale(E1,middle);
      Es2.Scale(E2,middle);
      double ll[2];
      if(Overlap(ll,Es1,Es2)) {over=middle;l[0]=ll[0];l[1]=ll[1];}
      else nover=middle;
      delta=over-nover;
    }
  return over;
}

bool Collision(double &dt,double l[2],Ellisse E1,Ellisse E2)//collision between ellipses
{
  Ellisse E1fin;
  Ellisse E2fin;
  Ellisse E1m;
  Ellisse E2m;
  E1fin.CopyAB(E1);
  E2fin.CopyAB(E2);  
  E1fin.Update(E1,dt);
  E2fin.Update(E2,dt);  
  bool over=Overlap(l,E1fin,E2fin);
  if(over)
    {
      E1m.CopyAB(E1);
      E2m.CopyAB(E2);       
      double delta=dt;
      double middle;
      double in=0;
      double fin=dt;
      while(delta>1e-12)
	{
	  middle=in+delta/2;
	  E1m.Update(E1,middle);
	  E2m.Update(E2,middle);      
	  if(Overlap(l,E1m,E2m)) fin=middle;
	  else in=middle;
	  delta=fin-in;
	}
      dt=middle;
    }
  return over;
}


bool Physical_evolution(int ne,Ellisse *E,double &phys_coll,int nw,Wall *P,double dtstep,double &t)
{
  double versor[2];
  double dp=0;
  double ps1=0,ps2=0;
  double pv1=0,pv2=0;
  double dr1[2];
  double dr2[2];
  double dt; 
  int collk;
  double loc_t=0;
  int conta=0;
  while(loc_t<dtstep)
    {
      collk=0;
      double l[2];
      dt=dtstep-loc_t;
      int coi,coj;
      for(int i=0;i<ne;i++)
	{
	  double ldt;	  
	  double ll[2];
	  for(int j=0;j<nw;j++)
	    {
	      ldt=dt;
	      E[i].Collision(ldt,ll,P[j]);
	      if((ldt<dt)&&(ldt>0))
		{
		  dt=ldt;
		  for(int k=0;k<2;k++) l[k]=ll[k];
		  coi=i;
		  coj=j;
		  collk=1;
		}
	    }
	  for(int j=i+1;j<ne;j++)
	    {
	      ldt=dt;
	      Collision(ldt,ll,E[i],E[j]);
	      if((ldt<dt)&&(ldt>0))
		{
		  dt=ldt;
		  for(int k=0;k<2;k++) l[k]=ll[k];
		  coi=i;
		  coj=j;
		  collk=2;
		}
	    }
	}
      if(dt<1e-4) {conta++;}//cout << "dt=" << dt << " conta=" << conta << endl;}
      if(conta>1000) return true;
      t+=dt;      
      loc_t+=dt;
      for(int i=0;i<ne;i++) E[i].Update(dt);
      if(collk)
	{
	  double en=0;
	  for(int i=0;i<ne;i++) en+=E[i].E();
	  if(collk==2)
	    {
	      E[coj].Versor(versor,l);	  
	      for(int i=0;i<2;i++)
		{
		  dr1[i]=l[i]-E[coi].x[i];
		  dr2[i]=l[i]-E[coj].x[i];
		}
	      ps1=Scalar(E[coi].v,versor);
	      ps2=Scalar(E[coj].v,versor);
	      pv1=Vprod(dr1,versor);
	      pv2=Vprod(dr2,versor);
	      double A,B;
	      A=1./E[coi].M+1./E[coj].M+pv1*pv1/E[coi].I+pv2*pv2/E[coj].I;
	      B=-2*(ps1-ps2-E[coi].omega*pv1+E[coj].omega*pv2);
	      dp=B/A;
	      for(int i=0;i<2;i++)
		{
		  E[coi].v[i]+=(versor[i]*dp)/E[coi].M;
		  E[coj].v[i]-=(versor[i]*dp)/E[coj].M;
		} 
	      E[coi].omega-=(dp*pv1)/E[coi].I;
	      E[coj].omega+=(dp*pv2)/E[coj].I;
	      phys_coll+=0.5*dp*dp*(1/E[coi].M+1/E[coj].M+pv1*pv1/E[coi].I+pv2*pv2/E[coj].I);
	    }
	  else
	    {
	      E[coi].Versor(versor,l);
	      versor[0]=-versor[0];
	      versor[1]=-versor[1];
	      for(int i=0;i<2;i++)
		{
		  dr1[i]=l[i]-E[coi].x[i];
		}
	      ps1=Scalar(E[coi].v,versor);
	      pv1=Vprod(dr1,versor);
	      double A,B;
	      A=1./E[coi].M+pv1*pv1/E[coi].I;
	      B=-2*(ps1-E[coi].omega*pv1);
	      dp=B/A;
	      double ep=(E[coi].v[0]*E[coi].v[0]+E[coi].v[1]*E[coi].v[1])*E[coi].M+E[coi].omega*E[coi].omega*E[coi].I;	      
	      for(int i=0;i<2;i++)
		{
		  E[coi].v[i]+=(versor[i]*dp)/E[coi].M;		  
		} 
	      E[coi].omega-=(dp*pv1)/E[coi].I;
	      phys_coll+=0.5*dp*dp*(1/E[coi].M+pv1*pv1/E[coi].I);
	    }      
	}
    }
  return false;
}

#endif  

