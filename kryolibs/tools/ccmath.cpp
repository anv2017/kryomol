/** Some of the subroutines from the ccmath 2.2.1 library*/

#include <math.h>
#include <stdlib.h>
#include "ccmath.h"

void ccmath::eigen(double *a,double *ev,int n)
{ double *dp;
  dp=(double *)calloc(n,sizeof(double));
  housev(a,ev,dp,n);
  qrevec(ev,a,dp,n); trnm(a,n);
  free(dp);
}

void ccmath::housev(double *a,double *d,double *dp,int n)
{ double sc,x,y,h;
  int i,j,k,m,e;
  double *qw,*qs,*pc,*p;
  qs=(double *)calloc(n,sizeof(double));
  for(j=0,pc=a; j<n-2 ;++j,pc+=n+1){
    m=n-j-1;
    for(i=1,sc=0.; i<=m ;++i) sc+=pc[i]*pc[i];
    if(sc>0.){ sc=sqrt(sc);
      if((x= *(pc+1))<0.){ y=x-sc; h=1./sqrt(-2.*sc*y);}
      else{ y=x+sc; h=1./sqrt(2.*sc*y); sc= -sc;}
      for(i=0,qw=pc+1; i<m ;++i){
        qs[i]=0.; if(i) qw[i]*=h; else qw[i]=y*h;
       }
      for(i=0,e=j+2,p=pc+n+1,h=0.; i<m ;++i,p+=e++){
        qs[i]+=(y=qw[i])* *p++;
	for(k=i+1; k<m ;++k){
          qs[i]+=qw[k]* *p; qs[k]+=y* *p++;
         }
        h+=y*qs[i];
       }
      for(i=0; i<m ;++i){
	qs[i]-=h*qw[i]; qs[i]+=qs[i];
       }
      for(i=0,e=j+2,p=pc+n+1; i<m ;++i,p+=e++){
        for(k=i; k<m ;++k) *p++ -=qw[i]*qs[k]+qs[i]*qw[k];
       }
     }
    d[j]= *pc; dp[j]=sc;
   }
  d[j]= *pc; dp[j]= *(pc+1); d[j+1]= *(pc+=n+1);
  free(qs);
  for(i=0,m=n+n,p=pc; i<m ;++i) *p-- =0.;
  *pc=1.; *(pc-=n+1)=1.; qw=pc-n;
  for(m=2; m<n ;++m,qw-=n+1){
    for(j=0,p=pc,*pc=1.; j<m ;++j,p+=n){
      for(i=0,qs=p,h=0.; i<m ;) h+=qw[i++]* *qs++;
      for(i=0,qs=p,h+=h; i<m ;) *qs++ -=h*qw[i++];
     }
    for(i=0,p=qw+m; i<n ;++i) *(--p)=0.;
    *(pc-=n+1)=1.;
   }
}

int ccmath::qrevec(double *ev,double *evec,double *dp,int n)
{ double cc,sc,d,x,y,h,tzr=1.e-15;
  int i,j,k,m,mqr=8*n;
  double *p;
  for(j=0,m=n-1;;++j){
    while(1){ if(m<1) return 0; k=m-1;
      if(fabs(dp[k])<=fabs(ev[m])*tzr) --m;
      else{ x=(ev[k]-ev[m])/2.; h=sqrt(x*x+dp[k]*dp[k]);
        if(m>1 && fabs(dp[m-2])>fabs(ev[k])*tzr) break;
	    if((cc=sqrt((1.+x/h)/2.))!=0.) sc=dp[k]/(2.*cc*h); else sc=1.;
        x+=ev[m]; ev[m--]=x-h; ev[m--]=x+h;
        for(i=0,p=evec+n*(m+1); i<n ;++i,++p){
	      h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
         }
       }
     }
    if(j>mqr) return -1;
    if(x>0.) d=ev[m]+x-h; else d=ev[m]+x+h;
    cc=1.; y=0.; ev[0]-=d;
    for(k=0; k<m ;++k){
      x=ev[k]*cc-y; y=dp[k]*cc; h=sqrt(x*x+dp[k]*dp[k]);
      if(k>0) dp[k-1]=sc*h;
      ev[k]=cc*h; cc=x/h; sc=dp[k]/h; ev[k+1]-=d; y*=sc;
      ev[k]=cc*(ev[k]+y)+ev[k+1]*sc*sc+d;
      for(i=0,p=evec+n*k; i<n ;++i,++p){
        h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
       }
     }
    ev[k]=ev[k]*cc-y; dp[k-1]=ev[k]*sc; ev[k]=ev[k]*cc+d;
   }
  return 0;
}

void ccmath::trnm(double *a,int n)
{ double s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p; *p++ = *q; *q=s; q+=n;
     }
   }
}

int ccmath::minv(double *a,int n)
{ int lc,*le; double s,t,tq=0.,zr=1.e-15;
  double *pa,*pd,*ps,*p,*q,*q0;
  int i,j,k,m;
  le=(int *)malloc(n*sizeof(int));
  q0=(double *)malloc(n*sizeof(double));
  for(j=0,pa=pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j>0){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
      	q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(le-j); free(q0); return -1;}
    *le++ =lc;
    if(lc!=j){
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
    *pd=t;
   }
  for(j=1,pd=ps=a; j<n ;++j){
    for(k=0,pd+=n+1,q= ++ps; k<j ;++k,q+=n) *q*= *pd;
   }
  for(j=1,pa=a; j<n ;++j){ ++pa;
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *q++ = *p;
    for(k=0; k<j ;++k){ t=0.;
      for(i=k,p=pa+k*n+k-j,q=q0+k; i<j ;++i) t-= *p++ * *q++;
      q0[k]=t;
     }
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,pd=pa=a+n*n-1; j>=0 ;--j){ --pa; pd-=n+1;
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *q++ = *p;
    for(k=n-1,ps=pa; k>j ;--k,ps-=n){ t= -(*ps);
      for(i=j+1,p=ps,q=q0; i<k ;++i) t-= *++p * *q++;
      q0[--m]=t;
     }
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *p= *q++;
   }
  for(k=0,pa=a; k<n-1 ;++k,++pa){
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
    for(j=0,ps=a; j<n ;++j,ps+=n){
      if(j>k){ t=0.; p=ps+j; i=j;}
      else{ t=q0[j]; p=ps+k+1; i=k+1;}
      for(; i<n ;) t+= *p++ *q0[i++];
      q0[j]=t;
     }
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,le--; j>=0 ;--j){
    for(k=0,p=a+j,q=a+ *(--le); k<n ;++k,p+=n,q+=n){
      t=*p; *p=*q; *q=t;
     }
   }
  free(le); free(q0);
  return 0;
}

void ccmath::vmul(double *vp,double *mat,double *v,int n)
{ double s,*q; int k,i;
  for(k=0; k<n ;++k){
    for(i=0,q=v,s=0.; i<n ;++i) s+= *mat++ * *q++;
    *vp++ =s;
   }
}

void ccmath::mmul(double *c,double *a,double *b,int n)
{ double *p,*q,s; int i,j,k;
  trnm(b,n);
  for(i=0; i<n ;++i,a+=n){
    for(j=0,q=b; j<n ;++j){
      for(k=0,p=a,s=0.; k<n ;++k) s+= *p++ * *q++;
      *c++ =s;
     }
   }
  trnm(b,n);
}
/** Resolve a symmetric linear system */
int ccmath::solvps(double *a,double *b,int n)
{ double *p,*q,*r,*s,t;
  int j,k;
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(q=a+j*n; q<p ;++q) *p-= *q* *q;
    if(*p<=0.) return -1;
    *p=sqrt(*p);
    for(k=j+1,q=p+n; k<n ;++k,q+=n){
      for(r=a+j*n,s=a+k*n,t=0.; r<p ;) t+= *r++ * *s++;
      *q-=t; *q/= *p;
     }
   }
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(k=0,q=a+j*n; k<j ;) b[j]-=b[k++]* *q++;
    b[j]/= *p;
   }
  for(j=n-1,p=a+n*n-1; j>=0 ;--j,p-=n+1){
    for(k=j+1,q=p+n; k<n ;q+=n) b[j]-=b[k++]* *q;
    b[j]/= *p;
   }
  return 0;
}

/** Resolve a general linear system*/
int ccmath::solv(double *a,double *b,int n)
{ int i,j,k,lc; double *ps,*p,*q,*pa,*pd;
  double *q0,s,t,tq=0.,zr=1.e-15;
  q0=(double *)calloc(n,sizeof(double));
  for(j=0,pa=a,pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
        q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(q0); return -1;}
    if(lc!=j){ t=b[j]; b[j]=b[lc]; b[lc]=t;
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
   }
  for(j=1,ps=b+1; j<n ;++j){
    for(k=0,p=a+n*j,q=b,t=0.; k<j ;++k) t+= *p++ * *q++;
    *ps++ -=t;
   }
  for(j=n-1,--ps,pd=a+n*n-1; j>=0 ;--j,pd-=n+1){
    for(k=j+1,p=pd,q=b+j,t=0.; k<n ;++k) t+= *++p * *++q;
    *ps-=t; *ps-- /= *pd;
   }
  free(q0); return 0;
}

void ccmath::mattr(double *a,double *b,int m,int n)
{ double *p; int i,j;
  for(i=0; i<n ;++i,++b)
    for(j=0,p=b; j<m ;++j,p+=n) *a++ = *p;
}




void ccmath::eigen(float *a,float *ev,int n)
{ float *dp;
  dp=(float *)calloc(n,sizeof(float));
  housev(a,ev,dp,n);
  qrevec(ev,a,dp,n); trnm(a,n);
  free(dp);
}

void ccmath::housev(float *a,float *d,float *dp,int n)
{ float sc,x,y,h;
  int i,j,k,m,e;
  float *qw,*qs,*pc,*p;
  qs=(float *)calloc(n,sizeof(float));
  for(j=0,pc=a; j<n-2 ;++j,pc+=n+1){
    m=n-j-1;
    for(i=1,sc=0.; i<=m ;++i) sc+=pc[i]*pc[i];
    if(sc>0.){ sc=sqrt(sc);
      if((x= *(pc+1))<0.){ y=x-sc; h=1./sqrt(-2.*sc*y);}
      else{ y=x+sc; h=1./sqrt(2.*sc*y); sc= -sc;}
      for(i=0,qw=pc+1; i<m ;++i){
        qs[i]=0.; if(i) qw[i]*=h; else qw[i]=y*h;
       }
      for(i=0,e=j+2,p=pc+n+1,h=0.; i<m ;++i,p+=e++){
        qs[i]+=(y=qw[i])* *p++;
	for(k=i+1; k<m ;++k){
          qs[i]+=qw[k]* *p; qs[k]+=y* *p++;
         }
        h+=y*qs[i];
       }
      for(i=0; i<m ;++i){
	qs[i]-=h*qw[i]; qs[i]+=qs[i];
       }
      for(i=0,e=j+2,p=pc+n+1; i<m ;++i,p+=e++){
        for(k=i; k<m ;++k) *p++ -=qw[i]*qs[k]+qs[i]*qw[k];
       }
     }
    d[j]= *pc; dp[j]=sc;
   }
  d[j]= *pc; dp[j]= *(pc+1); d[j+1]= *(pc+=n+1);
  free(qs);
  for(i=0,m=n+n,p=pc; i<m ;++i) *p-- =0.;
  *pc=1.; *(pc-=n+1)=1.; qw=pc-n;
  for(m=2; m<n ;++m,qw-=n+1){
    for(j=0,p=pc,*pc=1.; j<m ;++j,p+=n){
      for(i=0,qs=p,h=0.; i<m ;) h+=qw[i++]* *qs++;
      for(i=0,qs=p,h+=h; i<m ;) *qs++ -=h*qw[i++];
     }
    for(i=0,p=qw+m; i<n ;++i) *(--p)=0.;
    *(pc-=n+1)=1.;
   }
}

int ccmath::qrevec(float *ev,float *evec,float *dp,int n)
{ float cc,sc,d,x,y,h,tzr=1.e-15;
  int i,j,k,m,mqr=8*n;
  float *p;
  for(j=0,m=n-1;;++j){
    while(1){ if(m<1) return 0; k=m-1;
      if(fabs(dp[k])<=fabs(ev[m])*tzr) --m;
      else{ x=(ev[k]-ev[m])/2.; h=sqrt(x*x+dp[k]*dp[k]);
        if(m>1 && fabs(dp[m-2])>fabs(ev[k])*tzr) break;
	    if((cc=sqrt((1.+x/h)/2.))!=0.) sc=dp[k]/(2.*cc*h); else sc=1.;
        x+=ev[m]; ev[m--]=x-h; ev[m--]=x+h;
        for(i=0,p=evec+n*(m+1); i<n ;++i,++p){
	      h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
         }
       }
     }
    if(j>mqr) return -1;
    if(x>0.) d=ev[m]+x-h; else d=ev[m]+x+h;
    cc=1.; y=0.; ev[0]-=d;
    for(k=0; k<m ;++k){
      x=ev[k]*cc-y; y=dp[k]*cc; h=sqrt(x*x+dp[k]*dp[k]);
      if(k>0) dp[k-1]=sc*h;
      ev[k]=cc*h; cc=x/h; sc=dp[k]/h; ev[k+1]-=d; y*=sc;
      ev[k]=cc*(ev[k]+y)+ev[k+1]*sc*sc+d;
      for(i=0,p=evec+n*k; i<n ;++i,++p){
        h=p[0]; p[0]=cc*h+sc*p[n]; p[n]=cc*p[n]-sc*h;
       }
     }
    ev[k]=ev[k]*cc-y; dp[k-1]=ev[k]*sc; ev[k]=ev[k]*cc+d;
   }
  return 0;
}

void ccmath::trnm(float *a,int n)
{ float s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p; *p++ = *q; *q=s; q+=n;
     }
   }
}

int ccmath::minv(float *a,int n)
{ int lc,*le; float s,t,tq=0.,zr=1.e-15;
  float *pa,*pd,*ps,*p,*q,*q0;
  int i,j,k,m;
  le=(int *)malloc(n*sizeof(int));
  q0=(float *)malloc(n*sizeof(float));
  for(j=0,pa=pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j>0){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
      	q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(le-j); free(q0); return -1;}
    *le++ =lc;
    if(lc!=j){
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
    *pd=t;
   }
  for(j=1,pd=ps=a; j<n ;++j){
    for(k=0,pd+=n+1,q= ++ps; k<j ;++k,q+=n) *q*= *pd;
   }
  for(j=1,pa=a; j<n ;++j){ ++pa;
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *q++ = *p;
    for(k=0; k<j ;++k){ t=0.;
      for(i=k,p=pa+k*n+k-j,q=q0+k; i<j ;++i) t-= *p++ * *q++;
      q0[k]=t;
     }
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,pd=pa=a+n*n-1; j>=0 ;--j){ --pa; pd-=n+1;
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *q++ = *p;
    for(k=n-1,ps=pa; k>j ;--k,ps-=n){ t= -(*ps);
      for(i=j+1,p=ps,q=q0; i<k ;++i) t-= *++p * *q++;
      q0[--m]=t;
     }
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *p= *q++;
   }
  for(k=0,pa=a; k<n-1 ;++k,++pa){
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
    for(j=0,ps=a; j<n ;++j,ps+=n){
      if(j>k){ t=0.; p=ps+j; i=j;}
      else{ t=q0[j]; p=ps+k+1; i=k+1;}
      for(; i<n ;) t+= *p++ *q0[i++];
      q0[j]=t;
     }
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,le--; j>=0 ;--j){
    for(k=0,p=a+j,q=a+ *(--le); k<n ;++k,p+=n,q+=n){
      t=*p; *p=*q; *q=t;
     }
   }
  free(le); free(q0);
  return 0;
}

void ccmath::vmul(float *vp,float *mat,float *v,int n)
{ float s,*q; int k,i;
  for(k=0; k<n ;++k){
    for(i=0,q=v,s=0.; i<n ;++i) s+= *mat++ * *q++;
    *vp++ =s;
   }
}

void ccmath::mmul(float *c,float *a,float *b,int n)
{ float *p,*q,s; int i,j,k;
  trnm(b,n);
  for(i=0; i<n ;++i,a+=n){
    for(j=0,q=b; j<n ;++j){
      for(k=0,p=a,s=0.; k<n ;++k) s+= *p++ * *q++;
      *c++ =s;
     }
   }
  trnm(b,n);
}
/** Resolve a symmetric linear system */
int ccmath::solvps(float *a,float *b,int n)
{ float *p,*q,*r,*s,t;
  int j,k;
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(q=a+j*n; q<p ;++q) *p-= *q* *q;
    if(*p<=0.) return -1;
    *p=sqrt(*p);
    for(k=j+1,q=p+n; k<n ;++k,q+=n){
      for(r=a+j*n,s=a+k*n,t=0.; r<p ;) t+= *r++ * *s++;
      *q-=t; *q/= *p;
     }
   }
  for(j=0,p=a; j<n ;++j,p+=n+1){
    for(k=0,q=a+j*n; k<j ;) b[j]-=b[k++]* *q++;
    b[j]/= *p;
   }
  for(j=n-1,p=a+n*n-1; j>=0 ;--j,p-=n+1){
    for(k=j+1,q=p+n; k<n ;q+=n) b[j]-=b[k++]* *q;
    b[j]/= *p;
   }
  return 0;
}

/** Resolve a general linear system*/
int ccmath::solv(float *a,float *b,int n)
{ int i,j,k,lc; float *ps,*p,*q,*pa,*pd;
  float *q0,s,t,tq=0.,zr=1.e-15;
  q0=(float *)calloc(n,sizeof(float));
  for(j=0,pa=a,pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
        q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(q0); return -1;}
    if(lc!=j){ t=b[j]; b[j]=b[lc]; b[lc]=t;
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
   }
  for(j=1,ps=b+1; j<n ;++j){
    for(k=0,p=a+n*j,q=b,t=0.; k<j ;++k) t+= *p++ * *q++;
    *ps++ -=t;
   }
  for(j=n-1,--ps,pd=a+n*n-1; j>=0 ;--j,pd-=n+1){
    for(k=j+1,p=pd,q=b+j,t=0.; k<n ;++k) t+= *++p * *++q;
    *ps-=t; *ps-- /= *pd;
   }
  free(q0); return 0;
}

void ccmath::mattr(float *a,float *b,int m,int n)
{ float *p; int i,j;
  for(i=0; i<n ;++i,++b)
    for(j=0,p=b; j<m ;++j,p+=n) *a++ = *p;
}

double ccmath::lsqsv(double *x,int *pr,double *var,double *d,double *b,double *v,
		int m,int n,double th)
{ double ssq,sig,*y,*p;
  int i,k;
  y=(double *)calloc(n,sizeof(double));
  for(i=n,ssq=0.,p=b+n; i<m ;++i,++p) ssq+= *p* *p;
  for(i=k=0; i<n ;++i){
    if(d[i]<th){ y[i]=0.; ssq+=b[i]*b[i];}
    else{ y[i]=b[i]/d[i]; ++k;}
   }
  *pr=k;
  vmul(x,v,y,n);
  if(var!=NULL && m>n){
    sig=ssq/(double)(m-n);
    for(i=0; i<n ;++i){
      if(d[i]<th) y[i]=0.; else y[i]=sig/(d[i]*d[i]);
     }
    smgen(var,y,v,n);
   }
  free(y);
  return ssq;
}

void ccmath::smgen(double *a,double *eval,double *evec,int n)
{ double *p,*q,*ps,*r,*s,*t,*v=evec+n*n;
  for(ps=a,p=evec; p<v ;p+=n){
    for(q=evec; q<v ;q+=n,++ps){ *ps=0.;
      for(r=eval,s=p,t=q; r<eval+n ;)
        *ps+= *r++ * *s++ * *t++;
     }
   }
}

int ccmath::sv2lsq(double *d,double *a,double *b,int m,double *v,int n)
{ double *p,*p1,*q,*pp,*w,*e;
  double s,t,h,r,sv;
  int i,j,k,mm,nm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m+n,sizeof(double)); e=w+m;
  for(i=0,mm=m,p=a; i<n ;++i,--mm,p+=n+1){
    if(mm>1){ h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; w[0]+=h;
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r=r*s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(j=0,q=b+i,r=0.; j<mm ;) r+=w[j++]* *q++;
        for(j=0,q=b+i,r*=s; j<mm ;) *q++ -=w[j++]*r;
       }
      d[i]= -h;
     }
    if(mm==1) d[i]= *p;
   }
  for(i=0,p=a; i<n ;++i){
    for(j=0; j<n ;++j,++p){
      if(j<i) *p=0.;
      else if(j==i) *p=d[i];
     }
   }
  for(i=0,mm=n,nm=n-1,p=a; i<n ;++i,--mm,--nm,p+=n+1){
    if(i && mm>1){ sv=h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
        sv=1.+fabs(*p/h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  for(j=0,q=p+k,r*=s; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(j=0,q=b+i,r=0.; j<mm ;) r+=w[j++]* *q++;
        for(j=0,q=b+i,r*=s; j<mm ;) *q++ -=r*w[j++];
       }
      *p=sv; d[i]= -h;
     }
    if(mm==1) d[i]= *p;
    p1=p+1;
    if(nm>1){ sv=h=0.;
      for(j=0,q=p1,s=0.; j<nm ;++j,++q) s+= *q* *q;
      if(s>0.){
	h=sqrt(s); if(*p1<0.) h= -h;
        sv=1.+fabs(*p1/h);
	s+= *p1*h; s=1./s; t=1./(*p1+=h);
	for(k=n,ms=n*(n-i); k<ms ;k+=n){
	  for(j=0,q=p1,pp=p1+k,r=0.; j<nm ;++j) r+= *q++ * *pp++;
	  for(j=0,q=p1,pp=p1+k,r*=s; j<nm ;++j) *pp++ -=r* *q++;
	 }
	for(j=1,q=p1+1; j<nm ;++j) *q++ *=t;
       }
      *p1=sv; e[i]= -h;
     }
    if(nm==1) e[i]= *p1;
   }
  ldvmat(a,v,n);
  qrbdbv(d,e,b,v,n);
  for(i=0; i<n ;++i){
    if(d[i]<0.){ d[i]= -d[i];
      for(j=0,p=v+i; j<n ;++j,p+=n) *p= - *p;
     }
   }
  free(w);
  return 0;
} 

void ccmath::ldvmat(double *a,double *v,int n)
{ double *p0,*q0,*p,*q,*qq;
  double h,s;
  int i,j,k,mm;
  for(i=0,mm=n*n,q=v; i<mm ;++i) *q++ =0.;
  *v=1.; q0=v+n*n-1; *q0=1.; q0-=n+1;
  p0=a+n*n-n-n-1;
  for(i=n-2,mm=1; i>0 ;--i,p0-=n+1,q0-=n+1,++mm){
    if(*(p0-1)!=0.){
      for(j=0,p=p0,h=1.; j<mm ;++j,++p) h+= *p* *p;
      h= *(p0-1); *q0=1.-h;
      for(j=0,q=q0+n,p=p0; j<mm ;++j,q+=n) *q= -h* *p++; 
      for(k=i+1,q=q0+1; k<n ;++k){
        for(j=0,qq=q+n,p=p0,s=0.; j<mm ;++j,qq+=n) s+= *qq* *p++;
        s*=h;
        for(j=0,qq=q+n,p=p0; j<mm ;++j,qq+=n) *qq-=s* *p++;
        *q++ = -s;
       }
     }
    else{
      *q0=1.;
      for(j=0,p=q0+1,q=q0+n; j<mm ;++j,q+=n) *q= *p++ =0.;
     }
   }
}

int ccmath::qrbdbv(double *d,double *e,double *b,double *v,int n)
{ int i,j,k,nn,jj,nm;
  double u,x,y,f,g,c,s,t,w,*p,*q;
  for(j=1,t=fabs(d[0]); j<n ;++j)
    if((s=fabs(d[j])+fabs(e[j-1]))>t) t=s;
  t*=1.e-15; nn=100*n; nm=n;
  for(j=0; n>1 && j<nn ;++j){
    for(k=n-1; k>0 ;--k){
      if(fabs(e[k-1])<t) break;
      if(fabs(d[k-1])<t){
        for(i=k,s=1.,c=0.; i<n ;++i){
          f=s*e[i-1]; g=d[i]; e[i-1]*=c;
          d[i]=u=sqrt(f*f+g*g); s= -f/u; c=g/u;
          p=b+k-1; q=b+i;
          w=c* *p+s* *q; *q=c* *q-s* *p; *p=w;
	 }
        break;
       }
     }
    y=d[k]; x=d[n-1]; u=e[n-2];
    f=(y+x)*(y-x)-u*u; s=y*e[k]; g=s+s;
    u=sqrt(f*f+g*g);
    c=sqrt((u+f)/(u+u)); s/=(c*u);
    for(i=k; i<n-1 ;++i){
      g=e[i];
      if(i>k){
	f=s*e[i]; g*=c;
	e[i-1]=u=sqrt(x*x+f*f);
	c=x/u; s=f/u;
       }
      f=c*y+s*g; g=c*g-s*y;
      for(jj=0,p=v+i; jj<nm ;++jj,p+=nm){
        q=p+1;
        w=c* *p+s* *q; *q=c* *q-s* *p; *p=w;
       }
      s*=d[i+1]; d[i]=u=sqrt(f*f+s*s);
      y=c*d[i+1]; c=f/u; s/=u;
      x=c*g+s*y; y=c*y-s*g;
      p=b+i; q=p+1;
      w=c* *p+s* *q; *q=c* *q-s* *p; *p=w;
     }
    e[n-2]=x; d[n-1]=y;
    if(fabs(x)<t) --n;
    if(n==k+1) --n; 
   }
  return j;
}

int ccmath::svdlsq(double *d,double *a,double *b,int m,double *v,int n)
{ double *p,*p1,*q,*pp,*w,*e;
  double s,h,r,t,sv;
  int i,j,k,mm,nm,ms;
  if(m<n) return -1;
  w=(double *)calloc(m+n,sizeof(double)); e=w+m;
  for(i=0,mm=m,nm=n-1,p=a; i<n ;++i,--mm,--nm,p+=n+1){
    if(mm>1){ h=0.;
      for(j=0,q=p,s=0.; j<mm ;++j,q+=n){
	w[j]= *q; s+= *q* *q;
       }
      if(s>0.){
	h=sqrt(s); if(*p<0.) h= -h;
	s+= *p*h; s=1./s; t=1./(w[0]+=h);
	for(k=1,ms=n-i; k<ms ;++k){
	  for(j=0,q=p+k,r=0.; j<mm ;q+=n) r+=w[j++]* *q;
	  r*=s;
	  for(j=0,q=p+k; j<mm ;q+=n) *q-=r*w[j++];
	 }
        for(j=0,q=b+i,r=0.; j<mm ;) r+=w[j++]* *q++;
        for(j=0,q=b+i,r*=s; j<mm ;) *q++ -=r*w[j++];
       }
      d[i]= -h;
     }
    if(mm==1) d[i]= *p;
    p1=p+1; sv=h=0.;
    if(nm>1){
      for(j=0,q=p1,s=0.; j<nm ;++j,++q) s+= *q* *q;
      if(s>0.){
	h=sqrt(s); if(*p1<0.) h= -h;
        sv=1.+fabs(*p1/h);
	s+= *p1*h; s=1./s; t=1./(*p1+=h);
	for(k=n,ms=n*(m-i); k<ms ;k+=n){
	  for(j=0,q=p1,pp=p1+k,r=0.; j<nm ;++j) r+= *q++ * *pp++;
	  r*=s;
	  for(j=0,q=p1,pp=p1+k; j<nm ;++j) *pp++ -=r* *q++;
	 }
        for(j=1,q=p1+1; j<nm ;++j) *q++ *=t;
       }
      *p1=sv; e[i]= -h;
     }
    if(nm==1) e[i]= *p1;
   }
  ldvmat(a,v,n);
  qrbdbv(d,e,b,v,n);
  for(i=0; i<n ;++i){
    if(d[i]<0.){ d[i]= -d[i];
      for(j=0,p=v+i; j<n ;++j,p+=n) *p= - *p;
     }
   }
  free(w);
  return 0;
}
