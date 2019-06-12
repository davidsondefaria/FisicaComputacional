#define alfa 1
#define beta 1
#define k 0.5
#define npassos 3000
#include <stdio.h>
double f(double x,double a,double b){
if(x>=1)
	return b*x+a-b;
if(x<=1)
	return b*x -a +b;
else
	return a*x;


}

double dxdr(double x,double y,double a,double b){
	return alfa*(y-x-f(x,a,b));

}

double dydr(double x,double y,double z){
return x-y+z;

}

double dzdr(double y){

return -1*beta*y;
}

main(){
long double x1d1,x1d2,x1d3,y1d1,y1d2,y1d3,z1d1,z1d2,z1d3,x2d1,x2d2,x2d3,y2d1,y2d2,y2d3,z2d1,z2d2,z2d3;
int i;
double a=1;
double b=-1;
x1d1=0.7;
y1d1=z1d1=0.0;
x1d2=k;
y1d2=0;
z1d2=-k;
x1d3=-k;
y1d3=0;
z1d3=k;
//fim do setup
FILE *D1,*D2,*D3;
char za[50],zb[50],zc[50];
sprintf(za,"dimensao1passos%d.dat",npassos);
sprintf(zb,"dimensao2passos%d.dat",npassos);
sprintf(zc,"dimensao3passos%d.dat",npassos);


D1=fopen(za,"w");
D2=fopen(zb,"w");
D3=fopen(zc,"w");
for(i=0;i<npassos;i++){
	fprintf(D1,"%d\t%.10Lf\t%.10Lf\t%.10Lf\n",i,x1d1,y1d1,z1d1);
	fprintf(D2,"%d\t%Lf\t%Lf\t%Lf\n",i,x1d2,y1d2,z1d2);
	fprintf(D3,"%d\t%Lf\t%Lf\t%Lf\n",i,x1d3,y1d3,z1d3);
	x2d1+=dxdr(x1d1,y1d1,a,b);
	x2d2+=dxdr(x1d2,y2d2,a,b);
	x2d3+=dxdr(x1d3,y1d3,a,b);
	y2d1+=dydr(x1d1,y1d1,z1d1);
	y2d2+=dydr(x1d2,y2d2,z1d2);
	y2d3+=dydr(x1d3,y1d3,z1d3);
	z2d1+=dzdr(y1d1);
	z2d2+=dzdr(y2d2);
	z2d3+=dzdr(y1d3);
	x1d1=x2d1;
	x1d2=x2d2;
	x1d3=x2d3;
	y1d1=y2d1;
	y1d2=y2d2;
	y1d3=y2d3;
	z1d1=z2d1;
	z1d2=z2d2;
	z1d3=z2d3;
	
}
}
