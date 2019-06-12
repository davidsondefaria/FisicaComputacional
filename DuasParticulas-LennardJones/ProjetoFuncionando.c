#include <stdio.h>
#include <math.h>
int contorno1(int c, int TAM){
	if(abs(c)>TAM){
		c=c-2*TAM*(c/abs(c));
	}
	return c;
}

main(){
	FILE *energia,*posicao;
	double   m1, m2, ka, r0, dt;
	double   r1[2], r2[2], r1a[2], r2a[2], r1d[2], r2d[2], v1[2], v2[2];
	double   pi, def, mred, tau, ry, rz, r12x, r12y, r12, dpx[2], dpy[2], Lx, Ly;
	double   Ecin, Epot, Etot, for1[2], for2[2], somaecin, somaepot, somaetot;
	int  npassos, i, j;
	//manuseando arquivos
	energia=fopen("energia.dat","w");
	posicao=fopen("posicao.dat","w");
	//pedindo parametros
	puts("digite massa 1 e massa 2");
	scanf("%lf, %lf",&m1, &m2);  //  Massas das Particulas
	puts("digite a distancia de equilibrio");
	scanf("%lf",&r0);              // Distância de Equilibrio
	puts("digite o passo temporal");
	scanf("%lf",&dt);              // Passo temporal
	puts("digite o numero de passos");
	scanf("%d",&npassos);         // Numero de passos de integracao
	puts("digite as cordenadas x1 e y1");
	scanf("%lf, %lf", &r1[1], &r1[2]);    // Posicao inicial - Atomo 1
	puts("digite as cordenadas x2 e y2");
	scanf("%lf, %lf", &r2[1], &r2[2]);    // Posicao inicial - Atomo 2
	puts("digite a velocidade em x1 e y1");
	scanf("%lf, %lf", &v1[1], &v1[2]);    // Velocidades Inicial - Atomo 1
	puts("digite a velocidade em x2 e y2");
	scanf("%lf, %lf", &v2[1], &v2[2]);    // Velocidades Inicial - Atomo 2
	puts("digite o tamanho da rede em x e y");
	scanf("%lf, %lf",&Lx, &Ly);          // Dimensoes Laterais
	
	//definindo variaveis
	pi = 3.141592653590;
	ka = 12.0*((26.0/(pow(r0,14))) - (7.0/(pow(r0,8))));		//constante elastica
	def = sqrt(pow((r1[1] - r2[1]),2)+ pow(r1[2] - r2[2],2));	//distancia entre particulas
	mred = m1*m2/(m1 + m2);					//centro de massa
	tau = 2*pi*sqrt(mred/ka);
	
	//calculo de energias		
	Ecin = (m1*v1[1]*v1[1] + m2*v2[1]*v2[1] + m1*v1[2]*v1[2] + m2*v2[2]*v2[2])/2;
	Epot = 4.0*((1./(pow(def,12)))-(1./(pow(def,6))));
	Etot = Ecin + Epot;
	
// Primeira integracao utilizando o algoritmo de Euler:
	for (i=0;i<2;i++){
		r1a[i] = r1[i];
		r2a[i] = r2[i];
	}
	r12x = r1[1] - r2[1];
	r12y = r1[2] - r2[2];
	//calculo de forças
	for1[1] =  24.*((2./(pow(def,13)))-(1./(pow(def,7))))*r12x;
	for2[1] = -for1[1];
	for1[2] =  24.*((2./(pow(def,13)))-(1./(pow(def,7))))*r12y;
	for2[2] = -for1[2];
	
// Integrando as posicoes pelo metodo de Euler:
	
	r1[1] = r1[1] + v1[1]*dt + (for1[1]/(2*m1))*dt*dt;
	r2[1] = r2[1] + v2[1]*dt + (for2[1]/(2*m2))*dt*dt;
	r1[2] = r1[2] + v1[2]*dt + (for1[2]/(2*m1))*dt*dt;
	r2[2] = r2[2] + v2[2]*dt + (for2[2]/(2*m2))*dt*dt;
	
	
	somaecin = 0.00;
	somaepot = 0.00;
	somaetot = 0.00;
	
// Loop principal: vamos fazer aqui NPASSOS integracoes das equacoes de 1
// Movimento de cada particula:

	for (i=1;i<=npassos;i++){
	
		// Distancias entre as particulas (com condicoes de contorno periodicas):

		r12x = r1[1] - r2[1];
		r12y = r1[2] - r2[2];
		r12x = r12x - round(r12x/(2*Lx))*(2.*Lx);
		r12x = r12x - round(r12x/(2*Lx))*(2.*Lx);
		r12 = sqrt(pow(r12x,2) + pow(r12y,2));
		
		// Forca sobre as particulas:

		for1[1] =  24.*((2./(pow(r12,13)))-(1./(pow(r12,7))))*r12x;
		for2[1] = -for1[1];
		for1[2] =  24.*((2./(pow(r12,13)))-(1./(pow(r12,7))))*r12y;
		for2[2] = -for1[2];
		
		// Integracoes pelo algoritmo de Verlet:
		
		r1d[1] = 2.*r1[1] - r1a[1] + (for1[1]/m1)*dt*dt;
		r2d[1] = 2.*r2[1] - r2a[1] + (for2[1]/m2)*dt*dt;
		r1d[2] = 2.*r1[2] - r1a[2] + (for1[2]/m1)*dt*dt;
		r2d[2] = 2.*r2[2] - r2a[2] + (for2[2]/m2)*dt*dt;
		
		// Velocidades:
		
		dpx[1] = r1d[1] - r1a[1];
		dpx[2] = r2d[1] - r2a[1];
		dpy[1] = r1d[2] - r1a[2];
		dpy[2] = r2d[2] - r2a[2];
	
		// CCP:

		dpx[1] = dpx[1] - round(dpx[1]/(2.0*Lx))*(2.0*Lx);
		dpx[2] = dpx[2] - round(dpx[2]/(2.0*Lx))*(2.0*Ly);
		dpy[1] = dpy[1] - round(dpy[1]/(2.0*Ly))*(2.0*Lx);
		dpy[2] = dpy[2] - round(dpy[2]/(2.0*Ly))*(2.0*Ly);	
		
		v1[1] = dpx[1]/(2.0*dt);
		v1[2] = dpy[1]/(2.0*dt);
		v2[1] = dpx[2]/(2.0*dt);
		v2[2] = dpy[2]/(2.0*dt);
		
		// Emprego de condicoes de contorno periodicas nas posicoes integradas:
		
		r1d[1] = r1d[1] - round(r1d[1]/(2.0*Lx))*(2.0*Lx);
		r1d[2] = r1d[2] - round(r1d[2]/(2.0*Ly))*(2.0*Ly);
		r2d[1] = r2d[1] - round(r2d[1]/(2.0*Lx))*(2.0*Lx);
		r2d[2] = r2d[2] - round(r2d[2]/(2.0*Ly))*(2.0*Ly);

		//write(*,*) 'PASSO = ', i
		//write(*,'(4f12.5)') r1[1], r2[1], r12x, r12;     
		//write(*,'(4f12.5)') dpx[1], dpx[2], dpy[1], dpy[2];
		//write(*,'(4f12.5)') v1[1], v1[2], v2[1], v2[2];
		
		// Calculo das energias e guardando-as no arquivo de saida:
		
		Ecin = (m1*v1[1]*v1[1] + m2*v2[1]*v2[1] + m1*v1[2]*v1[2] + m2*v2[2]*v2[2])/2.0;
		Epot = ((4./(pow(r12,12)))-(4/(pow(r12,6))));
		Etot = Ecin + Epot;
	
		// Escrevendo as energias e posicoes:
		
		printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, Ecin, Epot, Etot,r1d[1],r1d[2],r2d[1],r2d[2]);
		fprintf(energia,"%d\t%f\t%f\t%f\n",i,Ecin,Epot,Etot);
		fprintf(posicao,"%f\t%f\t%f\t%f\n ",r1d[1],r1d[2],r2d[1],r2d[2]);
		
		
		// Acumulando as Energias:
		
		somaecin = somaecin + Ecin;
		somaepot = somaepot + Epot;
		somaetot = somaetot + Etot;
		
		// Atualizacoes
		
		r1a[1] = r1[1];
		r1a[2] = r1[2];
		r2a[1] = r2[1];
		r2a[2] = r2[2];
		r1[1] = r1d[1];
		r1[2] = r1d[2];
		r2[1] = r2d[1];
		r2[2] = r2d[2];

// fim do loop sobre npassos
}

// Calculando os valores medios das energias:
		fclose(energia);
		fclose(posicao);
		somaecin = somaecin/(float)(npassos);
		somaepot = somaepot/(float)(npassos);
		somaetot = somaetot/(float)(npassos);
/*
print*
print*, "-----------------------------------------------------------------"
print*, "                 VALORES MEDIOS DAS ENERGIAS                     "
*/
puts("-----------------------------------------------------------------");
printf("\nEnergia Cinetica Media = %lf\n", somaecin);
printf("\nEnergia Potencial Media = %lf\n", somaepot);
printf("\nEnergia Mecanica Media = %lf\n", somaetot);
puts("-----------------------------------------------------------------");

// Fim do programa

}

//-------------------------------------------------------------------------

