clear all;clc;close all;disp('program vibration fullerene ');
nat=60; n=3*nat; % nat: Nombres d atomes 
C=4.5;% constantes de forces en mdyn/A
Def=0.095;
teta0 = (120.7/180.0)*3.1415;
n=3*nat; % Nombres d atomes * nombres de degr√©s de libert√© 
unite=1303.23; % unite pour avoir des cm-1
id=10000; % nombre de point DOS
cof=1.5; % facteur d echelle pour les frequences
atom='C' ;
masse = 12;
zat = 6;
x=[3.451266498; 3.451266498;-3.451266498;-3.451266498; 0.685000000;-0.685000000; 0.685000000;-0.685000000; 0.000000000; 0.000000000; 0.000000000; 0.000000000; 3.003809890; 3.003809890; 3.003809890; 3.003809890;-3.003809890;-3.003809890;-3.003809890;-3.003809890; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 2.580456608; 2.580456608; 2.580456608; 2.580456608;-2.580456608;-2.580456608;-2.580456608;-2.580456608; 0.724000000; 0.724000000;-0.724000000;-0.724000000; 0.724000000; 0.724000000;-0.724000000;-0.724000000; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890]; 
y =[0.685000000;-0.685000000; 0.685000000;-0.685000000; 0.000000000; 0.000000000; 0.000000000; 0.000000000; 3.451266498; 3.451266498;-3.451266498;-3.451266498; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 3.003809890; 3.003809890; 3.003809890; 3.003809890;-3.003809890;-3.003809890;-3.003809890;-3.003809890; 0.724000000; 0.724000000;-0.724000000;-0.724000000; 0.724000000; 0.724000000;-0.724000000;-0.724000000; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.580456608; 2.580456608; 2.580456608; 2.580456608;-2.580456608;-2.580456608;-2.580456608;-2.580456608];
z=[ 0.000000000;0.000000000; 0.000000000; 0.000000000; 3.451266498; 3.451266498;-3.451266498;-3.451266498; 0.685000000;-0.685000000; 0.685000000;-0.685000000; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 1.171456608;-1.171456608; 3.003809890; 3.003809890; 3.003809890; 3.003809890;-3.003809890;-3.003809890;-3.003809890;-3.003809890; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 1.409000000; 1.409000000;-1.409000000;-1.409000000; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.279809890;-2.279809890; 2.580456608; 2.580456608; 2.580456608; 2.580456608;-2.580456608;-2.580456608;-2.580456608;-2.580456608; 0.724000000; 0.724000000;-0.724000000;-0.724000000; 0.724000000; 0.724000000;-0.724000000;-0.724000000];
disp(size(x));
disp(size(y));
disp(size(z));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phidef=zeros(n,n);phi=zeros(n,n);phiel=zeros(n,n);
%tij=zeros(3,3);tik=zeros(3,3);tjk=zeros(3,3);tki=zeros(3,3);tkj=zeros(3,3);tji=zeros(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%************************************************************************
%*     CONSTANTES DES ELONGATIONS
%************************************************************************

%c constantes de forces aij en mdyn/a

 aij=[4.0,2.35,1.21,1.21,-1.05,0,0];

%************************************************************************
%*   CONSTANTES DES DEFORMATIONS%c 
%constantes de forces aijk en mdyn.a.rad-2
%************************************************************************
aijk=zeros(7,7);      
       aijk(2,2)=1;
       aijk(1,2)=0.25;
       
       aijk(1,3)=0.095;
       aijk(1,4)=0.095;
       aijk(2,3)=0.095;
       aijk(2,4)=0.095;
       
       aijk(1,5)=0.325;
       aijk(2,5)=0.325;
       
%c symÈtrisation des constantes de force et des angles d'Èquilibre
       for i=1:7
          for j=i:7
                aijk(j,i)=aijk(i,j);
            end
            end
%************************************************************************
fid0 = fopen('structFullerene.xyz', 'wt');
    fprintf(fid0,'%d\n\n', nat);
for i=1:nat
    fprintf(fid0,'%c  \t %f \t %f \t %f\n', atom,x(i),y(i),z(i));
end
%---------------------------------------------------%
	 	 kk=0;
		 nx(1)=0;
	  for i=1:nat	   
	     for j=1:nat
	        xij=x(i)-x(j);
			yij=y(i)-y(j);
            zij=z(i)-z(j);
			rij=sqrt(xij*xij+yij*yij+zij*zij);
          if j~=i
		    if rij<1.42
			kk=kk+1;
			iac(kk)=j;
			ntype(kk)=1;
            else
                if rij>1.42 && rij <1.5
                    kk=kk+1;
                    iac(kk)=j;
                    ntype(kk)=2;
                else
                    if rij > 2.3 && rij <2.4
                        kk=kk+1;
                        iac(kk)=j;
                        ntype(kk)=3;
                    else
                        if rij >2.4 && rij<2.5
                         kk=kk+1;
                        iac(kk)=j;
                        ntype(kk)=4;
                        else
                            if rij > 2.5 && rij < 3.0
                             kk=kk+1;
                            iac(kk)=j;
                             ntype(kk)=5;
                            else
                                if rij>3.0 && rij < 3.6
                                 kk=kk+1;
                                 iac(kk)=j;
                                 ntype(kk)=6;
                                else
                                    if rij >3.6 && rij < 3.8
                                     kk=kk+1;
                                    iac(kk)=j;
                                    ntype(kk)=7;
                                    end
                                end
                            end
                        end
                    end
                end
            end
          end
         end
		 nx(i+1)=kk;
      end
      
%******************************************
for i=1:nat
   for k=nx(i)+1:nx(i+1)
        j=iac(k); mj=ntype(k);
        if j>i
        xij=x(i)-x(j);
        yij=y(i)-y(j);
        zij=z(i)-z(j);
        rij2=xij*xij+yij*yij+zij*zij;
        rij=sqrt(rij2);   
  %Derivees du second ordre du potentiel stretching pour avoir PHI
             ar(1)=xij;              
             ar(2)=yij;
             ar(3)=zij;             
           for k1=1:3
                kp1=k1+(i-1)*3;% pour la premi√®re atome i = 1 , kp1 = k1 pour la deuxieme atome i = 2, kp1 = k1 + 3, i = 3 kp = k1 + 6
              for k2=1:3
                   kp2=k2+(j-1)*3;
                       vij=-aij(mj)*ar(k1)*ar(k2)/rij2;
                   phiel(kp1,kp2)=vij;
                  
                end %k2            
           end % k1
        %end
   end %j 
   end
end
%------ Potentiel de Deformation ----------%
kkkk=0;
for i=1:nat
   for kk=nx(i)+1:nx(i+1)
        j=iac(kk); mj=ntype(kk);
        xji=x(j)-x(i);
        yji=y(j)-y(i);
        zji=z(j)-z(i);  
        rji2=xji*xji+yji*yji+zji*zji;
        rji=sqrt(rji2);
  
            for ll=nx(j)+1:nx(j+1) 
                k=iac(ll);
                if (k ~= j) && (k~=i) 
                   if k>i
                kkkk=kkkk+1;
                mk=ntype(ll);
                                       
         tij=zeros(3,3);tik=zeros(3,3);tjk=zeros(3,3);tki=zeros(3,3);tkj=zeros(3,3);tji=zeros(3,3);
           
            xjk=x(j)-x(k);
            yjk=y(j)-y(k);
            zjk=z(j)-z(k);            
            rjk2=xjk*xjk+yjk*yjk+zjk*zjk;
            rjk=sqrt(rjk2);
            rjirjk=xji*xjk+yji*yjk+zji*zjk;
            costet=rjirjk/(rji*rjk);
            costet0=costet;
            darc=sqrt(1-costet^2);
            
            dtetaxi=(((-xjk/(rji*rjk))+((xji*rjk*rjirjk)/(rji*(rji*rjk)^(2)))))/darc;
            dtetayi=(((-yjk/(rji*rjk))+((yji*rjk*rjirjk)/(rji*(rji*rjk)^(2)))))/darc;
            dtetazi=(((-zjk/(rji*rjk))+((zji*rjk*rjirjk)/(rji*(rji*rjk)^(2)))))/darc;
            dtetaxj=((xjk+xji)/(rji*rjk)-((rjirjk/(rji*rjk)^(2))*((xji*rjk/rji)+(xjk*rji/rjk))))/darc;
            dtetayj=((((yjk+yji)/(rji*rjk))-((rjirjk/(rji*rjk)^(2))*((yji*rjk/rji)+(yjk*rji/rjk)))))/darc;
            dtetazj=((((zjk+zji)/(rji*rjk))-((rjirjk/(rji*rjk)^(2))*((zji*rjk/rji)+(zjk*rji/rjk)))))/darc;
            dtetaxk=(((-xji/(rji*rjk))+((xjk*rji*rjirjk)/(rjk*(rji*rjk)^(2)))))/darc;
            dtetayk=(((-yji/(rji*rjk))+((yjk*rji*rjirjk)/(rjk*(rji*rjk)^(2)))))/darc;
            dtetazk=(((-zji/(rji*rjk))+((zjk*rji*rjirjk)/(rjk*(rji*rjk)^(2)))))/darc;
            %%%%%%%%% tij
            tij(1,1)=aijk(mj,mk)*dtetaxi*dtetaxj;
              tji(1,1)=tij(1,1);
              tij(1,2)=aijk(mj,mk)*dtetaxi*dtetayj;
              tji(2,1)=tij(1,2);
              tij(1,3)=aijk(mj,mk)*dtetaxi*dtetazj;
              tji(3,1)=tij(1,3);
              tij(2,2)=aijk(mj,mk)*dtetayi*dtetayj;
              tji(2,2)=tij(2,2);
              tij(2,3)=aijk(mj,mk)*dtetayi*dtetazj;
              tji(3,2)=tij(2,3);
              tij(3,3)=aijk(mj,mk)*dtetazi*dtetazj;
              tji(3,3)=tij(3,3);
              tij(2,1)=aijk(mj,mk)*dtetayi*dtetaxj;
              tji(1,2)=tij(2,1);
              tij(3,1)=aijk(mj,mk)*dtetazi*dtetaxj;
              tji(1,3)=tij(3,1);
              tij(3,2)=aijk(mj,mk)*dtetazi*dtetayj;
              tji(2,3)=tij(3,2);

              tik(1,1)=aijk(mj,mk)*dtetaxi*dtetaxk;
              tki(1,1)=tik(1,1);
              tik(1,2)=aijk(mj,mk)*dtetaxi*dtetayk;
              tki(2,1)=tik(1,2);
              tik(1,3)=aijk(mj,mk)*dtetaxi*dtetazk;
              tki(3,1)=tik(1,3);
              tik(2,2)=aijk(mj,mk)*dtetayi*dtetayk;
              tki(2,2)=tik(2,2);
              tik(2,3)=aijk(mj,mk)*dtetayi*dtetazk;
              tki(3,2)=tik(2,3);
              tik(3,3)=aijk(mj,mk)*dtetazi*dtetazk;
              tki(3,3)=tik(3,3);
              tik(2,1)=aijk(mj,mk)*dtetayi*dtetaxk;
              tki(1,2)=tik(2,1);
              tik(3,1)=aijk(mj,mk)*dtetazi*dtetaxk;
              tki(1,3)=tik(3,1);
              tik(3,2)=aijk(mj,mk)*dtetazi*dtetayk;
              tki(2,3)=tik(3,2);

              tjk(1,1)=aijk(mj,mk)*dtetaxj*dtetaxk;
              tkj(1,1)=tjk(1,1);
              tjk(1,2)=aijk(mj,mk)*dtetaxj*dtetayk;
              tkj(2,1)=tjk(1,2);
              tjk(1,3)=aijk(mj,mk)*dtetaxj*dtetazk;
              tkj(3,1)=tjk(1,3);
              tjk(2,2)=aijk(mj,mk)*dtetayj*dtetayk;
              tkj(2,2)=tjk(2,2);
              tjk(2,3)=aijk(mj,mk)*dtetayj*dtetazk;
              tkj(3,2)=tjk(2,3);
              tjk(3,3)=aijk(mj,mk)*dtetazj*dtetazk;
              tkj(3,3)=tjk(3,3);
              tjk(2,1)=aijk(mj,mk)*dtetayj*dtetaxk;
              tkj(1,2)=tjk(2,1);
              tjk(3,1)=aijk(mj,mk)*dtetazj*dtetaxk;
              tkj(1,3)=tjk(3,1);
              tjk(3,2)=aijk(mj,mk)*dtetazj*dtetayk;
              tkj(2,3)=tjk(3,2);
              
              
          for k1=1:3
                 kp1=k1+(i-1)*3;
                 for k2=1:3
                    kp2=k2+(j-1)*3;
                    if j > i
                       phidef(kp1,kp2)=phidef(kp1,kp2)+tij(k1,k2);
                    else
                       phidef(kp2,kp1)=phidef(kp2,kp1)+tji(k2,k1);
                    end
                 end
          end
           
            for k1=1:3
                 kp1=k1+(i-1)*3;
                 for k2=1:3
                    kp2=k2+(k-1)*3;
                    phidef(kp1,kp2)=phidef(kp1,kp2)+tik(k1,k2);
                 end
            end
            
            for k1=1:3
                 kp1=k1+(j-1)*3;
                 for k2=1:3
                    kp2=k2+(k-1)*3;
                    if k > j
                       phidef(kp1,kp2)=phidef(kp1,kp2)+tjk(k1,k2);
                    else
                       phidef(kp2,kp1)=phidef(kp2,kp1)+tkj(k2,k1);
                    end
                 end
            end
            
                   end
                end
                end
            end
end

for i=1:n
          for j=1:n
             phi(i,j)=phiel(i,j)+phidef(i,j);
          end
end

%************************************************************************
%*     SYMETRISATION DE LA MATRICE PHI2
%************************************************************************
      for i=1:n
        for j=1:n
            phi(j,i)=phi(i,j);
        end
      end

disp('CALCUL DES ELEMENTS DIAGONAUX DE LA MATRICE PHI');
%************************************************************************
       for k0=1:3
            for i=k0:3:n-3+k0
                for k=1:3
                s=0;
                   for j=k:3:n-3+k
                   s=s+phi(i,j);
                   end
                phi(i,i+k-k0)=-s;
                end
            end
       end
 disp('DIVISION DE PHI2 PAR LES MASSES AFIN OBTENIR LA MATRICE DYNAMIQUE D');

       for k=1:3
           for i=1:nat
           ik=(i-1)*3+k;           
             for l=1:3
                 for j=1:nat
                 jl=(j-1)*3+l;
                 D(ik,jl)=phi(ik,jl)/sqrt(masse*masse);
                 end
             end
           end
       end
%       
%**********************************************************************
%*     ECRITURE DES MATRICES
%**********************************************************************
fid1 = fopen('MatriceDyn.dat', 'wt');
       for i=1:nat
          ik=(i-1)*3;
          for j=1:nat
             jk=(j-1)*3;
             fprintf(fid1,'%5d %5d\n',i,j);
             for k=1:3
                fprintf(fid1,'%15.5f %15.5f %15.5f\n',D(ik+k,jk+1),D(ik+k,jk+2),D(ik+k,jk+3));
               % write (3,'(3(d15.8))')(d(ik+k,jk+l),l=1,3)
             end
          end
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Vp,omega2]= eig(D); % diagonalisation de la mat dyn D
omega=sqrt(abs(diag(omega2)));
F=sort(omega*unite);
for i=1:3*nat
	 for j=1:nat
	       for k=1:3
	          jk=(j-1)*3+k;
	          vt(k,i,j)=Vp(jk,i);
           end
     end
end
fid = fopen('FreqFlu.mkl', 'wt'); % Fichier data pour MOLEKEL (visualisation de freq propres)
fprintf(fid,'%s\n','$MKL');
fprintf(fid,'%s\n','$COORD');
      for i=1:nat	 
  fprintf(fid,'%10d  %15.6f  %15.6f  %15.6f\n',zat,x(i),y(i),z(i));
      end
	fprintf(fid,'%s\n','$END');
	fprintf(fid,'%s\n','$FREQ');
    
    for m=0:3:3*nat-1
 fprintf(fid,'%18s  %18s  %18s \n','X','X','X');
 fprintf(fid,'%20.3f %20.3f %20.3f \n',omega(m+1)*unite,omega(m+2)*unite,omega(m+3)*unite);
      for n=1:nat
        fprintf(fid,'%9.5f %9.5f %9.5f',vt(1,m+1,n),vt(2,m+1,n),vt(3,m+1,n));
        fprintf(fid,'%9.5f %9.5f %9.5f',vt(1,m+2,n),vt(2,m+2,n),vt(3,m+2,n));
        fprintf(fid,'%9.5f %9.5f %9.5f\n',vt(1,m+3,n),vt(2,m+3,n),vt(3,m+3,n));
      end
    end
fprintf(fid,'%s\n','$END');  
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la DOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 del=0.0004; %epsillon du Dirac
      fid=id;
      for  i=1:id
         fi=i;
         o=fi/fid*cof;
         um(i)=o;
         scq=0;
         for j=1:n
            scq=scq+1/((o-omega(j))-0.004j)/pi; % DOS
         end
        g(i)=imag(scq);
      end
   %Normalisation
   g=g./sum(g);
      for j=1:id
         um(j)=um(j)*unite;
      end
     plot(um,g)
     xlabel('Frequency (cm-1)');
     ylabel('VDOS (arbitrary unit)');