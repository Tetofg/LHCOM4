function [b,a]=CalculoCoef1(x,y)
clc;
%solución de sistema tipo AX=B
%Definición de matriz de coeficientes 
A=[x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   x(5)   x(4)   x(3)   x(2)   x(1)   -y(10)  -y(9)   -y(8)   -y(7)   -y(6)   -y(5)   -y(4)   -y(3)   -y(2)   -y(1);
   x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   x(5)   x(4)   x(3)   x(2)   -y(11)  -y(10)  -y(9)   -y(8)   -y(7)   -y(6)   -y(5)   -y(4)   -y(3)   -y(2);
   x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   x(5)   x(4)   x(3)   -y(12)  -y(11)  -y(10)  -y(9)   -y(8)   -y(7)   -y(6)   -y(5)   -y(4)   -y(3);
   x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   x(5)   x(4)   -y(13)  -y(12)  -y(11)  -y(10)  -y(9)   -y(8)   -y(7)   -y(6)   -y(5)   -y(4);
   x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   x(5)   -y(14)  -y(13)  -y(12)  -y(11)  -y(10)  -y(9)   -y(8)   -y(7)   -y(6)   -y(5);
   x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   x(6)   -y(15)  -y(14)  -y(13)  -y(12)  -y(11)  -y(10)  -y(9)   -y(8)   -y(7)   -y(6);
   x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   x(7)   -y(16)  -y(15)  -y(14)  -y(13)  -y(12)  -y(11)  -y(10)  -y(9)   -y(8)   -y(7);
   x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   x(8)   -y(17)  -y(16)  -y(15)  -y(14)  -y(13)  -y(12)  -y(11)  -y(10)  -y(9)   -y(8);
   x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  x(9)   -y(18)  -y(17)  -y(16)  -y(15)  -y(14)  -y(13)  -y(12)  -y(11)  -y(10)  -y(9); 
   x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  x(10)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15)  -y(14)  -y(13)  -y(12)  -y(11)  -y(10); 
   x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  x(11)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15)  -y(14)  -y(13)  -y(12)  -y(11); 
   x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  x(12)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15)  -y(14)  -y(13)  -y(12);
   x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  x(13)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15)  -y(14)  -y(13);
   x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  x(14)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15)  -y(14);
   x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  x(15)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16)  -y(15);
   x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  x(16)  -y(25)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17)  -y(16);
   x(27)  x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  x(17)  -y(26)  -y(25)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18)  -y(17); 
   x(28)  x(27)  x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  x(18)  -y(27)  -y(26)  -y(22)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19)  -y(18); 
   x(29)  x(28)  x(27)  x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  x(19)  -y(28)  -y(27)  -y(26)  -y(25)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20)  -y(19);
   x(30)  x(29)  x(28)  x(27)  x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  x(20)  -y(29)  -y(28)  -y(27)  -y(26)  -y(25)  -y(24)  -y(23)  -y(22)  -y(21)  -y(20);
   x(31)  x(30)  x(29)  x(28)  x(27)  x(26)  x(25)  x(24)  x(23)  x(22)  x(21)  -y(30)  -y(29)  -y(28)  -y(27)  -y(26)  -y(25)  -y(24)  -y(23)  -y(22)  -y(21)  ];
%definición de matriz de resultados
B=[y(11); y(12); y(13); y(14); y(15); y(16); y(17); y(18); y(19); y(20); y(21); y(22); y(23); y(24); y(25); y(26); y(27); y(28); y(29); y(30); y(31)];

%Solución de los coeficientes C
%C=inv(A)*B;

C=A\B;

fprintf('Los Coeficientes quedan de la siguiente forma:\n b_0=  %f' ,C(1,1));
fprintf(' \n b_1=  %f' ,C(2,1));
fprintf(' \n b_2=  %f' ,C(3,1));
fprintf(' \n b_3=  %f' ,C(4,1));
fprintf(' \n b_4=  %f' ,C(5,1));
fprintf(' \n b_5=  %f' ,C(6,1));
fprintf(' \n b_6=  %f' ,C(7,1));
fprintf(' \n b_7=  %f' ,C(8,1));
fprintf(' \n b_8=  %f' ,C(9,1));
fprintf(' \n b_9=  %f' ,C(10,1));
fprintf(' \n b_10= %f' ,C(11,1));
fprintf(' \n a_1=  %f' ,C(12,1));
fprintf(' \n b_2=  %f' ,C(13,1));
fprintf(' \n a_3=  %f' ,C(14,1));
fprintf(' \n a_4=  %f' ,C(15,1));
fprintf(' \n a_5=  %f' ,C(16,1));
fprintf(' \n a_6=  %f' ,C(17,1));
fprintf(' \n a_7=  %f' ,C(18,1));
fprintf(' \n a_8=  %f' ,C(19,1));
fprintf(' \n a_9=  %f' ,C(20,1));
fprintf(' \n a_10= %f \n' ,C(21,1));

%Comparación evaluando coeficientes --------------------------------------------------------------
        fprintf('Comparación de valores para valores entre [20] y [30] , \n Primera columna: Función encontrada  Segunda Fila: Función original')
        M=A*C;
        for i=1:21
        m(i)=y(i+10);
        end
        l=transpose(m);
        Comparacion=[M l]

%Encontrando H(Z)
        c=transpose(C);
       %Ordenando coeficientes poara el denominador y el numerador,
     
      
        for i=1:11
            b(i)=c(1,i);
        end
        a(1)=1;
        for j=12:21
           a(j-10)=c(j);
        end
       
        syms z;
            Z=[1 ;z^(-1);z^(-2);z^(-3);z^(-4);z^(-5);z^(-6);z^(-7);z^(-8);z^(-9);z^(-10)];
        h1=(b*Z);
        h2=(a*Z);
        fprintf('La función H (sin fracción parciales) para el sistema es:') 
        H=h1/h2   
        
 %con fracciones parciales
 
    fprintf('\n Coeficientes de fracciones parciales: \n')
    [r,p,k]=residuez(b,a)
    r=transpose(r);
    p=transpose(p);
    
    fprintf('La función H (con fracción parciales) para el sistema es:') 
    h=sum(r./(1-p*z^(-1)))+sum(k*[1 ; z^(-1)])
     for i=1:31
      z=i;
      HH(i)=sum(r./(1-p*z^(-1)))+sum(k*[1 ; z^(-1)]);
     end
     HH;
     
     
     %Reescribiendo fel vector a y b para el total de las posiciones
      %  for i=12:31
      %      b(i)=0;
      %  end
      %  for j=11:31
      %     a(j)=0;
      %  end
       
%ANALISIS ESPECTRAL
            n=0:30;
            f=1/10;%Frecuencia de 0.1Hz
            N1=30; %numero de muestras
            N2=120;
        %Analisis Espectral Para X
            %Transformadas:
            X1=abs(fft(x,N1));
            X2=abs(fft(x,N2));

            %Rango normalizado para transformadas:
            F1x=[(0:N1-1)/N1];
            F2x=[(0:N2-1)/N2];

        %Analisis Espectral Para Y
            %Transformadas:
            Y1=abs(fft(y,N1));
            Y2=abs(fft(y,N2));

            %Rango normalizado para transformadas:
            F1y=[(0:N1-1)/N1];
            F2y=[(0:N2-1)/N2];

 %Grafica de funciones
 
     %Funciones X y Y
            figure (1)
            subplot(4,2,1)
            stem(n,x)
            title('X[n]')

            subplot(4,2,2)
            stem(n,y)
            title('Y[n]')
    
     %Gráficas del espectro:   
            subplot(4,2,3)
            stem(F1x,X1,'.')
            title('Espectro en frecuencia para X[n] con N=30')

            subplot(4,2,4)
            stem(F2x,X2,'.')
            title('Espectro en frecuencia para X[n] con N=120')

            subplot(4,2,5)
            stem(F1y,Y1,'.')
            title('Espectro en frecuencia para Y[n] con N=30')

            subplot(4,2,6)
            stem(F2y,X2,'.')
            title('Espectro en frecuencia paraY[n] con N=120')

    %Parte Real e imaginaria
            subplot(4,2,7)
            stem(n,real(HH))
            title('Parte real H[Z]')

            subplot(4,2,8)
            stem(n,imag(HH))
            title('Parte imaginaria H[Z]')
    
    %respuesta en amplitud y frecuencia
            figure(2)
            freqz(b,a)
            title('Respuesta en Amplitud y Fase de H(z)')
    
    %Respuesta al impulso
             %l=1:21;
             figure(3)
             y= filter(b,a,n)
             plot(n,y)
             title('Respuesta al impulso y al escalón unitario')
        
    %Plolos  y Ceros
            figure(4)
            zplane(b,a)
            title('Diagrama de Polos y Ceros')
            axis([-3.5 3.5 -3.5 3.5]) % establece los límites del eje x e y

