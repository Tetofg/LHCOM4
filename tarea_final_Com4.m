%Para la escritura del vector---------------
clc
clear %limpiamos la ventana de comandos
pkg load symbolic
pkg load signal
%Definición del vector x[n]
n = -14:14;
x1 = 1+(n/9).^3; %vector generado

%Acondicionamiento de x para 31 posciciones
n = 1:31;
x = zeros(1,31);
x=x1;
x(30) = 0;
x(31) = 0;

%Vector y
y=cos(x)+2;

%Grafica vector [x[n]]
figure(5)
subplot(2,1,1)
stem(n, x)
title('X[n]')

%Grafica vector [y[n]]
subplot(2,1,2)
stem(n, y)
title('Y[n]')


[b,a] = CalculoCoef1(x,y);
