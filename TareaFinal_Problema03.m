%Para la escritura del vectorize
clc
clear %limpiamos la ventana de comandos
%Establecemos en Octave
pkg load signal;
pkg load symbolic;
pkg load optim;
%Se lee el archivo de audio x, almacenado en la actual carpeta de trabajo
[x, fs] = audioread('x2_U017.wav'); 

lx=length(x); % Cálculo del tamaño del vector de audio

%Se lee el archivo de audio y, almacenado en la actual carpeta de trabajo
[y, fs] = audioread('y2_U017.wav');

ly=length(y); % Cálculo del tamaño del vector de audio

%Grafica de los vectores
n=1:lx;
figure(5);
subplot(2,1,1);
stem(n,x);
title('X[n]');
subplot(2,1,2);
stem(n,y);
title('Y[n]');

%Para la solucion del problema
[b,a] = CalculoCoef3(x,y,lx);
