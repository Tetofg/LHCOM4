# Para la escritura del vector
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from sympy import Symbol

def CalculoCoef1(x, y):

    A = np.array([[x[10], x[9],  x[8],  x[7],  x[6], x[5], x[4], x[3], x[2], x[1], x[0],  -y[9],  -y[8], -y[7],  -y[6],  -y[5],  -y[4],  -y[3],  -y[2],  -y[1],  y[0]],
                  [x[11], x[10], x[9],  x[8],  x[7], x[6], x[5], x[4], x[3], x[2], x[1],  -y[10], -y[9], -y[8],  -y[7],  -y[6],  -y[5],  -y[4],  -y[3],  -y[2],  y[1]],
                  [x[12], x[11], x[10], x[9],  x[8], x[7], x[6], x[5], x[4], x[3], x[2],  -y[11], -y[10],-y[9],  -y[8],  -y[7],  -y[6],  -y[5],  -y[4],  -y[3],  y[2]],
                  [x[13], x[12], x[11], x[10], x[9], x[8], x[7], x[6], x[5], x[4], x[3],  -y[12], -y[11],-y[10], -y[9],  -y[8],  -y[7],  -y[6],  -y[5],  -y[4],  y[3]],
                  [x[14], x[13], x[12], x[11], x[10],x[9], x[8], x[7], x[6], x[5], x[4],  -y[13], -y[12],-y[11], -y[10], -y[9],  -y[8],  -y[7],  -y[6],  -y[5],  y[4]],
                  [x[15], x[14], x[13], x[12], x[11],x[10],x[9], x[8], x[7], x[6], x[5],  -y[14], -y[13],-y[12], -y[11], -y[10], -y[9],  -y[8],  -y[7],  -y[6],  y[5]],
                  [x[16], x[15], x[14], x[13], x[12],x[11],x[10],x[9], x[8], x[7], x[6],  -y[15], -y[14],-y[13], -y[12], -y[11], -y[10], -y[9],  -y[8],  -y[7],  y[6]],
                  [x[17], x[16], x[15], x[14], x[13],x[12],x[11],x[10],x[9], x[8], x[7],  -y[16], -y[15],-y[14], -y[13], -y[12], -y[11], -y[10], -y[9],  -y[8],  y[7]],
                  [x[18], x[17], x[16], x[15], x[14],x[13],x[12],x[11],x[10],x[9], x[8],  -y[17], -y[16],-y[15], -y[14], -y[13], -y[12], -y[11], -y[10], -y[9],  y[8]],
                  [x[19], x[18], x[17], x[16], x[15],x[14],x[13],x[12],x[11],x[10],x[9],  -y[18], -y[17],-y[16], -y[15], -y[14], -y[13], -y[12], -y[11], -y[10], y[9]],
                  [x[20], x[19], x[18], x[17], x[16],x[15],x[14],x[13],x[12],x[11],x[10], -y[19], -y[18],-y[17], -y[16], -y[15], -y[14], -y[13], -y[12], -y[11], y[10]],
                  [x[21], x[20], x[19], x[18], x[17],x[16],x[15],x[14],x[13],x[12],x[11], -y[20], -y[19],-y[18], -y[17], -y[16], -y[15], -y[14], -y[13], -y[12], y[11]],
                  [x[22], x[21], x[20], x[19], x[18],x[17],x[16],x[15],x[14],x[13],x[12], -y[21], -y[20],-y[19], -y[18], -y[17], -y[16], -y[15], -y[14], -y[13], y[12]],
                  [x[23], x[22], x[21], x[20], x[19],x[18],x[17],x[16],x[15],x[14],x[13], -y[22], -y[21],-y[20], -y[19], -y[18], -y[17], -y[16], -y[15], -y[14], y[13]],
                  [x[24], x[23], x[22], x[21], x[20],x[19],x[18],x[17],x[16],x[15],x[14], -y[23], -y[22],-y[21], -y[20], -y[19], -y[18], -y[17], -y[16], -y[15], y[14]],
                  [x[25], x[24], x[23], x[22], x[21],x[20],x[19],x[18],x[17],x[16],x[15], -y[24], -y[23],-y[22], -y[21], -y[20], -y[19], -y[18], -y[17], -y[16], y[15]],
                  [x[26], x[25], x[24], x[23], x[22],x[21],x[20],x[19],x[18],x[17],x[16], -y[25], -y[24],-y[23], -y[22], -y[21], -y[20], -y[19], -y[18], -y[17], y[16]],
                  [x[27], x[26], x[25], x[24], x[23],x[22],x[21],x[20],x[19],x[18],x[17], -y[26], -y[25],-y[24], -y[23], -y[22], -y[21], -y[20], -y[19], -y[18], y[17]],
                  [x[28], x[27], x[26], x[25], x[24],x[23],x[22],x[21],x[20],x[19],x[18], -y[27], -y[26],-y[25], -y[24], -y[23], -y[22], -y[21], -y[20], -y[19], y[18]],
                  [x[29], x[28], x[27], x[26], x[25],x[24],x[23],x[22],x[21],x[20],x[19], -y[28], -y[27],-y[26], -y[25], -y[24], -y[23], -y[22], -y[21], -y[20], y[19]],
                  [x[30], x[29], x[28], x[27], x[26],x[25],x[24],x[23],x[22],x[21],x[20], -y[29], -y[28],-y[27], -y[26], -y[25], -y[24], -y[23], -y[22], -y[21], y[20]]])

    B = [y[10], y[11], y[12], y[13], y[14], y[15], y[16], y[17], y[18], y[19], y[20], y[21], y[22], y[23], y[24], y[25], y[26], y[27], y[28], y[29], y[30]]
    # Solución de los coeficientes C
    # C=inv(A)*B;
    C = np.linalg.inv(A) @ B

    print('Los Coeficientes quedan de la siguiente forma:')
    print(f'b_0 = {C[0]}')
    print(f'b_1 = {C[1]}')
    print(f'b_2 = {C[2]}')
    print(f'b_3 = {C[3]}')
    print(f'b_4 = {C[4]}')
    print(f'b_5 = {C[5]}')
    print(f'b_6 = {C[6]}')
    print(f'b_7 = {C[7]}')
    print(f'b_8 = {C[8]}')
    print(f'b_9 = {C[9]}')
    print(f'b_10 = {C[10]}')
    print(f'a_1 = {C[11]}')
    print(f'a_2 = {C[12]}')
    print(f'a_3 = {C[13]}')
    print(f'a_4 = {C[14]}')
    print(f'a_5 = {C[15]}')
    print(f'a_6 = {C[16]}')
    print(f'a_7 = {C[17]}')
    print(f'a_8 = {C[18]}')
    print(f'a_9 = {C[19]}')
    print(f'a_10 = {C[20]}')


    # Comparación evaluando coeficientes --------------------------------------------------------------
    print('Comparación de valores para valores entre [20] y [30] , \n Primera columna: Función encontrada  Segunda Fila: Función original')
    M = A @ C
    Comparacion = M
    print(Comparacion)

    c = C.T
    b = np.zeros((11, 1))
    a = np.zeros((11, 1))

    # Ordenando coeficientes para el denominador y el numerador
    for i in range(11):
        b[i] = c[i]

    a[0] = 1
    for j in range(12, 22):
        a[j-11] = c[j-1]

    z = Symbol('z')
    Z = np.array([[1], [z**(-1)], [z**(-2)], [z**(-3)], [z**(-4)], [z**(-5)], [z**(-6)], [z**(-7)], [z**(-8)], [z**(-9)], [z**(-10)]])

    h1 = np.dot(b.T, Z)[0]
    h2 = np.dot(a.T, Z)[0]

    print('La función H (sin fracción parciales) para el sistema es:')
    H = h1 / h2

    # Coeficientes de fracciones parciales
    print('\n Coeficientes de fracciones parciales: \n')
    r, p, k = signal.residuez(np.squeeze(b), np.squeeze(a))
    r = np.transpose(r)
    p = np.transpose(p)

    # Función H (con fracciones parciales) para el sistema
    print('La función H (con fracción parciales) para el sistema es:')
    h = np.sum(r / (1 - p * z**(-1.0))) + np.sum(k * np.array([1, z**(-1.0)]))
    print(h)

    # Valores de H para z de 1 a 31
    HH = np.zeros(31)
    for i in range(31):
        z = i + 1
        HH[i] = np.sum(r / (1 - p * z**(-1.0))) + np.sum(k * np.array([1, z**(-1.0 )]))
    print(HH)

    n = np.arange(0, 31)
    f = 1/10 # Frecuencia de 0.1Hz
    N1 = 30 # número de muestras
    N2 = 120

    # Análisis Espectral Para X
    # Transformadas
    X1 = np.abs(np.fft.fft(x, N1))
    X2 = np.abs(np.fft.fft(x, N2))

    # Rango normalizado para transformadas
    F1x = np.arange(0, N1)/N1
    F2x = np.arange(0, N2)/N2

    # Análisis Espectral Para Y
    # Transformadas
    Y1 = np.abs(np.fft.fft(y, N1))
    Y2 = np.abs(np.fft.fft(y, N2))

    # Rango normalizado para transformadas
    F1y = np.arange(0, N1)/N1
    F2y = np.arange(0, N2)/N2

    # Grafica de funciones

    fig, axs = plt.subplots(4, 2)

    axs[0, 0].stem(n, x)
    axs[0, 0].set_title('X[n]')

    axs[0, 1].stem(n, y)
    axs[0, 1].set_title('Y[n]')

    # Gráficas del espectro:
    axs[1, 0].stem(F1x, X1, '.')
    axs[1, 0].set_title('Espectro en frecuencia para X[n] con N=30')

    axs[2, 0].stem(F2x, X2, '.')
    axs[2, 0].set_title('Espectro en frecuencia para X[n] con N=120')

    axs[1, 1].stem(F1y, Y1, '.')
    axs[1, 1].set_title('Espectro en frecuencia para Y[n] con N=30')

    axs[2, 1].stem(F2y, X2, '.')
    axs[2, 1].set_title('Espectro en frecuencia paraY[n] con N=120')

    # Ocultar ejes no utilizados
    axs[3, 0].axis('off')
    axs[3, 1].axis('off')

    plt.show()


    # Parte Real e Imaginaria
    fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    axs[0].stem(n, np.real(HH))
    axs[0].set_title('Parte real H[Z]')
    axs[1].stem(n, np.imag(HH))
    axs[1].set_title('Parte imaginaria H[Z]')

    # Respuesta en amplitud y frecuencia
    fig = plt.figure(figsize=(8, 6))
    w, h = signal.freqz(b, a)
    plt.plot(w, abs(h))
    plt.title('Respuesta en Amplitud y Fase de H(z)')
    plt.show()
    # Respuesta al impulso
    fig = plt.figure(figsize=(8, 6))
    y = signal.lfilter(np.squeeze(b), np.squeeze(a), n)
    plt.plot(n, y)
    plt.title('Respuesta al impulso y al escalón unitario')
    plt.show()
    # Polos y Ceros
    fig, ax = plt.subplots(figsize=(8, 8))
    z, p, k = signal.tf2zpk(np.squeeze(b), np.squeeze(a))
    plt.scatter(np.real(z), np.imag(z), marker='o', facecolors='none', edgecolors='b')
    plt.scatter(np.real(p), np.imag(p), marker='x', color='r')
    plt.legend(['Zeros', 'Polos'], loc='best')
    plt.title('Diagrama de Polos y Ceros')
    plt.axis([-3.5, 3.5, -3.5, 3.5])
    plt.show()




# Definición del vector x[n]
n = np.arange(-14, 14)
x1 = 1 + (n / 9) ** 3  # vector generado

# Acondicionamiento de x para 31 posiciones
n = np.arange(1, 32)
x = np.zeros(31)
x[:len(x1)] = x1
x[29:] = 0

# Vector y
y = np.cos(x) + 2

# Grafica vector [x[n]]
plt.figure(5)
plt.subplot(2, 1, 1)
plt.stem(n, x)
plt.title('X[n]')

# Grafica vector [y[n]]
plt.subplot(2, 1, 2)
plt.stem(n, y)
plt.title('Y[n]')
plt.show()

CalculoCoef1(x,y)