import numpy as np
from random import random
import matplotlib.pyplot as plt


#Entrada dos valores de alguns parâmetros
L = float(input("Entre com o valor do raio do anel externo: \n"))  # raio externo do anel
R0 = float(input("Entre com o valor do fator de escala: \n")) # fator de escala da função de distribuição
LC = float(input("Entre com o valor da posição de partida do caminhante: \n")) # Posição de partida do caminhante

#Definindo as constantes

PI = 3.14159265358979323846
RV = 1
TOTALDISTANCE = 100000000
LARGESTFLIGHT = 1000*L
travel = 0
# Definindo as funções

def SQR(x):
    if x != 0:
        x = x*x
        return x
    else: 
        return 0

def rng_levy48(alpha, rr):
    mu = alpha - 1
    mu1 = mu - 1
    xmu = 1/mu
    xmu1 = xmu-1
    phi = (random() - 0.5)*PI
    ee = -np.log(random())
    return rr*np.sin(mu*phi)/pow(np.cos(phi),xmu)*pow(np.cos(phi*mu1)/ee, xmu1)


    
