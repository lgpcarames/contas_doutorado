import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


teste_dados = pd.read_csv("concentric_levy.csv")


#plot inside-percent x log(distance)
title = 'Teste Power-law truncada: Eficiência x alpha\n com escala=1'
plt.title(title)

plt.ylabel("Eficiência")
plt.xlabel("alpha")
plt.plot(teste_dados['alpha'], teste_dados['eta'], 'bo-', label='Eficiência')
#plt.plot(dados_array[0]['alpha'], dados_array[0]['outside-percent'], 'r*-', label='outside-targets')
#plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

plt.show()




#plot inside-percent x log(distance)
title = 'Teste power-law truncada: proporção alvos inside e outside x alpha\n com escala=1'
plt.title(title)

plt.ylabel("Proporção")
plt.xlabel("alpha")
plt.plot(teste_dados['alpha'], teste_dados[' inside-percent'], 'bo-', label='inside-targets')
plt.plot(teste_dados['alpha'], teste_dados[' outside-percent'], 'r*-', label='outside-targets')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

plt.show()