import numpy as np
import matplotlib.pylab as plt
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score

from constantes import ESPECIES


def leer_archivo_fasta(ruta):
    secuencia = []
    with open(ruta, 'rb') as bufer:
        for linea in bufer:
            secuencia.append(linea.decode('utf8').strip().lower())
    secuencia.pop(0)
    secuencia = ''.join(secuencia)
    return secuencia

def entropia_shannon(observaciones, base=4):
    """
    Función que calcula la entropía
    de un conjunto de observaciones
    H(X) = -1/log(base) * (suma ( p(xi) *log(p(xi)) ) )
    para cada xi, de los valores que aparecen en `observaciones`
    Parameters:
        - observaciones de valores a ser contados
        - base=2  base del logaritmo, que determina la unidad en que se mide H
    Return:
        - float H(X) donde X es la v.a. aleatoria que fue observada
    """
    S = 0
    N = len(observaciones)
    f = {}  # Diccionario de frecuencias
    for x in observaciones:
        f[x] = f.get(x, 0) + 1      # calculamos la frecuencia del símbolo x
    F = np.array([f[x] / N for x in f])
    logF = np.log(F)
    S = -F.dot(logF) * (1 / np.log(base))
    return S

def crea_grafica(especie, ruta, secuencia, ventana=20):
    limite_final = len(secuencia) - (len(secuencia) % 3)
    resultados = []
    for i in range(0, limite_final, ventana):
        buffer = secuencia[i:i+ventana]
        resultados.append(entropia_shannon(buffer))
    promedio = np.mean(resultados)
    sigma = np.std(resultados)
    plt.figure(figsize=(18,6), dpi=100)
    plt.plot(resultados)
    plt.axhline(y=promedio, color='black', linestyle='--', linewidth=1, label=f'Promedio = {promedio:.3f}')
    plt.axhline(y=promedio+sigma, color='green', linestyle='--', linewidth=1, label=f'Sigma = {promedio+sigma:.3f}')
    plt.axhline(y=promedio-sigma, color='green', linestyle='--', linewidth=1, label=f'Sigma = {promedio-sigma:.3f}')
    plt.xlabel('posicion [bp]')
    plt.ylabel('Entropia [bits]')
    plt.title(f'Especie: {especie}')
    plt.legend(loc='lower right')
    plt.savefig(ruta+especie+".png")

def primera_comparacion():
    especies_HBA = {}
    especies_HBB = {}
    for sp in ESPECIES:
        especies_HBA[sp] = leer_archivo_fasta("./data/HBA/"+sp+".fasta")
        especies_HBB[sp] = leer_archivo_fasta("./data/HBB/"+sp+".fasta")
    # ####################################################
    for nombre, secuencia in especies_HBA.items():
        crea_grafica(nombre,"./graficas/Entropia_Shannon/HBA/", secuencia)
    for nombre, secuencia in especies_HBB.items():
        crea_grafica(nombre,"./graficas/Entropia_Shannon/HBB/", secuencia)

def crea_matriz_im(especies, tipo):
    nombres = []
    secuencias = []
    for nombre, secuencia in especies.items():
        nombres.append(nombre)
        secuencias.append(secuencia)
    N = len(nombres)
    M = np.zeros( (N,N) )
    i = 0
    while i < N:
        j = 0
        while j <= i:
            X = np.array(list(secuencias[i]))
            Y = np.array(list(secuencias[j]))
            # M[i,j] = mutual_info_score(X,Y)
            M[i,j] = normalized_mutual_info_score(X,Y)
            j += 1
        print(i)
        i += 1
    plt.figure(figsize=(10, 8))
    plt.imshow(M, cmap='hot')
    plt.colorbar()
    plt.xticks(range(len(nombres)), nombres, rotation=45, ha='right')
    plt.yticks(range(len(nombres)), nombres)
    plt.grid(color='white')
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    # plt.title(f'Información Mutua: {tipo}', pad=5)
    plt.title(f'Información Mutua: {tipo} (Normalizada)', pad=5)
    # plt.savefig('./graficas/IM/informacionMutual'+tipo+'.png', dpi=300, bbox_inches='tight')
    plt.savefig('./graficas/IM/informacionMutual'+tipo+'Normalizada.png', dpi=300, bbox_inches='tight')
    plt.close()

def segunda_comparacion():
    especies_HBA = {}
    especies_HBB = {}
    for sp in ESPECIES:
        especies_HBA[sp] = leer_archivo_fasta("./data_alineada/HBA/"+sp+".fasta")
        especies_HBB[sp] = leer_archivo_fasta("./data_alineada/HBB/"+sp+".fasta")
    crea_matriz_im(especies_HBA, "HBA")
    crea_matriz_im(especies_HBB, "HBB")

def im_pares(especie_uno, especie_dos, tipo):
    nombre_uno = especie_uno[0]
    nombre_dos = especie_dos[0]
    secuencia_uno = np.array(list(especie_uno[1]))
    secuencia_dos = np.array(list(especie_dos[1]))
    resultados = []
    for i in range(201):
        resultados.append(normalized_mutual_info_score(secuencia_uno[i:i+10], secuencia_dos[i:i+10]))
    plt.figure(figsize=(18,6), dpi=100)
    plt.plot(resultados, 'o-')
    plt.legend()
    plt.xlabel('posicion [bp]')
    plt.ylabel('Información mutua normalizada')
    plt.title(tipo+': Informacion mutua: '+nombre_uno+' vs '+nombre_dos)
    plt.savefig(f'./graficas/Mutua/{tipo}/{nombre_uno}vs{nombre_dos}.png')

def pares_graficas(especies, tipo):
    nombres = []
    secuencias = []
    for nombre, secuencia in especies.items():
        nombres.append(nombre)
        secuencias.append(secuencia)
    N = len(nombres)
    M = np.zeros( (N,N) )
    i = 0
    while i < N:
        j = 0
        while j <= i:
            if i != j:
                secuencia1 = [nombres[i], secuencias[i]]
                secuencia2 = [nombres[j], secuencias[j]]
                im_pares(secuencia1, secuencia2, tipo)
            j += 1
        print(i)
        i += 1

def tercera_comparacion():
    especies_HBA = {}
    especies_HBB = {}
    for sp in ESPECIES:
        especies_HBA[sp] = leer_archivo_fasta("./data/HBA/"+sp+".fasta")
        especies_HBB[sp] = leer_archivo_fasta("./data/HBB/"+sp+".fasta")
    pares_graficas(especies_HBA, 'HBA')
    pares_graficas(especies_HBB, 'HBB')

def main():
    # primera_comparacion()
    # segunda_comparacion()
    # tercera_comparacion()

if __name__ == "__main__":
    main()