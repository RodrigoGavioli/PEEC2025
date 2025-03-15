import os
import matplotlib.pyplot as plt
def saving_graphs(folder, name):

    i = 1

    graph = f"{name}_{i}.png"
    caminho = os.path.join(folder, graph)

    if os.path.exists(folder) == False:
        return print("this folder does not exist:", folder)
    
    while os.path.exists(caminho):
        i += 1
        graph = f"{name}_{i}.png"
        caminho = os.path.join(folder, graph)

    plt.savefig(caminho, dpi=300, bbox_inches='tight')
    print(f"Gráfico salvo em: {caminho}")

    return caminho