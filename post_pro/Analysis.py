import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def extract_simulation_data(filename="../data/results.txt"):
    with open(filename, "r") as file:
        lines = file.readlines()

    grid_config = {}
    grid_config["N"] = float(lines[1].split(":")[1].split(",")[0].strip())
    grid_config["L"] = float(lines[1].split(":")[2].split(",")[0].strip())
    grid_config["h"] = float(lines[1].split(":")[3].split(",")[0].strip())
    grid_config["CFL"] = float(lines[1].split(":")[4].split(",")[0].strip())
    grid_config["dt"] = float(lines[1].split(":")[5].split(",")[0].strip())
    grid_config["nt"] = float(lines[1].split(":")[6].strip())


    start_idx = lines.index("Initial condition (U points):\n") + 1
    U_initial = np.array([float(x) for x in lines[start_idx].strip().split()])


    time_steps = []
    U_steps = []
    i = start_idx + 2  
    while i < len(lines):
        time = float(lines[i].split("=")[1].strip())
        U = np.array([float(x) for x in lines[i+2].strip().split()])
        time_steps.append(time)
        U_steps.append(U)
        i += 4  
    
    return grid_config, U_initial, time_steps, U_steps

if __name__ == '__main__':
    grid_config, U_initial, time_steps, U_steps = extract_simulation_data()
    print("Grid Configuration:", grid_config)
    L=grid_config["L"]
    N=grid_config["N"]
    x=np.linspace(-L/2,3*L/2,int(N), endpoint=False)
    # Création de la figure pour l'animation
    fig, ax = plt.subplots()
    line, = ax.plot(x, U_steps[0], 'b-', lw=2)  # Courbe initiale
    ax.set_xlim(-L/2, 3*L/2)
    ax.set_ylim(min(map(min, U_steps)), max(map(max, U_steps)))
    ax.set_xlabel("x")
    ax.set_ylabel("U(x, t)")
    ax.set_title("Évolution de U en fonction de x et t")
    
    # Fonction de mise à jour de l'animation
    def update(frame):
        line.set_ydata(U_steps[frame])  # Mise à jour des valeurs U
        ax.set_title(f"Évolution de U - t = {time_steps[frame]:.4f}")
        return line,
    
    # Création de l'animation
    ani = animation.FuncAnimation(fig, update, frames=len(U_steps), interval=100, blit=False)
    ani.save("../plot/simulation.gif", writer="pillow", fps=5)

    # Affichage
    plt.show()
