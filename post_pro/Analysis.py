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
    grid_config["nt"] = float(lines[1].split(":")[6].split(",")[0].strip())
    grid_config["c"] = float(lines[1].split(":")[7].split(",")[0].strip())
    grid_config["sigma"] = float(lines[1].split(":")[8].strip())


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

def exact_gaussian(t,x,c,sigma):
    u = np.exp(-((x-t)**2)/(sigma**2))
    return u;

def diagnosis(U,h,A,sigma,t,x,c):
    I=np.zeros(len(U))
    E=np.zeros(len(U))
    R=np.zeros(len(U))
    for i in range(len(U)):
        for j in range(len(U[i])):
            I[i]+= U[i][j]
            E[i]+= U[i][j]**2
            R[i]+= (U[i][j] - exact_gaussian(t[i],x[j],c,sigma))**2
    I= I*(h/(sigma*A))
    E = E*(h/(sigma*(A**2)))
    R= R*(h/(sigma*(A**2)))
    return I,E,R
    
    
def plot_global_diagnostics(U, h, A, sigma, t, x, c, L):
    I, E, R = diagnosis(U, h, A, sigma, t, x, c)
    ctL = (c * t) / L

    plt.figure(figsize=(8, 6))

    plt.plot(ctL, I, label=r'$I_n^h$', linestyle='-', color='blue')
    plt.plot(ctL, E, label=r'$E_n^h$', linestyle='--', color='black')
    plt.plot(ctL, R, label=r'$R_n^h$', linestyle='-.', color="green")
    
    # Labels and title
    plt.xlabel(r'$ct / L$', fontsize=14)
    plt.ylabel('Diagnostics', fontsize=14)
    plt.title('Evolution of Global Diagnostics', fontsize=16)
    plt.legend(fontsize=12)
    
    # Grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Show the figure
    plt.show()
    

def anim(x,U_steps,L):
    ###### Animation ######
    fig, ax = plt.subplots()
    line, = ax.plot(x, U_steps[0], 'b-', lw=2)  # Courbe initiale
    ax.set_xlim(-L/2, 3*L/2)
    ax.set_ylim(min(map(min, U_steps)), max(map(max, U_steps)))
    ax.set_xlabel("x")
    ax.set_ylabel("U(x, t)")
    ax.set_title("Évolution de U en fonction de x et t")
    def update(frame):
        line.set_ydata(U_steps[frame])  # Mise à jour des valeurs U
        ax.set_title(f"Évolution de U - t = {time_steps[frame]:.4f}")
        return line,   
    ani = animation.FuncAnimation(fig, update, frames=len(U_steps), interval=100, blit=False)
    ani.save("../plot/simulation.gif", writer="pillow", fps=5)
    ######################## 
    




if __name__ == '__main__':
    grid_config, U_initial, time_steps, U_steps = extract_simulation_data()
    print("Grid Configuration:", grid_config)
    L=grid_config["L"]
    N=grid_config["N"]
    h=grid_config["h"]
    c=grid_config["c"]
    sigma=grid_config["sigma"]
    x=np.linspace(-L/2,L/2,int(N), endpoint=False)
    
    anim(x,U_steps,L)
    
    
    U = np.concatenate(([U_initial], U_steps))
    t = np.concatenate(([0], time_steps))
    A=1 ##attention si on modif le code
    plot_global_diagnostics(U, h, A, sigma, t, x, c, L)
    
    
