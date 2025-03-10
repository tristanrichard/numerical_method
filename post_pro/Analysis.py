import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def extract_simulation_data(filename):
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
        
    U_steps = np.concatenate(([U_initial], U_steps))
    time_steps = np.concatenate(([0], time_steps))
    
    return grid_config, time_steps, U_steps

def exact_gaussian(t,x,c,sigma,L):
    x_shifted = x - c * t
    x_periodic = (x_shifted + L/2) % L - L/2  
    u = np.exp(-((x_periodic) ** 2) / (sigma ** 2))
    
    return u

def diagnosis(U,h,A,sigma,t,x,c,L):
    I=np.zeros(len(U))
    E=np.zeros(len(U))
    R=np.zeros(len(U))
    for i in range(len(U)):
        for j in range(len(U[i])):
            I[i]+= U[i][j]
            E[i]+= U[i][j]**2
            R[i]+= (U[i][j] - exact_gaussian(t[i],x[j],c,sigma,L))**2
    I= I*(h/(sigma*A))
    E = E*(h/(sigma*(A**2)))
    R= R*(h/(sigma*(A**2)))
    return I,E,R
    
    
def plot_global_diagnostics(U, h, A, sigma, t, c, L,N,config):
    x = np.linspace(-L/2, L/2,int(N), endpoint=False)
    I, E, R = diagnosis(U, h, A, sigma, t, x, c,L)
    ctL = (c * t) / L

    plt.figure(figsize=(8, 6))

    plt.plot(ctL, I, label=r'$I_n^h$', linestyle='-', color='blue')
    plt.plot(ctL, E, label=r'$E_n^h$', linestyle='--', color='black')
    plt.plot(ctL, R, label=r'$R_n^h$', linestyle='-.', color="green")
    
    # Labels and title
    plt.xlabel(r'$ct / L$', fontsize=14)
    plt.ylabel('Diagnostics', fontsize=14)
    plt.title(f'Evolution of Global Diagnostics {config}', fontsize=16)
    plt.legend(fontsize=12)
    
    # Grid and styling
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Show the figure
    plt.show()
    

def anim(x,U_steps,L,time_steps):
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
    
def plot_versus(U,t,grid_config,scheme):
    L=grid_config["L"]
    N=int(grid_config["N"])
    dt=grid_config["dt"]
    c=grid_config["c"]
    sigma= grid_config["sigma"]
    
    x = np.linspace(-L/2, 3*L/2,2*int(N), endpoint=False)
    index_05 = int(0.5/dt) +1
    index_1 = int(1/dt) +1
    U_05 = np.concatenate((U[index_05],U[index_05]),axis=0)
    U_1 = np.concatenate((U[index_1],U[index_1]),axis=0)
    
    x1 = np.linspace(-L/2, L/2,1000, endpoint=False)
    x2 = np.linspace(L/2, 3*L/2,1000, endpoint=False)
    U_05_real = exact_gaussian(t[index_05-1], x1, c, sigma,L)
    U_1_real = exact_gaussian(t[index_1-1], x1, c, sigma,L)
    U_05_real = np.concatenate((U_05_real,U_05_real),axis=0)
    U_1_real = np.concatenate((U_1_real,U_1_real),axis=0)

    x1 = np.concatenate((x1,x2),axis=0)
    
    # Plotting the results
    plt.figure(figsize=(10, 8))
    plt.plot(x, U_05, label=r'Numerical solution at $\frac{ct}{L}=0.5$', color='royalblue', linestyle='--', linewidth=2)
    plt.plot(x1, U_05_real, label=r'Exact solution at $\frac{ct}{L}=0.5$', color='royalblue', linestyle='-', linewidth=2)
    plt.plot(x, U_1, label=r'Numerical solution at $\frac{ct}{L}=1$', color='orangered', linestyle='--', linewidth=2)
    plt.plot(x1, U_1_real, label=r'Exact solution at $\frac{ct}{L}=1$', color='orangered', linestyle='-', linewidth=2)


    plt.xlabel(r'$\frac{x}{L}$', fontsize=25)
    plt.ylabel(r'$\frac{u(x,t)}{U}$', fontsize=25)
    #plt.legend(fontsize=12,loc='lower right', labelspacing=0.1)
    #plt.ylim([-0.7,1])
    plt.grid(True)
    plt.tick_params(axis='both', labelsize=15)
    plt.savefig(f'../plot/versus_{scheme}_N{N}.png', dpi=300, bbox_inches='tight')
    plt.show()




if __name__ == '__main__':
    #### E2 32####
    grid_config_E2_N32, time_steps_E2_N32, U_steps_E2_N32 = extract_simulation_data("../data/results_E2_N32.txt")
    # anim(x, U_steps, L, time_steps)
    A = 1  # Attention si on modifie le code   
    #plot_global_diagnostics(U_steps_E2_N32, grid_config_E2_N32["h"], A, grid_config_E2_N32["sigma"], time_steps_E2_N32, grid_config_E2_N32["c"], grid_config_E2_N32["L"], grid_config_E2_N32["N"], "E2_N32")
    plot_versus(U_steps_E2_N32, time_steps_E2_N32, grid_config_E2_N32, "E2")
    
    #### E4 32####
    grid_config_E4_N32, time_steps_E4_N32, U_steps_E4_N32 = extract_simulation_data("../data/results_E4_N32.txt")
    A = 1   
    # plot_global_diagnostics(U_steps_E4_N32, grid_config_E4_N32["h"], A, grid_config_E4_N32["sigma"], time_steps_E4_N32, grid_config_E4_N32["c"], grid_config_E4_N32["L"], grid_config_E4_N32["N"], "E4_N32")
    plot_versus(U_steps_E4_N32, time_steps_E4_N32, grid_config_E4_N32, "E4")
    
    #### I4 32####
    grid_config_I4_N32, time_steps_I4_N32, U_steps_I4_N32 = extract_simulation_data("../data/results_I4_N32.txt")
    A = 1    
    # plot_global_diagnostics(U_steps_I4_N32, grid_config_I4_N32["h"], A, grid_config_I4_N32["sigma"], time_steps_I4_N32, grid_config_I4_N32["c"], grid_config_I4_N32["L"], grid_config_I4_N32["N"], "I4_N32")
    plot_versus(U_steps_I4_N32, time_steps_I4_N32, grid_config_I4_N32, "I4")
    
    #### ED 32####
    grid_config_ED_N32, time_steps_ED_N32, U_steps_ED_N32 = extract_simulation_data("../data/results_ED_N32.txt")
    A = 1     
    # plot_global_diagnostics(U_steps_ED_N32, grid_config_ED_N32["h"], A, grid_config_ED_N32["sigma"], time_steps_ED_N32, grid_config_ED_N32["c"], grid_config_ED_N32["L"], grid_config_ED_N32["N"], "ED_N32")
    plot_versus(U_steps_ED_N32, time_steps_ED_N32, grid_config_ED_N32, "ED")

    

    
    
