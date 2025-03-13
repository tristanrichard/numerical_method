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
    
    
def plot_global_diagnostics(U_32, U_64, U_128, t_32, t_64, t_128, A, config_32, config_64, config_128, scheme):
    # Calcul des coordonnées spatiales
    x_32 = np.linspace(-config_32["L"]/2, config_32["L"]/2, int(config_32["N"]), endpoint=False)
    x_64 = np.linspace(-config_64["L"]/2, config_64["L"]/2, int(config_64["N"]), endpoint=False)
    x_128 = np.linspace(-config_128["L"]/2, config_128["L"]/2, int(config_128["N"]), endpoint=False)

    # Calcul des diagnostics
    I_32, E_32, R_32 = diagnosis(U_32, config_32["h"], A, config_32["sigma"], t_32, x_32, config_32["c"], config_32["L"])
    I_64, E_64, R_64 = diagnosis(U_64, config_64["h"], A, config_64["sigma"], t_64, x_64, config_64["c"], config_64["L"])
    I_128, E_128, R_128 = diagnosis(U_128, config_128["h"], A, config_128["sigma"], t_128, x_128, config_128["c"], config_128["L"])

    # Normalisation du temps
    ctL_32 = (config_32["c"] * t_32) / config_32["L"]
    ctL_64 = (config_64["c"] * t_64) / config_64["L"]
    ctL_128 = (config_128["c"] * t_128) / config_128["L"]

    # Création de la figure avec 3 sous-graphiques
    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    # Liste des diagnostics
    diagnostics = [("I", I_32, I_64, I_128, 'black', [1.7, 1.8], "lower right"),
                   ("E", E_32, E_64, E_128, 'blue', [0.8,1.3], "lower left"),
                   ("R", R_32, R_64, R_128, 'brown', [0, 0.9], "upper left")]

    for i, (diag_name, diag_32, diag_64, diag_128, color, ylim, loc) in enumerate(diagnostics):
        ax[i].plot(ctL_32, diag_32, label=r'N = 32pt', linestyle='solid', color=color)
        ax[i].plot(ctL_64, diag_64, label=r'N = 64pt', linestyle='dashdot', color=color)
        ax[i].plot(ctL_128, diag_128, label=r'N = 128pt', linestyle='dotted', color=color)

        ax[i].set_ylabel(rf'${diag_name}^{{n}}_{{h}}$', fontsize=20)
        ax[i].legend(fontsize=13, loc=loc)
        ax[i].set_xlim([0, 1])
        ax[i].set_ylim(ylim)
        ax[i].grid(True)
        ax[i].tick_params(axis='both', labelsize=14)

    # Ajouter le label commun de l'axe x
    ax[2].set_xlabel(r'$\frac{ct}{L}$', fontsize=20)

    # Ajustement automatique pour éviter le chevauchement
    plt.tight_layout()

    # Sauvegarde et affichage
    plt.savefig(f'../plot/diagnosis_combined_{scheme}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
def plot_convergence(U_32, U_64, U_128, t_32, t_64, t_128, A, config_32, config_64, config_128, schemes):
    x_32 = []
    x_64 = []
    x_128 = []
    R_32 = []
    R_64 = []
    R_128 = []

    # Loop over the schemes
    for i in range(len(U_32)):
        # Calculate the spatial coordinates for each resolution
        x_32.append(np.linspace(-config_32[i]["L"]/2, config_32[i]["L"]/2, int(config_32[i]["N"]), endpoint=False))
        x_64.append(np.linspace(-config_64[i]["L"]/2, config_64[i]["L"]/2, int(config_64[i]["N"]), endpoint=False))
        x_128.append(np.linspace(-config_128[i]["L"]/2, config_128[i]["L"]/2, int(config_128[i]["N"]), endpoint=False))
        
        # Compute diagnostics
        I_32, E_32, R_32_all = diagnosis(U_32[i], config_32[i]["h"], A, config_32[i]["sigma"], t_32[i], x_32[i], config_32[i]["c"], config_32[i]["L"])
        I_64, E_64, R_64_all = diagnosis(U_64[i], config_64[i]["h"], A, config_64[i]["sigma"], t_64[i], x_64[i], config_64[i]["c"], config_64[i]["L"])
        I_128, E_128, R_128_all = diagnosis(U_128[i], config_128[i]["h"], A, config_128[i]["sigma"], t_128[i], x_128[i], config_128[i]["c"], config_128[i]["L"])
        
        # Find the relevant time step index for each resolution
        index1 = int(1/config_32[i]["dt"]) + 1
        R_32.append(R_32_all[index1])
        
        index2 = int(1/config_64[i]["dt"]) + 1
        R_64.append(R_64_all[index2])
        
        index3 = int(1/config_128[i]["dt"]) + 1
        R_128.append(R_128_all[index3])

    # Define the colors
    color = ['black', 'darkgreen', 'darkorange', 'darkblue']
    
    # Plotting the global error (R) for each scheme
    plt.figure(figsize=(8, 6))  # Set figure size
    
    ax = plt.gca()  # Get current axes
    
    for i, scheme in enumerate(schemes):
        ax.loglog([1/8, 1/4, 1/2], 
                  [R_128[i], R_64[i], R_32[i]], 
                  marker='o', 
                  label=f'{scheme}', 
                  linestyle='-', 
                  markersize=8,
                  color=color[i])
    
    # Set labels for x and y axes with appropriate font sizes
    ax.set_xlabel(r'$\frac{h}{\sigma}$', fontsize=25)
    ax.set_ylabel(r'$R^{n}_{h}$', fontsize=25)
    
    # Set legend with a smaller font size
    ax.legend(fontsize=12)
    
    # Set the ticks for the x-axis and y-axis and adjust their font size
    plt.tick_params(axis='both', labelsize=15)
    
    # Customize the grid: Make it grid lines more visible and style them nicely
    ax.grid(True, which='both', linestyle='--', linewidth=0.7, color='grey')
    
    # Set x-axis ticks and labels (customizing the fractions for readability)
    ax.set_xticks([1/8, 1/4, 1/2])  # Define x-ticks
    ax.set_xticklabels([r'$\frac{1}{8}$', r'$\frac{1}{4}$', r'$\frac{1}{2}$'], fontsize=15)
    
    # Save the plot as a PNG image with tight bounding box for better layout
    plt.tight_layout()
    plt.savefig('../plot/convergence.png', dpi=300, bbox_inches='tight')
    
    # Display the plot
    plt.show()
    
    order_E2 = np.log(np.sqrt(R_32[0]/R_128[0]))/np.log(4)
    print(f"order of E2 is {order_E2}")
        
    order_E4 = np.log(np.sqrt(R_32[1]/R_128[1]))/np.log(4)
    print(f"order of E4 is {order_E4}")
    
    order_I4 = np.log(np.sqrt(R_32[2]/R_128[2]))/np.log(4)
    print(f"order of I4 is {order_I4}")
    
    order_ED = np.log(np.sqrt(R_32[3]/R_128[3]))/np.log(4)
    print(f"order of ED is {order_ED}")

      

    
        

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
    #plot_versus(U_steps_E2_N32, time_steps_E2_N32, grid_config_E2_N32, "E2")
    
    #### E4 32####
    grid_config_E4_N32, time_steps_E4_N32, U_steps_E4_N32 = extract_simulation_data("../data/results_E4_N32.txt")
    A = 1   
    #plot_versus(U_steps_E4_N32, time_steps_E4_N32, grid_config_E4_N32, "E4")
    
    #### I4 32####
    grid_config_I4_N32, time_steps_I4_N32, U_steps_I4_N32 = extract_simulation_data("../data/results_I4_N32.txt")
    A = 1    
    #plot_versus(U_steps_I4_N32, time_steps_I4_N32, grid_config_I4_N32, "I4")
    
    #### ED 32####
    grid_config_ED_N32, time_steps_ED_N32, U_steps_ED_N32 = extract_simulation_data("../data/results_ED_N32.txt")
    A = 1     
    #plot_versus(U_steps_ED_N32, time_steps_ED_N32, grid_config_ED_N32, "ED")
    
    #### E2 64 ####
    grid_config_E2_N64, time_steps_E2_N64, U_steps_E2_N64 = extract_simulation_data("../data/results_E2_N64.txt")
    # anim(x, U_steps, L, time_steps)
    A = 1  # Attention si on modifie le code   
    #plot_versus(U_steps_E2_N64, time_steps_E2_N64, grid_config_E2_N64, "E2")
    
    #### E4 64 ####
    grid_config_E4_N64, time_steps_E4_N64, U_steps_E4_N64 = extract_simulation_data("../data/results_E4_N64.txt")
    A = 1   
    #plot_versus(U_steps_E4_N64, time_steps_E4_N64, grid_config_E4_N64, "E4")
    
    #### I4 64 ####
    grid_config_I4_N64, time_steps_I4_N64, U_steps_I4_N64 = extract_simulation_data("../data/results_I4_N64.txt")
    A = 1    
    #plot_versus(U_steps_I4_N64, time_steps_I4_N64, grid_config_I4_N64, "I4")
    
    #### ED 64 ####
    grid_config_ED_N64, time_steps_ED_N64, U_steps_ED_N64 = extract_simulation_data("../data/results_ED_N64.txt")
    A = 1     
    #plot_versus(U_steps_ED_N64, time_steps_ED_N64, grid_config_ED_N64, "ED")
    
    #### E2 128 ####
    grid_config_E2_N128, time_steps_E2_N128, U_steps_E2_N128 = extract_simulation_data("../data/results_E2_N128.txt")
    # anim(x, U_steps, L, time_steps)
    A = 1  # Attention si on modifie le code   
    #plot_versus(U_steps_E2_N128, time_steps_E2_N128, grid_config_E2_N128, "E2")
    
    #### E4 128 ####
    grid_config_E4_N128, time_steps_E4_N128, U_steps_E4_N128 = extract_simulation_data("../data/results_E4_N128.txt")
    A = 1   
    #plot_versus(U_steps_E4_N128, time_steps_E4_N128, grid_config_E4_N128, "E4")
    
    #### I4 128 ####
    grid_config_I4_N128, time_steps_I4_N128, U_steps_I4_N128 = extract_simulation_data("../data/results_I4_N128.txt")
    A = 1    
    #plot_versus(U_steps_I4_N128, time_steps_I4_N128, grid_config_I4_N128, "I4")
    
    #### ED 128 ####
    grid_config_ED_N128, time_steps_ED_N128, U_steps_ED_N128 = extract_simulation_data("../data/results_ED_N128.txt")
    A = 1     
    #plot_versus(U_steps_ED_N128, time_steps_ED_N128, grid_config_ED_N128, "ED")
    
    
    A=1
    #plot_global_diagnostics(U_steps_E2_N32,U_steps_E2_N64,U_steps_E2_N128,time_steps_E2_N32,time_steps_E2_N64,time_steps_E2_N128, A,grid_config_E2_N32,grid_config_E2_N64,grid_config_E2_N128,"E2")
    #plot_global_diagnostics(U_steps_E4_N32,U_steps_E4_N64,U_steps_E4_N128,time_steps_E4_N32,time_steps_E4_N64,time_steps_E4_N128, A,grid_config_E4_N32,grid_config_E4_N64,grid_config_E4_N128,"E4")
    #plot_global_diagnostics(U_steps_I4_N32,U_steps_I4_N64,U_steps_I4_N128,time_steps_I4_N32,time_steps_I4_N64,time_steps_I4_N128, A,grid_config_I4_N32,grid_config_I4_N64,grid_config_I4_N128,"I4")
    #plot_global_diagnostics(U_steps_ED_N32,U_steps_ED_N64,U_steps_ED_N128,time_steps_ED_N32,time_steps_ED_N64,time_steps_ED_N128, A,grid_config_ED_N32,grid_config_ED_N64,grid_config_ED_N128,"ED")
    
    plot_convergence([U_steps_E2_N32,U_steps_E4_N32,U_steps_I4_N32,U_steps_ED_N32],[U_steps_E2_N64,U_steps_E4_N64,U_steps_I4_N64,U_steps_ED_N64],[U_steps_E2_N128,U_steps_E4_N128,U_steps_I4_N128,U_steps_ED_N128],[time_steps_E2_N32,time_steps_E4_N32,time_steps_I4_N32,time_steps_ED_N32],[time_steps_E2_N64,time_steps_E4_N64,time_steps_I4_N64,time_steps_ED_N64],[time_steps_E2_N128,time_steps_E4_N128,time_steps_I4_N128,time_steps_ED_N128], A,[grid_config_E2_N32,grid_config_E4_N32,grid_config_I4_N32,grid_config_ED_N32],[grid_config_E2_N64,grid_config_E4_N64,grid_config_I4_N64,grid_config_ED_N64],[grid_config_E2_N128,grid_config_E4_N128,grid_config_I4_N128,grid_config_ED_N128],["E2","E4","I4","ED"])


    

    
    
