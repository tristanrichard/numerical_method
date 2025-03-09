import numpy as np

def extract_simulation_data(filename="data/results.txt"):
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
        if "#### End of iteration" in lines[i]:
            time = float(lines[i+1].split("=")[1].strip())
            U = np.array([float(x) for x in lines[i+2].strip().split()])
            time_steps.append(time)
            U_steps.append(U)
        i += 3  
    
    return grid_config, U_initial, time_steps, U_steps


grid_config, U_initial, time_steps, U_steps = extract_simulation_data()

# Affichage des résultats
print("Grid Configuration:", grid_config)
print("Initial U:", U_initial)
print("Time steps:", time_steps)
print("U at each time step:", U_steps)