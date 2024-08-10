import numpy as np

# Function to calculate cross-section area of reinforcement in each row
def rebar_area(array):
    arr_As = np.array([])
    z = np.array([])
    for i, x in enumerate(array):
        As_row = (array[i][0])**2 * np.pi/4 * array[i][1]
        arr_As = np.append(arr_As, As_row)
        z = np.append(z, array[i][2])      
    return arr_As, z

# Function to calculate cross-section area of strands in each row
def strand_area(array, Ap):
    arr_Ap = np.array([])
    z = np.array([])
    for i, x in enumerate(array):
        Ap_row = (array[i][1]) * Ap
        arr_Ap = np.append(arr_Ap, Ap_row)
        z = np.append(z, array[i][0]) 
    return arr_Ap, z

# Function to calculate effective depths - tensile reinforcement
def dtens(H, z):
    d = H - z
    return d

# Function to calculate effective depths - compressive reinforcement
def dcomp(z):
    d = z
    return d