import os
import scipy.io as sio
import numpy as np
from scipy import signal
from scipy.io import loadmat
import matplotlib.pyplot as plt
import math  # needed for the floor function
# port control

def testReadInAndVisualize():
    # Set parameters
    filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'
    filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat'

    noStates = 6  # corresponds to 'n' in MATLAB code
    n = noStates
    ny = 2
    poly = 0  # switch used for 12 states
    NTMstr = ''
    Vinf1 = 8  # Vinf1 = utemp(1,1,t0) in MATLAB
    Vstr = ''
    yawmode = 0
    percentTrain = 0.7
    t0 = 240 -1
    # sys_red = {}  # 'sys_red' is supposed to be a dictionary
    folderName = f"Vin_Vinf{Vinf1:.1f}{Vstr}{NTMstr}_states{noStates:02d}".replace('.', 'dot')
    fileName = f"stateName_K{noStates:02d}_P{poly}.mat"  # stateName

    # Define turbine coordinates and mesh specifications
    Wp = {
        'turbine': {
            'Crx': [400., 1032.1],  # X-coordinates of turbines (multiplied by 1e3 in the original code)
            'Cry': [400., 398.5747]  # Y-coordinates of turbines
        },
        'mesh': {
            'Lx': 1882.1,  # Domain length in x-direction
            'Ly': 800.0,  # Domain length in y-direction
            'Nx': 50,  # Number of cells in x-direction
            'Ny': 25  # Number of cells in y-direction
        }
    }

    # Set path to data directories
    codedir = os.path.realpath(__file__)
    parentdir = os.path.dirname(os.path.dirname(codedir))
    dataOutDir = os.path.join(parentdir, 'DataInOutWfSim')

    # Conditional directory naming based on 'poly' and 'yawmode'
    if poly == 1 and yawmode == 0:
        dirdmdName = 'eDMDresults_UasOutput_poly'
    elif poly == 0 and yawmode == 0:
        dirdmdName = 'eDMDresults_UasOutput'
    elif poly == 1:
        dirdmdName = 'eDMDresults_UasOutput_MIMO_poly'
    else:
        dirdmdName = 'eDMDresults_UasOutput_MIMO'

    dmddir = os.path.join(dataOutDir, dirdmdName)
    dirFig = os.path.join(dmddir, folderName)  # this depends on vinf
    dirFile = os.path.join(dirFig, fileName)

    # Load Koopman matrix and simulation run
    K = loadKoopmanSystem(dirFile)
    tmpId = loadSimulationRun(parentdir, filenameId)
    Inputs, Inputs_val, Outputs, Outputs_val, Deterministic, Deterministic_val, Input_u, Inputvalid_u = assignVariableNames(tmpId, yawmode, percentTrain, t0, Wp)
    # K_psi(Deterministic, PolyLiftingFunction=0, noStates=6, structPrev=None, yawmode=0)

    # Generate prediction data
    InputsAll = np.vstack((Inputs, Input_u))
    ysim, FITje = eDMDsimulate(Deterministic, InputsAll, K, poly, n, yawmode=0)

    InputsAll = np.vstack((Inputs_val, Inputvalid_u))
    ysim_val, FITje_val = eDMDsimulate(Deterministic_val, InputsAll, K, poly, n, yawmode=0)

    print('Fitje')
    print(FITje)
    print('Fitje_val')
    print(FITje_val)

    # Create a new figure
    plt.figure()

    # Loop over the first two indices (1 and 2 in MATLAB, 0 and 1 in Python)
    for idx in range(2):
        # Add subplot (2 rows, 1 column, current plot)
        plt.subplot(2, 1, idx + 1)

        # Plot the line from "Deterministic"
        plt.plot(Deterministic_val[idx, :], '-', label='real')

        # Plot 'ysim' on the same graph. t:
        plt.plot(ysim_val[idx, :], '--', label='estimate')

        # Add a legend to the plot
        plt.legend()

    # Show the figure
    plt.show()
    # Set the desired file name
    filename = filenameId + "Py.png"

    # Save the figure as a PNG file with the specified name
    plt.savefig(filename)


def loadKoopmanSystem(dirFile):
    data = loadmat(dirFile)  # Load data
    # not used: FITje = data['FITje'] FITje_val = data['FITje_val'] dirName = data['dirName'] xo = data['xo']
    Koop = data['Koop']
    K = Koop['K'][0, 0].T
    # not used:  A = K[0:n, 0:n].T B = K[n:, 0:n].T  C = np.concatenate((np.eye(2), np.zeros((2, 4))), axis=1) D =
    # np.zeros((2, 3))
    # Adjust dimensions based on your specific system sys_red0 = signal.StateSpace(A, B, C, D, dt=1)
    return K

def loadSimulationRun(parentdir, filenameId):
    # Load data: Simulation data
    DataIn = os.path.join(parentdir, 'DataT2OLWFSim')
    sepStr = filenameId.find('_')
    if not os.path.isdir(DataIn):
        print('Warning: DataIn directory does not exist')
        return

    # load identification data
    dirName = os.path.join(DataIn, filenameId[:sepStr] + '_OL_Ct', filenameId)
    tmpId = sio.loadmat(dirName)  # This loads the .mat file as a dictionary
    # 'Ct1', 'Ct2', 'FT1', 'FT2', 'PT1', 'PT2', 'Ur1', 'Ur2', 'p', 'phi1', 'phi2', 'u', 'v': keys in the dictionary
    return tmpId

def assignVariableNames(tmpId, yawmode, percentTrain, t0, Wp):
    # Assuming tmpId is a dictionary where keys are the field names
    # field_names = list(tmpId.keys())
    # Creating new global variables based on the field names
    # for field in field_names:
    # The following line creates a new global variable with the name stored in 'field'
    # and sets its value to the corresponding value from the 'tmpId' dictionary
    # globals()[field] = tmpId[field]
    # Assuming 'Ur1', 'percentTrain', and 't0' are already defined
    # as well as 'Wp' dictionary from the previous example

    Ur1 = tmpId['Ur1'][0,:]
    Ur2 = tmpId['Ur2'][0,:]
    Ct1 = tmpId['Ct1'][0,:]
    Ct2 = tmpId['Ct2'][0,:]
    phi1 = tmpId['phi1'][0,:]
    PT1 = tmpId['PT1'][0,:]
    PT2 = tmpId['PT2'][0,:]
    FT1 = tmpId['FT1'][0,:]
    FT2 = tmpId['FT2'][0,:]
    u = tmpId['u'][:,:]

    tident = math.floor(Ur1.shape[0] * percentTrain)   # floor is used for rounding down
    tval = tident + t0

    # Calculate the grid indices for the turbine's position
    nTx = round(Wp['turbine']['Crx'][0] / Wp['mesh']['Lx'] * Wp['mesh']['Nx']) - 1 -1
    nTy = round(Wp['turbine']['Cry'][0] / Wp['mesh']['Ly'] * Wp['mesh']['Ny']+0.01) -1

    # Reshape 'u' into a 3-dimensional array based on the dimensions of 'Wp.mesh'
    uPython = u.T
    kk = u.shape[1] // Wp['mesh']['Ny']  # Using floor division to get an integer result
    utempPython = np.reshape(uPython,(kk, Wp['mesh']['Ny'],Wp['mesh']['Nx']))
    # plt.figure()
    # plt.imshow(utempPython[1, :, :], aspect='auto', cmap='viridis')
    # plt.show()

    # Extract a specific series from 'utemp' at position [nTx, nTy, :]
    utempT = utempPython[:,nTy,nTx]
    # plt.figure()
    # plt.plot(utempT)
    # plt.show()

    # Slicing the array for different phases (training and validation)
    Input_u = utempT[t0:tident]  # Extracting a portion from 't0' to 'tident'
    Inputvalid_u = utempT[tval:]  # Extracting from 'tval' to the end

    # For VinfAll and Vinf, assuming it's a kind of velocity or reference measurement
    VinfAll = utempPython [:,0, 0]
    Vinf = VinfAll[tval:]  # Extracting the validation phase data points

    if yawmode == 0:
        Inputs = np.vstack((Ct1[t0:tident], Ct2[t0:tident]))
        Inputs_val = np.vstack((Ct1[tval:], Ct2[tval:]))
    else:
        Inputs = np.vstack((Ct1[t0:tident], Ct2[t0:tident], phi1[t0:tident]))
        Inputs_val = np.vstack((Ct1[tval:], Ct2[tval:], phi1[tval:]))

    # Prepare outputs
    Outputs = np.vstack((PT1[t0:tident], PT2[t0:tident], FT1[t0:tident], FT2[t0:tident]))
    Outputs_val = np.vstack((PT1[tval:], PT2[tval:], FT1[tval:], FT2[tval:]))

    # Prepare 'deterministic' states (effective wind speeds)
    Deterministic = np.vstack((Ur1[t0:tident], Ur2[t0:tident]))
    Deterministic_val = np.vstack((Ur1[tval:], Ur2[tval:]))
    return Inputs, Inputs_val, Outputs, Outputs_val, Deterministic, Deterministic_val, Input_u, Inputvalid_u

def K_psi(in1, PolyLiftingFunction=0, noStates=6, structPrev=None, yawmode=0):
    # Initialize structPrev if it's None
    if structPrev is None:
        structPrev = {
            'Ur1_prev1': in1[0,1],
            'Ur2_prev1': in1[0,2],
            'dUr1_prev1': 0,
            'dUr2_prev1': 0,
            'M1': in1[0,1],
            'M2': in1[0,2],
            'k': 1
        }

    # Non-linear lifting function
    Ur1_prev1 = structPrev['Ur1_prev1']
    Ur2_prev1 = structPrev['Ur2_prev1']
    dUr1_prev1 = structPrev['dUr1_prev1']
    dUr2_prev1 = structPrev['dUr2_prev1']
    M1 = structPrev['M1']
    M2 = structPrev['M2']
    k = structPrev['k']

    # Assign values in in1-vector to variables
    if in1.ndim > 1:
        Ur1 = in1[0, :]
        Ur2 = in1[1, :]
    else:
        Ur1 = in1[0]
        Ur2 = in1[1]

    ct1 = ct2 = phi1 = U1 = None
    if in1.shape[0] > 2 and in1.ndim > 1:  # if there are more rows, meaning more variables
        ct1 = in1[2, :]
        ct2 = in1[3, :]
    elif in1.shape[0] > 2:
        ct1 = in1[2]
        ct2 = in1[3]

    if in1.shape[0] > 5 and in1.ndim > 1:
        phi1 = in1[4, :]
        U1 = in1[5, :]
    elif in1.shape[0] > 5:
        phi1 = in1[4]
        U1 = in1[5]
    elif in1.shape[0] > 4 and in1.ndim > 1 and yawmode:
        phi1 = in1[4, :]
    elif in1.shape[0] > 4 and yawmode:
        phi1 = in1[4]
    elif in1.shape[0] > 4 and in1.ndim > 1:
        U1 = in1[4, :]
    elif in1.shape[0] > 4:
        U1 = in1[4]
        # else: do nothing

    # Create lifting functions and other calculations
    if in1.ndim > 1:
        Ur1_prev = np.hstack(([Ur1_prev1], Ur1[:-1]))  # previous states
        Ur2_prev = np.hstack(([Ur2_prev1], Ur2[:-1]))
    else:
        Ur1_prev = Ur1_prev1 # previous states
        Ur2_prev = Ur2_prev1

    dUr1 = Ur1 - Ur1_prev  # Difference to previous state
    dUr2 = Ur2 - Ur2_prev

    if in1.ndim > 1:
        dUr1_prev = np.hstack(([dUr1_prev1], dUr1[:-1]))
        dUr2_prev = np.hstack(([dUr2_prev1], dUr2[:-1]))
    else:
        dUr1_prev = dUr1_prev1
        dUr2_prev = dUr2_prev1

    ddUr1 = dUr1 - dUr1_prev
    ddUr2 = dUr2 - dUr2_prev

    dUr1sqr = Ur1 ** 2 - Ur1_prev ** 2
    dUr2sqr = Ur2 ** 2 - Ur2_prev ** 2

    n = 25
    lVec = min(n, k)
    M1 = (lVec - 1) / lVec * M1 + 1 / lVec * Ur1_prev
    M2 = (lVec - 1) / lVec * M2 + 1 / lVec * Ur2_prev

    DUr1 = Ur1 - M1
    DUr2 = Ur2 - M2

    # More calculations depending on PolyLiftingFunction
    if PolyLiftingFunction:
        Xaug2 = np.vstack(
            [Ur1, Ur2, Ur1 ** 2, Ur2 ** 2, Ur1 ** 3, Ur2 ** 3, Ur1 ** 4, Ur2 ** 4, Ur1 ** 5, Ur2 ** 5, Ur1 ** 6,
             Ur2 ** 6])
        stateName = 'Ur1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;Ur1.^4;Ur2.^4;Ur1.^5;Ur2.^5;Ur1.^6;Ur2.^6'
    else:
        # The operations are analogous to MATLAB but notice that indexing starts from 0
        XaugAll = np.vstack(
            [Ur1, Ur2, Ur1 ** 2, Ur2 ** 2, Ur1 ** 3, Ur2 ** 3, DUr1, DUr2, DUr1 ** 2, DUr2 ** 2, DUr1 ** 3,
             DUr2 ** 3,
             dUr1, dUr2, dUr1sqr, dUr2sqr, ddUr1, ddUr2, dUr1_prev, dUr2_prev])
        Xaug2 = XaugAll[0:noStates, :]

    # Construct state names
    allNames = ['Ur1', 'Ur2', 'Ur1^2', 'Ur2^2', 'Ur1^3', 'Ur2^3', 'DUr1', 'DUr2', 'DUr1^2', 'DUr2^2', 'DUr1^3',
                'DUr2^3',
                'dUr1', 'dUr2', 'dUr1sqr', 'dUr2sqr', 'ddUr1', 'ddUr2', 'dUr1_prev', 'dUr2_prev']
    stateName = ';'.join(allNames[:noStates])

    # Update structPrev for next iteration
    if in1.ndim > 1:
        structPrev['Ur1_prev1'] = Ur1[-1]
        structPrev['Ur2_prev1'] = Ur2[-1]
        structPrev['dUr1_prev1'] = dUr1[-1]
        structPrev['dUr2_prev1'] = dUr2[-1]
    else:
        structPrev['Ur1_prev1'] = Ur1
        structPrev['Ur2_prev1'] = Ur2
        structPrev['dUr1_prev1'] = dUr1
        structPrev['dUr2_prev1'] = dUr2

    structPrev['M1'] = M1
    structPrev['M2'] = M2
    structPrev['k'] = k + 1

    # Output
    return Xaug2, stateName, structPrev

def eDMDsimulate(Deterministic, InputsAll, K, poly, n, yawmode=0):
    """
    eDMDsimulate generates simulated output from the states
    Parameters:
    Deterministic (np.array): Array representing deterministic states
    InputsAll (np.array): Array of inputs
    K (np.array): Koopman matrix
    poly (TYPE): Polynomial order (translation requires knowledge of original 'poly' function or type)
    n (int): Dimension or other parameter (specific interpretation requires context from original code)
    yawmode (int, optional): Mode related to 'yaw' (default is 0 if not provided)
    Returns:
    ysim (np.array): Simulated outputs
    FITje (np.array): Fitness or error metric between the deterministic states and simulated outputs
    """

    # Size states, inputs, outputs
    ny = 2  # number outputs:  Ur1, Ur2 (no lifted states)
    Ky = K[:ny, :]  # Get the relevant rows of the Koopman matrix

    # Initialize variables and structure similar to 'structPrev' in MATLAB
    structPrev = {
        'Ur1_prev1': Deterministic[0, 0],
        'Ur2_prev1': Deterministic[1, 0],
        'dUr1_prev1': 0,
        'dUr2_prev1': 0,
        'M1': Deterministic[0, 0],  # Moving mean
        'M2': Deterministic[1, 0],
        'k': 1
    }

    # Use non-linear Koopman functions based on deterministic states
    xsim = np.empty((Ky.shape[1], Deterministic.shape[1]), dtype=np.float64)
    xsim.fill(np.nan)  # Initialize with NaN
    tmp, _ , _ = K_psi(np.concatenate((Deterministic[:, 0], InputsAll[:, 0])), poly, n, structPrev, yawmode)
    tmpIn = np.array(InputsAll[:, 0])
    xsim[:, 0] = np.concatenate((np.squeeze(tmp), tmpIn))

    for idx in range(Deterministic.shape[1] - 1):
        xsim[:ny, idx + 1] = np.dot(Ky, xsim[:, idx])
        tmp, _, structPrev = K_psi(
            np.concatenate((xsim[:ny, idx + 1], InputsAll[:, idx + 1])), poly, n, structPrev, yawmode
        )
        tmpIn = np.array(InputsAll[:, idx + 1])
        tmpIn_reshaped = tmpIn[np.newaxis, :]
        xsim[:, idx + 1] =  np.concatenate((np.squeeze(tmp),tmpIn))

    ysim = xsim[:ny, :]

    # Calculate the covariance matrices. Note that we need to transpose the matrices,
    # because np.cov assumes rows as variables and columns as observations.
    covariance_matrix = np.cov((Deterministic - ysim))
    total_covariance = np.cov(Deterministic)

    # Calculate fitness metric (FITje). We're avoiding division by zero by adding a small number in the denominator.
    FITje = np.maximum(np.diag(100 * (np.eye(Deterministic.shape[0]) -
                                      covariance_matrix / (total_covariance + 1e-9))), 0)

    return ysim, FITje

# NOTE: 'K_psi' and 'cov' are assumed to be user-defined functions from the MATLAB code.
# You should replace these with the appropriate Python translations of these functions.

if __name__ == "__main__":
    testReadInAndVisualize()

""" 
num_time_steps = 15000
time_vector = np.arange(num_time_steps)  # This creates an array of time points [0, 1, 2, ..., num_time_steps-1].

# Calculate the step response using the signal.step function.
t_out, responses = signal.dstep(sys_red0, t=time_vector)
# dc_gain = signal.dcgain(sys)

# 'response' is a list of arrays, one for each output. We need to extract the array from the list.
# If your system has more than one output, you'll need to handle each element of the list.
# y_out = response[0]
responses = np.array(responses)
print(responses[-1])

# We will get 6 series for a 2x3 system.
num_series = responses.shape[0] * responses.shape[2]

# Reshape the data to be a 2D array (num_series, num_time_steps)
responses_reshaped = responses.transpose(0, 2, 1).reshape(num_series, -1)

# Create subplots for each time series
fig, axs = plt.subplots(responses.shape[0], responses.shape[2], figsize=(10, 5 * num_series))

# Plot each time series
for i in range(num_series):
    y_out = responses_reshaped[i, :]
    axs[i].plot(t_out, y_out, basefmt=" ")
    axs[i].set_title(f'Step Response - Series {i+1}')
    axs[i].set_xlabel('Time steps')
    axs[i].set_ylabel('Amplitude')
    axs[i].grid(True)

# Adjust layout and show plots
plt.tight_layout()
plt.show()
"""
