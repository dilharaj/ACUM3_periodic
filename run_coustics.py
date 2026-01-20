import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
P_REF = 20e-6  # Reference pressure in Pa

case_folder = r"D:\Ansys\ACUM\cases\Tempest-2\case29"
surf_offset = 1917  # Adjust based on your file naming convention

omega = 367.0  # Angular velocity in rad/s
R = 0.3  # Reference radius in meters
NB = 11  # Number of blades
rhoRef = 1.225
a0 = 340
nTime = 360
    
def periodic_interp_fft(theta, p, refine):
    n = len(theta)
    P = np.fft.rfft(p)
    P_pad = np.zeros(refine * P.size, dtype=complex)
    P_pad[:P.size] = P
    p_fine = np.fft.irfft(P_pad, n=refine*n)
    return p_fine



# --- A-weighting function ---
def a_weighting_dB(f):
    f1, f2, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
    f_sq = f**2
    num = (f4**2) * (f_sq**2)
    den = (f_sq + f1**2) * (f_sq + f4**2) * np.sqrt((f_sq + f2**2) * (f_sq + f3**2))
    RA = num / den
    A = 20 * np.log10(RA) + 2.0
    A = np.where(f == 0, -np.inf, A)
    return A

# --- PSD from time signal ---
def compute_psd(x, fs):
    N = len(x)
    x = x - np.mean(x)
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(N, 1/fs)
    Sxx = (1.0 / (fs * N)) * np.abs(X)**2
    if N % 2 == 0:
        Sxx[1:-1] *= 2
    else:
        Sxx[1:] *= 2
    return f, Sxx

# --- OASPL functions ---
def oaspl_time(x):
    rms = np.sqrt(np.mean((x - np.mean(x))**2))
    return 20 * np.log10(rms / P_REF)

def oaspl_psd(f, Sxx):
    df = np.mean(np.diff(f))
    mean_sq = np.sum(Sxx * df)
    return 10 * np.log10(mean_sq / (P_REF**2))

def oaspl_A_weighted(f, Sxx):
    A_dB = a_weighting_dB(f)
    A_lin = 10 ** (A_dB / 10)
    A_lin = np.where(np.isfinite(A_dB), A_lin, 0.0)
    df = np.mean(np.diff(f))
    mean_sq = np.sum(Sxx * A_lin * df)
    return 10 * np.log10(mean_sq / (P_REF**2))



def make_surfaces(case_folder,omega,R):
    
    print("Running Acoustics Post-Processing...")
    #f_input = "input"
    f_input = rf"{case_folder}/input"
    Mref = omega*R/340
    # with open(f_input, 'w') as f:
    #     f.write("***Environmental\n")
    #     f.write(f"a0             	{a0}\n")
    #     f.write(f"rhoRef          {rhoRef}\n")
    #     f.write("pRef            101325\n")
    #     f.write("\n***Geometric\n")
    #     f.write(f"oM              {omega}\n")
    #     f.write("nSurf           1\n")
    #     f.write(f"nTime           {nTime}\n")
    #     f.write("dPsi            1\n")
    #     f.write("Minf            0\n")
    #     f.write(f"Mref            {Mref}\n")
    #     f.write("xy_angle        0\n")
    #     f.write("xz_angle        0       \n")
    #     f.write("CFDscale        1.0\n")
    #     f.write("OBSscale        1.0\n")
    #     f.write("\n")
    #     f.write("***Data\n")
    #     f.write("impermeable     1\n")
    #     f.write("periodic        1\n")
    #     f.write("lowpass         0\n")
    #     f.write("BB_noise        0\n")

    
    nel_per_file = 25000

    for ii in range(nTime):
        data = []
        print(f"Processing surface data for time step {ii+1}/{nTime}...")
        with open(rf"{case_folder}/surf_data/surf-{ii+surf_offset}", 'r') as f:
            next(f)  # skip line 0
            for line in f:
                data += [line.split()]
        n = len(data)
        nfiles = n // nel_per_file + 1

        if ii==0:
            print(f"Total number of surface elements: {n}")
            print(f"Number of surface files to be created: {nfiles}")

        cp = np.zeros(n,dtype=float)
        area = np.zeros(n,dtype=float)
        x = np.zeros(n,dtype=float)
        y = np.zeros(n,dtype=float)
        z = np.zeros(n,dtype=float)
        nx = np.zeros(n,dtype=float)
        ny = np.zeros(n,dtype=float)
        nz = np.zeros(n,dtype=float)

        for i in range(n):
            # for j in range(len(data[i+1])):
            #     print(data[i+1][j],", ")
            cp[i] = float(data[i][4])/(0.5*rhoRef*(Mref)**2*a0**2)
            area[i] = float(data[i][5])
            x[i] = float(data[i][1])
            y[i] = float(data[i][2])
            z[i] = float(data[i][3])
            nx[i] = -float(data[i][6])/float(data[i][5])
            ny[i] = -float(data[i][7])/float(data[i][5])
            nz[i] = -float(data[i][8])/float(data[i][5])

        # angles = 2*np.pi * np.arange(nTime) / nTime
        # cos_t = np.cos(angles)
        # sin_t = np.sin(angles)

        for k in range(nfiles):
            start = k*nel_per_file
            end = min((k+1)*nel_per_file, n)
            surf_file = rf"{case_folder}/surfaces/surface_{k+1}.dat"
            if ii==0:
                with open(surf_file, 'w', buffering=1024*1024) as f:
                    f.write(f"{end-start} {nTime}\n")
                 
                    lines = [] 
                    for j in range(start, end):
                        
                        lines.append(
                            "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n"
                            % (x[j], y[j], z[j], cp[j], area[j], nx[j], ny[j], nz[j])
                        )  

                    f.write("".join(lines))
            else:
                with open(surf_file, 'a', buffering=1024*1024) as f:
                   
                    lines = [] 
                    for j in range(start, end):
                        
                        lines.append(
                            "%.7e %.7e %.7e %.7e %.7e %.7e %.7e %.7e\n"
                            % (x[j], y[j], z[j], cp[j], area[j], nx[j], ny[j], nz[j])
                        )   

                    f.write("".join(lines))
        
    print("Surface data prepared for ACUM.")
    
    # fpA = r"Z:\UFX_dilhara\Ansys\3D_CFD_Optimizer\OptiSlang\src\Variation\pA.out"
    # fpL = r"Z:\UFX_dilhara\Ansys\3D_CFD_Optimizer\OptiSlang\src\Variation\pL.out"
    # fpT = r"Z:\UFX_dilhara\Ansys\3D_CFD_Optimizer\OptiSlang\src\Variation\pT.out"
    return
    
def run_acoustics():
    cmd = [
        "wsl",
        "/mnt/d/Ansys/ACUM/ACUM3_periodic-main/ACUM3_periodic-main/acum_cuda"
    ]

    with open(rf"case_folder/acum_log.txt", 'w') as log_file:
        log_file.write("Running ACUM with command:\n")
        log_file.write(" ".join(cmd) + "\n")
        result = subprocess.run(
            cmd,
            cwd=case_folder,
            stdout=log_file,
            stderr=log_file,
            text=True
        )

    print(result.stdout)
    print(result.stderr)
    if result.returncode != 0:
        print("Error: ACUM execution failed")
        return
    return

def post_process():
        

    fpA = rf"{case_folder}\pA.out"
    fpL = rf"{case_folder}\pL.out"
    fpT = rf"{case_folder}\pT.out"

    # Read results
    lst = []
    with open(fpA, 'r') as f:
        for line in f:
            lst += [line.split()]
    
    nObs = int(lst[1][1])
    theta = np.zeros(nTime,dtype=float)
    pA = np.zeros((nObs,nTime),dtype=float)

    for i in range(nTime):
        theta[i] = float(lst[2 + i][0])
        for j in range(nObs):
            pA[j][i] = float(lst[2 + i][1 + j])

    
    lst = []
    with open(fpL, 'r') as f:
        for line in f:
            lst += [line.split()]
    
    pL = np.zeros((nObs,nTime),dtype=float)
    for i in range(nTime):
        for j in range(nObs):
            pL[j][i] = float(lst[2 + i][1 + j])

    lst = []
    with open(fpT, 'r') as f:
        for line in f:
            lst += [line.split()]
    
    pT = np.zeros((nObs,nTime),dtype=float)

    for i in range(nTime):
        for j in range(nObs):
            pT[j][i] = float(lst[2 + i][1 + j])

    
    p_total = pA

    

    plt.figure()
    iobs = 1
    # plt.plot(theta*180/np.pi, pA[iobs, :])
    # plt.plot(theta*180/np.pi, pL[iobs, :], color='blue')
    # plt.plot(theta*180/np.pi, pT[iobs, :], color='red')

    plt.plot(theta*180/np.pi, p_total[iobs, :], linestyle='-')
    plt.xlabel("Theta (deg)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Total Pressure at Observer 1")
    plt.grid()
    plt.show()
    plt.savefig(rf"{case_folder}/acoustics_pressure_time_obs{iobs+1}.png")

    # SPL/OASPL Calculation

    rpm = omega*60/(2*np.pi)        
    f_rot = rpm / 60.0

    dtheta = theta[1] - theta[0]
    dt = dtheta / omega

    t = theta / (omega)
    fs = 1.0 / np.mean(np.diff(t))

    SPL = np.zeros((nObs, nTime//2 + 1))
    OASPL = np.zeros((nObs))
    OASPL_Aw = np.zeros((nObs))

    for i in range(nObs):
        freq, Sxx = compute_psd(p_total[i, :], fs)
        df = np.mean(np.diff(freq))
        SPL[i, :] = 10 * np.log10(Sxx * df/ (P_REF**2))
        OASPL_Aw[i] = oaspl_A_weighted(freq, Sxx)
        OASPL[i] = oaspl_time(p_total[i, :])





    # Plot SPL for first observer
    #plt.figure()
    plt.semilogx(freq, SPL[iobs])
    plt.ylim(0, np.max(SPL[iobs])+10)
    #plt.axvline(NB * f_rot, color='r', linestyle='--')
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("SPL [dB]")
    plt.show()
    plt.savefig(rf"{case_folder}/acoustics_SPL_obs{iobs+1}.png")

    # Save results
    with open(rf"{case_folder}/acoustics_SPL.csv", 'w') as f:
        f.write("Frequency(Hz)," + ",".join([f"Obs_{j+1}_SPL(dB)" for j in range(nObs)]) + "\n")
        for i in range(len(freq)):
            f.write(f"{freq[i]}," + ",".join([f"{SPL[j][i]}" for j in range(nObs)]) + "\n")
    with open(rf"{case_folder}/acoustics_OASPL.csv", 'w') as f:
        f.write("Observer,OASPL(dB),OASPL_A(dB)\n")
        for j in range(nObs):
            f.write(f"{j+1},{OASPL[j]},{OASPL_Aw[j]}\n") 

    mean_OASPL = np.mean(OASPL)
    mean_OASPL_A = np.mean(OASPL_Aw)
    max_OASPL = np.max(OASPL)
    max_OASPL_A = np.max(OASPL_Aw)
    print(f"Mean OASPL: {mean_OASPL:.2f} dB, Mean A-weighted OASPL: {mean_OASPL_A:.2f} dB")
    print(f"Max OASPL: {max_OASPL:.2f} dB, Max A-weighted OASPL: {max_OASPL_A:.2f} dB")

    return mean_OASPL, mean_OASPL_A, max_OASPL, max_OASPL_A


# make_surfaces(case_folder,omega,R)
# run_acoustics()
mean_OASPL, mean_OASPL_A, max_OASPL, max_OASPL_A = post_process()

