import os
import numpy as np

def grep(pattern,list_,):
    r=[]
    for row in list_:
        #print(line)
        if pattern in row:
            r.append(row)
    return r

def grep_until_end_pattern(pattern, end_pattern, list_,ignore_n_init_rows=0,ignore_n_final_rows = 0):
    r = []
    
    for i, row in enumerate(list_):
        if pattern in row:
            newres = []
            for j, subrow in enumerate(list_[i:]):
                #print(subrow)
                newres.append(subrow)
                if end_pattern in subrow:
                    break
            newres = newres[ignore_n_init_rows:-ignore_n_final_rows]
            r.append(newres)
    return r


def orca2json(path_to_output_file, system, opt=True, freq=True):

    calcdir = "/".(path_to_output_file.split("/")[:-1])

    # Try to open output file
    try:
        with open(path_to_output_file) as outf:
            outf = outf.readlines()   
    except: 
        # If no output file is found, assume that the calculation did not start
        print("Did not start", calcdir)
        return {"System": system, "Done":"Nostart"}

    # Get atoms, xyz
    try:
        coord = grep_until_end_pattern("CARTESIAN COORDINATES (ANGSTROEM)", "CARTESIAN COORDINATES (A.U.)", orca_output, 2,3 )[-1]
    except:
    # Get total energy
        coord = None
    try:
        total_energy = grep("Total Energy       :", outf)[-1]
        total_energy = float(total_energy.split()[3])
    except:
        print("No final energy", calcdir)
        total_energy = None
    # Orbital_energies
    try:
        orbital_energies = grep_until_end_pattern("ORBITAL ENERGIES",
                                                 "MULLIKEN POPULATION ANALYSIS",
                                                  outf, ignore_n_init_rows=4,
                                                  ignore_n_final_rows=3)[-1]
        occupation_energy = np.asfarray([[row.split()[1], row.split()[2]] for row in orbital_energies])
        homos = occupation_energy[occupation_energy[:,0]>0][:,1]
        lumos = occupation_energy[occupation_energy[:,0]==0][:,1]
    except:
        print("Problem with HOMOLUMO", calcdir)
        homos = None
        lumos = None

    # Mulliken atomic charges
    try:
        mulliken_charges = grep_until_end_pattern("MULLIKEN ATOMIC CHARGES",
                                                 "MULLIKEN REDUCED ORBITAL CHARGES",
                                                 outf,
                                                 ignore_n_init_rows=2,
                                                 ignore_n_final_rows=4)[-1]
        mulliken_charges = np.asfarray([row.rstrip().split()[-1] for row in mulliken_charges])
    except:
        print("No Mulliken", calcdir)
        mulliken_charges = None
    # Loewdin atomic charges
    try:
        loewdin_charges = grep_until_end_pattern("LOEWDIN ATOMIC CHARGES",
                                                 "LOEWDIN REDUCED ORBITAL CHARGES",
                                                 outf,
                                                 ignore_n_init_rows=2,
                                                 ignore_n_final_rows=4)[-1]
        loewdin_charges = np.asfarray([row.rstrip().split()[-1] for row in loewdin_charges])
    except:
        print("No loewdin", calcdir)
        loewdin_charges = None
    # Mayer population analysis
    try:
        mayer_valences = grep_until_end_pattern("MAYER POPULATION ANALYSIS",
                                      "Mayer bond orders larger than 0.1", outf,
                                      ignore_n_init_rows=11,
                                      ignore_n_final_rows=2)[-1]

        mayer_valences = np.asfarray([row.strip().split()[-3:]for row in mayer_valences])
        mayer_total_valence = mayer_valences[:,0]
        mayer_bonded_valence = mayer_valences[:,1]
        mayer_free_valence  = mayer_valences[:,2]
    except:
        print("No mayer ", calcdir)
        mayer_total_valence = None
        mayer_bonded_valence = None
        mayer_free_valence = None
    # Mayer bond orders
    try:
        mayer_bo = grep_until_end_pattern("Mayer bond orders larger than 0.1",
                                         "TIMINGS", outf,
                                         ignore_n_init_rows=1,
                                         ignore_n_final_rows=4)[0]
        mayer_bo = " ".join([row.rstrip() for row in mayer_bo]).split("B")[1:]
        n_atoms = mayer_bonded_valence.shape[0]
        mayer_bo_matrix = np.zeros((n_atoms,n_atoms))
        for row in mayer_bo:
            row = row.split()
            at1 = int(row[1].split("-")[0])
            at2 = int(row[3].split("-")[0])
            bo = float(row[-1])
            mayer_bo_matrix[at1,at2] = mayer_bo_matrix[at2,at1] = bo
    except:
        print("No mayer bo", calcdir)
        mayer_bo_matrix = None
    # D3 energy
    try:
        d3bj_energy = float(grep("Dispersion correction",outf)[-1].split()[-1])
    except:
        print("No d3", calcdir)
        d3bj_energy = None
    # Dipole moment
    try:
        dipole_debye = float(grep("Magnitude (Debye)", outf)[-1].split()[-1])
    except:
        print("No dipole", calcdir)
        dipole_debye = None
    result_dict =  {"System": system, "Coordinates": coord, "TotalEnergy" : total_energy,
           "OccupiedOrbitals": homos,
           "VirtualOrbitals": lumos,
           "MullikenCharges": mulliken_charges,
           "LoewdinCharges": loewdin_charges,
           "MayerTotValence": mayer_total_valence,
            "MayerFreeValence": mayer_free_valence,
            "MayerBondOrderMatrix": mayer_bo_matrix,
            "D3BJEnergy": d3bj_energy,
            "DipoleMoment":dipole_debye, "Done": "SP"}

    
    # Check if geometry has converged
    if opt:
        geoopt = grep("THE OPTIMIZATION HAS CONVERGED", outf)
        if len(geoopt)>0:
            geo_converged = True
            try: 
                with open(os.path.join(calcdir,"orca.xyz")) as optxyz:
                    optxyz = "".join(optxyz.readlines()[2:])
            except:
                optxyz = None
                print("No geo", calcdir)
        else:
            geo_converged = False
        # Check for imaginary frequencies
    if freq and geo_converged:
        
        # Get IR
        try:
            ir_spectrum = grep_until_end_pattern("VIBRATIONAL FREQUENCIES", 
                                                 "NORMAL MODES", outf,
                                                 ignore_n_init_rows=3,
                                                 ignore_n_final_rows=4)[0]
            ir_spectrum = [float(row.split()[1]) for row in ir_spectrum]
            ir_spectrum = np.asarray(ir_spectrum)
        
            n_imaginary_freq = ir_spectrum[ir_spectrum<-1].shape[0]
        except:
            print("No IR, calcdir")
            ir_spectrum = None
            n_imaginary_freq = None
        # Get dispersion energy
                
        
        #Thermodyn
        try:
            thermal_correction = float(grep("Total thermal correction", outf)[0].split()[-4])
            zpe = float(grep("Non-thermal (ZPE) correction", outf)[0].split()[-4])

            enthalpy = float(grep("Total Enthalpy", outf)[-1].split()[-2])
            entropy = float(grep("Final entropy term", outf)[-1].split()[-4])
            gibbs_free_enthalpy = float(grep("Final Gibbs free enthalpy", outf)[-1].split()[-2])
        except:
            print("No TD", calcdir)
            thermal_correction = None
            zpe = None
            enthalpy = None
            entropy = None
            gibbs_free_enthalpy = None
        result_dict =  {"System": system, "Coordinates": coord, "TotalEnergy" : total_energy,
               "OccupiedOrbitals": homos,
               "VirtualOrbitals": lumos,
               "MullikenCharges": mulliken_charges,
               "LoewdinCharges": loewdin_charges,
               "MayerTotValence": mayer_total_valence,
                "MayerFreeValence": mayer_free_valence,
                "MayerBondOrderMatrix": mayer_bo_matrix,
                "OptCoordinates" : optxyz,
                "GeoConverged" : geo_converged,
                "VibrationalFreq": ir_spectrum,
                "NImaginaryFreq": n_imaginary_freq,
                "D3BJEnergy": d3bj_energy,
                "DipoleMoment":dipole_debye,
                "ThermalCorrection": thermal_correction,
                "ZPE":zpe,
                "Enthalpy":enthalpy,
                "Entropy":entropy,
                "GibbsFreeEnergy":gibbs_free_enthalpy, "Done": "Optfreq"}
    return result_dict
