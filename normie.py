import numpy as np
import os

def grep_until_end_pattern(pattern, end_pattern, list_,ignore_n_init_rows=0,ignore_n_final_rows = 0):
    r = []
    
    for i, row in enumerate(list_):
        if pattern in row:
            newres = []
            for j, subrow in enumerate(list_[i:]):
                newres.append(subrow)
                if end_pattern in subrow:
                    break
            newres = newres[ignore_n_init_rows:-ignore_n_final_rows]
            r.append(newres)
    return r



def get_normal_modes(orca_output): 
    """
    Input:  readlines() version of orca 4 output file
    Output: 3Nx3N numpy array, normal coordinate matrix
    """
    normal_ugly = grep_until_end_pattern("NORMAL MODES", "IR SPECTRUM", orca_output, 7,4)[0]
    rs = []
    r = []
    for row in normal_ugly:
        row = row.strip().split()
        if row[0] == "0" and len(r) > 0:
            del r[-1]
            rs.append(r)
            r = []
        r.append(row[1:])
    rs.append(r)
    rs = [np.asfarray(r) for r in rs[1:]]
    normal_coord = np.concatenate(rs, axis=1)
    return normal_coord
def get_atoms_xyz(orca_output):
    """
    Input:  readlines() version of orca 4 output file
    Output: list of atoms, Nx3 numpy array of cartesian cordinates in Å
    """
    coord = grep_until_end_pattern("CARTESIAN COORDINATES (ANGSTROEM)", "CARTESIAN COORDINATES (A.U.)", orca_output, 2,3 )[-1]
    atoms = [row.split()[0] for row in coord] 
    xyz = np.asfarray([row.split()[1:] for row in coord])
    return atoms, xyz
        
def get_freqs(orca_output):
    """
    Input:  readlines() version of orca 4 output file
    Output: 3N numpy array of imaginary frequencies
    """
    vibfreq = grep_until_end_pattern("VIBRATIONAL FREQUENCIES", "NORMAL MODES", orca_output, 3, 4)[0]
    freqs = np.asfarray([row.split()[1] for row in vibfreq])
    return freqs
def displace_coordinates(normal_coord, xyz, index, scale):
    """
    Displaces coordinates with the Cartesian displacement of normal mode, scaled by scalar scaling factor
        
    """
    N = xyz.shape[0]
    disp = scale*np.reshape(normal_coord[:,index], (N,3))
    new_xyz = xyz + disp
    return new_xyz

def displace_imaginary_frequencies(path, freqs, normal_coord, at, xyz, scale, write=True,\
                                   imaginary_threshold = 0, only_smallest = True,
                                  write_original=True):
    
    """
    Displaces coordinates along the normal modes of imaginary frequencies.
    
    Input:
    
    freqs : np array of 3N imaginary frequencies
    normal_coord : np array of 3Nx3N cartesian displacements
    at, xyz = coords
    scale = scale factor
    write = writes files with "i53.12.xyz" -format
    imaginary_threshold : frequencies smaller than this value will be taken as imaginary
    only_smallest : operate only on the smallest imaginary frequency
    
    Output:
        returns list of displaced coordinates
    
    """
    
    imaginary_freq_indices = np.where(freqs<imaginary_threshold)[0]
        
    disp_coords = []
    
    if write_original:
        write_coordinates(at,xyz, "0", path)
        
    for i, imaginary_freq_idx in enumerate(imaginary_freq_indices):
        if only_smallest and i>0:
            break
        new_xyz = displace_coordinates(normal_coord, xyz,imaginary_freq_idx, scale)
        disp_coords.append(new_xyz)
        if write:
            fname = str(freqs[imaginary_freq_idx]).replace("-","i")
            write_coordinates(at,new_xyz, "{}".format(fname), path)
    return disp_coords

def atoms_cartesian2xyz(at, xyz, xyz_format=True):
    """Writes proper xyz files, if xyz_format=True the header is written as well
    Input: 
    at : list of atoms
    xyz : numpy array of cartesian coordinates
    xyz_format : If true write the "N\n\n" header
    """
    
    N = len(at)
    if xyz_format:
        xyzout = "{}\n\n".format(N)
    else:
        xyzout =""
    for a,x in zip(at,xyz):
        xyzout += "{:5}{:10.6f}{:10.6f}{:10.6f}\n".format(a,x[0],x[1],x[2])
    return xyzout

def write_coordinates(at, xyz, outname, outpath = ""):
    "Writes coordinates to outname.xyz"
    xyzf = atoms_cartesian2xyz(at,xyz, xyz_format=True)
    with open(os.path.join(outpath, "{}.xyz".format(outname)), "w") as outf:
        for row in xyzf:
            outf.write(row)
            
            
def read_orca_file(path_to_outfile):
    try:
        with open(path_to_outfile) as orca_output:
            orca_output = orca_output.readlines()    
        return orca_output
    except:
        print("No output file found")
