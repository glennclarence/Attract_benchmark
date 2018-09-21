import numpy as np
import os
import Bio
import subprocess
import pandas as pd

#compares two arrays of positions of atoms and calculates the optimal mode defomation by calculating the scalar product of the positional difference


def getDelta(pos_bound, pos_unbound,num_atoms):
    assert len(pos_bound) == len(pos_unbound)
    delta = np.zeros((3,num_atoms))
    for k in range(3):
        delta[k] = pos_bound[k] - pos_unbound[k]
    return delta


def compare(pos_bound, pos_unbound, modes, num_atoms, num_modes):
    assert len(pos_bound) == len(pos_unbound)
    delta = np.zeros((3,num_atoms))
    scalar = np.zeros(num_modes)
    #build the dor product of the modes with the delta vector
    for k in range(3):
        delta[k] = pos_bound[k] - pos_unbound[k]
        for m in range(num_modes):
            scalar[m] += np.sum(np.inner(delta[k], modes[m][k]))
    return scalar


def norm_modeVec( modes, num_modes):
    amplitude = np.zeros(num_modes)
    for m in range(num_modes):
        amplitude[m] = np.linalg.norm(modes[m])
    return amplitude

def read_pdb(f):
  x,y,z, res, resn, atom,id_atom = [], [],[], [], [], [], []
  typ, charge = [],[]
  curr_resid = None
  resindex = 0
  for l in open(f).readlines():
    if not l.startswith("ATOM"): continue
    x.append(float(l[30:38]))
    y.append(float(l[38:46]))
    z.append(float(l[46:54]))
    id_atom.append(int(l[6:11]))
    atom.append(l[12:16])
    cresn = l[17:20]
    resid = l[21:26]

    if resid != curr_resid:
      curr_resid = resid
      resindex += 1
    res.append(int(resindex))
    resn.append(cresn)

    typ.append(int(l[57:59]))
    charge.append(l[61:67])
  return id_atom, atom, resn, res, x,y,z, typ, charge

def read_modes(filename_modes):
    with open(filename_modes, 'r') as f:
        lines = f.readlines()
        d_modes = {}
        mode_idx = 1
	amplitude = np.zeros(25)
	numModes = 0
        for i, line in enumerate(lines):
            if len(line.split()) == 2 and int(line.split()[0]) == mode_idx:
                x = []
                y = []
                z = []
		numModes += 1
                d_modes[mode_idx ] = ( float(line.split()[1]), x,y,z )
                mode_idx += 1
                #eigenvalues[mode_idx - 2] = float(line.split()[1])
            elif len(line.split()) == 3:
                d_modes[mode_idx - 1][1].append( float(line.split()[0]) )
                d_modes[mode_idx - 1][2].append( float(line.split()[1]) )
                d_modes[mode_idx - 1][3].append( float(line.split()[2]) )
		amplitude[mode_idx - 1] += float(line.split()[0]) **2 + float(line.split()[1]) **2 + float(line.split()[2]) **2
	amplitude = np.sqrt(amplitude)
	for mode in range(numModes):
		for i in range(len(d_modes[mode + 1][1])):
			d_modes[mode + 1][1][i] /= amplitude[mode + 1] if  amplitude[mode + 1] > 0 else 1
			d_modes[mode + 1][2][i] /= amplitude[mode + 1] if  amplitude[mode + 1] > 0 else 1
			d_modes[mode + 1][3][i] /= amplitude[mode + 1] if  amplitude[mode + 1] > 0 else 1
	d_modes

    return d_modes, numModes


def writePDB(num_atoms,id_atom, atom, resn, res, x,y,z, typ, charge, filename_pdb, chain = "" ):
    elem = 1

    file = open(filename_pdb, "w")
    for i in range(num_atoms):
        file.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}  {:3d}  {:5s} 0 1.00\n".format("ATOM", id_atom[i], atom[i], " ",resn[i],   chain,res[i]," ", x[i],  y[i], z[i], typ[i], charge[i]))
    file.close()


def createPDB(num_atoms, num_modes, x,y,z, modes, dof_modes):
    assert num_atoms == len(x)
    x1 = np.asarray(x)
    y1 = np.asarray(y)
    z1 = np.asarray(z)
    for mode in range(num_modes):
        for i in range(num_atoms):
            x1[i] += dof_modes[mode] * modes[mode][0][i]
            y1[i] += dof_modes[mode] * modes[mode][1][i]
            z1[i] += dof_modes[mode] * modes[mode][2][i]
    return x1, y1, z1


def getOverlap(id_mode, num_atoms, delta,modes):
    m = id_mode
    sum = 0
    for i in range(num_atoms):
        scalar= delta[0][i]* modes[m][0][i] + delta[1][i]* modes[m][1][i] + delta[2][i]* modes[m][2][i]
        abs = np.linalg.norm(np.asarray((delta[0][i],delta[1][i],delta[2][i]))) * np.linalg.norm(np.asarray((modes[m][0][i],modes[m][1][i],modes[m][2][i])))
        oi = np.fabs( scalar / abs)
        sum += oi**2

    return np.sqrt(sum/num_atoms)


def countResidues(num_atoms, name_residues, id_residues):
    assert num_atoms == len(name_residues)
    dict_residues = {}
    curr_res = 1000
    for id, res in zip(id_residues,name_residues):
        found = False
        if not res in dict_residues:
            curr_res = id
            dict_residues[res] = 1
        else:
            for key in dict_residues:
                if key == res and curr_res != id:
                    curr_res = id
                    dict_residues[key] += 1
                    found = True

    return dict_residues


def getPivot( num_atoms , pos_x,pos_y,pos_z):
    pivot = np.zeros(3)
    for x,y,z in zip(pos_x, pos_y, pos_z):
        pivot += np.asarray((x, y, z))
    pivot /= num_atoms
    return pivot




#id_atom, atom, cresn, resid, x_def,y_def,z_def, typ, charge= read_pdb("./test_compare/xyrr_modified_working.pdb")
#id_atom, atom, cresn, resid, x_def,y_def,z_def, typ, charge= read_pdb("./test_compare/oldattract/xyrca.pdb")

#num_atoms = len(coords)
#coord_deformed = np.zeros((3,num_atoms))
#for dim in range(3):
#    for i in range(num_atoms):
#        coord_deformed[dim][i] = coords[i][dim]
#num_atoms = len(x_def)



#id_atom, atom, cresn, resid, x,y,z, typ, charge = read_pdb("./test_compare/xyur.pdb")
#id_atom, atom, cresn, resid, x,y,z, typ, charge = read_pdb("./test_compare/oldattract/xyuca.pdb")

#coord_undeformed = np.zeros((3,num_atoms))
#for dim in range(3):
#    for i in range(num_atoms):
#        coord_undeformed[dim][i] = coords[i][dim]

#d_modes , num_modes= read_modes("./test_compare/modes_5m.dat")
#d_modes , num_modes= read_modes("./test_compare/oldattract/modes_5m.dat")



#modes = np.zeros((num_modes, 3, num_atoms))
#for mode in range(num_modes):
#    for dim in range(3):
  #      modes[mode][dim][:] = np.asarray(d_modes[mode +1 ][dim+1])

#dof_modes = compare( np.asarray((x_def,y_def,z_def)), np.asarray((x,y,z)), modes, num_atoms, num_modes )
#print dof_modes
#print dof_modes
#delta = getDelta( np.asarray((x_def,y_def,z_def)), np.asarray((x,y,z)),num_atoms)
#for i in range(num_modes):
    #print getOverlap(i, num_atoms, delta, modes)

#x_modes, y_modes, z_modes =  createPDB(num_atoms, num_modes, x,y,z, modes, dof_modes)
#writePDB(num_atoms,id_atom, atom, cresn, resid, x_modes,y_modes,z_modes, typ, charge,"test.pdb", chain = "" )


#dict_res =  countResidues(num_atoms, cresn, resid)

#print dict_res
#print sum(dict_res.values())
#print getPivot(num_atoms, x,y,z)
#print norm_modeVec(np.array([[4,5,6],[0]]), 1)

def evaluateModes(num_modes, filename_modes, filename_boundPdb, filename_unboundPdb, filename_dof_in,filenames_outPdb, filename_result):
    id_atom_b, atom_b, cresn_b, resid_b, x_b, y_b, z_b, typ_b, charge_b = read_pdb(filename_boundPdb)
    id_atom_u, atom_u, cresn_u, resid_u, x_u, y_u, z_u, typ_u, charge_u= read_pdb(filename_unboundPdb)#, typ_u, charge_u
    num_atoms = len(x_u)

    d_modes,num_modes_file = read_modes(filename_modes)
    num_modes_max = max(num_modes)
    overlap = np.zeros(num_modes_max+1)
    rmsd= np.zeros(num_modes_max+1)
    dof_modes = np.zeros(num_modes_max + 1)
    moderange= np.arange(0,num_modes_max+1)
    ev = np.zeros(num_modes_max + 1)
    modes = np.zeros((max(num_modes), 3, num_atoms))
    for mode in range(max(num_modes)):
        ev[mode+1] = d_modes[mode + 1][0]
        for dim in range(3):
            modes[mode][dim][:] = np.asarray(d_modes[mode + 1][dim + 1])

    dof_modes[1:] = compare(np.asarray((x_b, y_b, z_b)), np.asarray((x_u, y_u, z_u)), modes, num_atoms, num_modes_max)

    delta = getDelta(np.asarray((x_b, y_b, z_b)), np.asarray((x_u, y_u, z_u)), num_atoms)
    for i in range(num_modes_max):
        overlap[i+1] = getOverlap(i, num_atoms, delta, modes)
    rmsd[0] = os.popen("python2 $ATTRACTDIR/lrmsd.py {} {} {}".format(filename_dof_in, filename_unboundPdb,
                                                                          filename_boundPdb)).readlines()[0].split()[1]

    for i,mode in enumerate(num_modes):
        x_modes, y_modes, z_modes = createPDB(num_atoms, mode, x_u, y_u, z_u, modes, dof_modes[1:])
        writePDB(num_atoms, id_atom_b, atom_b, cresn_b, resid_b, x_modes, y_modes, z_modes, typ_b, charge_b,filenames_outPdb[i], chain="")
        rmsd[mode]= os.popen( "python2 $ATTRACTDIR/lrmsd.py {} {} {}".format(filename_dof_in, filenames_outPdb[i], filename_boundPdb)).readlines()[0].split()[1]


    df = pd.DataFrame()
    df["mode"] = moderange
    df["overlap"] = overlap
    df["rmsd"] = rmsd
    df["dof"] = dof_modes
    df["eigenvalue"] = ev
    df.to_csv(filename_result)



