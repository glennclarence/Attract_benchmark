

import numpy as np
def read_pdb(f):
  coor, res, resn, atom = [], [], [], []
  curr_resid = None
  resindex = 0
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    coor.append((x,y,z))
    resid = l[21:26]
    cresn = l[17:20]
    if resid != curr_resid:
      curr_resid = resid
      resindex += 1
    res.append(resindex)
    resn.append(cresn)
    atom.append(l[12:16])
  return coor, res, resn, atom




def read_modes(filename_modes):

    with open(filename_modes, 'r') as f:
        lines = f.readlines()
        d_modes = {}
        mode_idx = 1
        for i, line in enumerate(lines):
            if len(line.split()) == 2 and int(line.split()[0]) == mode_idx:
                x = []
                y = []
                z = []
                d_modes[mode_idx ] = ( float(line.split()[1]), x,y,z )
                mode_idx += 1
                #eigenvalues[mode_idx - 2] = float(line.split()[1])
            elif len(line.split()) == 3:
                d_modes[mode_idx - 1][1].append( float(line.split()[0]) )
                d_modes[mode_idx - 1][2].append( float(line.split()[1]) )
                d_modes[mode_idx - 1][3].append( float(line.split()[2]) )

    return d_modes

def safeModeFile( filename_mode, d_modes):
    with open(filename_mode, 'w') as f:

        for key, value in d_modes.iteritems():
            line = "  {} {}\n".format(str(key), str(value[0]))
            f.write(line)
            for idx in range(len(value[1])):
                line = "{} {} {}\n".format(value[1][idx],value[2][idx],value[3][idx])
                f.write(line)
            line="\n"
            f.write(line)
        f.close()

def mode_transpose(file_inputModes, file_refPdb, file_targetPdb, file_outputModes):
    residuesTarget  = read_pdb(file_targetPdb)[1]
    residuesRef = read_pdb(file_refPdb)[1]
    d_modesRef = read_modes(file_inputModes)
    d_modesTarget = {}
    for key, value in d_modesRef.iteritems():
        x = np.zeros(len(residuesTarget))
        y = np.zeros(len(residuesTarget))
        z = np.zeros(len(residuesTarget))
        d_modesTarget[key] = ( value[0], x,y,z )

    for i, resTar in enumerate(residuesTarget):
        modeIdx = residuesRef.index(resTar)
        for key, value in d_modesTarget.iteritems():
            d_modesTarget[key][1][i] = d_modesRef[key][1][modeIdx]
            d_modesTarget[key][2][i] = d_modesRef[key][2][modeIdx]
            d_modesTarget[key][3][i] = d_modesRef[key][3][modeIdx]

    safeModeFile(file_outputModes,  d_modesTarget)



def test( file_testModeFile, file_referenceModeFile):
    test= read_modes(file_testModeFile)
    ref = read_modes(file_referenceModeFile)
    for (keyTest, valTest), (keyRef, valRef) in zip(test.iteritems(), ref.iteritems()):
        if keyTest != keyRef:
            print "wrong key ", keyTest, keyRef
        if valTest[0] != valRef[0]:
            print "wrong eigenvalue", valTest[0], valRef[0]
        if len(valRef[1]) != len(valTest[1]):
            print "wrong value length", len(valRef[1]), len(valTest[1])

        for idx in range(len(valTest[1])):
            if valTest[1][idx] != valRef[1][idx]:
                print "wrong x value", valTest[1][idx], valRef[1][idx]

            if valTest[2][idx] != valRef[2][idx]:
                print "wrong y value", valTest[2][idx], valRef[2][idx]

            if valTest[3][idx] != valRef[3][idx]:
                print "wrong z value", valTest[3][idx], valRef[3][idx]




#path = "/home/glenn/test_Mode_transpose/"
#mode_transpose(path + "1AVX-ligand-for-docking-5-modes.dat", path + "1AVX-ligand-for-docking-reduce.pdb", path + "1AVX-ligand-for-docking-heavy.pdb", path + "test")
#test( path + "test", path + "referenceModes.dat")