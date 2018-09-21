import prody as dy
import matplotlib.pyplot as plt
import os


path = "/home/glenn/Documents/WebNma/best"

names = []
file =open(os.path.join(path, "dir.txt"))

for line in file.readlines():
    names.append(line.strip("\n"))
    print line
file.close()
for name in names:
    protPath = os.path.join(path, name)
    proteinName = ["{}A-unbound".format(name),"{}B-unbound".format(name) ]
    for prot in proteinName:
        inPDB = os.path.join(protPath, prot+".pdb")
        print inPDB
        outSTRIDE = os.path.join(protPath, prot)
        outSecondary = os.path.join(protPath,prot+".sec")
        dy.execSTRIDE(inPDB, outSTRIDE)
        pdb = dy.parsePDB(inPDB)
        dy.parseSTRIDE(outSTRIDE+".stride", pdb)
        secondaryStruc = pdb.getSecstrs()
        resIdx = pdb.getResindices()
        secFile = open(outSecondary, 'w+')
        curres = -1
        for sec,res in zip(secondaryStruc, resIdx):
            if res != curres:
                secFile.write("{}\n".format(sec))
                curres = res
        secFile.close()




names=['3MXW']
for name in names:
    proteinName = ["{}A-unbound".format(name), "{}B-unbound".format(name)]
    for prot in proteinName:
        fileSecStruc = "{}/{}/{}.sec".format(path,name,prot)
        fileDeformation = "{}/{}/webnma/{}_fluctuationsplot.dat".format(path,name,prot)
        secFile = open(fileSecStruc)
        sec =[]
        for line in secFile.readlines():
            sec.append(line.split()[0])
        secFile.close()

        defFile = open(fileDeformation)
        defo = []
        for line in defFile.readlines():
            defo.append(float(line.split()[1]))
        plt.clf()
        plt.plot(sec)
        plt.plot(defo)
        plt.show()