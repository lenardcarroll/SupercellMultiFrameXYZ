import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser()
#Cheat sheet of arguments to be used with the script
parser.add_argument("-inp1", "--input1", dest = "input1", default = "VECTORS", help="Name of input file that contains only the scaling factor and the cell vectors")
parser.add_argument("-inp2", "--input2", dest = "input2", default = "0.xyz", help="Name of .xyz file that contains the cartesian coordinates")
parser.add_argument("-unit", "--unit", dest = "unit", default = "2x2x1", help="Size of desired unit cell along abc axes. Example, 2x2x1")
parser.add_argument("-out", "--output", dest = "output", default = "POSCAR", help="Name of your output file")
args = parser.parse_args()

#Determines the supercell dimensions
xpos = []
for i in range(len(args.unit)):
    if args.unit[i]=='x':
        xpos.append(i)
unit1 = int(args.unit[:xpos[0]])
unit2 = int(args.unit[xpos[0]+1:xpos[1]])
unit3 = int(args.unit[xpos[1]+1:])

#Read in vector file
f = open(args.input1,"r")
content = f.readlines()
#Get scaling factor
ScaleFactor = float(content[0])
#Get cell vector
Vectors = []
for i in range(1,4):
    Vectors.append(content[i].split())
VectorsAdj = []
#Multiply scaling factor/lattice constant with vector
for i in Vectors:
    VectorsAdj.append([float(i[0])*ScaleFactor,float(i[1])*ScaleFactor,float(i[2])*ScaleFactor])
#Find which atoms are fixed
FixedVals = content[4]
commapos = []
for i in range(len(FixedVals)):
    if FixedVals[i]==',':
        commapos.append(i)
FRange = []
if len(commapos)==0:
    x = FixedVals[:len(FixedVals)]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            FRange.append(i)
    else:
        if x!='None':
             FRange.append(int(x)-1)
else:
    for i in range(len(commapos)):
        if i==0:
            x = FixedVals[:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    FRange.append(i)
            else:
                FRange.append(int(x)-1)
        else:
            x = FixedVals[commapos[i-1]+1:commapos[i]]
            if '-' in x:
                xindex = x.index('-')
                firstval = int(x[:xindex])
                secondval = int(x[xindex+1:])
                for i in range(firstval-1,secondval):
                    FRange.append(i)
            else:
                FRange.append(int(x)-1)
if len(commapos)>0:
    x =  FixedVals[commapos[len(commapos)-1]+1:]
    if '-' in x:
        xindex = x.index('-')
        firstval = int(x[:xindex])
        secondval = int(x[xindex+1:])
        for i in range(firstval-1,secondval):
            FRange.append(i)
    else:
        FRange.append(int(x)-1)

#Read in coordinates of multi-frame XYZ file
with open(args.input2) as f:
    numberofatoms = int(f.readline().rstrip())
g = open(args.input2,"r")
content = g.readlines()
#determines how many frames the file has
numberofrows = int(len(content)/numberofatoms)

#Writes output file
h = open(args.output,"w")
h.close()

#Goes through each frame
for j in range(numberofrows):
    AtomCoords = []
    Atoms = []
    for i in range(j*(numberofatoms+2)+2,j*(numberofatoms+2)+numberofatoms+2):
        Atoms.append(content[i].split()[0])
        AtomCoords.append([float(content[i].split()[1]),float(content[i].split()[2]),float(content[i].split()[3])])
    
    #converts frame to POSCAR file
    h = open("POSCAR","w")
    print("XYZ to POSCAR",file=h)
    print("%0.16f" % ScaleFactor,file=h)
    for i in Vectors:
        print("%0.16f" % float(i[0]),"%0.16f" % float(i[1]),"%0.16f" % float(i[2]),file=h)
    AtomList = []
    AtomNum = []
    for i in range(len(Atoms)):
        if i==0:
            AtomList.append(Atoms[i])
            AtomNum.append(1)
        else:
            if Atoms[i] != AtomList[len(AtomList)-1]:
                AtomList.append(Atoms[i])
                AtomNum.append(1)
            else:
                AtomNum[len(AtomNum)-1] += 1
    for i in AtomList:
        print(i, end=" ",file=h)
    print("",file=h)
    for i in AtomNum:
        print(i, end=" ",file=h)
    print("",file=h)
    print("Selective dynamics",file=h)
    print("Cartesian",file=h)

    AtomTot = 0
    for i in AtomNum:
        AtomTot+=i

    for i in range(AtomTot):
        MatrixMult = AtomCoords[i]
        if i in FRange:
            print('%0.16f' % MatrixMult[0], '%0.16f' % MatrixMult[1],'%0.16f' % MatrixMult[2],'F','F','F',file=h)  
        else:
            print('%0.16f' % MatrixMult[0], '%0.16f' % MatrixMult[1],'%0.16f' % MatrixMult[2],'T','T','T',file=h)  
    h.close()
    #Convert POSCAR file to .cif file
    import os
    os.system('vasp2cif.py POSCAR')
    #Read in .cif file and make a supercell
    from pymatgen.core.structure import Structure
    structure = Structure.from_file('POSCAR.cif') 
    structure.make_supercell((unit1,unit2,unit3))
    #Write the supercell to a file called test
    h = open(args.output,"a")
    print(len(structure),file=h)
    print("Unit cell %dx%dx%d" % (unit1,unit2,unit3),file=h)
    g = open("test","w")
    print(structure,file=g)
    #Read in test file and get the coordinates
    with open('test') as f:
        lines_after_7 = f.readlines()[7:]
    atoms = []
    data_x = []
    data_y = []
    data_z = []
    for line in lines_after_7:
        line=line.strip()
        columns=line.split()
        data_x.append(float(columns[2]))
        data_y.append(float(columns[3]))
        data_z.append(float(columns[4]))
        atoms.append((columns[1]))
    #Save supercell coordinates to the main file
    for i in range(len(structure)):
        AtomCoordinates = []
        AtomCoordinates = [data_x[i],data_y[i],data_z[i]]
        newVectorsAdj = [[VectorsAdj[0][0]*unit1,VectorsAdj[0][1]*unit2,VectorsAdj[0][2]]*unit3,[VectorsAdj[1][0]*unit1,VectorsAdj[1][1]*unit2,VectorsAdj[1][2]]*unit3,[VectorsAdj[2][0]*unit1,VectorsAdj[2][1]*unit2,unit3*VectorsAdj[2][2]]]
        MatrixMulti = np.matmul(AtomCoordinates,(newVectorsAdj))
        print(atoms[i],MatrixMulti[0],MatrixMulti[1],MatrixMulti[2],file=h)
    h.close()
    g.close()
