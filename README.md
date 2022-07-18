This script makes a supercell of the structures which contains multiple frames.

Each XYZ file frame is converted into a POSCAR file, then the POSCAR file converted to a .cif file. Pymatgen is then used to make the supercell from the .cif file, coordinates read from the supercell and saved to a .xyz file.

The python script works with:

```
python -i SCMFXYZ.py -inp1 <VECTORFILE> -inp2 <XYZ FILE> -u <SUPERCELL INFO> -out <OUTPUT XYZ FILE>
```

for example:

```
python -i SCMFXYZ.py -inp1 Vectors.txt -inp Structure.xyz -u 2x2x1 -out Output.xyz
```

The vector file must look like:

```
 1.00000000000000     
 24.6914512445513985  0.0000000000000000  0.0000000000000000
  7.4065517955292197 12.8319727856189001  0.0000000000000000
  0.0000000001313560  0.0000000000558372 10.0000000000000000
1-20,24,25,29,35-74
```

The first line should be the scaling factor/lattice constant, lines 2-4 the cell vectors and the last line all the atoms that are frozen. If none are frozen, add 'None' at the end.

