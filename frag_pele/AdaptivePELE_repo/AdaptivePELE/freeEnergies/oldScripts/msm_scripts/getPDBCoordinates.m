%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple matlab script to extract PDB coordinates      %
%% Written by D. Lecina                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, Rlig] = getPDBCoordinates(pdbfilename, ligandName)

pdb = pdbread(pdbfilename);

proteinAtoms = pdb.Model.Atom;
numberOfProteinAtoms = size(pdb.Model.Atom,2);

x = [];
y = [];
z = [];

for i = 1:numberOfProteinAtoms
    x = [x; proteinAtoms(i).X];
    y = [y; proteinAtoms(i).Y];
    z = [z; proteinAtoms(i).Z];
end

heteroAtoms = pdb.Model.HeterogenAtom;
numberOfHeteroAtoms = size(pdb.Model.HeterogenAtom,2);

xlig = [];
ylig = [];
zlig = [];

for i = 1:numberOfHeteroAtoms
    if ~ strcmpi(heteroAtoms(i).resName, ligandName)
        x = [x; heteroAtoms(i).X];
        y = [y; heteroAtoms(i).Y];
        z = [z; heteroAtoms(i).Z];
    else
        xlig = [xlig; heteroAtoms(i).X];
        ylig = [ylig; heteroAtoms(i).Y];
        zlig = [zlig; heteroAtoms(i).Z];
    end
end

R = [x y z];
Rlig = [xlig ylig zlig];
