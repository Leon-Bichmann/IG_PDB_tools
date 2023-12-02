#!/usr/bin/env python
"""
This script allows to retrieve epitope residues of IGs from the CoV-AbDab DB
"""
__author__="Leon Bichmann"


#import python libraries
from Bio.PDB import *
import pandas as pd
from collections import Counter 


#read and filter CoV-AbDab database for single entry pdbs
df=pd.read_csv("../../Downloads/CoV-AbDab_130623.csv",sep=",")
df=df[(df["Structures"].str!="ND") & (df["Structures"].str.contains("http"))]
df=df[~df["Structures"].str.contains(";")]
df=df[["VHorVHH","CDRH3","VL","CDRL3","Structures"]].dropna()


#initialize Bio.PDB classes
pdbl = PDBList()
parser = PDBParser()
ppb = CaPPBuilder()


#iterate through CoV-AbDab DB
for row in df.iterrows():
    cdrl3=row[1]["CDRL3"]
    print(cdrl3)
    cdrh3=row[1]["CDRH3"]
    pdbid=str(row[1]["Structures"]).split("/")[-1]
    pdb_struct=pdbl.retrieve_pdb_file(pdbid,file_format="pdb")
    structure = parser.get_structure(pdbid,pdb_struct)[0] #keep first structure in case of ensemble
    atoms = Selection.unfold_entities(structure, "A")

    #define VH, VL and CDR chains and residue positions for each structure
    #keep last in case of multiple symmetry structures
    chains=structure.get_chains()
    for c in chains:
    	#retrieve sequence using polypeptide builder
        pps=ppb.build_peptides(structure[c.id])
        seq="".join([str(pp.get_sequence()) for pp in pps])
        if cdrl3 in seq:
            print("VL chain: " + str(c.id))
            vl=c.id
            vl_seq=seq
            vl_cdr3_idx=vl_seq.find(cdrl3) # estimate vl_cdr3 start index position
        if cdrh3 in seq:
            print("VH chain: " + str(c.id))
            vh=c.id
            vh_seq=seq
            vh_cdr3_idx=vh_seq.find(cdrh3) # estimate vh_cdr3 start index position

    #retrieve neighbourhood residues
    ns = NeighborSearch(atoms)
    neighbour_rs=[]
    radius=8 # angstrom
    
    #VH CDR3 neighbourhood
    for i in range(vh_cdr3_idx-1,vh_cdr3_idx+1+len(cdrh3)):
        vh_as=structure[vh][i].get_atoms()
        for a in vh_as:
            for na in ns.search(a.coord,radius):
                if na.get_parent().get_parent().id!=vl and na.get_parent().get_parent().id!=vh:
                    neighbour_rs.append(na.get_parent())

    #VL CDR3 neighbourhood    
    for i in range(vl_cdr3_idx-1,vl_cdr3_idx+1+len(cdrl3)):
        vl_as=structure[vl][i].get_atoms()
        for a in vl_as:
            for na in ns.search(a.coord,radius):
                if na.get_parent().get_parent().id!=vl and na.get_parent().get_parent().id!=vh:
                    neighbour_rs.append(na.get_parent())

    #Identify antigen chain based on majority vote                
    antigen_chain = Counter([r.get_parent().id for r in neighbour_rs]).most_common(1)[0][0]
    antigen_epitope_reslist=[]
    for r in neighbour_rs:
        if r.get_parent().id==antigen_chain:
            antigen_epitope_reslist.append(r.get_resname()+"_"+str(r.id[1])+antigen_chain)
 

    #print neighbourhood residue list
    print(list(set(antigen_epitope_reslist)))
    break # remove stop for parsing the whole db
