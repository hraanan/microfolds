import Bio
import os
import sys
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
import math
import numpy
from collections import Counter
import random 
from Bio.PDB import *
import gzip
def get_center(res_list):
    coord = []
    
    for atom in residue:
        #   print(atom.coord)
           at=atom.coord
           x=at[0]
           y=at[1]
           z=at[2]
           atcord=[x,y,z]
           coord.append(atcord)
    x=0
    y=0
    z=0
    i=0
    for point in coord:
        i=i+1
        x=x+point[0]
        y=y+point[1]
        z=z+point[2]
    x=x/i
    y=y/i
    z=z/i
    center=numpy.array([x,y,z])    
    return center;









pdbl=PDB.PDBList()
Error_out=open("microfolds_out.txt","w")

cng=0
fd=['1fdn']
AA=['PHE','TRP','TYR','ALA','CYS','ASP','GLU','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL']
CF=[' DA',' DC',' DG',' DT','  A','  C','  G','  U','HOH','UNK','UNX']
Metals=['FE','MN','CU','CD','OS','CO','NI','W','PT','MO','U','TA','V','AU','IR','Y','GD','RU','YB','SM','PD','AG','EU','RH','PR','RE','LU','TB','HF','HO','DY','ZR','CR','LA','CE','ER','AM','CM','TH','PU','SC','PA']
cofactor=['BCB','CLA','CHL','BCL','CL0','PMR','PHO']
idxfile='entries.idx.cont.txt'
idx=open(idxfile,"r")
idx.readline()
idx.readline()
idx.readline()
EC=""
for line in idx:
    
    try:
        
        protein=line[0:4]
#        print ('pdb_code:'+protein)
        #protein='3kzi' 
        #print ('pdb_code:'+protein)
        protein=protein.lower()
        Error_out.write('pdb_code:'+protein+'\n')

        parser = PDB.PDBParser(PERMISSIVE=1,get_header=1,QUIET=1)
        curdir=os.getcwd()
        filename='/home/hraanan/pdb_download/'+protein[1:3]+'/pdb'+protein+'.ent.gz'
        final_file='/home/hraanan/pdb_download/'+protein[1:3]+'/pdb'+protein+'.ent'
        #print (protein,'/home/hraanan/pdb_download/'+protein[1:3]+'/pdb'+protein+'.ent.gz')
        #print ('unziping')
        gz = gzip.open(filename, 'rb') 
        with open(final_file, 'wb') as out: 
            out.writelines(gz) 
        gz.close()
        #print ('unziping done')
        #os.remove(filename)
        structure = parser.get_structure(protein,'/home/hraanan/pdb_download/'+protein[1:3]+'/pdb'+protein+'.ent')
        head= structure.header['head']
        comp = structure.header['compound']
        EC==""
        try:
            comp=comp['1']
#        except KeyError:
#            try:
#                EC=comp['ec_number']
#            except KeyError:
#                try:
#                    EC=comp['ec']
        except KeyError:
                    EC='-.-.-.-'
        try:
            EC=comp['ec']
        except KeyError:
            pass
        try:
            EC=comp['ec_number']
        except KeyError:
            pass
        if EC=="": 
            EC='-.-.-.-'
        #print(EC)
##
##        
               
        sf4ID=[]
        sf4coord=[]
        for model in structure:
            if model.id==0:
             atom_list = Selection.unfold_entities(model, 'A') # A for atoms
             ns = NeighborSearch(atom_list)
             lig=[]
             for chain in model:
                 for residue in chain:
                     #if residue.resname not in AA and residue.resname not in CF :
                      #print(chain.id,residue.resname)
                      if residue.resname in cofactor:                       
                       # print(chain.id,residue.resname)
                        atom_in_res=[]
                        for atom in residue:
                            atom_in_res.append(atom.element)
                           
                        #if any(x in Metals for x in atom_in_res)==False:
                            #print ('not metal')
                         #   continue
                            
                        center = get_center(residue)
                        #print ('center',center)
                        lig=protein,chain.id,residue.id[1],residue.resname,center
                        #print(lig)
                        all_neighbors = ns.search(center, 15.0,"R") # 15.0 for distance in angstrom
                        microfold_name=protein+'.'+residue.resname+'_'+ chain.id +'_'+str(residue.id[1])+'_'+head+'_'+EC
                        microfold_name=microfold_name.replace(' ','')
                        microfold_dir=residue.resname
                        microfold_dir=microfold_dir.replace(' ','')
                        print(microfold_name)
                        if not os.path.exists('/home/hraanan/MicrofoldsPDBs/'+microfold_dir):
                            os.makedirs('/home/hraanan/MicrofoldsPDBs/'+microfold_dir)
                        Select = Bio.PDB.Select
                        class MicroSelect(Select):
                           def accept_residue(self, residue):
                               if residue in all_neighbors and residue.resname!='HOH':
                                   return 1
                               else:
                                   return 0
                        io=PDBIO()
                        io.set_structure(structure)
                        #print('/home/hraanan/MicrofoldsPDBs/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())                        
                        io.save('/home/hraanan/MicrofoldsPDBs/'+microfold_dir+'/'+microfold_name+'.pdb', MicroSelect())
    except:
        e = sys.exc_info()[0]
        Error_out.write('xxx\n')
        Error_out.write( "<p>Error: %s</p>" % e )
        Error_out.write('xxx\n')
       
        continue
                               
    
Error_out.close()
#prot.close()
print("end")
