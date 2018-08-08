# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 12:24:39 2017

@author: hraanan
"""

#
#%pylab inline
from IPython.display import HTML

from pypdb.pypdb import *

import pprint



out_file=open('list_of_pdbs_for_organic_new_center.txt','w')

organic_cofactors_list=[]#['BPB','BPH','68G','CL7','RCC']#['BCB','CLA','CHL','BCL','CL0','PMR','PHO','BPB','BPH','68G','CL7','RCC']#['4YP','BHQ','BNT','DCQ','DMW','DQN','EMO','FLV','KIA','NQ','PL9','PLQ','PQ9','PQN','PQQ','TPQ']
organic_cofactors_pdb_list=[]
organic_cofactors_pdb_file=open('manual_cofactor_list_with_quinone.txt','r')
for line in organic_cofactors_pdb_file:
    line=line.split('\t')
    organic_cofactors_list.append(line[1])

for lig in organic_cofactors_list:
    print(lig)
    search_dict = make_query(lig,'ChemCompIdQuery')
    found_pdbs=[]
    found_pdbs = do_search(search_dict)
    print(lig+str(len(found_pdbs)))     
    organic_cofactors_pdb_list= organic_cofactors_pdb_list + list(set(found_pdbs) - set(organic_cofactors_pdb_list))    
       
    #print(len(organic_cofactors_pdb_list))
out_file.writelines(["%s\n" % item  for item in organic_cofactors_pdb_list])
out_file.close()