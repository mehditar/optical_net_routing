import numpy as np
import pandas as pd
import csv
import time

from trans_store_clean import *
from demand_update import one_demand_update



# reading the files
def file_reader(file_name):
    reader = csv.reader(open(file_name, "r"), delimiter=",")
    x = list(reader)
    return np.array(x).astype("float")

[PATHS, INFO, PATHORDER, LINKNAME] = ADR()
NETWORKDIS = file_reader("csvND.dat") 
NETWORKLINK = file_reader("csvNL.dat")
DEMANDSEQ =  file_reader("QN.dat")
TOTALDEMANDS=np.size(DEMANDSEQ,1)
par = file_reader("Par.dat")


# NETWORk config
PARAM = [0] * 12
PARAM[1] = par[0, 1]
PARAM[9] = 2 # 2:Heuristic, 1:ILP
PARAM[8] = 1 # algorith---original:1, no reduction:2, mimimumbitrate:3
PARAM[7] = 50 # of transcievers per fibler link
PARAM[3] = 1 # Ac
PARAM[0] = 50 # of slices per fible
PARAM[4] = 50 # ms Total
PARAM[5] = 8 # ms 
PARAM[6] = 1 # ms replace
PARAM[2] = 50 # total trx.(ILP Exlusive PARAMeters)
PARAM[10] = 0 # Zeta: a number from [0,1], to apply the 
              # worth of the routes. (Heuristic Specific)
PARAM[11] = 7 # k first paths-considered.(Heuristic Specific)


#creation of network_state_storage_arrays/variables
active = np.zeros((1000, 11))
passive = np.zeros((1000, 8))
served = {}
served[1] = 0
served[2] = set()
served[3] = {}
key_dict = {}
dem_id_dict ={}
state = [0, 0, 0, -1, 0, 0, 0, 0]

#flags
mxws = min(PARAM[2], PARAM[7], PARAM[0]) - 1
last_demand = 0
last_time = 0
flag_1 = 1
flag_2 = 1



#-------Network_operation_starts

while (state[6] + flag_1 + flag_2) != 0:
    
    [active, passive] = DBR(active, passive, state, last_time, mxws)
    reduction = 0
  
   

   
    if state[1] == 0:
        Id = state[2]
     
        [active, passive,reduction] = one_demand_update(active, passive, DEMANDSEQ[:, Id],
                                          state, served, key_dict,dem_id_dict, 1, 0, 1,
                                          DEMANDSEQ[9, Id], NETWORKDIS, NETWORKLINK,PARAM,
                                          PATHS, INFO, PATHORDER, LINKNAME)


        
    #endiftype0
    if state[1] == 1 :
        Sid_ash = Finder(active[0: state[6], 0: 2], 0, 1, state[2])
        active = Finder3(active, state[6], Sid_ash)
        state[6] = state[6] - 1
        served[1] = served[1] - 1
        served[2] = served[2] - {state[2]}
        reduction = 1
        del served[3][state[2]]
        DEMANDSEQ[4, state[2]] = 1
        DEMANDSEQ[5, state[2]] = state[0]
        
  
    if reduction == 1:
        if state[6] != 0:
            for Sid in range(state[6]) :
                state[2] = int(active[Sid, 1])
                Id = state[2]
                [active, passive,reduction] = one_demand_update(active, passive, DEMANDSEQ[:, Id],
                                                  state, served, key_dict, dem_id_dict, 0,
                                                  Sid, 0, active[Sid, 10], NETWORKDIS, NETWORKLINK,
                                                  PARAM, PATHS, INFO, PATHORDER, LINKNAME)

           
        ind = 0
        q_size = state[7]
        flag_2 = 1
        while True:
            if state[7] == 0:
                break
            state[2] = int(passive[ind, 1])
            Id = state[2] 
            [active, passive, reduction] = one_demand_update(active, passive, DEMANDSEQ[:, Id],
                                               state, served, key_dict, dem_id_dict, 1,
                                               0, 0, passive[ind, 7], NETWORKDIS, NETWORKLINK, 
                                               PARAM, PATHS, INFO, PATHORDER, LINKNAME)
            if Id == passive[ind, 1]:
                ind = ind + 1
            if ind == state[7]:
                break
        
        
        if q_size == state[7]:
             flag_2 = 0
        
    last_time = state[0]          
    if last_demand == TOTALDEMANDS - 1:
        [state, b] = NDF(active, state, "haha", 1)
        flag_1 = 0
    else:
        [state, b] = NDF(active, state, DEMANDSEQ[:, last_demand + 1], 0)
        last_demand = last_demand + b
        
        
#---end of network operation
# calculating the total served demands ratio
a = sum(DEMANDSEQ[4,:]) 
b = sum(DEMANDSEQ[4,:]*DEMANDSEQ[3,:])/sum(DEMANDSEQ[3,:])
c = np.mean((DEMANDSEQ[5,:]-DEMANDSEQ[1,:])*DEMANDSEQ[3,:])
d = np.mean((DEMANDSEQ[2,:]-DEMANDSEQ[1,:])*DEMANDSEQ[3,:])




    
  
                 


    
    

