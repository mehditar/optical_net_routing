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

NETWORKDIS = file_reader("csvND.dat") # this is a matrix that defines the 
# length between each two nodes if they are connected dirrectly through fiber cables.
NETWORKLINK = file_reader("csvNL.dat") # the same matrix as in the previous part
# here the the files indicates the existenxe of a 
# direct path between each source and destination

[PATHS, INFO, PATHORDER, LINKNAME] = ADR() # this function imports the four files.
# the files combined contain the 50 shortest paths between each pair 
#of the nodes in the network. I used MATLAB to calculte these paths

DEMANDSEQ =  file_reader("QN.dat")
TOTALDEMANDS=np.size(DEMANDSEQ,1)
par = file_reader("Par.dat")


# NETWORk configuration values
PARAM = [0] * 12
PARAM[1] = par[0, 1]
PARAM[9] = 2 # 2:Heuristic, that is if we want to use the heuristic algorithm
# to conduct the routing or 1 if we want the ILP to do so
PARAM[8] = 1 # algorith---original:1, no reduction:2, mimimumbitrate:3 this 
# number changes the beahviour of the heuristic algorithm and basically 
# provides three different algorithm
PARAM[7] = 50 # maximum number of transceivers that can be allocated  to 
#any particiliar port of fiber
PARAM[3] = 1 # Ac: a number to repeat the search if the search is nut successful in a try
PARAM[0] = 50 # number of frequency slices per fible. it should be noted that if a link
# exists between two nodes we assume there are two fiber between them, each in one direction
PARAM[4] = 50 # ms Total
PARAM[5] = 8 # ms 
PARAM[6] = 1 # ms replace
PARAM[2] = 50 # total trx.(ILP Exlusive PARAMeters). total number of the transceivers that
# are can exist in any nodes of the network.
PARAM[10] = 0 # Zeta: a number from [0,1], to apply the 
              # worth of the routes. (Heuristic Specific)
PARAM[11] = 7 # k first paths-considered.(Heuristic Specific)


#creation of network_state_storage_arrays/variables
active = np.zeros((1000, 11)) # this matrix will keep the information about the demands that are
# active (i.e,being served) including the data rates arrical time,
#the number of lighpaths that are allocated to each
#demand and the current data rate that they have 
passive = np.zeros((1000, 8)) # this matrix does a similiar job for the
# demands that have arrived to the network but have
# not received any recources and have been put in line for resource allocation
served = {} # this dictionary will contains all the detials about any 
# demand resources. we need these information to be able to release the 
# resources when the demand is served completely and the connection is terminated
served[1] = 0
served[2] = set()
served[3] = {}
key_dict = {} # this is another torage similiar to served
dem_id_dict ={} # this is another torage similiar to served
state = [0, 0, 0, -1, 0, 0, 0, 0] # this contains some general information about the networks
# such as the number of demeand that are receiving recources, or are in the line for
# for resource allocation, the curren time, the curren demand id that is being considered

#flags
mxws = min(PARAM[2], PARAM[7], PARAM[0]) - 1 # the parameter makes sure the 
# that the different parameters in the algorithm about the
# number of the transcieevrs mach. that is the minimum 
last_demand = 0
last_time = 0
flag_1 = 1 # becomes zero when the there is no active demands
flag_2 = 1 # becomes zero when there is no demands in line (passive demands)



#-------Network_operation_starts

while (state[6] + flag_1 + flag_2) != 0:# this is the beggining of the interatin. the iteration 
    # cintinues untill all the demands are taken into account and are all either served or rejected
    
    [active, passive] = DBR(active, passive, state, last_time, mxws) # this function updates 
    # active and passive after each iteration and depeneding on weather
    # the search for resource was successful or not. the demands that have been in line for resource
    # allocaton and their end time has arrive are simple rejected and removed from the lists
    reduction = 0 # this variable indicates wether any resource allocation resulted in some resource 
    # being release (reduced) or not. if yes the search to fill those released 
    # functions is conducted by trying to find additional paths for the 
    # active paths or first path for the demands that are in line for resource allocation
  
   

   
    if state[1] == 0: # thi indicates that the demand that is being considered
        # is a new demand that has just recently arrived to the network
        Id = state[2]
     
        [active, passive,reduction] = one_demand_update(active, passive, DEMANDSEQ[:, Id],
                                          state, served, key_dict,dem_id_dict, 1, 0, 1,
                                          DEMANDSEQ[9, Id], NETWORKDIS, NETWORKLINK,PARAM,
                                          PATHS, INFO, PATHORDER, LINKNAME)
        # here one compplete reouting search is conducted using one_demand_update


        
    #endiftype0
    if state[1] == 1 : # execution of this part indicates that the demand was served completely
        # and here the recources that were allocated to the demand are released 
        Sid_ash = Finder(active[0: state[6], 0: 2], 0, 1, state[2]) # this simpy find the id if the demand 
        # in the active table
        active = Finder3(active, state[6], Sid_ash) # this function removes the demand from active table
        state[6] = state[6] - 1 # update and served are also updated
        served[1] = served[1] - 1
        served[2] = served[2] - {state[2]}
        reduction = 1
        del served[3][state[2]]
        DEMANDSEQ[4, state[2]] = 1 # here this value, previously zero, receives 1 as the new value to
        # indicates that the demand completely served
        DEMANDSEQ[5, state[2]] = state[0]
        
  
    if reduction == 1:
        if state[6] != 0:
            for Sid in range(state[6]) : # in this section for any active demand, as
                # earliar the algorithms seaches to find additional paths and dedicates
                # in other fords for each demand from top to bottom of the active demand table
                # there will be a routing search to fill the reources that became free (reduction=0)
                state[2] = int(active[Sid, 1])
                Id = state[2]
                [active, passive,reduction] = one_demand_update(active, passive, DEMANDSEQ[:, Id],
                                                  state, served, key_dict, dem_id_dict, 0,
                                                  Sid, 0, active[Sid, 10], NETWORKDIS, NETWORKLINK,
                                                  PARAM, PATHS, INFO, PATHORDER, LINKNAME)

           
        ind = 0
        q_size = state[7]
        flag_2 = 1
        while True: # a similiar search is conducted for each in-line demand in the passive list
            # to find for new orutes
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
            if ind == state[7]: # at the end of the table the search  is complete. 
                break
        
        
        if q_size == state[7]: # if all the demand are given at least one route 
             # and passive list become empty.
             flag_2 = 0
        
    last_time = state[0]          
    if last_demand == TOTALDEMANDS - 1: # this is just to keep trck of the last and the next demand 
        # that has been and will be considered based on the arrival time
        [state, b] = NDF(active, state, "_", 1)
        flag_1 = 0
    else:
        [state, b] = NDF(active, state, DEMANDSEQ[:, last_demand + 1], 0)
        last_demand = last_demand + b
        
        
#---end of network operation
# calculating the total served demands ratio
a = sum(DEMANDSEQ[4,:]) # this calculated the tial demand that were served completely.
b = sum(DEMANDSEQ[4,:]*DEMANDSEQ[3,:])/sum(DEMANDSEQ[3,:]) # this calcuates the ratio of demands
# that were served to all demands.
c = np.mean((DEMANDSEQ[5,:]-DEMANDSEQ[1,:])*DEMANDSEQ[3,:])
d = np.mean((DEMANDSEQ[2,:]-DEMANDSEQ[1,:])*DEMANDSEQ[3,:])




    
  
                 


    
    

