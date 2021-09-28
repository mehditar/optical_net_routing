import copy
from integer_linear_programming import router_ILP
from heuristic import router_heuristic
from trans_store_clean import *
def one_demand_update(active, passive, demand, state, served, key_dict,
        dem_id_dict, resrclss_dem, Sid, is_arrvd, bitrate_min,
        NETWORKDIS, NETWORKLINK, PARAM,
        PATHS, INFO, PATHORDER, LINKNAME):
    

    
    repeat = True
    reduction = 0
    Id = state[2]
    change = 0
    attempts = 0
    while (repeat == True) & (attempts < int(PARAM[3])):
  
            attempts += 1
            input_table = TP(served, key_dict, active, state,
                             PARAM, dem_id_dict)
            state[4] = input_table[5]
            
            if PARAM[9] == 1:
                
               search_res = router_ILP([input_table[0], input_table[4]], input_table[1],
                         input_table[2], input_table[3], served, key_dict,
                         demand[6], demand[7], 1-resrclss_dem, Sid, bitrate_min, state,
                         dem_id_dict, NETWORKDIS, NETWORKLINK, PARAM, PATHS,
                         INFO, PATHORDER, LINKNAME)
            else:
                
               search_res = router_heuristic([input_table[0], input_table[4]], input_table[1], 
                          input_table[2], input_table[3], served, key_dict,
                          demand[6], demand[7], 1 - resrclss_dem, Sid, bitrate_min, state,
                          dem_id_dict,NETWORKDIS, NETWORKLINK, PARAM, PATHS, INFO,
                          PATHORDER, LINKNAME)
            
            repeat = search_res[0]

            
            
            if search_res[0] == True:
               state[4] = search_res[1]
               state[3] = state[3] + len(search_res[2])
               
            
               if resrclss_dem == 1:
                   served = USC(served, Id, search_res)
                   for S in search_res[2]:
                       key_dict[S] = [Id, 0]
               if resrclss_dem == 0:
                   [served, active] = USA(served, active, state, Id, search_res)
                   for S in search_res[2]:
                       key_dict[S] = [Id, served[3][Id][4]]
           
        
        
        
               if search_res[7] == 1:
                  [served, active] = USAR(served, key_dict, active, state, search_res)
                  reduction = 1

               if resrclss_dem == 1:
                   Info = [Id, demand[8], change, demand[1], demand[2], demand[3],
                         search_res[6], demand[3] / search_res[6] + state[0], len(search_res[2]),
                         demand[3] / (demand[2] - state[0])]
                   active = Finder2_2(active, state[6], Info)
                   state[6] = state[6] + 1
                   Sid = Finder(active[0: state[6], 0: 2], 0, 1, Id)
                   resrclss_dem = 0
                   if is_arrvd==0:
                      Sid_ash = Finder(passive[0: state[7], 0: 2], 0, 1, Id)
                      passive = Finder3(passive, state[7], Sid_ash)
                      state[7] -= 1
                   
           
           
           
            if (search_res[0] == False ) & resrclss_dem == 1 & is_arrvd == 1:
              
                Info = [Id, demand[8], change, demand[1], demand[2], demand[3], demand[9]]
                passive = Finder2_2(passive, state[7], Info)
                state[7] += 1     
    
    return [active, passive, reduction]