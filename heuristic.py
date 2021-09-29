import numpy as np
from aLightPathSelector import LPS
from bandwidth import max_band_selec
from trans_store_clean import *
def router_heuristic(t_1andt_5, t_2, t_3, t_4, served, key_dict, src_node, dest_nodes,
        is_active, Sid, bitrate_min, state, dem_id_dict, NETWORKDIS,
        NETWORKLINK, PARAM, PATHS, INFO, PATHORDER, LINKNAME):
   # this function performs one complete search to find a path and modulationd and spectrum for a demand
    

    

    
   n_of_nodes = np.size(NETWORKLINK, 0) # number of nodes of the network
   t_1 = t_1andt_5[0]
   signal_cnt = len(t_3) # number fo signals that are dedicated to different demands.
   active_cnt = np.size(t_1, 0) # the number of active demands
   trx_mx = int(PARAM[7]) 
   f_slice_cnt = int(PARAM[0])
   k_shrt_pths = int(PARAM[11])
   my_signals_cnt_mx = int(PARAM[6])
   goal = PARAM[4]
   Zeta = PARAM[10]
   non_redctn = int(PARAM[8])
   trx_tot_mx = int(PARAM[2])
   tot_pths_per_pair = int(np.size(PATHS, 2) / (n_of_nodes * (n_of_nodes - 1)))
   ind_strt = (PATHORDER[int(src_node), int(dest_nodes)] - 1) * tot_pths_per_pair
   paths = PATHS[:, :, ind_strt: ind_strt + k_shrt_pths] # here the k shortest paths nbetween the source
   # and the destinations are retirved
   path_rates = INFO[3, ind_strt: ind_strt + k_shrt_pths] # the data rates for each path is retrived
   path_mods = INFO[2, ind_strt: ind_strt + k_shrt_pths]# the modulation level that can be used 
   #for each path, calcu;ated based on the length of the path, is retrived
   path_worths = INFO[4, ind_strt: ind_strt + k_shrt_pths]# and a value that is a combination of length 
   # and rate of each path s retrieved
   ##
   if is_active == 1:
       path_worths = path_worths / (t_1[Sid][6])
   else:
       path_worths = np.ones((k_shrt_pths, 1))


   if is_active == 1:
       my_signal_s = list(t_2[Sid])
       my_signals_cnt = len(my_signal_s)
       my_signals_ever = (served[3][state[2]][4]) + 1
   else:
       my_signal_s = []
       my_signals_cnt = 0
       my_signals_ever = 0

  
   #######################condition checking########################
   if (my_signals_cnt == PARAM[5]) or (PARAM[4] == my_signals_ever):
        success = False; new_goal = 0; to_add_set = {}; modulation = 0; links_in_path = 0;
        bitrate = 0; bitrate_tot = 0; reduction = False; to_relsc_set = set()
        return [success, new_goal, to_add_set, modulation, links_in_path, bitrate
                ,bitrate_tot, reduction, to_relsc_set, dem_id_dict]
  

   #----------------------------------------------------------------
   #----------------------------------------------------------------
   net_mtrx = np.zeros((n_of_nodes, n_of_nodes, f_slice_cnt))
   # for this heuristic algorithm, matrix operation are used extensively
   # the information about the network and recources that are aither unused or dedictaed 
   # will be transfered to these matrices.
   net_mtrx_sub = np.zeros((n_of_nodes, n_of_nodes, f_slice_cnt))
   my_net = np.zeros((n_of_nodes, n_of_nodes, f_slice_cnt)) # this matrix will contain the information
   # about the k paths that are supposed to be examined indivudually.
   #
   _relsc_cap = t_1[:, 5] - t_1[:, 2]
   #
   for i in range(signal_cnt): # 
     mn = min(t_3[i][1])
     mn_f = dem_id_dict[mn][1][3]
     path_id = dem_id_dict[mn][3]
     mx = max(t_3[i][1])
     mx_f = dem_id_dict[mx][1][3]
     for j in range(mx_f - mn_f + 1):
         # the information are loaded here
         f = mn_f + j
         net_mtrx_sub[:, :, f] = net_mtrx_sub[:, :, f] + (PATHS[:, :, path_id] * (mn + j + 1))
         net_mtrx[:, :, f] = net_mtrx[:, :, f] + (PATHS[:, :, path_id] * (i + 1))
         if i in my_signal_s:
             my_net[:, :, f] = my_net[:, :, f] + (PATHS[:, :, path_id] * (i + 1))
   #
   t_3_map = np.zeros((1, signal_cnt))
   for i in range(active_cnt):
        for j in list(t_2[i]):
            t_3_map[0, j] = i
   #--------------------------------------------------------------------
       
   ttx = 0
   if is_active == 0:
       ttx = 0
   else:
       # the information about the maximun number of transciver that
       # can be allocated are calculated.
       TTX = np.zeros((len(t_2[Sid])))
       for i in range(my_signals_cnt):
           TTX[i] = len(t_3[my_signal_s[i]][1]) - 1
           ttx += TTX[i]

   #derivatives
   if is_active == 1:
       my_relsc_cap = _relsc_cap[Sid]
       new_my_relsc_cap = my_relsc_cap.copy()
       _relsc_cap[Sid] = _relsc_cap[Sid] + 1000000
   else:
       my_relsc_cap = -1 * bitrate_min
   output = [0, 0, 0, 0, 0]
   loss = 0
   

   for i in range(k_shrt_pths):
    # for each one of the k shortest path are tested. for each the maximum band that can be
    # allocated is searched. at the end all path and the allocated results are compared and 
    # the one that has overall either high band or low ength or both is selected and returned. 
    path = paths[:, :, i]
    path_rate = path_rates[i]
    path_worth = path_worths[i]
    if my_signals_cnt == my_signals_cnt_mx:
        # in this loop if the demand is active and the number of paths allocated to it hits the cap
        # then the individual paths are, one at a time, removed and the best possible substitute path
        # is searched.
        for j in range(my_signals_cnt):
            sig = my_net * (my_net == (my_signal_s[j] + 1))
            new_net = net_mtrx - sig
            new_my_net = my_net - sig
            new_my_relsc_cap = my_relsc_cap - (len(t_3[my_signal_s[j]][1]) - 1
                                              )* t_3[my_signal_s[j]][2]
            ttx_remnng = trx_tot_mx - (ttx - TTX[j])
            res = max_band_selec(new_net, new_my_net, _relsc_cap, new_my_relsc_cap, t_3,t_3_map,
                     path, path_rate, path_worth, Zeta, trx_mx, f_slice_cnt,
                    active_cnt, is_active, non_redctn, 9000, ttx_remnng)
            if res[0] > output[0]:
                output[0] = res[0]
                output[1] = res[1]
                output[2] = res[2]
                output[3] = i
                output[4] = j
                loss = (len(t_3[my_signal_s[j]][1]) - 1) * t_3[my_signal_s[j]][2]*(
                    1 + 10 * Zeta * t_3[my_signal_s[j]][3]) * (1 + (non_redctn - 1) * 100000000)
                
        
    else:
            # if this is new demand, it means no paths are already stablished for the demand
            # therefore, the search is conducted only by considering the candidate paths
            ttx_remnng = trx_tot_mx - ttx
            res = max_band_selec(net_mtrx, my_net, _relsc_cap, my_relsc_cap, t_3,
                                 t_3_map, path, path_rate, path_worth,Zeta, trx_mx,
                                 f_slice_cnt, active_cnt, is_active, non_redctn, i, ttx_remnng)           
            if res[0] > output[0]:
                output[0] = res[0]
                output[1] = res[1]
                output[2] = res[2]
                output[3] = i
                output[4] = -1
                
 
   #aftercalculations:
   # here the search is complete and the results are loaded and returned 
   # in the next sextiona

   success = False
   new_goal = 0
   to_add_set = {}
   modulation = 0
   links_in_path = 0
   bitrate = 0
   bitrate_tot = 0
   reduction = False
   to_relsc_set = set()
   freqncy_s = set()
   if (output[0] - loss) > 0:
        success = True
        new_goal = output[0] - loss
        #
        temp = np.unique(paths[:, :, output[3]] * LINKNAME)
        link_s = set((temp[1: ]).tolist()) # the links in the selected path are saved
        path_id = int(ind_strt + output[3]) # the id of the path is retieved
        worth = INFO[4, path_id]
        #
        freqncy_s = set(list(range(output[1], output[1] + output[2]))) # the spectrum range that
        # is avaialbe alongth this path is retrived
        modulation = path_mods[output[3]] # modulation technique that can be used for the path
        # is returned
        bitrate = path_rates[output[3]] # the data rate per frequency slice is returned
        #
        targetpath = paths[:, :, output[3]]
        set_1 = set()
        for i in range(f_slice_cnt):
            if i in freqncy_s:
               temp = np.unique(net_mtrx_sub[:, :, i] * targetpath)
               set_1 = set_1.union((temp[1:] - 1).tolist())
        set_2 = set()
        if output[4] > -1:
           set_2 = t_3[my_signal_s[output[4]]][1]
        to_relsc_set = set_1.union(set_2) # all the resources that have 
        # conflict with this path and spectrum are determined
        if len(to_relsc_set) > 0:
          reduction = True
        else:
            reduction = False
  
        to_add_set = set()
        Rid = state[3] + 1
        bitrate_tot = bitrate * (len(freqncy_s) - 1)
        links_in_path = len(link_s)
        while len(freqncy_s) != 0: # the above information about the path are stored  in dem_id_dic
            a = min(freqncy_s)
            freqncy_s = freqncy_s - {a}
            dem_id_dict[Rid] = {}
            dem_id_dict[Rid][1] = [len(link_s), bitrate, modulation, a]
            dem_id_dict[Rid][2] = link_s
            dem_id_dict[Rid][3] = path_id
            dem_id_dict[Rid][4] = worth
            to_add_set.add(Rid)
            if len(freqncy_s) != 0:
                Rid += 1
      
        
      
        
  
   return [success, new_goal, to_add_set, modulation, links_in_path, bitrate,
           bitrate_tot, reduction, to_relsc_set, dem_id_dict]
