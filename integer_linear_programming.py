#from gurobipy import *
import time as time
import xlrd
import numpy as np
import pandas as pd
import csv
from trans_store_clean import *
def router_ILP(t_1andt_5, t_2, t_3, t_4, served, key_dic,
       src_node, dest_nodes, is_active,
       Sid, bitrate_min, state, dem_id_dict,
       NETWORKDIS, NETWORKLINK, PARAM, PATHS, INFO, PATHORDER, LINKNAME):


 

 signl_s = set()
 rev_active = set()
 demand = np.zeros((6, 1))
 demand[4] = src_node
 demand[5] = dest_nodes
 t_1 = t_1andt_5[0]
 active_cnt = np.size(t_1, 0)
 t_5 = t_1andt_5[1]
 light_cnt = np.size(t_4, 1)
 signl_cnt = len(t_3)
 for i in range(signl_cnt):
    signl_s.add(i)
 if is_active == 1:
     my_signls = list(t_2[Sid])
     signl_s = signl_s-t_2[Sid]
     my_signls_1 = len(my_signls)
     my_signls_2 = (served[3][state[2]][4]) + 1
 else:
     my_signls = []
     my_signls_1 = 0
     my_signls_2 = 0
 if is_active == 0:
     t_2replace = {}
     t_2replace[0] = set()
 else:
     t_2replace = t_2      
 for i in range(active_cnt):
     if i == Sid:
        if is_active == 0:
            rev_active.add(i)
     else:
        rev_active.add(i)
 
        
 try:

   
    m = Model("TS")
    
    link_mtrx = NETWORKLINK
    dist_mtrx = NETWORKDIS
    n_of_nodes = len(NETWORKLINK)
    n_of_links = int(sum(sum(link_mtrx)) / 2)
    n_of_dests = (demand.shape[0]) - 5
    n_of_demands = 1
    f_slice_cnt = int(PARAM[0])
    MAXDIST = sum(sum(dist_mtrx)) / 2
    trx_mx_tot = int(PARAM[2])
    A = 1#!!!!
    B = -0.0967
    
    
    
    W = range(f_slice_cnt) 
    M = [0, 1, 2, 3, 4]
    D = range(n_of_dests)
    link_range = range(int(2 * n_of_links))
    MOD_REACH = [0, 625, 1250, 2500, 5000]
    Data_Rate_PER_SLICE = [0, 50, 37, 25, 12]
    K = {0}
    Kdic = {}
    
    
    V = range(n_of_nodes)
    V_set = set()
    V_dict = {}
    link_lngth = {}
    paired_links = []
    for i in V:
        V_set.add(i)
    for v in V:
        V_dict[v] = {}
        V_dict[v]["in_links"] = set()
        V_dict[v]["out_links"] = set()
    link_ind = 0
    for i in V:
        for j in V:
                if link_mtrx[i][j] == 1:
                    V_dict[i]["in_links"].add(link_ind)
                    V_dict[j]["out_links"].add(link_ind)
                    link_lngth[link_ind] = dist_mtrx[i][j]
                    link_ind += 1
    for i in V:
        for j in V:
            if NETWORKLINK[i][j] == 1:
                paired_links.append([V_dict[i]["in_links"].
                                     intersection(V_dict[j]["out_links"]).pop(),
                                 V_dict[i]["out_links"].
                                     intersection(V_dict[j]["in_links"]).pop()])
    
    Zandset = np.zeros((n_of_demands, n_of_dests, 2)) 
    for k in K:
        for d in D:
            Zandset[k][d][0] = demand[4][k]
            Zandset[k][d][1] = demand[5 + d][k]
    
    Sigmaset = set()       
    for i in range(int(n_of_links)):
        Sigmaset.add((2 * i))
    
          
    

   
    H = m.addVars(2 * n_of_links, n_of_demands, n_of_dests, lb = 0,
                ub = 1, vtype = GRB.BINARY, name= "H" )
    Sub = m.addVars(n_of_demands, n_of_dests, n_of_dests,
                  lb = 0, ub = 1, vtype = GRB.BINARY, name = "Sub")
    cr = m.addVars(n_of_demands, n_of_dests, n_of_dests, lb = 0,
                 ub = 1,vtype = GRB.BINARY, name = "cr")
    Subr = m.addVars(n_of_demands, n_of_dests, n_of_dests, lb = 0,
                   ub = 1,vtype = GRB.BINARY, name = "Subr")
    R = m.addVars(2 * n_of_links, n_of_demands, n_of_dests, lb = 0, ub = 1,
                vtype = GRB.BINARY, name = "R")
    Reg = m.addVars(n_of_nodes, n_of_demands, n_of_dests, lb = 0, ub = 1,
                  vtype = GRB.BINARY, name = "Reg")
    Z = m.addVars(2 * n_of_links, n_of_demands, n_of_dests, f_slice_cnt, lb = 0,
                ub = 1, vtype = GRB.BINARY, name = "Z")
    Zand = m.addVars(2, n_of_demands, n_of_dests, f_slice_cnt, lb = 0,
                   ub = 1, vtype = GRB.BINARY, name = "Zand")
    Mod = m.addVars(n_of_demands, n_of_dests, 5, lb = 0, ub = 1,
                  vtype = GRB.BINARY, name= "Mod")
    Modr = m.addVars(n_of_demands, n_of_dests, 5, lb = 0, ub = 1,
                   vtype = GRB.BINARY, name = "Modr")
    ft = m.addVars(n_of_demands, n_of_dests, lb = 0, ub=f_slice_cnt,
                 vtype = GRB.INTEGER, name = "ft")
    fr = m.addVars(n_of_demands, n_of_dests, lb = 0, ub = f_slice_cnt,
                 vtype = GRB.INTEGER, name = "fr")
    Sigma = m.addVars(2 * n_of_links, n_of_demands, f_slice_cnt, lb=0, ub=1,
                    vtype = GRB.BINARY, name = "Sigma")
    Com = m.addVars(n_of_demands, n_of_dests, n_of_dests,
                  lb = 0, ub = 1, vtype = GRB.BINARY, name = "Com")
    Alpha = m.addVars(n_of_demands, n_of_dests, n_of_dests,
                    lb = 0, ub = 1, vtype = GRB.BINARY, name = "Alpha")
    s1 = m.addVars(n_of_demands, n_of_dests, lb = 0, ub = 1,
                 vtype = GRB.BINARY, name = "s1")
    s2 = m.addVars(n_of_demands, n_of_dests, lb = 0, ub = 1,
                 vtype = GRB.BINARY, name = "s2")
    s3 = m.addVars(n_of_demands, n_of_dests, lb = 0, ub = 1,
                 vtype = GRB.BINARY, name = "s3")
    Lights = m.addVars(light_cnt, lb = 0, ub = 1, vtype = GRB.BINARY,
                     name = "lights")
    Mos = m.addVars(n_of_demands, n_of_dests, f_slice_cnt, 5, lb = 0, ub = 1,
                  vtype = GRB.BINARY, name = "mos")
    datarate_goal = m.addVar(vtype = GRB.INTEGER, name = "datarate_goal")
    new_sig_datarate = m.addVar(vtype = GRB.INTEGER, name = "new_sig_datarate")
    Signal = m.addVars(signl_cnt, lb = 0, ub = 1, vtype = GRB.BINARY, name = "signal")
    Dofail = m.addVar(lb = 0, ub = 1, vtype = GRB.BINARY,name = "Dofail") 
    Datarates = m.addVars(active_cnt, vtype = GRB.INTEGER, name = "datrates")
    


    #creating constrs 
    m.addConstrs(((quicksum(H[j, k, d] for j in V_dict[i]["in_links"]) - 1 <= 0)
                  for i in V for k in K for d in D), "aa")
    m.addConstrs(((quicksum(H[j, k, d] for j in V_dict[i]["out_links"]) - 1 <= 0)
                  for i in V for k in K for d in D), "ab")
    m.addConstrs(((quicksum(H[j, k, d] for j in V_dict[demand[5 + d][k]]["in_links"])
                   - quicksum(H[j, k, d] for j in V_dict[demand[5 + d][k]]["out_links"]) == 1)
                  for k in K for d in D), "ac")
   
    m.addConstrs(((quicksum(H[j, k, d] for j in V_dict[demand[4][k]]["out_links"])
                   - quicksum(H[j, k, d] for j in V_dict[demand[4][k]]["in_links"]) == 1)
                  for k in K for d in D), "acc")
    m.addConstrs(((quicksum(H[j, k, d]for j in V_dict[i]["in_links"])
                   - quicksum(H[j, k, d]for j in V_dict[i]["out_links"]) == 0)
                  for k in K for d in D for i in 
                  (V_set - set(demand[[4], [k]]) - set(demand[[5 + d], [k]]))), "ae")
    m.addConstrs((quicksum(Z[j, k, d, w]for w in W) - H[j, k, d] * 400 <= 0
                  for j in link_range for k in K for d in D), "ah")
    
    m.addConstrs((Reg[ia, k, da] + Reg[ib, k, db] + Sub[k, da, db] -2 <= 0 
                  for ia in V_set for ib in V_set if ia != ib 
                  for k in K for da in D for db in D), "ap")
    m.addConstrs((quicksum(Reg[i, k, d]for i in V_set) == 0 
                  for k in K for d in D), "aq")

    m.addConstrs((R[j, k, d] - H[j, k, d] <= 0 for j in link_range 
                  for k in K for d in D), "az")
    m.addConstrs((quicksum(R[j, k, d]for j in V_dict[demand[5 + d][k]]["in_links"])
                  - quicksum(Reg[i, k, d]for i in V_set) == 0 
                  for k in K for d in D), "ba")
    m.addConstrs((quicksum(R[j, k, d]for j in V_dict[i]["out_links"])
                  - quicksum(R[j, k, d]for j in V_dict[i]["in_links"]) - Reg[i,k,d] == 0 
                  for k in K for d in D for i in (V_set - set(demand[[4, 5 + d], [k, k]]))), "bb")
  
    
   
    m.addConstrs((cr[k, da, db] * 2 - Sub[k, da, db] -
                  quicksum(Reg[i, k, db]for i in V_set) <= 1  
                  for k in K for da in D for db in D), "awer")
    
    m.addConstrs((Mod[k, d, m]
                  * MAXDIST + quicksum((link_lngth.get(j) * (H[j, k, d] - R[j, k, d])for j in link_range)) -
                  (MOD_REACH[m] * (A + B * quicksum(Sub[k, d, db] - cr[k, d, db]for db in D if db != d)))
                  - (1 / 2) <= MAXDIST for k in K for d in D for m in M), "bh")
    m.addConstrs((Modr[k, d, m]
                  * MAXDIST + quicksum((link_lngth.get(j) * (R[j, k, d])for j in link_range)) -
                  (MOD_REACH[m] * (A + B * quicksum(Subr[k, d, db]for db in D if db != d)))
                  - (1 / 2) <= MAXDIST for k in K for d in D for m in M), "bh2")

  
    m.addConstrs((Mod[k, da, m] - Mod[k, db, m] - Sub[k, da, db] + 1 >= 0
                  for k in K for da in D for db in D for m in M), "plg")
    m.addConstrs((Modr[k, da, m] - Modr[k, db, m] - Subr[k, da, db] + 1 >= 0
                  for k in K for da in D for db in D for m in M), "plg2")
   
    m.addConstr(quicksum(Mod[0, 0, m] for m in M) == 1, "bk")
    
    m.addConstrs(ft[k, d] == quicksum
                 (quicksum(Z[j, k, d, w]for w in W)
                  for j in V_dict[demand[4][k]]["out_links"])for d in D for k in K)
    
    m.addConstrs(fr[k,d]==quicksum
                 (quicksum(Z[j, k, d, w]for w in W)
                  for j in V_dict[demand[5 + d][k]]["in_links"])for d in D for k in K)
    
    m.addConstrs((2 * Zand[b, k, d, w] - quicksum(Z[j, k, d, w]for 
                                           j in V_dict[Zandset[k][d][b]]["in_links"] | V_dict[Zandset[k][d][b]]["out_links"] )
                  - ( -1 * quicksum(Z[j, k, d, w - 1]for j in V_dict[Zandset[k][d][b]]["in_links"]
                                | V_dict[Zandset[k][d][b]]["out_links"]) + 1) <= 0 
                  for d in D for k in K for b in range(2) 
                  for w in W if w != 0), "kjhd")
    
    m.addConstrs((2 * Zand[b, k, d, w] - quicksum(Z[j, k, d, w]for
                                           j in V_dict[Zandset[k][d][b]]["in_links"]
                                                  |V_dict[Zandset[k][d][b]]["out_links"] )
                  - ( -1 * quicksum(Z[j, k, d, w - 1]for j in V_dict[Zandset[k][d][b]]["in_links"]
                                | V_dict[Zandset[k][d][b]]["out_links"]) + 1) >= (-3 / 2) 
                  for d in D for k in K for b in range(2) for w in W if w != 0), "kjfd2")
    

    
    m.addConstrs((Zand[b, k, d, 0] == quicksum
                  (Z[j, k, d, 0]for j in V_dict[Zandset[k][d][b]]["in_links"]
                   | V_dict[Zandset[k][d][b]]["out_links"]) for d in D 
                  for k in K for b in range(2)), "kjhd_prim")
  
    m.addConstrs((quicksum(Zand[b, k, d, w]for w in W)
                  <= 1 for k in K for d in D for b in range(2)), "haha")
   
    m.addConstrs((quicksum(Z[j, k, d, w]for j in V_dict[i]["out_links"])
                  - quicksum(Z[j, k, d, w]for j in V_dict[i]["in_links"]) + Reg[i, k, d] >= 0
                  for w in W for k in K for d in D for i in
                  (V_set - set(demand[[4], [k]]) - set(demand[[5 + d], [k]]))), "bm")
    
    m.addConstrs((quicksum(Z[j, k, d, w]for j in V_dict[i]["out_links"])
                  - quicksum(Z[j, k, d, w]for j in V_dict[i]["in_links"]) - Reg[i, k, d] <= 0
                  for w in W for k in K for d in D for i in
                  (V_set - set(demand[[4, 5 + d], [k, k]]))), "bm2")
  
    m.addConstrs((quicksum(w * (Zand[0, k, da, w] - Zand[0, k, db, w])
                           for w in W) - f_slice_cnt * (1 - Sub[k, da, db]) <= 0 
                  for k in K for da in D for db in D),"juih")
    
    m.addConstrs((quicksum(w * (Zand[1, k, da, w] - Zand[1, k, db, w])
                           for w in W) - f_slice_cnt * (2 - Subr[k, da, db] - Com[k, da, db]) <= 0 
                  for k in K for da in D for db in D), "juih2")
  
    m.addConstrs((Z[j, k, da, w] + Z[j, k, db, w] - Sub[k, da, db] <= 1 
                  for w in W for k in K for da in D for db in D
                  for j in link_range), "ttf")
    
    m.addConstrs((Subr[k, da, db] - Sub[k, da, db] <= 0 
                  for k in K for da in D for db in D), "uvh")
    
    m.addConstrs((Com[k, da, db] - Sub[k, da, db] >= 0 
                  for k in K for da in D for db in D), "uvh2")

    m.addConstrs((R[j, k, da] + R[j, k, db] <= (2 + Com[k, da, db] - Subr[k, da, db])
                  for j in link_range for k in K for da in D for db in D), "azpl")
    
    m.addConstrs((f_slice_cnt * Alpha[k, da, db] -
                  quicksum(w * (Zand[1, k, da, w] - Zand[1, k, db, w])for w in W)
                  - (1 / 2) >= 0 for k in K for da in D for db in D), "cop1")
    
    m.addConstrs((f_slice_cnt * Alpha[k, da, db] -
                  quicksum(w * (Zand[1, k, da, w] - Zand[1, k, db, w])for w in W)
                  - (1 / 2) <= f_slice_cnt for k in K for da in D for db in D), "cop2")

    m.addConstrs((quicksum(w * (Zand[1, k, da, w] - Zand[1, k, db, w])for w in W)
                  + f_slice_cnt * 2 * (Alpha[k, db, da] - 1) +
                  fr[k, da] - f_slice_cnt * 2 * (1 + Subr[k, da, db] - Com[k, da, db]) <= 0 
                  for k in K for da in D for db in D),"akh")
    
    m.addConstrs((Sub[k, da, db] - Sub[k, db, da] == 0 
                  for k in K for da in D for db in D), "ghn")
    
    m.addConstrs((Subr[k, da, db] - Subr[k, db, da] == 0 
                  for k in K for da in D for db in D), "ghn2")
    
    m.addConstrs((Sub[k, da, db] == 1 for k in K for da in D
                  for db in D), "khn")
  
    
    
    m.addConstrs((H[paired_links[a][0], 0, 0] + H[paired_links[a][1], 0, 0] <= 1 
                  for a in range(len(paired_links))), "loop")
   
    m.addConstrs((1 - Lights[g]) * len(dem_id_dict[t_4[0, g]][2]) >=
                 quicksum(Z[j, 0, 0, dem_id_dict[t_4[0, g]][1][3]] for j in dem_id_dict[t_4[0, g]][2])
                 for g in range(light_cnt))
    
    if PARAM[8] == 2:
      m.addConstrs(Signal[g] == 1 for g in signl_s)
    
       
    
   
    m.addConstrs((Lights[t_5[g]] - Lights[t_5[g-1]] + 1) * 300 >=
                 quicksum(Lights[t_5[k]]for k in greater_membrs(t_3[j][1], g)) 
                 for i in range(active_cnt) for j in t_2[i] 
                 for g in (t_3[j][1] - {min(t_3[j][1])} - {max(t_3[j][1])}))
    
    
    m.addConstrs(Datarates[i] == quicksum(((quicksum(Lights[t_5[kk]]
                                                    for kk in t_3[j][1]) - Signal[j]) * t_3[j][2])for
                                         j in t_2[i])for i in rev_active)
    m.addConstrs(Datarates[i] >= t_1[i,2] for i in rev_active)
   
    
    
    m.addConstrs(2 * Mos[k,d,w,m] <=
                 quicksum(Z[j, k, d, w]for j in V_dict[demand[5 + d][k]]["in_links"])
                 + Mod[k, d, m] for k in K for d in D for m in M for w in W)
    m.addConstr(new_sig_datarate ==
                is_active * quicksum(((quicksum(Lights[t_5[kk]]for kk in t_3[j][1])
                                       - Signal[j]) * t_3[j][2])for j in t_2replace[Sid])
                + quicksum(Mos[0, 0, w, m] * Data_Rate_PER_SLICE[m] for w in W for m in M)
                - quicksum(Mod[0, 0, m] * Data_Rate_PER_SLICE[m]for m in M))
    m.addConstr(is_active *
                quicksum(((quicksum(Lights[t_5[kk]]for kk in t_3[j][1])
                           - Signal[j]))for j in t_2replace[Sid]) +
                quicksum(Z[j, 0, 0, w]
                         for j in V_dict[demand[5 + 0][0]]["in_links"]for w in W) - 1 <= trx_mx_tot)
   
    m.addConstr(quicksum(Z[j, 0, 0, w]
                         for j in V_dict[demand[5 + 0][0]]["in_links"]for w in W) - 1 <= PARAM[7])
    m.addConstrs(Signal[i] * 300 >=
                 quicksum(Lights[t_5[j]]
                          for j in t_3[i][1])for i in range(signl_cnt))
    m.addConstrs(Signal[i] <=
                 quicksum(Lights[t_5[j]]
                          for j in t_3[i][1])for i in range(signl_cnt))
    m.addConstr(quicksum(Signal[i] for i in my_signls) <= PARAM[6] - 1)
    
    if my_signls_2 == PARAM[4]:
         m.addConstr(Dofail == 1 / 2)
            
    if my_signls_1 == PARAM[5]:

         m.addConstr(Dofail == 1/2)

    
    m.addConstr(new_sig_datarate >= bitrate_min)
    m.addConstr(datarate_goal == PARAM[1] * (new_sig_datarate + quicksum(Datarates[i]
                                                    for i in rev_active)))
    m.addConstr(datarate_goal >= state[4] + 12)
    
  
    m.setPARAM('TimeLimit', n_of_demands * 130)
    m.setPARAM('OutputFlag', False)
  
    m.optimize()


  

 except AttributeError:
    print('Encountered an attribute error')
 m.write("myexample.mps")






 
 link_s = set()
 frqncy_s = set()
 to_relsc_set = set()
 to_add_set = set()
 reduction = 0   
 new_goal = 0
 modulation = 0
 links_in_path = 0
 bitrate = 0
 bitrate_tot = 0
 reduction = 0
 success = False
 
 if m.status == 2:
        success = True
        new_goal = datarate_goal.x
        valid_modulations = set()
        for j in link_range:
            for w in W:
                b = 0
                b = (Z[j, 0, 0, w].x)
                if b > .5:
                    link_s.add(j)
                    frqncy_s.add(w) 
        links_in_path = len(link_s)
        for i in range(light_cnt):
            if Lights[i].x <= 0.5:
                reduction = 1
                to_relsc_set.add(t_4[0, i])
  
                
        for m in M:
            if Mod[0, 0, m].x > .5:
                valid_modulations.add(m)
        modulation = min(valid_modulations)
        Rid = state[3] + 1
        bitrate = sum([Mod[0, 0, m].x * Data_Rate_PER_SLICE[m] for m in M])
        bitrate_tot = bitrate * (len(frqncy_s) - 1)
        
        while len(frqncy_s) != 0:
            a = min(frqncy_s)
            frqncy_s = frqncy_s - {a}
            dem_id_dict[Rid] = {}
            dem_id_dict[Rid][1] = [len(link_s), Drateperslice, min(valid_modulations), a]
            dem_id_dict[Rid][2] = link_s
            to_add_set.add(Rid)
            if len(frqncy_s) != 0:
                Rid += 1
        
       
        
       
        
 
        
     
                
 return [success, new_goal, to_add_set, modulation, links_in_path, bitrate,
         bitrate_tot, reduction, to_relsc_set]

