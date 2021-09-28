import numpy as np
import csv
import random
def USA(served, active, state, id, output):
    Iid = served[3][id][4] + 1
    served[3][id][1] = served[3][id][1] + 1
    served[3][id][2].add(Iid)
    served[3][id][3][Iid] = {}
    served[3][id][3][Iid]["Gen"] = [len(output[2]), output[5], (
        len(output[2]) - 1) * output[5], output[4]]
    served[3][id][3][Iid]["Main"] = output[2]
    served[3][id][4] = Iid
    
    Sid_ash = Finder(active[0:served[1], 0:2], 0, 1, id)
    active[Sid_ash, 7] = active[Sid_ash, 7] + (len(output[2]) - 1) * output[5]
    active[Sid_ash, 9] = active[Sid_ash, 9] + len(output[2])
    active[Sid_ash, 8] = active[Sid_ash, 6] / active[Sid_ash, 7] + state[0]
    
    return [served, active]
#

def USC(served, id, output):
    served[1] = served[1] + 1
    served[2].add(id)
    served[3][id] = {}
    served[3][id][1] = 1
    served[3][id][2] = set()
    served[3][id][2].add(0)
    served[3][id][3] = {}
    served[3][id][3][0] = {}
    served[3][id][3][0]["Gen"] = [len(output[2]), output[5], (
        len(output[2]) - 1) * output[5], output[4]]
    served[3][id][3][0]["Main"] = output[2]
    served[3][id][4] = 0
    return served
#


def USAR(served,key_dict,active,state,output):
    reduced_already = set()
    for S in output[8]:
         guard = 0
         if S not in reduced_already:
                   served[3][key_dict[S][0]][3][key_dict[S][1]][
                       "Main"] = served[3][key_dict[S][0]][3][key_dict[S][1]]["Main"] - {S}
                   served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][
                       0] = served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][0] - 1
                   served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][
                       2] = served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][
                       2] - served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][1]

                   b = served[3][key_dict[S][0]][3][key_dict[S][1]]["Gen"][1]
                   if len(served[3][key_dict[S][0]][3][key_dict[S][1]]["Main"]) == 1:
                       reduced_already = reduced_already.union(served[3][
                           key_dict[S][0]][3][key_dict[S][1]]["Main"])
                       served[3][key_dict[S][0]][1] = served[3][key_dict[S][0]][1] - 1
                       served[3][key_dict[S][0]][2] = served[
                           3][key_dict[S][0]][2] - {key_dict[S][1]}
                       del served[3][key_dict[S][0]][3][key_dict[S][1]]
                       guard = 1
                   #updatecorrespondingid_in active
                   Sid_ash = Finder(active[0:served[1], 0: 2], 0, 1, key_dict[S][0])
                   active[Sid_ash, 7] = active[Sid_ash, 7] - b
                   active[Sid_ash, 9] = active[Sid_ash, 9] - 1 - guard
                   active[Sid_ash, 8] = active[Sid_ash, 6] / active[Sid_ash, 7] + state[0]

    return[served,active]
#


def TP(served, key_dict, active, state, PARAM, dem_id_dict):
   
    a_size = state[6]
    t_1 = np.zeros((a_size, 6))
    t_1[:, 0] = active[0:a_size, 0]
    t_1[:, 1] = active[0:a_size, 1]
    t_1[:, 2] = active[0:a_size, 10]
    t_1[:, 3] = active[0:a_size, 9]
    t_1[:, 5] = active[0:a_size, 7]
    for i in range(a_size): 
      t_1[i, 4] = served[3][active[i, 1]][1]
   
    temp = sum(t_1[:, 3])

    n = 0
    m = 0
    t_4 = np.zeros((1, int(temp)))
    t_2 = {}
    t_3 = {}
    t_5 = {}
    ind_t = 0
    ind_t_2 = 0
    for i in range(a_size):
        Id = active[i, 1]
        signals = served[3][Id][2]
        ind_t = ind_t + len(signals)

        t_2[i] = set()
        for j in signals:
            t_2[i].add(n)
            t_3[n] = {}
            t_3[n][1] = set()
            t_3[n][2] = served[3][Id][3][j]["Gen"][1]
            for u in served[3][Id][3][j]["Main"]:
                t_3[n][1].add(u) 
                t_4[0, m] = u
                t_5[u] = m
                m += 1
            t_3[n][3] = dem_id_dict[u][4] 
            ind_t_2 += dem_id_dict[u][4]
            n += 1
    #------------------

    old_goal=0
    if state[6] != 0:
        old_goal = PARAM[1] * sum(active[0:state[6], 7])

    return [t_1, t_2, t_3, t_4, t_5, old_goal]
#


def ADR():
    
    def file_reader(file_name,data_type):
        reader = csv.reader(open(file_name, "r"), delimiter=",")
        x = list(reader)
        return np.array(x).astype(data_type)
    
    
 

    P = file_reader("a-PATHS.dat","int")
    Nodes = int(np.size(P, 0))
    Frame = int(np.size(P, 1) / (Nodes * Nodes * (Nodes - 1)))
    PATHS = np.zeros((Nodes, Nodes, (Nodes * (Nodes - 1) * Frame)))
    for i in range(Nodes * (Nodes - 1) * Frame):
        PATHS[:, :, i] = P[:, i * Nodes:i * Nodes + Nodes]
    PATHS = PATHS.astype(int)
    
    
    INFO = file_reader("a-INFO.dat","float")
    PathOrder= file_reader("a-PathOrder.dat","int")
    LinkName = file_reader("a-LinkName.dat","int")

    return [PATHS,INFO,PathOrder,LinkName]
#

def NDF(active, state, q_next, no_new_demand):
    if no_new_demand == 0:
        if state[6] == 0:
            state[2] = int(q_next[0])
            state[1] = 0
            state[0] = q_next[1]
            b = 1
        else:
            a = min(active[0:state[6], 8])
            if a <= q_next[1]:
                for i in range(state[6]):
                    if active[i, 8] == a:
                        Sid_ash = i
                state[2] = int(active[Sid_ash, 1])
                state[1] = 1
                state[0] = a
                b = 0
            else:
                state[2] = int(q_next[0])
                state[1] = 0
                state[0] = q_next[1]
                b = 1
    else:
         a = 0 
         if state[6] != 0:
             a = min(active[0:state[6], 8])
             for i in range(state[6]):
                 if active[i, 8] == a:
                     Sid_ash = i
             c = int(active[Sid_ash, 1])
         else:
             c = int(0)   
         state[2] = c
         state[1] = 1
         state[0] = a
         b = 0
    return [state, b]
#

def Finder3(Array, Size, Id):
    if Size == 1:
        Array[Id, 1:] = Array[Id, 1:] * 0
        return Array
    elif Id == 0:
        Array[0:Size - 1, 1:] = Array[1:Size, 1:]
        Array[Size - 1, 1:] = 0
        return Array
    elif Id == Size - 1:
        Array[Size - 1, 1:] = Array[Size - 1, 1:] * 0
        return Array
    else:
        Array[Id:Size - 1, 1:] = Array[Id + 1:Size, 1:]
        Array[Size - 1, 1:] = 0
        return Array
    
    
                
#

def Finder2(Array, Size, info_array):
    if Size >= 1:
        if info_array[1] >= Array[0, 2]:
            Array[1: Size + 1, 1:] = Array[0:Size, 1:]
            Array[0, 1:] = info_array
            
        elif info_array[1] < Array[Size - 1,2]:
            Array[Size, 1:] = info_array
        else:
            for s in range(1, Size):
                if (info_array[1] >= Array[s, 2]) & (info_array[1] < Array[s - 1, 2]):
                    Array[s + 1:Size + 1, 1:] = Array[s:Size, 1:]
                    Array[s, 1:] = info_array
    else:
         Array[0, 1:] = info_array
    return Array
#

def Finder(Array, Dimention, Dvalue, Value):
    if Dimention == 0:
        for ind in range(len(Array[:, Dvalue])):
            if Array[ind, Dvalue] == Value:
                return ind
    if Dimention == 1:
        for ind in range(len(Array[Dvalue, :])):
            if Array[Dvalue, ind] == Value:
                return ind
#


def DBR(active,passive,state,last_time,mxws):
    a = last_time
    b = state[0]
    if a == b:  
      return [active,passive]

    else:
        if state[6] == 0:
          active = active
        else:
          for i in range(state[6]):
              active[i, 6] = active[i, 6] - (active[i, 7] * (b - a))
              active[i, 8] = active[i, 6] / active[i, 7] + b
              active[i, 10] = active[i, 6] / (active[i, 5] - b)
        if state[7] == 0:
          passive = passive
        else:
          i = 0
          while i != "tamam" :
              if passive[i, 5] <= b or ((passive[i, 6] / (passive[i, 5] - b)) > mxws * 50):
                  passive = Finder3(passive, state[7], i)
                  state[7] = state[7] - 1
              else:
                  i = i + 1
              if i == state[7]:
                  i = "tamam"
        
          for i in range(state[7]):
             passive[i, 7] = passive[i, 6] / (passive[i, 5] - b)
    return [active, passive]
#
def greater_membrs(my_set,value):
    out=set()
    for i in my_set:
        if i>value:
            out.add(i)
    return out
#


def Finder2_2(arr, size, info_arr):
    if size >= 1:
        efficiency = 1
        raandom = 0
        fcfs = 0
        inefficiency = 0
        if efficiency == 1:
                if info_arr[1] >= arr[0, 2]:
                    arr[1:size + 1, 1:] = arr[0:size, 1:]
                    arr[0, 1:] = info_arr
            
                elif info_arr[1] < arr[size - 1, 2]:
                    arr[size, 1:] = info_arr
                else:
                    for s in range(1, size):
                        if (info_arr[1] >= arr[s, 2]) & (info_arr[1] < arr[s - 1, 2]):
                            arr[s + 1: size + 1, 1:] = arr[s:size, 1:]
                            arr[s, 1:] = info_arr
        #--
        ###
        if inefficiency == 1:
                if info_arr[1] <= arr[0, 2]:
                    arr[1: size + 1, 1:] = arr[0:size, 1:]
                    arr[0, 1:] = info_arr
            
                elif info_arr[1] > arr[size - 1, 2]:
                    arr[size, 1:] = info_arr
                else:
                    for s in range(1, size):
                        if (info_arr[1] <= arr[s, 2]) & (info_arr[1] > arr[s - 1, 2]):
                            arr[s + 1:size + 1, 1:] = arr[s:size, 1:]
                            arr[s, 1:] = info_arr
        #-----
        ######
        if fcfs == 1:
            arr[1:size + 1, 1:] = arr[0:size, 1:]
            arr[0, 1:] = info_arr
        #----
        if raandom == 1:
            omi = random.randint(0, size)
            if omi == 0:
                    arr[1:size + 1, 1:] = arr[0:size, 1:]
                    arr[0, 1:] = info_arr
            
            elif omi == size:
                    arr[size, 1:] = info_arr
            else:
                    
                   arr[omi + 1:size + 1, 1:] = arr[omi:size, 1:]
                   arr[omi, 1:] = info_arr
            
        #----   
        
    else:
         arr[0, 1:] = info_arr
    return arr
                