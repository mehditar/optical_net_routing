import numpy as np
from trans_store_clean import *
def max_band_selec(newNet, newmyNet, _relsc_cap, new_my_relsc_cap,
        t_3, t_3_map, path, path_rate, path_worth,
        Zeta, trx_mx, f_slice_cnt, active_cnt, is_active, 
        non_redctn, ith_path, ttx_remnng):
    
    

    
    #min_bit_-start
    delta = 10000000
    temp = _relsc_cap.reshape(active_cnt, 1)
    _relsc_cap = np.tile(temp, (1, f_slice_cnt))
    min_bit_bef = np.zeros((active_cnt, f_slice_cnt))
    min_bit_aft = np.zeros((active_cnt, f_slice_cnt))
    #--min_bit_-end
    #---
    #---
    
    #midf-start
    path = path[:, :, np.newaxis]
    ints = newNet * np.tile(path, (1, 1, f_slice_cnt))
    signal = np.unique(ints)
    signal = signal[1:]
    #--
    #-condition:NoSIGNAL
    if len(signal) == 0:
        thickness = min([trx_mx+1, f_slice_cnt])
        objective = (thickness - 1) * path_rate * (1 + 10 * Zeta * path_worth)
        ind = 0
        if (thickness - 1) * path_rate < - new_my_relsc_cap:
          return [-1, ind,thickness]  
            
        return [objective, ind, thickness]
    #---
    #------
    signals = len(signal)
    sheet = np.zeros((signals, f_slice_cnt))
    sheet_bit = np.zeros((signals, f_slice_cnt))
    GIGA = (1+ (non_redctn - 1) * delta)
    for j in range(signals):
        sheet[j, :] = np.amax(np.amax((
            ints == signal[j]), axis = 0), axis = 0)
        sheet_bit[j, :] = sheet[j, :] * t_3[signal[j] - 1][2] * (
            1 + 10 * Zeta * t_3[signal[j] - 1][3])
    mid = sheet.copy()
    temp_1 = mid.copy()
    temp_1[:, 0:-1] = mid[:, 1:]
    temp_1[:, -1] = np.zeros((signals))
    temp_2 = mid.copy()
    temp_2[:, 1:] = mid[:, 0:-1]
    temp_1[:, 0] = np.zeros((signals))
    mid = np.logical_and(mid, np.logical_and(temp_1, temp_2))
    bef_mid = mid.copy()
    aft_mid = mid.copy()
    #midf-end
    
    #gaurd-start
    guard = sheet.copy()
    guard_bit = sheet_bit.copy()
    bef_guard = guard.copy()
    after_guard = guard.copy()
    #gaurd-end
    
    #general
    Pth = np.tile(path, (1, 1, f_slice_cnt))
    cost = np.zeros((1, f_slice_cnt))
    benefit = np.ones((1, f_slice_cnt)) * path_rate * (
        1 + 10 * Zeta * path_worth)
    mx = np.zeros((f_slice_cnt, 2))
    #general
    
    #mybit-start
    my_bit = newmyNet.copy()
    my_cost = np.zeros((1, f_slice_cnt))
    my_new_bit = np.ones((1, f_slice_cnt)) * path_rate
    #mybit-end
    
    #-----
    #------
    for i in range(f_slice_cnt):
        st = np.unique(newNet[:, :, i] * Pth[:, :, i])
        my = np.unique(my_bit[:, :, i] * Pth[:, :, i])
        for j in st:
            if j != 0:
                cost[0, i] = cost[0, i] + t_3[j - 1][2] * (
                    1 + 10 * Zeta * path_worth) * GIGA
                min_bit_bef[int(t_3_map[0, int(j - 1)]), i] = min_bit_bef[
                    int(t_3_map[0, int(j - 1)]), i] + t_3[j - 1][2]
        #mybit
        for j in my:
            if j != 0:
                my_cost[0, i] = my_cost[0, i] + t_3[j - 1][2] * GIGA
        #my_bit
    #---
    #---
        
    bef = benefit - cost
    aft = np.zeros((1, f_slice_cnt))
    
    #mybit
    my_new_bit = my_new_bit - my_cost
    mybef = my_new_bit.copy()
    myaft = np.zeros((1, f_slice_cnt))
    #mybit
    
    #---
    #---
    for i in range(min([trx_mx + 1, f_slice_cnt, int(ttx_remnng + 1)])):
        
        #--min_bit_-start
        min_bit_aft[:, i:] = min_bit_aft[:, i:] + min_bit_bef[
            :, 0:f_slice_cnt - i]
        aa = np.amax((min_bit_aft - _relsc_cap) > 0, axis = 0) * delta
        #--min_bit_-end
        
        #-midf-start
        aft_mid[:, i:] = np.logical_and(aft_mid[:, i:], bef_mid[
            :, 0:f_slice_cnt - i])
        bb = np.amax(aft_mid, axis = 0) * delta
        #-midf-end
        
        #gaurd-start
        temp = np.zeros((signals, f_slice_cnt))
        temp[:, i:] = bef_guard[:, 0:f_slice_cnt - i]
        after_guard = np.logical_and(after_guard, temp)
        #
        b_1 = after_guard.copy()
        b_1[:, 0:-1] = after_guard[:, 1:]
        b_1[:, -1] = np.zeros((signals))
        b_1 = np.logical_not(b_1)
        b_2 = after_guard.copy()
        b_2[:, 1:] = after_guard[:, 0:-1]
        b_2[:, 0] = np.zeros((signals))
        b_2 = np.logical_not(b_2)
        #
        no_shift = np.dot(np.ones((1, signals)), (
            np.logical_and(after_guard, np.logical_and(b_1, b_2)) * guard_bit))
        no_shift = no_shift.reshape(1, f_slice_cnt)
        w_shift = np.zeros((1, f_slice_cnt))
        w_shift[0, 0:f_slice_cnt - i] = no_shift[0, i:f_slice_cnt]
        cc = w_shift
        
        #
        bef = bef + cc
        #gaurd-end
        
        #mybit
        myaft[:, i:] = myaft[:, i:] + mybef[:, 0:f_slice_cnt - i]
            
        dd = np.amax((myaft - path_rate) < - new_my_relsc_cap , axis = 0) * delta
        #mybit
        
        aft[0, i:] = aft[0, i:] + bef[0, 0:f_slice_cnt - i]
        aa = aa.reshape(1, f_slice_cnt)
        bb = bb.reshape(1, f_slice_cnt)
        dd = dd.reshape(1, f_slice_cnt)
        mm = aft - aa - bb - dd
        mx[i, 0] = np.amax(mm[0, i:])
        mx[i, 1] = np.argmax(mm[0, i:])
        mx[i, 1] = mx[i, 1] + i
    #--
    #---
    
    objective = np.amax(mx[:, 0], axis = 0)
    thickness = np.argmax(mx[:, 0], axis = 0)
    thickness = int(thickness + 1)
    objective = objective - path_rate * (1 + 10 * Zeta * path_worth)
    ind = int(mx[int(thickness - 1), 1] - thickness + 1)
    return [objective, ind, thickness]
   