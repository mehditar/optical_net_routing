import numpy as np
from trans_store_clean import *
def max_band_selec(newNet, newmyNet, _relsc_cap, new_my_relsc_cap,
        t_3, t_3_map, path, path_rate, path_worth,
        Zeta, trx_mx, f_slice_cnt, active_cnt, is_active, 
        non_redctn, ith_path, ttx_remnng):
    
    # this function receives two matrices. for each of these three dinemtional matrices
    # the first two demntions represents the links between nodes and the third dimentions
    # the the frequency slices for link. for each the zero means there is no link and
    # in the third dimention any other number represents the existenxe of the link and rate of 
    # data transfer if that link-frequency slice is already in use, or if it is selected in this
    # round. it should be noted that this rate is subject to change as it depends on the distance
    # between the source and destination not the lenght of the link.
    # to find the best spectrum (or frquency slice range) in this particular situation,
    # we can apply brute force method which will create three nested for loops and,
    # given the high number of the nodes and frequency slice, ithe algorithm looses its agility.
    # here we use an iterative rolling strategy combind with getting maximum to find the best spectrum.

    
    #min_bit_-start
    delta = 10000000
    temp = _relsc_cap.reshape(active_cnt, 1)
    _relsc_cap = np.tile(temp, (1, f_slice_cnt))
    min_bit_bef = np.zeros((active_cnt, f_slice_cnt))# this matrix is just to make sure that whatever 
    # spectrum comes out metts the minimum bitrate qualification needed for the demand.
    min_bit_aft = np.zeros((active_cnt, f_slice_cnt))
    #--min_bit_-end
    #---
    #---
    
    #midf-start
    path = path[:, :, np.newaxis]# the path is replicated into a three dimentional
    # matrix.
    ints = newNet * np.tile(path, (1, 1, f_slice_cnt))# here the signals that are in 
    # use with for other demand and have spectrum conflit with the current path are determined
    signal = np.unique(ints)
    signal = signal[1:]
    
    #--
    #-condition:NoSIGNAL
    if len(signal) == 0:
        # if there is no signal in the network them the bandwith is euivalent to the minimum of all
        # frequency slice and the number of allowed transceivers
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
    # another limitation on the spectrum is that it cannot be the middle spectrum of another signal
    # because this causes the guardband issue. the purpose of the mid, mid_bef, and mid_after
    # matices is to filter out the choices that dont meet this condition
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
    #general. this defines the total spectrum use increment that will be resulted
    
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
        # this is wherre the rolling is perfomed . for each roll, all the matrices are rolled once
        # and only the valid results are returned. 
        
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
        temp = aft - aa - bb - dd # here only the valid results will have positiv results.
        # and the one with the maximum value is stored as the maximum gain that we can get for this 
        # particular spectrum range
        mx[i, 0] = np.amax(temp[0, i:])
        mx[i, 1] = np.argmax(temp[0, i:])
        mx[i, 1] = mx[i, 1] + i
    #--
    #---
    
    objective = np.amax(mx[:, 0], axis = 0) # here the range the provides the maximum gain among all 
    # other spectrum is selected. Please note that a bigger spectrum may not necessirirly result in 
    # a bigger gain as use of it may require emiminating some spectrum of other  
    #signal that are in use to avoid spectrum conflict.
    thickness = np.argmax(mx[:, 0], axis = 0)
    thickness = int(thickness + 1)
    objective = objective - path_rate * (1 + 10 * Zeta * path_worth)
    ind = int(mx[int(thickness - 1), 1] - thickness + 1)
    return [objective, ind, thickness]
   
