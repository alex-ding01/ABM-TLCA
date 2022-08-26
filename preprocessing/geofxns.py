import numpy as np
from scipy import spatial


def min_dist(lc, code):

    fodder_bool = np.logical_or(lc == 1, lc == 2)
    fodder = np.where(fodder_bool)

    fodder_arr = np.array((fodder[0], fodder[1])).transpose()

    
    if np.isin(code, lc):
        # array contains location rows and columns of the energy crops
        code_loc = np.where(lc == code)

        # print(code_loc[0])
        # print(code_loc[1])

        # convert the array format for calculate
        code_loc_arr = np.array((code_loc[0], code_loc[1])).transpose()

        # calculate the distance and id
        tree = spatial.cKDTree(code_loc_arr)
        mindist, minid = tree.query(fodder_arr)
        # print(mindist, minid)

        # reconstruct 2D np array with distance values
        code_loc_arr_val = np.zeros(code_loc_arr.shape[0])
        idx = np.vstack((code_loc_arr, fodder_arr))
        # print(idx)
        dist = np.vstack((code_loc_arr_val[:, None], mindist[:, None]))

        # convert current unit to km
        dist = np.true_divide(dist, 10)

        dist = np.around(dist, decimals = 4)

        out = np.zeros(lc.shape)
        # out.fill(np.nan) # adding this sentence will cause out nan

        for i in np.arange(dist.size):
            out[idx[i, 0]][idx[i, 1]] = dist[i]

    else:
        # the initial landscape has no misc2chp, thus distance is so large, assume to be 1000000000
        # the initial landscape has no misc2bm, thus distance is so large, assume to be 1000000000
        out = np.zeros(lc.shape)

        for i in np.arange(len(fodder_arr)):
            out[fodder_arr[i, 0]][fodder_arr[i, 1]] = 1000000000 

        
        # engy_code = np.unique(lc)

        all_code = np.array([1, 2, 3, 4])
        # # print(all_code)
        engy_code = np.delete(all_code, np.where(np.logical_or(all_code == 1, all_code == 2)))
        engy_code = np.delete(engy_code, np.where(engy_code == code))
        # print(engy_code)
        engy = np.where(lc == engy_code)

        engy_arr = np.array((engy[0], engy[1])).transpose()
        
        for i in np.arange(len(engy_arr)):
            out[engy_arr[i, 0]][engy_arr[i, 1]] = 1000000000 



    return out

