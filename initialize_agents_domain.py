"""
Created on Mon Aug 12 11:15:12 2019

@author: kek25
"""
import numpy as np

import agents.farmer as farmer
import agents.d_cell as cell
import agents.company as company
import preprocessing.get_agent_data as getAgent
from numpy.lib.stride_tricks import as_strided
import random
from functools import reduce 


def initialize_domain(ny, nx):
    """ Create empty domain array

    :param ny:  Number of columns in domain
    :type ny:   Int

    :param nx:  Number of rows in domain
    :type nx:   Int

    :return:    Empty numpy array filled with class Dcell at each pixel
    :type:      Numpy Array
    """
    domain = np.empty((ny, nx), dtype=object)

    for i in np.arange(ny):

        for j in np.arange(nx):

            domain[i][j] = cell.Dcell()

    return domain


def make_windows(x, winSizeL, winSizeW):

    # this function devide the landscape (x) into windows
    # winSizeL: window size Length
    # winSizeW: window size Width

    # return window indexers and windows

    indexer = []
    windows = []           
   
    end_l = x.shape[0]+1
    end_w = x.shape[1]+1
    
    for i in range(0, end_l, winSizeL):
        for j in range(0, end_w, winSizeW):
            window = x[i:i+winSizeL, j:j+winSizeW] # each individual window
            
            if (i < int(x.shape[0])) & (j < int(x.shape[1])): # make sure no empty window is attached
                
                indexer.append([i, j])
                windows.append(window)
                
    return indexer, windows

def change_value_random(x, crop_id_potential, n_random_cells):

    random.seed(30)

    # this function converse random cells into new crop 
    
    conversed_x = x.copy()
    
    crop_number = np.count_nonzero(x==2)

    arr = np.where(x == 2) # assuming only previous crop is annual crop, i.e., 2
    cell_arr = np.array((arr[0], arr[1])).transpose()

    selected_cell_indexs = random.sample(range(len(cell_arr)), n_random_cells)

    for ind in selected_cell_indexs:

        conversed_x[cell_arr[ind][0], cell_arr[ind][1]] = crop_id_potential

    return conversed_x

def change_blocks_random(batches, batches_indexer, n_random_blocks, n_random_cells, crop_id_potential):

    # within the selected region, we further select n_random_blocks to converse cells
    # these n_random_blocks need to have enough area for conversion
    # batches: list of batches, calculated from make_windows
    # batches_indexer: list of batches indexers, calculated from make_windows
    # n_random_blocks: user defined to get certain batches for conversion
    # n_random_cells: user defined to get certain cells within selected batch for conversion
    # crop_id_potential: conversed energy crop name

    random.seed(30)
    
    conversed_batches = []
    conversed_batches_indexers = []
    
    potential_conversed_batches = []
    potential_conversed_batches_indexers = []
    
    for i in range(len(batches)):
        
        conversed_batch = batches[i].copy()
        conversed_batch_indexer = batches_indexer[i]
        
        crop_number = np.count_nonzero(batches[i]==2)

        if crop_number >= n_random_cells:
            
            
            potential_conversed_batches.append(batches[i])
            potential_conversed_batches_indexers.append(batches_indexer[i])
            
    selected_batches_indexers = random.sample(range(len(potential_conversed_batches_indexers)), n_random_blocks)

    for i in range(len(selected_batches_indexers)):

        selected_batch_id = selected_batches_indexers[i]

        conversed_batch= change_value_random(potential_conversed_batches[selected_batch_id], crop_id_potential, n_random_cells)
        conversed_batch_indexer = potential_conversed_batches_indexers[selected_batch_id]

        conversed_batches.append(conversed_batch)
        conversed_batches_indexers.append(conversed_batch_indexer)

    return conversed_batches_indexers, conversed_batches
        
def make_region_h_data(region, batch_length, batch_width, n_random_blocks, n_random_cells, crop_id_potential, selected=True):

    # this function stack conversed blocks with non-conversed blocks into regions on the horizontal direction
    # later the resulted all regions in the horizontal direction can be vertical stacked into land scape

    # region: one of the regions that are  divided from land scape. 
    # batch_length, batch_width: the region is devided by into batches with batch_length, batch_width

    
    n_batchs_row = int(region.shape[0]/batch_length)
    n_batchs_col = int(region.shape[1]/batch_width)

    batches_indexer, batches = make_windows(region, batch_length, batch_width)

    if selected:

        new_batches = batches.copy()

        conversed_batches_indexers, conversed_batches = change_blocks_random(batches, batches_indexer, n_random_blocks, n_random_cells, crop_id_potential)

        for i in range(len(batches_indexer)):
            for j in range(len(conversed_batches_indexers)):
                if batches_indexer[i] == conversed_batches_indexers[j]:
                    new_batches[i] = conversed_batches[j]

    else:
        new_batches =  batches.copy()

    region_chunks = [new_batches[x:x+n_batchs_col] for x in range(0, len(new_batches), n_batchs_col)]
    

    h_data_lst = []
    for chunk in region_chunks:
        h_data = reduce(lambda a, b: np.hstack([a, b]), chunk)
        h_data_lst.append(h_data)
        
    new_region = reduce(lambda a, b: np.vstack([a, b]), h_data_lst)

    return new_region

def initial_demonstration_area(x, n_regions_row=2, n_regions_col=1, batch_length=2, batch_width=2, n_random_blocks=1, n_random_cells=1, crop_id_potential=3, region_selected=[0, 0]):
    
    # batch_length and batch_width requirement: x.shape[0]/batch_length and x.shape[1]/batch_width can be divided by 2.
    regions_indexer, regions = make_windows(x, int(x.shape[0]/n_regions_row), int(x.shape[1]/n_regions_col))
    new_regions = regions.copy()    
    
    region_indexers = []
    for i in range(n_regions_row):
        for j in range(n_regions_col):
            region_indexers.append([i, j])
            
    new_regions_v = []
    for i in range(n_regions_row):
        
        new_regions_h = []
        for j in range(n_regions_col):
            
            if [i, j] == region_selected:
                selected = True        
            else:
                selected = False
                
            region_number = region_indexers.index([i, j])
            
            new_region = make_region_h_data(regions[region_number], batch_length, batch_width, n_random_blocks, n_random_cells, crop_id_potential, selected)
            new_regions_h.append(new_region)
        
        new_regions_v.append(reduce(lambda a, b: np.hstack([a, b]), new_regions_h))
        

    new_lc = reduce(lambda a, b: np.vstack([a, b]), new_regions_v)
    
    return new_lc

def place_agents(ny, nx, lc, key_file):

    """ Place agents on the landscape based on land cover and associated categorization
    """
    agent_array = np.empty((ny, nx), dtype='U10')

    agent_cat = key_file['GROUPE_CUL_E']
    code = key_file['GROUPE_COD']
    
    # array([1])
    forage_grass = np.array(code[agent_cat == 'forage production grassland']).astype(int)

    # array([2])
    forage_annual = np.array(code[agent_cat == 'forage production annual crop']).astype(int)
    # energy_maiz = np.array(code[agent_cat == 'maize']).astype(int)

    # array([4])
    energy_misc = np.array(code[agent_cat == 'miscanthus']).astype(int)


    for i in forage_grass:
        agent_array[lc == i] = 'f_grass' # fodder grass

    for i in forage_annual:
        agent_array[lc == i] = 'f_annual' # fodder annual

    # for i in energy_maiz:
    #     agent_array[lc == i] = 'e_maize' # energy maize

    for i in energy_misc:
        agent_array[lc == i] = 'e_misc' # energy miscanthus 

    if 3 in np.unique(lc):
        agent_array[lc == 3] = 'e_misc' # 3 represent biogas, which is new 
    

    agent_array[lc == 6] = 'c_bm' # company for biogas
    agent_array[lc == 7] = 'c_chp'


    return agent_array

def agents(agent_array, domain, size_cdf, switch_size, ny, nx, lc):

    for i in np.arange(ny):

        for j in np.arange(nx):

            if (agent_array[i][j] == 'f_grass')|(agent_array[i][j] == 'f_annual')|(agent_array[i][j] == 'e_misc'):


                # a dictionary of farmer attributes
                agent_data = getAgent.farmer_data(size_cdf, switch_size)
            
                # a structure of agent type that will be added to the Dcell
                new_agent = farmer.Farmer(CropType = agent_array[i][j], 
                                          FarmSize=agent_data["FarmSize"], 

                                          nFields=agent_data['nFields'],

                                          SizeAlpha=agent_data['SizeAlpha'],
                                          SizeBeta=agent_data['SizeBeta']
                                          )  # this is passing actual agent data


                domain[i][j].add_agent(new_agent)

            elif (agent_array[i][j] == 'c_bm')|(agent_array[i][j] == 'c_chp'):


                agent_data = getAgent.company_data(lc[i][j])

                # print(agent_data)

                new_agent = company.Company(Product=agent_data["Product"])

                domain[i][j].add_agent(new_agent)



    return domain

def init_profits(profit_signals, nt, ny, nx, crop_id_all, crop_ids):
    """Initialize np array of profits

    """

    profits_actual = np.zeros((nt, ny, nx))

    # print(profits_actual)

    for i in np.arange(ny):

        for j in np.arange(nx):

            crop_ind = crop_id_all[0, i, j]
            crop_ix = np.where(crop_ids == crop_ind)

            # print(crop_ind)
            # print(crop_ix)

            if crop_ind in crop_ids:

                if crop_ind == 4:  # initial miscanthus is used as combustion
                    profits_actual[0, i, j] = profit_signals[-1, 0]

                else:
                    profits_actual[0, i, j] = profit_signals[crop_ix[0][0], 0]

            else:
                profits_actual[0, i, j] = 0

    # print(profits_actual)

    return profits_actual
