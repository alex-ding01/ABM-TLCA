
import argparse
from datetime import datetime
import os
import sys
import subprocess

import netCDF4 as netcdf
import pandas as pd
import numpy as np
import gdal, osr
import matplotlib.pyplot as plt
import scipy.special as sp

import crop_functions.crop_decider as crpdec
import initialize_agents_domain as init_agent

import preprocessing.geofxns as gf
import preprocessing.get_agent_data as get_agent
import preprocessing.score_cult as sc

import postprocessing.lca_profiling as tlca
import postprocessing.save_tlca as saving

import agents.d_cell as cell

import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

import brightway2 as bw
from config_reader import ConfigReader

try:
    import pkg_resources
except ImportError:
    pass


class Janus:

    def __init__(self, config_file=None, args=None, save_result=True, plot_results=True):

        self.c = ConfigReader(config_file)

        # initialize landscape and domain
        self.scenario, self.lc, self.yd, self.domain, self.dist2misc4bm, self.dist2misc4chp, self.Ny, self.Nx = self.initialize_landscape_domain()

        # # initialize crops
        self.num_crops, self.crop_ids, self.crop_id_all = self.initialize_crops()

        # # initialize profits
        self.profit_signals, self.profits_actual, self.crop_use_ids, self.num_crops_use = self.initialize_profit()

        # # initialize agents
        self.agent_domain, self.agent_array = self.initialize_agents()

        for r in range(0, self.c.nr):

            print(f'running {r} time')

            # make agent decisions
            self.decisions()

            # tlca profiling
            self.tlca_all, self.luf_all = self.tlca_profiling()
            
            # plot the land use change results
            self.plot_results(running_time=r)

            # aggregate
            self.env_ind_year, self.luf_ind_year, self.area_id_year = self.aggregate()

	        # save results
            np.savetxt(os.path.join(self.c.output_dir, f'envi_r{r}.txt'), self.env_ind_year)
            np.savetxt(os.path.join(self.c.output_dir, f'lanf_r{r}.txt'), self.luf_ind_year)
            np.savetxt(os.path.join(self.c.output_dir, f'area_r{r}.txt'), self.area_id_year)

    def initialize_landscape_domain(self):

        scenario = int(self.c.profits_file[-5]) 

        # import the initial land cover data
        lc_raster = gdal.Open(self.c.f_init_lc_file)
        lc = lc_raster.GetRasterBand(1).ReadAsArray()

        yd_raster = gdal.Open(self.c.f_yd_file)
        yd = yd_raster.GetRasterBand(1).ReadAsArray()

        ny, nx = lc.shape

        dist2bgmisc = gf.min_dist(lc, 3) # the initial landscape has no misc4bm, thus distance is so large
        dist2chpmisc = gf.min_dist(lc, 4) # the initial landscape has misc4chp, thus distance

        # rename(factory number to 6, 7) for simplicity and plotting reason. For some reason, cannot plot 10 digit value
        lc[lc == 21] = 6
        lc[lc == 22] = 7
        lc[lc == -9999] = 0
        lc[lc == 23] = 0
        lc[np.isnan(lc)] = 0

       	if scenario == 1: 
            selected_region = self.c.regions[0]

            lc = init_agent.initial_demonstration_area(lc, self.c.n_regions_row, self.c.n_regions_col, 
                self.c.batch_length, self.c.batch_width, self.c.n_random_blocks, self.c.n_random_cells, 3, selected_region)
	        
        elif scenario == 2: 
            selected_region = self.c.regions[1]

            lc = init_agent.initial_demonstration_area(lc, self.c.n_regions_row, self.c.n_regions_col, 
                self.c.batch_length, self.c.batch_width, self.c.n_random_blocks, self.c.n_random_cells, 3, selected_region)
	        

        # initialize domain: create empty domain array
        domain = init_agent.initialize_domain(ny, nx)


        return scenario, lc, yd, domain, dist2bgmisc, dist2chpmisc, ny, nx

    def initialize_crops(self):


        forage_energy = np.where(np.logical_or(self.c.key_file['LABEL'] == 'forage', self.c.key_file['LABEL'] == 'energy-crop'))

        crop_ids_load = np.unique(np.int64(self.c.key_file['GROUPE_COD'][forage_energy[0]]))

        num_crops = len(np.unique(crop_ids_load))

        crop_ids = crop_ids_load.reshape(num_crops, 1)

        crop_id_all = np.zeros((self.c.Nt, self.Ny, self.Nx))

        crop_id_all[0, :, :] = self.lc # initial year is current land cover

        crop_id_all[crop_id_all == 6] = 0  # 6 and 7 are the lcoations of factory
        crop_id_all[crop_id_all == 7] = 0


        return num_crops, crop_ids, crop_id_all

    def initialize_profit(self):


        profits = pd.read_csv(self.c.profits_file, header=0)
        profit_signals = np.transpose(profits.values)

        crop_use = pd.read_csv(self.c.profits_file)

        crop_use_ids_lst = [int(c) for c in crop_use.columns]

        num_crops_use = len(crop_use_ids_lst)

        crop_use_ids = np.reshape(crop_use_ids_lst, (num_crops_use, 1))

        profits_actual = init_agent.init_profits(profit_signals, self.c.Nt, self.Ny, self.Nx, self.crop_id_all, self.crop_ids)


        return profit_signals, profits_actual, crop_use_ids, num_crops_use

    def initialize_agents(self):

        # size is a dataframe indicating size category and number of farms in that categories
        sizes = get_agent.sizes()

        # size cdf array, [[0, 0.02], [1, 0.3], [2, 0.9], [3, 1]]
        size_cdf = get_agent.make_size_cdf(sizes)

        # define agent type at each location: 'f_grass' fodder grass, 'f_annual' fodder annual, 'e_misc' energy miscanthus 
        agent_array = init_agent.place_agents(self.Ny, self.Nx, self.lc, self.c.key_file)

        # define agent attribute for each agent in the agent_array
        agent_domain = init_agent.agents(agent_array, self.domain,
                                         size_cdf, self.c.switch_size,
                                         self.Ny, self.Nx, self.lc)


        return agent_domain, agent_array

    def decisions(self):

        print('finish initialization, now decide')
        print('initial land use is: ')
        print(self.crop_id_all[0, :, :])

        for i in np.arange(1, self.c.Nt):

            print('#########################')
            print(f'for the {i} year')

            self.dist2bgmisc = gf.min_dist(self.crop_id_all[i-1, :, :], 3)
            self.dist2chpmisc = gf.min_dist(self.crop_id_all[i-1, :, :], 4)

            scenario = int(self.c.profits_file[-5])

            for j in np.arange(self.Ny):

                for k in np.arange(self.Nx):

                    
                    # scenario 1: subsidy zone north
                    # scenario 2: subsidy zone south

                    if self.agent_domain[j, k].FarmerAgents:
                    
                        # assess profit
                        profit_last, profit_pred = crpdec.assess_profit_potential(self.crop_id_all[i-1, j, k],
                                                                       self.profits_actual[i-1, j, k],
                                                                       self.profit_signals[:, i],
                                                                       self.c.y_average,
                                                                       self.yd[j, k],
                                                                       self.num_crops_use,
                                                                       self.crop_use_ids)

                        profit_neighbor_bgmisc = crpdec.assess_profit_neighbor(
                                                                        j, k, 
                                                                        self.crop_id_all[i-1, :, :], 
                                                                        self.profits_actual[i-1, :, :], 
                                                                        self.c.close_bar, 
                                                                        3)
                        profit_neighbor_chpmisc = crpdec.assess_profit_neighbor(
                                                                        j, k, 
                                                                        self.crop_id_all[i-1, :, :], 
                                                                        self.profits_actual[i-1, :, :], 
                                                                        self.c.close_bar, 
                                                                        4)

                        alpha_f_bgmisc, beta_f_bgmisc = crpdec.define_familiarity(self.dist2bgmisc[j, k], self.c.switch_size, self.c.close_bar, self.c.far_bar)
                        alpha_f_chpmisc, beta_f_chpmisc = crpdec.define_familiarity(self.dist2chpmisc[j, k], self.c.switch_size, self.c.close_bar, self.c.far_bar)

                        # identify the most profitable crop

                        crop_choice, profit_choice = crpdec.profit_maximizer(
                                                                    self.agent_domain[j, k].FarmerAgents[0].SizeAlpha, 
                                                                    self.agent_domain[j, k].FarmerAgents[0].SizeBeta,
                                                                    alpha_f_bgmisc, beta_f_bgmisc, 
                                                                    alpha_f_chpmisc, beta_f_chpmisc,
                                                                    self.c.fmin,
                                                                    self.c.fmax,
                                                                    self.c.n,
                                                                    profit_last,
                                                                    self.crop_use_ids,
                                                                    profit_pred,profit_neighbor_bgmisc, profit_neighbor_chpmisc
                                                                    )

                        crop_id_new = self.crop_id_all[i-1, j, k]

                        loc = np.where(self.crop_use_ids == crop_id_new)

                        # decide whether to switch and add random variation to actual profit
                        self.crop_id_all[i, j, k], self.profits_actual[i, j, k] = crpdec.make_choice(self.crop_id_all[i-1, j, k],
                                                                                                    self.profit_signals[loc[0][0], i],
                                                                                                    crop_choice,profit_choice, self.c.success_rate,
                                                                                                    self.c.fail_cost,
                                                                                                    seed = False)

            print('land use decisions:')        
            print(self.crop_id_all[i, :, :])
            
    def tlca_profiling(self):
        print('tlca profiling')

        [bidb, eidb] = sc.ecoinvent_set_up()
        licor_database = bw.Database('LiCOR TLCA')

        tlca_all = np.zeros((self.c.Ni, self.c.Nt, self.Ny, self.Nx))
        luf_all = np.zeros((len(self.c.product), self.c.Nt, self.Ny, self.Nx))

        dcf = sc.emis_cf(self.c.dcf_file)

        indicator_lst = []
        indicator_abb_lst = []

        for m in range(self.c.Ni):

            indicators = dcf['lcia_name'].unique()[m]

            # print(indicators)
            indicator_lst.append(indicators)

            indicator_abb = sc.get_initials(indicators)
            indicator_abb_lst.append(indicator_abb)


        for t in np.arange(0, self.c.Nt):

            for m in np.arange(0, self.c.Ni):

                # indicator
                tlca_all[m, t, :, :] = tlca.tlca_current(self.crop_id_all, t, self.c.boundary, indicator_abb_lst[m])
                

            for n in np.arange(0, len(self.c.product)-1):

                luf_all[n, t, :, :] = tlca.energy_current(self.crop_id_all, t, self.c.product[n])

            luf_all[-1, t, :, :] = tlca.lsu_current(self.crop_id_all, t)

        return tlca_all, luf_all

    def aggregate(self):
        # aggregate at the territorial level
        env_ind_year = np.zeros((self.c.Ni, self.c.Nt))
        luf_ind_year = np.zeros((len(self.c.product), self.c.Nt))
        area_id_year = np.zeros((self.num_crops_use, self.c.Nt))

        for t in range(self.c.Nt):

            for i in range(self.c.Ni):
                env_ind_year[i, t] = self.tlca_all[i, t].sum()

            for i in range(len(self.c.product)):
                luf_ind_year[i, t] = self.luf_all[i, t].sum()

            cur_crop = self.crop_id_all[t, :, :]

            for c in np.arange(self.num_crops_use):
                bools = (cur_crop == self.crop_use_ids[c])
                area_id_year[c, t] = np.sum(bools)

        return env_ind_year, luf_ind_year, area_id_year


if __name__ == '__main__':


    res = Janus(config_file = '/Users/tianran/Desktop/dtlca/data/configtd.yml')
