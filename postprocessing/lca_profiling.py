

import gdal
import sys
import os
import numpy as np
sys.path.append('/Users/tianran/Desktop/dtlca')

from config_reader import ConfigReader

config_file = '/Users/tianran/Desktop/dtlca/data/configtd.yml'
c = ConfigReader(config_file)

splca_folder = c.lca_folder


def energy_current(crop_id_all, nt,  energy_product):

	splca_new_folder = os.path.join(splca_folder, 'TLCA_new')

	crop_id_nt = crop_id_all[nt, :, :]

	energy_tot = np.zeros(crop_id_nt.shape)

	if energy_product == 'biog':
		array = np.where(crop_id_nt==3)
		array_len = len(array[0])
	else:
		array = np.where(crop_id_nt==4)
		array_len = len(array[0])

	if array_len != 0:
		# the LUF from misc4biog(3) and misc4comb(4)
		for filename in os.listdir(splca_new_folder):

			# (filename[4]==4) & (filename[6:-4]==energy_product) miscanthus converse to energy product
			if (filename[-3::] == 'txt') & (filename[6:-4]==energy_product) & (filename[4]=='4'):

				luf_raster = gdal.Open(f'{splca_new_folder}/{filename}')
				luf = luf_raster.GetRasterBand(1).ReadAsArray()

		for p in range(array_len):
			energy_tot[array[0][p], array[1][p]] = luf[array[0][p], array[1][p]]

	return energy_tot

def lsu_current(crop_id_all, nt):

	folder = os.path.join(splca_folder, 'TLCA_pre')

	crop_id_nt = crop_id_all[nt, :, :]

	lsu_tot = np.zeros(crop_id_nt.shape)

	array = np.where(np.logical_or(crop_id_nt==1, crop_id_nt==2))
	array_len = len(array[0])

	if array_len != 0:
		# the LUF from misc4biog(3) and misc4comb(4)
		for filename in os.listdir(folder):

			if (filename[-3::] == 'txt') & (filename[0:3]=='luf'):

				luf_raster = gdal.Open(f'{folder}/{filename}')
				luf = luf_raster.GetRasterBand(1).ReadAsArray()

		for p in range(array_len):
			lsu_tot[array[0][p], array[1][p]] = luf[array[0][p], array[1][p]]

	return lsu_tot




def tlca_current(crop_id_all, nt, in_or_of, indicator):

	"""
	indicator: abbreviated indicator name
	in_or_of: in, of, tot
	energy_product: biogas_m3pha, electricity_kwhpha, heat_mjpha

	in this work, only consider 1 crop, 2 pathways. 
	- miscanthus2combustion, miscanthus2biogas (4_elec, 4_heat, 4_biog)


	return: spatialized e-tlca, spatialized energy product
	"""

	splca_pre_folder = os.path.join(splca_folder, 'TLCA_pre')
	splca_new_folder = os.path.join(splca_folder, 'TLCA_new')

	crop_id_nt = crop_id_all[nt, :, :]


	# the basic land use map is grassland (1), annual (2) and misc4comb (4)
	# the new conversion could be misc4biog (3)

	def tlca(boundary):

		# the basic impact from 1 and 2
		for filename in os.listdir(splca_pre_folder):
			# find the lc_lca for the needed indicator
			if ('.txt' in filename) & (filename[3:-4]==indicator) & (filename[0:2]==boundary):

				lca_raster = gdal.Open(f'{splca_pre_folder}/{filename}')
				lca = lca_raster.GetRasterBand(1).ReadAsArray()

		tlca_tot = lca.copy()

		# the impact from misc4comb(4)
		array = np.where(crop_id_nt==4)
		array_len = len(array[0])

		if array_len != 0:  # lc map has that crop
			for filename in os.listdir(splca_new_folder):

				# (filename[3]==4) & (filename[12]==0) miscanthus converse to direct combustion
				if (filename[-3::] == 'txt') & (filename[14:-4]==indicator) & (filename[0:2]==boundary) & (filename[3]=='4'):

					if filename[12] == '0':

						lca_raster = gdal.Open(f'{splca_new_folder}/{filename}')
						lca = lca_raster.GetRasterBand(1).ReadAsArray()

			for p in range(array_len):
				tlca_tot[array[0][p], array[1][p]] = lca[array[0][p], array[1][p]]


		# the impact from misc4biog(3)
		array = np.where(crop_id_nt==3)
		array_len = len(array[0])

		if array_len != 0:


			for filename in os.listdir(splca_new_folder):

				# (filename[3]==4) & (filename[12]==1) miscanthus converse to biogas
				if (filename[-3::] == 'txt') & (filename[14:-4]==indicator) & (filename[0:2]==boundary) & (filename[3]=='4'):

					if (filename[12]=='1'):

						lca_raster = gdal.Open(f'{splca_new_folder}/{filename}')
						lca = lca_raster.GetRasterBand(1).ReadAsArray()


			for p in range(array_len):
				tlca_tot[array[0][p], array[1][p]] = lca[array[0][p], array[1][p]]

		return tlca_tot


	if in_or_of == 'in':
		tlca_result = tlca('in')

	elif in_or_of == 'tot':
		tlca_result = tlca('in') + tlca('of')


	return tlca_result

