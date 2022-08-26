import gdal
import sys
import os
from osgeo import ogr, gdal, osr
import rasterio as rio  
sys.path.append('/Users/tianran/Desktop/dtlca')
import preprocessing.geofxns as gf

os.chdir("/Users/tianran/OneDrive - Université Libre de Bruxelles/LiCOR/opt2/Script/paper4_v4")
# print(os.getcwd())

import score_cult as c
# print(c.__file__)

import numpy as np
import pandas as pd
import geopandas as gpd
import brightway2 as bw

pd.set_option("display.max_columns", 500)

crs_licor = 31370
misc_yield_wet = 16459/0.42 # 16459/0.42 kg/ha


shp_folder = '/Users/tianran/Desktop/janustd/data_td/Input4stlca/shp'
rast_folder = '/Users/tianran/Desktop/janustd/data_td/Input4stlca/rasters'

txt_folder = '/Users/tianran/Desktop/janustd/data_td/Input4stlca/rasters/txt'
# txt_folder = '/Users/tianran/Desktop/janustd/data_td/Input4stlca/rasters/sample'

output_folder = '/Users/tianran/Desktop/janustd/data_td/STLCA/full'
# output_folder = '/Users/tianran/Desktop/janustd/data_td/STLCA/sample'

nuts1 = gpd.read_file('data/wl_n1/l1_wl.shp')

lc_fp = os.path.join(txt_folder, 'txt_lc_fac.txt')

lc_raster = gdal.Open(lc_fp)
lc = lc_raster.GetRasterBand(1).ReadAsArray()
lc[lc == 21] = 6
lc[lc == 22] = 7
lc[lc == -9999] = 0
lc[lc == 23] = 0

print(np.unique(lc))

# potential locations is 1, otherwise is 0
lc_new = lc

lc_new[lc_new == 1] = 1
lc_new[lc_new == 2] = 1
lc_new[lc_new == 4] = 1
lc_new[lc_new == 6] = 0
lc_new[lc_new == 7] = 0
lc_new[lc_new == 21] = 0
lc_new[lc_new == 22] = 0
lc_new[lc_new == 23] = 0

lc_new[lc_new == -9999] = 0

[bidb, eidb] = c.ecoinvent_set_up()
licor_database = bw.Database('LiCOR TLCA')

iw_method = []
for m in bw.methods:
	if ('IMPACTWorld' in m[0]) & ('Endpoint' in m[0]):
		iw_method.append(m)

dcf = c.emis_cf('data/Brightway_IW_damage_1_41.bw2package')


imp = bw.ExcelImporter("/Users/tianran/Desktop/miscspa/inventory.xlsx")
imp.apply_strategies()
imp.match_database("ecoinvent 3.5 cutoff", fields=('name', 'unit', 'filename','location', 'reference product'))
imp.match_database("ecoinvent 3.5 cutoff", fields=('name', 'unit', 'location', 'reference product'))
imp.match_database("ecoinvent 3.5 cutoff", fields=('name', 'unit', 'location'))

imp.match_database(fields=['name', 'location'])
imp.statistics()
imp.write_database()


inv_misc_cult = [x for x in licor_database if ('miscanthus cultivation' in x['name'])][0]
inv_maiz = [x for x in licor_database if ('land use maize' in x['name'])][0]
inv_gras = [x for x in licor_database if ('land use grass' in x['name'])][0]

ef_misc = c.cal_tef_misc(inv_misc_cult, iw_method)
ef_maiz = c.cal_tef_pre(inv_maiz, iw_method)
ef_gras = c.cal_tef_pre(inv_gras, iw_method)

trans_name = 'transport, tractor and trailer, agricultural'
inv_org = [act for act in eidb if trans_name in act['name'] and act['location'] == 'RoW'][0]  # FU=1tkm

inv_comb = [x for x in licor_database if ('combustion' in x['name'])][0]
inv_ad = [x for x in licor_database if ('anaerobic digestion' in x['name'])][0]

ef = c.cal_ef_truck(inv_org, trans_name, iw_method)
tef_comb = c.cal_tef3(inv_comb, iw_method) # combusting 2024kg miscanthus 42% DM
tef_ad = c.cal_tef3(inv_ad, iw_method) # combusting 2024kg miscanthus 42% DM

# tools and preparing data to generate spatial tlca
###################################################
###################################################
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj][0]

def get_initials(name):

	names = name.split(' ')

	initials = ''.join([f"{n[0]}" for n in names])
	
	return initials.lower()

def save_txt(meta_src_fp, dst_fp, dst_data):
		with rio.open(meta_src_fp) as src:
			ras_meta = src.profile

		ras_meta['dtype'] = "float32"

		with rio.open(dst_fp, 'w', **ras_meta) as dst:
			dst.write(dst_data, 1)

def shp_to_tiff(InputShape, OutputRasterPath, OutputName):

	"""
	OutputName: name.tiff
	"""

	# 1) get shapefile information

	# making the shapefile as an object
	sources_ds = ogr.Open(InputShape)

	# getting layer information of shapefile
	source_layer = sources_ds.GetLayer()

	# 2) creating the destination raster data source

    # pixel size is proportional to size of shapefile
	pixelWidth = pixelHeight = 100 # depending how fine you want your raster

	# get extent values to set size of output raster
	x_min, x_max, y_min, y_max = source_layer.GetExtent()

	# calculate the resolution of the raster
	cols = int((x_max - x_min) / pixelHeight)
	rows = int((y_max - y_min) / pixelWidth)

	# define the raster file path including the name and extension (tiff)
	OutRaster = os.path.join(OutputRasterPath, OutputName)

	# get GeoTiff dryver, create the raster dataset
	target_ds = gdal.GetDriverByName('GTiff').Create(OutRaster, cols, rows, 1, gdal.GDT_Float32)

	# transforms between pixel raster space to projection coordinate space
	# x_upperlimit, x_size, x_skew, y_upperlimit, y_skew, y_size
	target_ds.SetGeoTransform((x_min, pixelWidth, 0, y_max, 0, -pixelHeight))

	# 3) adding the spatial reference
	target_dsSRS = osr.SpatialReference()
	target_dsSRS.ImportFromEPSG(31370)
	target_ds.SetProjection(target_dsSRS.ExportToWkt())

	# get required raster band
	band = target_ds.GetRasterBand(1)
	band.SetNoDataValue(-9999)

	# main conversion method
	gdal.RasterizeLayer(target_ds, [1], source_layer, options=["ATTRIBUTE=Val"])

	targets_ds = None

	return OutRaster

def get_splca_misc_cult_inoff(indicator):

	"""Cultivation uses the provided one from ecoinvent database. 
	To link with the conversion step provided by (Nguyen and Hermansen 2015), 
	I assume that the produced miscanthus chopped in ecoinvent represent 1 kg 42% DM and 
	the yield is 16459/0.42 kg/ha.
	the avoided emissions from fertilizers for anaerobic digestration process is included in the conversion process"""

	dcf_ind = dcf[dcf.lcia_name == indicator]

	# ef_misc.to_excel("/Users/tianran/Desktop/janustd/data_td/Input4stlca/misc_emis_pkg.xlsx")
	
	# ef_misc = pd.read_csv("/Users/tianran/Desktop/janustd/data_td/Input4stlca/misc_emis_pkg.csv")
	# ef_misc.drop(columns=['Unnamed: 0', 'crop_name', 'inpt_name', 'inpt_type', 'inpt_amnt'], inplace=True)

	ef_misc_in = ef_misc[ef_misc.emis_sorc == 'in-territory emission']
	ef_misc_off = ef_misc[ef_misc.emis_sorc == 'off-territory emission']

	tlca_dcf = ef_misc.merge(dcf_ind, on=['emis_name', 'emis_cate'])
	tlca_dcf_in = tlca_dcf[tlca_dcf.emis_sorc == 'in-territory emission']
	tlca_dcf_off = tlca_dcf[tlca_dcf.emis_sorc == 'off-territory emission']


	if indicator == 'Freshwater acidification':

		NH3 = 'Ammonia'
		HNO3 = 'Nitrate'
		NOX = 'Nitrogen oxides'
		SO2 = 'Sulfur dioxide'
		SO4 = 'Sulfate'

		NH3_in = ef_misc_in[ef_misc_in.emis_name == NH3].emis_amnt.sum()
		HNO3_in = ef_misc_in[ef_misc_in.emis_name == HNO3].emis_amnt.sum()
		NOX_in = ef_misc_in[ef_misc_in.emis_name == NOX].emis_amnt.sum()
		SO2_in = ef_misc_in[ef_misc_in.emis_name == SO2].emis_amnt.sum()
		SO4_in = ef_misc_in[ef_misc_in.emis_name == SO4].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_AcidFW_NH3') and ('.txt') in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			if ('txt_AcidFW_HNO3') and ('.txt') in filename:

				rcf_HNO3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_HNO3 = rcf_HNO3_raster.GetRasterBand(1).ReadAsArray()

				rcf_HNO3[rcf_HNO3 == -9999] = 0

			if ('txt_AcidFW_NOX') and ('.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0

			if ('txt_AcidFW_SO2') and ('.txt') in filename:

				rcf_SO2_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO2 = rcf_SO2_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO2[rcf_SO2 == -9999] = 0

			if ('txt_AcidFW_SO4') and ('.txt') in filename:

				rcf_SO4_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO4 = rcf_SO4_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO4[rcf_SO4 == -9999] = 0

		tlca_in = (NH3_in*rcf_NH3 + HNO3_in*rcf_HNO3 + NOX_in*rcf_NOX + SO2_in*rcf_SO2 + SO4_in*rcf_SO4)

	elif indicator == 'Terrestrial acidification':

		NH3 = 'Ammonia'
		NOX = 'Nitrogen oxides'
		SO2 = 'Sulfur dioxide'

		NH3_in = ef_misc_in[ef_misc_in.emis_name == NH3].emis_amnt.sum()

		NOX_in = ef_misc_in[ef_misc_in.emis_name == NOX].emis_amnt.sum()

		SO2_in = ef_misc_in[ef_misc_in.emis_name == SO2].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_AcidTerr_NH3') and ('.txt') in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			if ('txt_AcidTerr_NOX') and ('.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0

			if ('txt_AcidTerr_SO2') and ('.txt') in filename:

				rcf_SO2_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO2 = rcf_SO2_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO2[rcf_SO2 == -9999] = 0


		tlca_in = (NH3_in*rcf_NH3 + NOX_in*rcf_NOX + SO2_in*rcf_SO2)

	elif indicator == 'Freshwater eutrophication':

		BOD = 'BOD5, Biological Oxygen Demand'
		COD = 'COD, Chemical Oxygen Demand'
		Phosphate = 'Phosphate'
		Pacid = 'Phosphoric acid'
		Phosphorus = 'Phosphorus'

		BOD_in = ef_misc_in[ef_misc_in.emis_name == BOD].emis_amnt.sum()

		COD_in = ef_misc_in[ef_misc_in.emis_name == COD].emis_amnt.sum()

		Phosphate_in = ef_misc_in[ef_misc_in.emis_name == Phosphate].emis_amnt.sum()

		Pacid_in = ef_misc_in[ef_misc_in.emis_name == Pacid].emis_amnt.sum()

		Phosphorus_in = ef_misc_in[ef_misc_in.emis_name == Phosphorus].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_EutroFW_BOD.txt') in filename:

				rcf_BOD_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_BOD = rcf_BOD_raster.GetRasterBand(1).ReadAsArray()
				rcf_BOD[rcf_BOD == -9999] = 0

			elif ('txt_EutroFW_COD.txt') in filename:

				rcf_COD_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_COD = rcf_COD_raster.GetRasterBand(1).ReadAsArray()

				rcf_COD[rcf_COD == -9999] = 0

			elif ('txt_EutroFW_PHOSACID.txt') in filename:

				rcf_Pacid_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Pacid = rcf_Pacid_raster.GetRasterBand(1).ReadAsArray()

				rcf_Pacid[rcf_Pacid == -9999] = 0

			elif ('txt_EutroFW_PHOSPHORUS.txt') in filename:

				rcf_Phosphorus_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Phosphorus = rcf_Phosphorus_raster.GetRasterBand(1).ReadAsArray()

				rcf_Phosphorus[rcf_Phosphorus == -9999] = 0

			elif ('txt_EutroFW_PHOSPHATE.txt') in filename:

				rcf_Phosphate_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Phosphate = rcf_Phosphate_raster.GetRasterBand(1).ReadAsArray()

				rcf_Phosphate[rcf_Phosphate == -9999] = 0

		tlca_in = (BOD_in*rcf_BOD + COD_in*rcf_COD + Pacid_in*rcf_Pacid + Phosphorus_in*rcf_Phosphorus + Phosphate_in*rcf_Phosphate)

	elif indicator == 'Marine eutrophication':

		NH3 = 'Ammonia'
		NOX = 'Nitrogen oxides'

		NH3_in = ef_misc_in[ef_misc_in.emis_name == NH3].emis_amnt.sum()

		NOX_in = ef_misc_in[ef_misc_in.emis_name == NOX].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_EutroMar_NH3.txt')in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			elif ('txt_EutroMar_NOX.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0


		tlca_in = (NH3_in*rcf_NH3 + NOX_in*rcf_NOX)

	elif 'Climate change' in indicator:

		dcf_co2 = dcf_ind[dcf.emis_name=='Carbon dioxide, from soil or biomass stock']['dcf'].values[0]

		luc_co2_raster =  gdal.Open(os.path.join(txt_folder, 'txt_luc_co2_kgpha.txt'))
		luc_co2 = luc_co2_raster.GetRasterBand(1).ReadAsArray()

		luc_co2[luc_co2 == -9999] = 0


		tlca_in = lc_new * sum(tlca_dcf_in['emis_amnt'] * tlca_dcf_in['dcf']) + luc_co2 * dcf_co2

	else:

		tlca_in = lc_new * sum(tlca_dcf_in['emis_amnt'] * tlca_dcf_in['dcf']) 

	tlca_off = lc_new * sum(tlca_dcf_off['emis_amnt'] * tlca_dcf_off['dcf'])

	# print(tlca_in)
	# print(tlca_off)
	tlca_in = np.around(tlca_in, decimals = 10)
	tlca_off = np.around(tlca_off, decimals = 10)

	# the returned result represent impact/ha
	return tlca_in, tlca_off

def get_transport_inoff(indicator):



	ef_lcia = ef.merge(dcf, on=['emis_name', 'emis_cate'])
	ef_lcia['lcia_score'] = ef_lcia.emis_amnt * ef_lcia.dcf / 1000


	ef_lcia.drop(columns=['inve_name'], inplace=True)
	ef_lcia = ef_lcia.groupby(['lcia_name', 'emis_sorc'])['lcia_score'].sum().reset_index()

	ef_lcia_ind = ef_lcia[ef_lcia.lcia_name == indicator]


	val_ser_in = ef_lcia_ind[ef_lcia_ind['emis_sorc'] == 'in-territory emission']['lcia_score'].values

	if val_ser_in.size == 0:
		val_in= 0
	else:
		val_in = val_ser_in[0]


	val_ser_off = ef_lcia_ind[ef_lcia_ind['emis_sorc'] == 'off-territory emission']['lcia_score'].values

	if val_ser_off.size == 0:
		val_off= 0
	else:
		val_off = val_ser_off[0]

	return val_in, val_off # transporting 1 kg 1 km

def get_conversion_inoff(indicator, pathway):

	"""
	pathway: "combustion" = 4 or "anaerobic digestion" = 3
	"""

	if pathway == 'combustion':

		tef_path = tef_comb

	else:
		tef_path = tef_ad

	lca_path = tef_path.merge(dcf, on=['emis_name', 'emis_cate'])

	lca_path['lcia_score'] = lca_path.emis_amnt * lca_path.dcf

	lca_path = lca_path.groupby(['emis_sorc', 'lcia_name'])['lcia_score'].sum().reset_index()

	lca_path_ind = lca_path[lca_path.lcia_name == indicator]

	val_ser_in = lca_path_ind[lca_path_ind.emis_sorc == 'in-territory emission']['lcia_score'].values

	if val_ser_in.size == 0:
		val_in= 0
	else:
		val_in = val_ser_in[0]


	val_ser_off = lca_path_ind[lca_path_ind['emis_sorc'] == 'off-territory emission']['lcia_score'].values

	if val_ser_off.size == 0:
		val_off= 0
	else:
		val_off = val_ser_off[0]

	# print(val_in/2024)


	return val_in/2024, val_off/2024  #impact conversing 1 kg of 42% DM

def grid():

	soil = gpd.read_file("data/luc/soil_type.gpkg")
	clim = gpd.read_file("data/luc/climate_zone.gpkg")

	soil.to_crs(epsg=crs_licor, inplace=True)
	clim.to_crs(epsg=crs_licor, inplace=True)
	wl_nuts1 = nuts1.to_crs(crs=crs_licor)

	# clim and soil has attribute = 0, means no data
	clim = clim[clim.DN != 0]
	soil = soil[soil.DN != 0]

	# only calculate soil, clim in all municipalities
	soil_wl = gpd.overlay(soil, wl_nuts1, how='intersection').rename(columns={'DN': 'SOIL_CODE'})
	soil_clim_wl = gpd.overlay(soil_wl, clim).rename(columns={'DN': 'CLIM_CODE'})

	soil_clim_wl.to_file(os.path.join(shp_folder, 'soil_clim.shp'))

	return soil_clim_wl

# soil_clim_wl = grid()

def get_prev_lu_type():

	fp_lc = '/Users/tianran/Desktop/janustd/data_td/LandCover/SIGEC_PARC_AGRI_ANON__2019.shp'
	fp_factory = "/Users/tianran/Desktop/janustd/data_td/Facilities/facility_shp.shp"
	categories_csv = '/Users/tianran/Desktop/janustd/data_td/crop_categories.csv'
	
	lc_shp = gpd.read_file(fp_lc).to_crs(epsg=crs_licor)

	lc_shp[['CULT_COD']] = lc_shp[['CULT_COD']].astype(int)

	key = pd.read_csv(categories_csv).loc[:, ['GROUPE_COD', 'CULT_COD']]

	ind = key['CULT_COD'].astype(int)
	grouped_ind = key['GROUPE_COD'].astype(int)

	key_dict = dict(zip(ind, grouped_ind))

	lc_shp['CULT_COD_G'] = lc_shp['CULT_COD'].map(key_dict)
	# lc_shp.rename(columns={"CULT_COD_G": 'Val'}, inplace=True)
	print(lc_shp)


	fac_shp = gpd.read_file(fp_factory).to_crs(epsg=crs_licor).loc[:, ['CULT_COD_G', 'geometry']]
	# fac_shp.rename(columns={"CULT_COD_G": 'Val'}, inplace=True)

	lc_shp.to_crs(epsg=31370, inplace=True)
	fac_shp.to_crs(epsg=31370, inplace=True)

	lc_fac_shp = pd.concat([lc_shp, fac_shp])
	lc_fac_shp = gpd.GeoDataFrame(lc_fac_shp, geometry='geometry')

	lc_fac_shp.rename(columns={"CULT_COD_G": 'Val'}, inplace=True)

	lc_clim_soil = gpd.overlay(soil_clim_wl, lc_fac_shp)

	lc_fac_shp = lc_clim_soil.loc[:, ['Val', 'geometry']]


	output_lc_fac = os.path.join(shp_folder, 'lc_fac.shp')

	lc_fac_shp.to_file(output_lc_fac)

	shp_to_tiff(output_lc_fac, rast_folder, 'lc_fac.tiff')


	# calculate the LUC score for each soil, climate and prev type.
	luc_co2_misc = pd.read_excel('data/LUC/clim_soil_prelu_misc.xlsx').dropna(axis=0)
	# unit = tonnes/ha
	luc_co2_misc = luc_co2_misc.loc[:, ['Val', 'CLIM_NAME', 'CLIM_CODE', 'SOIL_NAME', 'SOIL_CODE', 'LUC_CO2']]
	luc_co2_misc['LUC_CO2'] = luc_co2_misc['LUC_CO2'] * 1000

	luc_co2 = lc_clim_soil.merge(luc_co2_misc, on=['Val', 'CLIM_CODE', 'SOIL_CODE'], how='outer')[['geometry', 'LUC_CO2']]
	
	luc_co2.rename(columns={"LUC_CO2": 'Val'}, inplace=True)

	output_luc_co2 = os.path.join(shp_folder, 'luc_co2_kgpha.shp')

	luc_co2.to_file(output_luc_co2)
	
	shp_to_tiff(output_luc_co2, rast_folder, 'luc_co2_kgpha.tiff')



	return lc_fac_shp

def get_misc_yield():


	yild_misc_fp = 'data/Bioenergy_crop_yields.nc'
	yild_misc = c.read_yieldnc(yild_misc_fp, 'Miscanthus').to_crs(epsg=crs_licor)
	yild_misc.val = yild_misc.val/0.42 # wet matter yield
	yild_misc = gpd.overlay(nuts1, yild_misc)
	lc = gpd.read_file(os.path.join(shp_folder, 'lc_fac.shp'))

	yild_misc = gpd.overlay(lc, yild_misc).loc[:, ['val', 'geometry']]

	yild_misc.rename(columns={'val': 'Val'}, inplace=True)

	output_yield = os.path.join(shp_folder, 'misc_yield_wet_tpha.shp')

	yild_misc.to_file(output_yield)

	shp_to_tiff(output_yield, rast_folder, 'misc_yield_wet_tpha.tiff')

	return yild_misc

def get_rcf():

	cf_AcidFW_fp = 'data/rCF/SHP/AcidFW_Damage_native.shp'
	cf_AcidTerr_fp = 'data/rCF/SHP/AcidTerr_Damage_native.shp'
	cf_EutroFW_fp = 'data/rCF/SHP/EutroFW_Damage_native.shp'
	cf_EutroMar_fp = 'data/rCF/SHP/EutroMar_Damage_native.shp'

	wl_nuts1 = nuts1.to_crs(epsg=4326)

	cf_AcidFW = gpd.read_file(cf_AcidFW_fp).loc[:, ['HNO3', 'NH3', 'NOX', 'SO2', 'SO4', 'geometry']]
	cf_AcidTerr = gpd.read_file(cf_AcidTerr_fp).loc[:, ['NH3', 'NOX', 'SO2', 'geometry']]
	cf_EutroFW = gpd.read_file(cf_EutroFW_fp).loc[:, ['BOD', 'COD', 'PHOSPHATE', 'PHOSACID', 'PHOSPHORUS', 'PHOSPENTOX', 'geometry']]
	cf_EutroMar = gpd.read_file(cf_EutroMar_fp).loc[:, ['HNO3', 'NH3', 'NOX', 'geometry']]

	lc = gpd.read_file(os.path.join(shp_folder, 'lc_fac.shp'))

	rcf_lst = [cf_AcidFW, cf_AcidTerr, cf_EutroFW, cf_EutroMar]
	
	for rcf in rcf_lst:

		rcf_name = namestr(rcf, locals())[3:]
		
		rcf_wl = gpd.overlay(wl_nuts1, rcf, how='intersection').to_crs(epsg=31370)
		
		data = gpd.overlay(lc, rcf_wl)

		print(rcf_name)
		# data.to_file(os.path.join(shp_folder, f'{rcf_name}_wl.shp'))

		# for each emission's cf in rcf file, not geometry
		for col in data.columns[2:-1]:


			shp = data.loc[:, [col, 'geometry']]
			shp.rename(columns={col: 'Val'}, inplace=True)

			output = os.path.join(shp_folder, f'{rcf_name}_{col}.shp')
			
			shp.to_file(output)
			
			shp_to_tiff(output, rast_folder, f'{rcf_name}_{col}.tiff')

# execution
###################################################
###################################################

# # shp of previous crop names
# get_prev_lu_type()

# shp of luc kg/ha co2
# get_luc_co2()

# # # shp of miscanthus wet yield
# get_misc_yield()

# shp of rcf
# get_rcf()

###################################################
###################################################

# use qgis to converse tiff to txt,
# then 

# functions to generate spatial tlca
###################################################
###################################################
def get_splca_pre_inoff(indicator, crop_number):

	# this function calculates for each land use (annual crop, or grassland)
	# the emissions, calculated from previous paper
	dcf_ind = dcf[dcf.lcia_name == indicator]

	tef_maiz = ef_maiz.copy()
	tef_gras = ef_gras.copy()

	tef_maiz['PRE_TYPE'] = 2
	tef_gras['PRE_TYPE'] = 1

	tef_pre = pd.concat([tef_gras, tef_maiz])

	tef_pre_grouped = tef_pre.groupby(['PRE_TYPE', 'emis_name', 'emis_cate', 'emis_sorc'])['emis_amnt'].sum().reset_index()

	tef_crop = tef_pre_grouped[tef_pre_grouped.PRE_TYPE == crop_number]

	tef_in = tef_crop[tef_crop.emis_sorc == 'in-territory emission']


	tlca_dcf = tef_crop.merge(dcf_ind, on=['emis_name', 'emis_cate'])
	tlca_dcf_in = tlca_dcf[tlca_dcf.emis_sorc == 'in-territory emission']
	tlca_dcf_of = tlca_dcf[tlca_dcf.emis_sorc == 'off-territory emission']


	if indicator == 'Freshwater acidification':

		NH3 = 'Ammonia'
		HNO3 = 'Nitrate'
		NOX = 'Nitrogen oxides'
		SO2 = 'Sulfur dioxide'
		SO4 = 'Sulfate'


		NH3_in = tef_in[tef_in.emis_name == NH3].emis_amnt.sum()

		HNO3_in = tef_in[tef_in.emis_name == HNO3].emis_amnt.sum()

		NOX_in = tef_in[tef_in.emis_name == NOX].emis_amnt.sum()

		SO2_in = tef_in[tef_in.emis_name == SO2].emis_amnt.sum()

		SO4_in = tef_in[tef_in.emis_name == SO4].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_AcidFW_NH3') and ('.txt') in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			if ('txt_AcidFW_HNO3') and ('.txt') in filename:

				rcf_HNO3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_HNO3 = rcf_HNO3_raster.GetRasterBand(1).ReadAsArray()

				rcf_HNO3[rcf_HNO3 == -9999] = 0


			if ('txt_AcidFW_NOX') and ('.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0

			if ('txt_AcidFW_SO2') and ('.txt') in filename:

				rcf_SO2_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO2 = rcf_SO2_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO2[rcf_SO2 == -9999] = 0

			if ('txt_AcidFW_SO4') and ('.txt') in filename:

				rcf_SO4_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO4 = rcf_SO4_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO4[rcf_SO4 == -9999] = 0

		tlca_in = NH3_in*rcf_NH3 + HNO3_in*rcf_HNO3 + NOX_in*rcf_NOX + SO2_in*rcf_SO2 + SO4_in*rcf_SO4

	elif indicator == 'Terrestrial acidification':

		NH3 = 'Ammonia'
		NOX = 'Nitrogen oxides'
		SO2 = 'Sulfur dioxide'

		NH3_in = tef_in[tef_in.emis_name == NH3].emis_amnt.sum()

		NOX_in = tef_in[tef_in.emis_name == NOX].emis_amnt.sum()

		SO2_in = tef_in[tef_in.emis_name == SO2].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_AcidTerr_NH3') and ('.txt') in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			if ('txt_AcidTerr_NOX') and ('.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0

			if ('txt_AcidTerr_SO2') and ('.txt') in filename:

				rcf_SO2_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_SO2 = rcf_SO2_raster.GetRasterBand(1).ReadAsArray()

				rcf_SO2[rcf_SO2 == -9999] = 0


		tlca_in = NH3_in*rcf_NH3 + NOX_in*rcf_NOX + SO2_in*rcf_SO2

	elif indicator == 'Freshwater eutrophication':

		BOD = 'BOD5, Biological Oxygen Demand'
		COD = 'COD, Chemical Oxygen Demand'
		Phosphate = 'Phosphate'
		Pacid = 'Phosphoric acid'
		Phosphorus = 'Phosphorus'

		BOD_in = tef_in[tef_in.emis_name == BOD].emis_amnt.sum()

		COD_in = tef_in[tef_in.emis_name == COD].emis_amnt.sum()

		Phosphate_in = tef_in[tef_in.emis_name == Phosphate].emis_amnt.sum()

		Pacid_in = tef_in[tef_in.emis_name == Pacid].emis_amnt.sum()

		Phosphorus_in = tef_in[tef_in.emis_name == Phosphorus].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_EutroFW_BOD.txt') in filename:

				rcf_BOD_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_BOD = rcf_BOD_raster.GetRasterBand(1).ReadAsArray()
				rcf_BOD[rcf_BOD == -9999] = 0

			if ('txt_EutroFW_COD.txt') in filename:

				rcf_COD_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_COD = rcf_COD_raster.GetRasterBand(1).ReadAsArray()

				rcf_COD[rcf_COD == -9999] = 0

			if ('txt_EutroFW_PHOSACID.txt') in filename:

				rcf_Pacid_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Pacid = rcf_Pacid_raster.GetRasterBand(1).ReadAsArray()

				rcf_Pacid[rcf_Pacid == -9999] = 0

			if ('txt_EutroFW_PHOSPHORUS.txt') in filename:

				rcf_Phosphorus_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Phosphorus = rcf_Phosphorus_raster.GetRasterBand(1).ReadAsArray()

				rcf_Phosphorus[rcf_Phosphorus == -9999] = 0

			if ('txt_EutroFW_PHOSPHATE.txt') in filename:

				rcf_Phosphate_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_Phosphate = rcf_Phosphate_raster.GetRasterBand(1).ReadAsArray()

				rcf_Phosphate[rcf_Phosphate == -9999] = 0

		tlca_in = BOD_in*rcf_BOD + COD_in*rcf_COD + Phosphate_in*rcf_Phosphate + Pacid_in*rcf_Pacid + Phosphorus_in*rcf_Phosphorus

	elif indicator == 'Marine eutrophication':

		NH3 = 'Ammonia'
		NOX = 'Nitrogen oxides'

		NH3_in = tef_in[tef_in.emis_name == NH3].emis_amnt.sum()

		NOX_in = tef_in[tef_in.emis_name == NOX].emis_amnt.sum()

		for filename in os.listdir(txt_folder):

			if ('txt_EutroMar_NH3.txt')in filename:

				rcf_NH3_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NH3 = rcf_NH3_raster.GetRasterBand(1).ReadAsArray()

				rcf_NH3[rcf_NH3 == -9999] = 0

			if ('txt_EutroMar_NOX.txt') in filename:

				rcf_NOX_raster = gdal.Open(f'{txt_folder}/{filename}')
				rcf_NOX = rcf_NOX_raster.GetRasterBand(1).ReadAsArray()

				rcf_NOX[rcf_NOX == -9999] = 0


		tlca_in = NH3_in*rcf_NH3 + NOX_in*rcf_NOX

	else:

		tlca_in = lc_new * sum(tlca_dcf_in['emis_amnt'] * tlca_dcf_in['dcf'])

	tlca_off = lc_new * sum(tlca_dcf_of['emis_amnt'] * tlca_dcf_of['dcf'])

	tlca_in = np.float32(np.around(tlca_in, decimals = 10))
	tlca_off = np.float32(np.around(tlca_off, decimals = 10))

	output_tlca_in = os.path.join(output_folder, f'crop{crop_number}_in_{get_initials(indicator)}.txt')
	output_tlca_off = os.path.join(output_folder, f'crop{crop_number}_of_{get_initials(indicator)}.txt')

    
	save_txt(lc_fp, output_tlca_in, tlca_in)
	save_txt(lc_fp, output_tlca_off, tlca_off)

	return tlca_in, tlca_off
	
def get_splca_misc_inoff(indicator, crop_number):

	lc_fp = os.path.join(txt_folder, 'txt_lc_fac.txt')

	lc_raster = gdal.Open(lc_fp)
	lc = lc_raster.GetRasterBand(1).ReadAsArray()
	lc[lc == 21] = 6
	lc[lc == 22] = 7
	lc[lc == -9999] = 0
	lc[lc == 23] = 0

	print(lc[lc==7])

	# crop_number = 3 or 4

	# cultivation stage is the same for 3 and 4
	cultivation_in, cultivation_of = get_splca_misc_cult_inoff(indicator)
	print(cultivation_in)

	# transportation and conversion step needs to know the spatial yield
	misc_yield_rs = gdal.Open(os.path.join(txt_folder, 'txt_misc_yield_wet_tpha.txt'))
	misc_yield = misc_yield_rs.GetRasterBand(1).ReadAsArray()

	misc_yield[misc_yield == -9999] = 0

	misc_yield = misc_yield * 1000

	transport_in, transport_of = get_transport_inoff(indicator)

	conversion2bm_in, conversion2bm_of = get_conversion_inoff(indicator, 'anaerobic digestion')
	conversion2cm_in, conversion2cm_of = get_conversion_inoff(indicator, 'combustion')

	# miscanthus for biomethane
	if crop_number == 3:
		
		# cultivation + transportation + conversion2bm

		# in

		score_in = cultivation_in \
					+ misc_yield * gf.min_dist(lc, 6) * transport_in \
					+ misc_yield * conversion2bm_in

		# off
		score_of = cultivation_of \
					+ misc_yield * gf.min_dist(lc, 6) * transport_of\
					+ misc_yield * conversion2bm_of

	# miscanthus for cumbustion
	else: #(crop_number == 4) 

		# cultivation + transportation + conversion2cm

		# in
		score_in = cultivation_in \
					+ misc_yield * gf.min_dist(lc, 7) * transport_in \
					+ misc_yield * conversion2cm_in

		# off
		score_of = cultivation_of \
					+ misc_yield * gf.min_dist(lc, 7) * transport_of \
					+ misc_yield * conversion2cm_of

	score_in = np.float32(np.around(score_in, decimals = 10))
	score_of = np.float32(np.around(score_of, decimals = 10))

	output_score_in = os.path.join(output_folder, f'crop{crop_number}_in_{get_initials(indicator)}.txt')
	output_score_of = os.path.join(output_folder, f'crop{crop_number}_of_{get_initials(indicator)}.txt')

	# print(cultivation_in[0, 111])
	# print(misc_yield[0, 111])
	# print(gf.min_dist(lc, 7)[0,111])
	# print(transport_in)
	# print(conversion2cm_in)

	# output_test1 = os.path.join(output_folder, f'test_cult_crop{crop_number}_in_{get_initials(indicator)}.txt')
	# output_test2 = os.path.join(output_folder, f'test_dist_crop{crop_number}_in_{get_initials(indicator)}.txt')
	# output_test3 = os.path.join(output_folder, f'test_conv_crop{crop_number}_in_{get_initials(indicator)}.txt')
	
	# print(transport_in)
	# save_txt(lc_fp, output_test1, cultivation_in)
	# save_txt(lc_fp, output_test2, np.float32(gf.min_dist(lc, 7)))
	# save_txt(lc_fp, output_test3, conversion2cm_in)



	save_txt(lc_fp, output_score_in, score_in)
	save_txt(lc_fp, output_score_of, score_of)

def get_splca_energy(crop_number):

	misc_yield_rs = gdal.Open(os.path.join(txt_folder, 'txt_misc_yield_wet_tpha.txt'))
	misc_yield = misc_yield_rs.GetRasterBand(1).ReadAsArray()

	misc_yield[misc_yield == -9999] = 0

	misc_yield = misc_yield * 1000

	if crop_number == 3:
		# producing biogas
		biogas_m3pkg = 564/2024
		electricity_kwhpkg = 0
		heat_mjpkg = 0

	if crop_number == 4:
		# producing electricity & heat
		biogas_m3pkg = 0
		electricity_kwhpkg = 955/2024
		heat_mjpkg = 9769/2024

	biogas_m3pha = misc_yield * biogas_m3pkg
	electricity_kwhpha = misc_yield * electricity_kwhpkg
	heat_mjpha = misc_yield * heat_mjpkg

	biogas_m3pha = np.around(biogas_m3pha, decimals = 1	)
	electricity_kwhpha = np.around(electricity_kwhpha, decimals = 1)
	heat_mjpha = np.around(heat_mjpha, decimals = 1)

	output_biogas_m3pha = os.path.join(output_folder, f'crop{crop_number}_biogas_m3pha.txt')
	output_electricity_kwhpha = os.path.join(output_folder, f'crop{crop_number}_electricity_kwhpha.txt')
	output_heat_mjpha = os.path.join(output_folder, f'crop{crop_number}_heat_mjpha.txt')

	

	save_txt(lc_fp, output_biogas_m3pha, biogas_m3pha)
	save_txt(lc_fp, output_electricity_kwhpha, electricity_kwhpha)
	save_txt(lc_fp, output_heat_mjpha, heat_mjpha)



###################################################
###################################################
# indicator = 'Climate change, ecosystem quality, long term'
# get_splca_pre_inoff(indicator, 1)


# get_transport_inoff('Climate change, ecosystem quality, long term')
# get_conversion_inoff('Climate change, ecosystem quality, long term', 'combustion')
# output_test2 = os.path.join(output_folder, f'test_dist.txt')
# save_txt(lc_fp, output_test2, np.float32(gf.min_dist(lc, 7)))

# test = '/Users/tianran/Desktop/janustd/data_td/STLCA/full/crop4_in_cceqlt.txt'
# a_a = gdal.Open(test)
# a = a_a.GetRasterBand(1).ReadAsArray()
# np.where(a != 0)


# get_splca_misc_inoff('Climate change, ecosystem quality, long term', 4)

# save all misc spatial impacts
for j in range(3, 5):
	get_splca_energy(j)

for n in dcf.lcia_name.unique():
	for i in range(1, 3):
		get_splca_pre_inoff(n, i)

	for j in range(3, 5):
		get_splca_misc_inoff(n, j)

# print(get_splca_misc_cult_inoff('Climate change, ecosystem quality, long term'))[0,111]
# print(get_transport_inoff('Freshwater ecotoxicity, long term'))


