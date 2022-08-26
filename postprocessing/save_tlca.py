

import gdal
import sys
import os
import numpy as np
import rasterio as rio  
sys.path.append('/Users/tianran/Desktop/dtlca')


def save_txt(meta_src_fp, dst_fp, dst_data):
	with rio.open(meta_src_fp) as src:
		ras_meta = src.profile

	ras_meta['dtype'] = "float32"

	with rio.open(dst_fp, 'w', **ras_meta) as dst:
		dst.write(dst_data, 1)


def converse_3d_to_2d(arr_3d):


	# reshaping the array from 3D matrice to 2D matrice.
	arr_reshaped = arr_3d.reshape(arr_3d.shape[0], -1)

	return arr_reshaped

	
  



