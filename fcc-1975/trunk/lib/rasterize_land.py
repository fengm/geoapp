'''
File: rasterize_shp.py
Author: Min Feng
Version: 0.1
Create: 2013-06-14 11:29:56
Description: rasterize shapefile into a image
'''

import logging

def percent_land(f_img):
	import geo_raster_c as ge
	_bnd = ge.geo_raster.open(f_img).get_band().cache()
	_sum = _bnd.data.sum()

	return float(_sum) / (_bnd.width * _bnd.height)

def land(obj, img, fzip):
	import config
	# print 'rasterize continental boundary'

	_tile = obj.tile

	if config.cfg.has_option('path', 'land_mask'):
		_f_mak = config.cfg.get('path', 'land_mask')

		if _f_mak.endswith('.shp'):
			logging.info('load land mask shapefile %s' % _f_mak)

			import geo_raster_ex_c as gx
			return gx.geo_band_stack_zip.from_shapefile(_f_mak, file_unzip=fzip).read_block(img)

		_f_land = config.cfg.get('path', 'land_mask') % _tile

		logging.info('use land mask %s (%s)' % (_tile, _f_land))
		import os
		if os.path.exists(_f_land):
			logging.info('load land mask %s' % _f_land)
			import geo_raster_c as ge
			return ge.open(_f_land).get_band().read_block(img)

	logging.info('rasterize continental boundary')

	_f_shp1 = config.get_at('path', 'adm_1')
	_f_out1 = fzip.generate_file('', '.tif')

	rasterize_shp_bnd(_f_shp1, img, _f_out1, fzip)
	_per = percent_land(_f_out1)
	logging.info('percent of land: %s' % _per)

	if _per >= 0.99:
		return None

	_f_shp2 = config.get_at('path', 'adm_2')
	_f_out2 = fzip.generate_file('', '.tif')
	rasterize_shp_bnd(_f_shp2, img, _f_out2, fzip)

	import geo_raster_c as ge
	return ge.geo_raster.open(_f_out2).get_band().cache()

def rasterize_shp_bnd(f_shp, img, f_out, fzip):
	logging.info('rasterize img')

	import geo_raster_ex_c as gx
	import geo_raster_c as ge

	_ext = gx.geo_polygon.from_raster(img)
	_ext = _ext.project_to(ge.proj_from_epsg())
	_ext = _ext.extent()
	_prj = img.proj.ExportToProj4().strip()

	_f_out = fzip.generate_file('', '.shp')

	# reproject the shapefile to the same projection of the target image
	_cmd = 'ogr2ogr -skipfailures -t_srs "%s" -clipsrc %f %f %f %f %s %s' % (_prj, _ext.minx, _ext.miny, _ext.maxx, _ext.maxy, _f_out, f_shp)
	logging.info('reprojecting to ' + _prj)

	import run_commands
	run_commands.run_exe(_cmd)

	rasterize_shp_bnd_same_prj(_f_out, img, f_out)

	logging.info('rasterized')

def rasterize_shp_img(f_shp, f_img, f_out):
	logging.info('rasterize shp ' + f_shp)

	import geo_raster_c as ge
	_img = ge.geo_raster.open(f_img)

	rasterize_shp_bnd(f_shp, _img, f_out)

def rasterize_shp_bnd_same_prj(f_shp, img, f_out):
	logging.info('rasterizing ' + f_shp + ' to ' + f_out)
	from osgeo import ogr, gdal

	_shp = ogr.Open(f_shp)
	_lyr = _shp.GetLayer()

	_img = gdal.GetDriverByName('GTiff').Create(f_out, img.width,\
			img.height, 1, gdal.GDT_Byte)
	_img.SetGeoTransform(img.geo_transform)
	_img.SetProjection(img.proj.ExportToWkt())

	gdal.RasterizeLayer(_img, [1], _lyr, burn_values=[1])

	del _img
	del _lyr
	del _shp

def usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('-i', '--input', dest='input', required=True)
	_p.add_argument('-o', '--output', dest='output', required=True)
	_p.add_argument('-r', '--refer', dest='refer', required=True)

	return _p.parse_args()

def main():
	_opts = usage()

	rasterize_shp_bnd(_opts.input, _opts.refer, _opts.output)

def init_env():
	import os, sys
	_d_in = os.path.join(sys.path[0], 'lib')
	if os.path.exists(_d_in):
		sys.path.append(_d_in)

	import logging_util
	logging_util.init()

if __name__ == '__main__':
	init_env()
	main()

