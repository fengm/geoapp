'''
File: build_tiles.py
Author: Min Feng
Version: 0.1
Create: 2015-07-21 17:32:47
Description:
'''

def tiles(level, ext=None):
	import math
	import geo_base_c as gb

	_b = 6378137.0
	_s = 256
	_p = _b * math.pi
	_r = (2 * _p) / (2 ** level)
	_c = _r / _s

	import geo_raster_c as ge
	_prj = ge.proj_from_epsg(3857)

	_rows = 2 ** level
	_cols = 2 ** level

	_num = -1
	for _row in xrange(_rows):
		for _col in xrange(_cols):
			_num += 1

			_x = -_p + (_col * _r)
			_y = -_p + (_row * _r)

			_ext = gb.geo_extent(_x, _y, _x + _r, _y + _r, _prj)
			if ext == None or _ext.is_intersect(ext):
				_geo = [_x, _c, 0, _y + _r, 0, -_c]
				yield level, _num, _col, _row, ge.geo_raster_info(_geo, _s, _s, _prj)

def load_shp(f):
	from osgeo import ogr
	import geo_base_c as gb

	_shp = ogr.Open(f)
	if _shp == None:
		raise Exception('Failed to load shapefile ' + f)

	_lyr = _shp.GetLayer()
	_objs = []
	_area = None

	for _f in _lyr:
		_obj = gb.geo_polygon(_f.geometry().Clone())
		_ext = _obj.extent()

		_objs.append(_obj)
		if _area == None:
			_area = _ext
		else:
			_area = _area.union(_ext)

	import geo_raster_c as ge
	_prj = ge.proj_from_epsg(3857)

	_reg = _area.to_polygon().segment_ratio(30).project_to(_prj)
	return _reg.extent()

def load_img(f, fzip):
	import geo_raster_c as ge

	_prj = ge.proj_from_epsg(3857)
	_reg = ge.open(fzip.unzip(f)).extent().to_polygon().segment_ratio(30).project_to(_prj)

	return _reg.extent()

class band:

	def __init__(self, f, fzip):
		if f.endswith('.shp'):
			import geo_raster_ex_c as gx
			self.bnd = [gx.geo_band_stack_zip(f, file_unzip=fzip)]
		else:
			import geo_raster_c as ge
			_img = ge.open(fzip.unzip(f))
			self.bnd = [_img.get_band(_b + 1) for _b in xrange(_img.band_num)]

		self.color = self.bnd[0].color_table if len(self.bnd) == 1 else None

	def _color(self, c):
		_cs = {}
		if c == None:
			return _cs

		for i in xrange(256):
			_v = c.GetColorEntry(i)
			if len(_v) == 3:
				_v = list(_v) + [255]

			_cs[i] = _v

		return _cs

	def _save(self, bnd, cs, f):
		import mod_image
		_dat = mod_image.convert(bnd, cs)

		import png
		png.from_array(_dat, 'RGBA').save(f)

	def read(self, bnd, f_out):
		if len(self.bnd) == 1:
			_bnd = self.bnd[0].read_block(bnd)
			self._save(_bnd, self._color(self.color), f_out)
			# _bnd.save(f_out, driver='GTIFF', color_table=self.color)
		else:
			from osgeo import gdal
			_img = gdal.GetDriverByName('PNG').Create(f_out, bnd.width,\
					bnd.height, len(self.bnd), gdal.GDT_Byte)

			for _b in xrange(len(self.bnd)):
				_img.get_band(_b+1).write(self.bnd[_b].read_block(bnd).data, 0, 0)

			_img.flush()

def main():
	_opts = _init_env()

	_level = 10

	import file_unzip
	import os

	with file_unzip.file_unzip() as _zip:
		_ext = load_shp(_opts.input) if _opts.input.endswith('.shp') else load_img(_opts.input, _zip)
		_img = band(_opts.input, _zip)

		for _lev, _num, _col, _row, _bnd in tiles(_level, _ext):
			_d = os.path.join(_opts.output, str(_lev), str(_col))
			try:
				os.path.exists(_d) or os.makedirs(_d)
			except Exception:
				pass

			_f = os.path.join(_d, '%s.png' % _row)
			_img.read(_bnd, _f)

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

	_p.add_argument('-i', '--input', dest='input', required=True)
	_p.add_argument('-o', '--output', dest='output', required=True)

	return _p.parse_args()

def _init_env():
	import os, sys

	_dirs = ['lib', 'libs']
	_d_ins = [os.path.join(sys.path[0], _d) for _d in _dirs if \
			os.path.exists(os.path.join(sys.path[0], _d))]
	sys.path = [sys.path[0]] + _d_ins + sys.path[1:]

	_opts = _usage()

	import logging_util
	logging_util.init(_opts.logging)

	import config
	config.load(_opts.config)

	import file_unzip as fz
	fz.clean(fz.default_dir(_opts.temp))

	return _opts

if __name__ == '__main__':
	main()

