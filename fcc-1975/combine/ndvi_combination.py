
def ndvi(bnd1, bnd2):
	import numpy as np

	_dat0 = (bnd1.data < 0) | (bnd2.data < 0)

	_dat1 = bnd1.data.astype(np.float32)
	_dat2 = bnd2.data.astype(np.float32)

	_dat3 = (_dat2 - _dat1) / (_dat2 + _dat1)
	_dat3[_dat3 > 1] = 1
	_dat3[_dat3 < -1] = -1
	_dat3 *= 1000
	_dat3[_dat0] = -9999

	return bnd1.from_grid(_dat3.astype(np.int16), nodata=-9999)

def combine(bs1, bs2):
	import mod_grid

	_b1_ndvi = bs1[0]
	_b1_band = bs1[1]

	_b2_ndvi = bs2[0].read_block(_b1_ndvi)
	_dat_ndvi = (_b1_ndvi.data > -9999) & (_b2_ndvi.data > -9999) & (_b2_ndvi.data > _b1_ndvi.data)
	for _b in xrange(len(_b1_band)):
		_b2_band = bs2[1][_b].read_block(_b1_ndvi)
		# print _dat_ndvi.shape
		# print _b1_band[_b].data.shape, _b2_band.data.shape
		# _b1_band[_b].data[_dat_ndvi] = _b2_band.data
		mod_grid.update(_b1_band[_b].data, _dat_ndvi, _b2_band.data)

	# _b1_ndvi.data[_dat_ndvi] = _b2_ndvi.data
	mod_grid.update(_b1_ndvi.data, _dat_ndvi, _b2_ndvi.data)
	return _b1_ndvi, _b1_band

def load_layers(f, d_out, fzip):
	import geo_raster_c as ge

	_img = ge.open(fzip.unzip(f))

	_bnds = []
	for _b in xrange(4):
		_bnds.append(_img.get_band(_b+1).cache())

	_bnd_ndvi = ndvi(_bnds[1], _bnds[2])

	import os
	_bnd_ndvi.save(os.path.join(d_out, os.path.basename(f)[:-7] + '_ndvi.tif'))

	return _bnd_ndvi, _bnds

def create_raster(f_out, bnd, band_num, nodata, pixel_type):
	import geo_raster_c as ge

	ge.geo_raster.create(f_out, [band_num, bnd.height, bnd.width], bnd.geo_transform,
			bnd.proj, driver='HFA', nodata=nodata, pixel_type=pixel_type)

	return ge.open(f_out, True)

def main():
	_opts = _init_env()

	import file_unzip

# 	_fs = '''/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM20260391976193GMD03.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM20260391975180AAA04.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391976202GMD05.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM20260391975198AAA04.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391974194AAA05.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391975189GMD04.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391974176GDS03.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391973199GDS03.img.gz
# /data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391973181AAA01.img.gz'''.splitlines()
	_fs = '''/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391974194AAA05.img.gz
/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/continents/LM10260391974176GDS03.img.gz'''.splitlines()

	_d_out = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi'
	_b_ndvi = None
	_b_band = None
	with file_unzip.file_unzip('/export2/data/mfeng/tmp/mss_ndvi') as _zip:
		import progress_percentage
		_ppp = progress_percentage.progress_percentage(len(_fs))

		for _f in _fs:
			_ppp.next()

			_bs = load_layers(_f, _d_out, _zip)
			if _b_ndvi == None:
				_b_ndvi = _bs[0]
				_b_band = _bs[1]
			else:
				combine((_b_ndvi, _b_band), _bs)

		_ppp.done()

		import os
		_b_ndvi.save(os.path.join(_d_out, 'combine_ndvi.tif'))

		import geo_raster_c as ge
		_img = create_raster(os.path.join(_d_out, 'combine_bands.img'), _b_ndvi, len(_b_band), -9999, ge.pixel_type('short'))

		for _b in xrange(len(_b_band)):
			_img.get_band(_b+1).write(_b_band[_b].read_block(_b_ndvi).data, 0, 0)

def _usage():
	import argparse

	_p = argparse.ArgumentParser()
	_p.add_argument('--logging', dest='logging')
	_p.add_argument('--config', dest='config')
	_p.add_argument('--temp', dest='temp')

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

