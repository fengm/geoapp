'''
File: mod_mss_fcc.py
Author: Min Feng
Version: 0.1
Create: 2015-07-06 13:23:10
Description:
'''

import logging

def metadata(p):
	import mongo_con

	with mongo_con.connect() as _con:
		_r = _con.tasks.mss_meta.find_one({'_code': str(p)})
		assert(_r != None)
		return {'SUN_AZIMUTH': _r['SUN_AZIMUTH'], 'SUN_ELEVATION': _r['SUN_ELEVATION']}

class landsat_obj:

	def __init__(self, f_inp, d_out, fzip):
		import landsat
		import os

		_inf = landsat.parse(os.path.basename(f_inp))
		if _inf == None:
			import re
			_m = re.search('p\d{3}r\d{3}', os.path.basename(f_inp))
			_inf = landsat.landsat_info('M', _m.group(), '19750601')

		assert(_inf != None)

		self._info = _inf
		self._finp = f_inp
		self._dout = d_out
		self._fout = os.path.join(d_out, '%s_%%s.%%s' % (str(_inf),))
		self._fzip = fzip

		self._open_bands()

	def band(self, tag):
		if tag not in self._bnds:
			raise Exception('faild to find band %s' % tag)

		return self._bnds[tag]

	def file_name(self, t, ext='tif'):
		return self._fout % (t, ext)

	def _open_bands(self):
		import geo_raster_c as ge
		import numpy as np

		if self._info.sensor == None:
			self._info.sensor = 'M'

		_bnd_fix = 1 if self._info.sensor == 'M' else 0
		logging.info('band fix number: %s' % _bnd_fix)

		_img = ge.open(self._fzip.unzip(self._finp))

		# QA band
		import params
		_dat_qa = np.empty((_img.height, _img.width), dtype=np.uint8)
		_dat_qa.fill(params.MSS_LAND)

		_bnd = {}
		for _b in xrange(1, _img.band_num + 1):
			_bbb = _img.get_band(_b).cache()
			_bnd['b%s' % (_b + _bnd_fix)] = _bbb
			_dat_qa[(_dat_qa == 1) & (_bbb.data == _bbb.nodata)] = 0

		self._temp = _img
		self._bnds = _bnd
		self._bnds['qa'] = _bnd.values()[0].from_grid(_dat_qa, nodata=0)

	def _cal_index(self, bnd1, bnd2):
		import numpy as np
		import config

		_dat1 = bnd1.data_ma.astype(np.float32)
		_dat2 = bnd2.data_ma.astype(np.float32)

		_dat = (_dat1 - _dat2) / (_dat2 + _dat1)
		_dat[(_dat < -1.0) & (_dat.mask == False)] = -1
		_dat[(_dat > 1.0) & (_dat.mask == False)] = 1
		_dat = _dat * config.cfg.getint('param', 'scale_value')

		_bnd = bnd1.from_ma_grid(_dat.astype(np.int16), nodata=-9999)

		return _bnd

	def _load_color_table(self, tag):
		import os
		import config

		_f_clr = os.path.join(config.cfg.get('conf', 'colors'), 'colors_%s.txt' % tag)
		logging.info('checking color table %s: %s' % (tag, _f_clr))
		if not os.path.exists(_f_clr):
			return None

		logging.info('loading color table %s: %s' % (tag, _f_clr))
		import geo_raster_c as ge
		return ge.load_colortable(_f_clr)

	def _mask_water(self):
		import params
		self._bnds['qa'].data[\
				(self._bnds['qa'].data == params.MSS_LAND) & \
				((self._bnds['mndwi'].data > 0) | \
				(self._bnds['ndwi'].data > 0)) \
				] = params.MSS_WATER

		import rasterize_land
		_bnd = rasterize_land.land(self._info, self._bnds['qa'], self._fzip)
		if _bnd == None: return
		self._bnds['qa'].data[(self._bnds['qa'].data != params.MSS_NODATA) & (_bnd.data != 1)] = params.MSS_WATER

	def _detect_thresholds(self):
		import mod_stat_band
		import detect_shreshold
		import config
		import params

		_ss_b3 = detect_shreshold.series(mod_stat_band.stat( \
			self._bnds['b3'], self._bnds['qa'], params.MSS_LAND), \
				0, 1000, 50)
		if config.cfg.getboolean('conf', 'debug'):
			_ss_b3.save(self.file_name('b3', 'csv'))

		_ps_b3 = _ss_b3.search_peak()

		# update the QA band with the thresholds
		_dat_qa = self._bnds['qa'].data
		_dat_b3 = self._bnds['b3'].data
		_dat_ndvi = self._bnds['ndvi'].data
		_ndvi_div = config.cfg.getfloat('param', 'ndvi_threshold') * config.cfg.getfloat('param', 'scale_value')

		_dat_qa[(_dat_qa == 1) & (_dat_b3 <= _ps_b3[0])] = params.MSS_FOREST1
		_dat_qa[(_dat_qa == 1) & (_dat_b3 <= _ps_b3[1])] = params.MSS_FOREST2
		_dat_qa[(_dat_qa == 1) & (_dat_ndvi >= _ndvi_div)] = params.MSS_FOREST3

		if config.cfg.getboolean('param', 'enable_ndvi'):
			_ss_ndvi = detect_shreshold.series(mod_stat_band.stat( \
				self._bnds['ndvi'], self._bnds['qa'], 1), \
					0, 1000, 50)
			if config.cfg.getboolean('conf', 'debug'):
				_ss_ndvi.save(self.file_name('ndvi', 'csv'))

	def _output_layer(self, tag, bnd):
		_bnd = bnd
		_tag = tag

		# if tag.startswith('i_'):
		# 	logging.info('convert fload raster to int16')

		# 	_tag = tag[2:]
			# import numpy as np
			# _dat = (_bnd.data_ma * 1000.0).astype(np.int16)
			# _bnd = bnd.from_grid(_bnd.data, nodata=-9999)

		logging.info('write %s' % _tag)
		_bnd.save(self.file_name(_tag, 'tif'), \
				color_table=self._load_color_table(_tag), \
				opts=['compress=lzw'])

	def output_results(self):
		import config

		if self._fout != None:
			logging.info('write results')

			for _k, _b in self._bnds.items():
				if config.cfg.getboolean('conf', 'debug'):
					self._output_layer(_k, _b)

	def _load_refer(self, f, bnd):
		if not f: return None

		import geo_raster_ex_c as gx

		logging.info('loading data: %s' % f)
		_bnd = gx.geo_band_stack_zip.from_shapefile(f, file_unzip=self._fzip)

		return _bnd.read_block(bnd)

	def _load_lyr_fcc(self):
		import config

		_bnd = self._load_refer(config.cfg.get('path', 'fcc'), self._bnds['qa'])
		_dat = _bnd.data

		_dat[(_dat == 11) | (_dat == 19)] = 1
		_dat[(_dat == 99) | (_dat == 91)] = 9

		return _bnd

	def _load_refers(self):
		import config
		# import agg_band

		self._bnds['dem'] = self._load_refer(config.cfg.get('path', 'dem'), self._bnds['qa'])
		self._bnds['fcc'] = self._load_lyr_fcc()

	def classify(self):
		logging.info('cal indexes')
		self._bnds['ndwi'] = self._cal_index(self._bnds['b2'], self._bnds['b4'])
		self._bnds['mndwi'] = self._cal_index(self._bnds['b2'], self._bnds['b5'])
		self._bnds['ndvi'] = self._cal_index(self._bnds['b4'], self._bnds['b3'])

		self._mask_water()
		self._detect_thresholds()
		self._load_refers()

		# logging.info('classify water')
		# import combine_grids
		# self._bnds['mak'] = combine_grids.classify_water(self._bnds)

def mss_fcc(f_inp, d_out, fzip):
	logging.info('process MSS scene: %s' % f_inp)

	import os
	(lambda x: os.path.exists(x) or os.makedirs(x))(d_out)

	_obj = landsat_obj(f_inp, d_out, fzip)
	_obj.classify()

	import mod_classify
	mod_classify.train(_obj)
	mod_classify.build(_obj)

	_obj.output_results()


