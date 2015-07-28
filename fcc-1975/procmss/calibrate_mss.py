
def cal_index(bnd1, bnd2):
	import numpy as np

	_dat1 = bnd1.data_ma.astype(np.float32)
	_dat2 = bnd2.data_ma.astype(np.float32)

	_dat = (_dat1 - _dat2) / (_dat1 + _dat2)
	_dat[_dat < -1.0] = -1.0
	_dat[_dat >  1.0] =  1.0

	return bnd1.from_grid((_dat * 1000.0).astype(np.int16).filled(-9999), nodata=-9999)

def cal_metrics(f):
	import geo_raster_c as ge

	_img = ge.open(f)
	_bnds = [_img.get_band(i+1).cache() for i in xrange(_img.band_num)]

	_out = {'ndvi': cal_index(_bnds[2], _bnds[1])}

	for _k, _v in _out.items():
		_f_out = (f[:-4] + '_%s.tif') % _k
		print 'output', _f_out
		_v.save(_f_out, opts=['compress=lzw'])

def distance(d):
	import config
	_ps = {}
	for _l in open(config.cfg.get('conf', 'sun_distance')).read().splitlines():
		_vs = _l.split('\t')
		if len(_vs) != 2:
			continue

		_ps[int(_vs[0])] = float(_vs[1])

	return _ps[int(d.strftime('%j'))]

	# import math

	# _d = int(d.strftime('%j'))
	# _om=(.9856*float(_d-4)) * math.pi / 180.0
	# _dsol=1.0 / ((1.0 - 0.01673 * math.cos(_om)) ** 2.0)

	# print _d, _dsol
	# return _dsol

def calibrate(f, f_out, fzip):
	import mod_mss

	_fs, _fm, _fo = mod_mss.load_file(f, fzip)
	if _fm == None:
		raise Exception('failed to find the metadata')

	_bs = [4, 5, 6, 7]
	_ms = mod_mss.parse_metadata(_fm)

	import geo_raster_c as ge
	import numpy as np
	import datetime
	import math

	_scale = True
	_date = datetime.datetime.strptime(_ms['DATE_ACQUIRED'], '%Y-%m-%d')

	_tmp = ge.open(_fs[0]).get_band()
	# _esun = [1805.42,1511.42,1241.4,908.018]
	_esun = [1823.0,1559.0,1276.0,880.1] # chander et al. 2009
	_zenith = math.radians(90.0 - _ms['SUN_ELEVATION'])
	_dis = distance(_date)

	_p_mult = [0.909, 0.667, 0.575, 0.480]
	_p_offs = [6.29134, -2.9748, -1.26654, 3.42008]

	_out = np.empty([len(_fs), _tmp.height, _tmp.width], dtype=np.int16)
	for _i in xrange(len(_fs)):
		_b = _bs[_i]

		_gain = _ms['RADIANCE_MULT_BAND_%s' % _b]
		_offs = _ms['RADIANCE_ADD_BAND_%s' % _b]

		_bnd = ge.open(_fs[_i]).get_band().read_block(_tmp)
		_dat = _bnd.data

		if _scale:
			_ddd = np.ma.array(_dat, mask=(_dat == 0)).astype(np.float32)

			_fac = (math.pi * (_dis ** 2)) / (_esun[_i] * math.cos(_zenith))
			print _b, _fac

			# print 'min:', (_ddd * _gain).min()
			# print 'max:', (_ddd * _gain).max()

			# _ddd = _ddd * _gain + _p_offs[_i]
			_ddd = _ddd * _p_mult[_i] + _p_offs[_i]
			# if _offs > _ddd.min():
			# 	print '****'
			# 	_ddd += _offs

			_ddd *= _fac
			_out[_i, :, :] = (_ddd * 10000.0).filled(-9999).astype(np.int16)

	mod_mss.create_output_raster(_out, _tmp, f_out, -9999)

def main():
	_opts = _init_env()
	del _opts

	import file_unzip
	with file_unzip.file_unzip('/export2/data/mfeng/tmp/mss_cal') as _zip:
		# _f_in = '/data/glcf-js-10-1/data/kcollin7-1/MSS/continents/LM20100241975236PAC00.tar.gz'
		_f_in = '/data/glcf-js-10-1/data/kcollin7-1/MSS/continents/LM10090241976239PAC08.tar.gz'
		_d_ot = '/data/glcf-nx-002/data/PALSAR/fcc_1975/sr/test/ndvi_8'

		import os
		os.path.exists(_d_ot) or os.makedirs(_d_ot)
		_f_out = os.path.join(_d_ot, os.path.basename(_f_in)[:-6] + 'img')

		calibrate(_f_in, _f_out, _zip)
		cal_metrics(_f_out)

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

