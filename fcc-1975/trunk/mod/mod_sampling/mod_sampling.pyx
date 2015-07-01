import numpy as np
import logging

cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

class sample:

	def __init__(self, val, row, col, dif):
		self.val = val
		self.row = row
		self.col = col
		self.dif = dif

def stat_band(bnd):
	cdef np.ndarray[np.uint8_t, ndim=2] _dat = bnd.data
	cdef _row, _col, _v, _nodata = bnd.nodata

	_ns = {}
	for _row in xrange(bnd.height):
		for _col in xrange(bnd.width):
			_v = _dat[_row, _col]
			if _v == _nodata:
				continue

			_ns[_v] = _ns.get(_v, 0.0) + 1

	return _ns

def sampling(bnd_qa, bnd_fcc, bnd_fcp, bnd_dif):
	import config
	import random

	cdef int _row, _col

	cdef int _nd_fcc = bnd_fcc.nodata, _v_fcc, _v_fcp, _v_qa
	cdef float _nd_dif = bnd_dif.nodata, _v_dif

	cdef np.ndarray[np.uint8_t, ndim=2] _dat_qa = bnd_qa.data
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_fcp = bnd_fcp.data
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_fcc = bnd_fcc.data
	cdef np.ndarray[np.float32_t, ndim=2] _dat_dif = bnd_dif.data

	_ss_f = []
	_ss_n = []

	logging.info('stat FCC band')
	_stat = stat_band(bnd_fcc)
	logging.info(str(_stat))

	_p_area = sum(_stat.values())
	logging.info('total pixels: %s' % _p_area)

	_p_change = 1 - ((_stat.get(91, 0.0) + _stat.get(19, 0.0)) / _p_area)
	logging.info('change rate: %s' % _p_change)

	# _inclusion_prop = config.cfg.getfloat('train', 'max_inclusion_prop')
	_max_samples = config.cfg.getfloat('train', 'max_samples')

	_valid_nc = [int(_v) for _v in config.cfg.get('train', 'valid_nonforest_types').replace(' ', '').split(',')]
	_valid_fc = [int(_v) for _v in config.cfg.get('train', 'valid_forest_types').replace(' ', '').split(',')]
	_valid_lc = _valid_nc + _valid_fc

	_inclusion_prop = _max_samples / sum([_stat.get(_lc, 1.0) for _lc in _valid_lc])

	_forest = config.cfg.getboolean('config', 'percent_forest') == False
	if not _forest and bnd_fcp == None:
		raise Exception('percent of forest layer is required for predicting of percent forest')

	logging.info('inclusion_prop: %s' % _inclusion_prop)

	import progress_percentage
	_ppp = progress_percentage.progress_percentage(bnd_qa.height, ' - collect samples')

	for _row in xrange(bnd_qa.height):
		_ppp.next()

		for _col in xrange(bnd_qa.width):
			_v_qa = _dat_qa[_row, _col]
			_v_fcc = _dat_fcc[_row, _col]
			_v_fcp = _dat_fcp[_row, _col]
			_v_dif = _dat_dif[_row, _col]

			if _v_qa != 1:
				continue

			if _v_fcc not in _valid_lc:
				continue

			if _inclusion_prop < 1.0 and random.random() > _inclusion_prop:
				continue

			if _v_fcc in _valid_fc:
				_ss_f.append(sample(1 if _forest else _v_fcp, _row, _col, _v_dif))
			else:
				_ss_n.append(sample(9 if _forest else _v_fcp, _row, _col, _v_dif))

			# if not (_v_fcc == 11 or _v_fcc == 99):
			# 	continue

			# if _inclusion_prop < 1.0 and random.random() > _inclusion_prop:
			# 	continue

			# if _v_fcc == 11:
			# 	_ss_f.append(sample(_row, _col, _v_dif))
			# else:
			# 	_ss_n.append(sample(_row, _col, _v_dif))

	_ppp.done()
	logging.info('init samples F: %s, N: %s' % (len(_ss_f), len(_ss_n)))

	_cmp = lambda x, y: cmp(x.dif, y.dif)
	_ss_f.sort(_cmp)
	_ss_n.sort(_cmp)

	print ' - filter samples'
	if config.cfg.getboolean('train', 'change_filter'):
		_p_max_change = config.cfg.getfloat('train', 'max_change')
		_pe_f = max(_p_change, _p_max_change)
		_pe_n = max(_p_change, _p_max_change)
		logging.info('change filter threshold: %s, %s' % (_pe_f, _pe_n))

		if len(_ss_f) > 10:
			_ss_f = _ss_f[:max(10, int(len(_ss_f) * _pe_f))]
		if len(_ss_n) > 10:
			_ss_n = _ss_n[:max(10, int(len(_ss_n) * _pe_n))]

	if config.cfg.getboolean('train', 'consistency_filter'):
		_p_forest_f = _stat.get(11) / max(_stat.get(99), 1.0)
		_p_forest_p = len(_ss_f) / float(max(len(_ss_n), 1))
		logging.info('forest/nonforest ratio: %s (FCC), %s (sampling)' % (_p_forest_f, _p_forest_p))

		if (0.01 < _p_forest_f < 0.99) and (0.01 < _p_forest_p < 0.99):
			logging.info('ratio filtering')

			if _p_forest_f < _p_forest_p:
				logging.info('ratio filter applied to forest samples')
				_ss_f = random.sample(_ss_f, int(_p_forest_f * len(_ss_n)))
			else:
				logging.info('ratio filter applied to nonforest samples')
				_ss_n = random.sample(_ss_n, int((1 / _p_forest_f) * len(_ss_f)))

			_p_forest_p = len(_ss_f) / float(max(len(_ss_n), 1))
			logging.info('updated forest/nonforest ratio: %s (sampling)' % (_p_forest_p))

	logging.info('final samples F: %s, N: %s' % (len(_ss_f), len(_ss_n)))
	print ' - samples F: %s, N: %s' % (len(_ss_f), len(_ss_n))

	return _ss_f, _ss_n

def write_names_file(name, vals, bs, f_out):
	_bs = bs
	_bs.sort()

	with open(f_out, 'w') as _fo:
		_fo.write(name + '\n')
		if vals:
			_fo.write('%s: %s\n' % (name, ','.join(['%s' % _v for _v in vals])))
		else:
			_fo.write('%s: continuous\n' % name)
		_fo.write('case weight: continuous\n')
		_fo.write('\n'.join(['b%s: continuous' % _b for _b in list(_bs)]))

def write_data_file(ss, bnd_qa, bnds, bnd_weight, f_out):
	_ks = bnds.keys()
	_ks.sort()

	_ext = bnd_qa.extent()
	_bnds = {}
	for _b in _ks:
		_bnds[_b] = bnds[_b].read_ext(_ext)

	cdef np.ndarray[np.uint8_t, ndim=2] _dat
	cdef np.ndarray[np.float32_t, ndim=2] _dat_weight = None
	if bnd_weight:
		_dat_weight = bnd_weight.data

	_nu = 0.0

	with open(f_out, 'w') as _fo:
		for _st in ss:
			for _s in ss[_st]:
				_vs = ['%s' % _s.val]

				# apply weight value here
				if bnd_weight == None:
					_vs.append('1.0')
				else:
					_v_weight = _dat_weight[_s.row, _s.col]
					_vs.append('%s' % _v_weight)

				_sk = False
				for _b in _ks:
					_dat = _bnds[_b].data
					_v = _dat[_s.row, _s.col]

					if _v == bnds[_b].nodata:
						_sk = True
						break

					_vs.append('%s' % _v)

				if _sk:
					continue

				_fo.write(','.join(_vs) + '\n')
				_nu += 1

		logging.info('written samples %s' % _nu)

	return _nu

