'''
File: detect_shreshold.py
Author: Min Feng
Version: 0.1
Create: 2015-07-07 13:00:36
Description:
'''

import logging

class series:

	def __init__(self, ss, v_min, v_max, v_div):
		self._ss = ss
		self._v_min = v_min
		self._v_max = v_max
		self._v_div = v_div

		_vn = []
		_vv = []

		import math
		for _v in xrange(v_min, v_max + 1, v_div):
			_vn.append(_v)
			_vv.append(math.log(self._agg(ss, _v, _v + v_div) + 1))

		self.tags = _vn
		self.vals = _vv
		self.num = len(_vn)

		self._slope()

	def _agg(self, ss, v_min, v_max):
		if v_min == v_max:
			return ss[v_min]

		_t = 0
		for _v in xrange(v_min, v_max):
			_t += ss[_v]
		return _t

	def _slope_val(self, idx):
		assert(not (idx < 0 or idx >= self.num))

		if idx == 0:
			return 0

		import math
		return math.atan((self.vals[idx] - self.vals[idx-1]) / 1.0)

	def _slope(self):
		_vv = []
		for i in xrange(self.num):
			_vv.append(self._slope_val(i))

		self.slps = _vv

		_vd = []
		for i in xrange(self.num - 1):
			_vd.append(_vv[i+1] - _vv[i])
		_vd.append(0)

		self.difs = _vd

	def _search_peak_end(self, idx):
		for i in xrange(idx, self.num-1):
			if self.difs[i] > 0 and self.difs[i] > self.difs[i+1]:
				return i

		raise Exception('failed to find the end of the peak')

	def _peak_test1(self, ps):
		_ps = filter(lambda x: x < 450, [self.tags[self.difs.index(_p)] for _p in ps])

		if len(_ps) == 0:
			return -1

		return self.tags.index(max(_ps))

	def _peak_test2(self, ps):
		_peak_idxs = [self.difs.index(_p) for _p in ps]
		_peak_bin = min(_peak_idxs)

		return _peak_bin

	def search_peak(self):
		# identify the peak
		_peaks = sorted([_v for _v in self.difs if _v < 0])
		logging.info('found %s initial peaks' % len(_peaks))
		if len(_peaks) == 0:
			return -1000, -1000

		_peaks = _peaks[: min(3, len(_peaks))]
		logging.info('identified peaks %s' % len(_peaks))

		_peak_bin = self._peak_test1(_peaks)
		logging.info('shreshold test1: %s' % _peak_bin)
		if _peak_bin < 0:
			_peak_bin = self._peak_test2(_peaks)
			logging.info('shreshold test2: %s' % _peak_bin)

		_peak_end = self._search_peak_end(_peak_bin)

		_peak_bin += 1
		_peak_end += 1

		logging.info('found thresholds: %s (%s), %s (%s)' % (self.tags[_peak_bin], _peak_bin,
			self.tags[_peak_end], _peak_end))
		return self.tags[_peak_bin], self.tags[_peak_end]

	def save(self, f):
		_ls = ['name,value,slope,diff']

		for i in xrange(len(self.tags)):
			_ls.append('%s,%s,%s,%s' % (self.tags[i], self.vals[i], self.slps[i], self.difs[i]))

		with open(f, 'w') as _fo:
			_fo.write('\n'.join(_ls))

