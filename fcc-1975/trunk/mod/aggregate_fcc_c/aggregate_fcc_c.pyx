'''
File: aggregate_fcc_c.pyx
Author: Min Feng
Version: 0.1
Create: 2013-12-02 13:12:19
Description:
'''

import numpy as np
import geo_raster_c as ge
import math
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def average(bnd_in, bnd_ot):
	# only support float32 data type
	assert(bnd_in.data.dtype == np.uint8)

	_geo_in = list(bnd_in.geo_transform)
	_geo_ot = list(bnd_ot.geo_transform)

	cdef float _cell_in = _geo_in[1]
	cdef float _cell_ot = _geo_ot[1]

	cdef float _dive = _cell_ot / _cell_in
	_size = [bnd_ot.height, bnd_ot.width]
	_offs = [(_geo_ot[3] - _geo_in[3]) / _geo_in[5],
				(_geo_ot[0] - _geo_in[0]) / _geo_in[1]]

	_nodata = bnd_in.nodata if bnd_in.nodata != None else -9999
	cdef np.ndarray[np.uint8_t, ndim=2] _dat = average_fcc_pixels(bnd_in.data,
			_offs[0], _offs[1], _dive,
			_nodata, _size[0], _size[1])

	return ge.geo_band_cache(_dat, _geo_ot, bnd_ot.proj,
				_nodata, ge.pixel_type())

cdef np.ndarray[np.uint8_t, ndim=2] average_fcc_pixels(np.ndarray[np.uint8_t, ndim=2] dat,
		float off_y, float off_x, float scale,
		int nodata, unsigned int rows, unsigned int cols):

	cdef unsigned int _rows_o, _cols_o
	cdef unsigned int _rows_n, _cols_n

	_rows_o = dat.shape[0]
	_cols_o = dat.shape[1]

	_rows_n = rows
	_cols_n = cols

	cdef unsigned int _row_o, _col_o
	cdef unsigned int _row_n, _col_n

	cdef int _row_min, _row_max
	cdef int _col_min, _col_max

	cdef float _row_min_f, _row_max_f
	cdef float _col_min_f, _col_max_f

	cdef double _ns
	cdef float _a
	cdef int _tp

	cdef int _nodata
	_nodata = nodata

	_dat = np.empty([_rows_n, _cols_n], np.uint8)
	_dat.fill(_nodata)

	_valid_vals = (11, 99, 19, 91)
	_chang_vals = (19, 91)

	_row_min_f = off_y - scale
	for _row_n from 0<=_row_n<_rows_n:
		_row_min_f = _row_min_f + scale
		_row_max_f = _row_min_f + scale

		_col_min_f = off_x - scale
		for _col_n from 0<=_col_n<_cols_n:
			_col_min_f = _col_min_f + scale
			_col_max_f = _col_min_f + scale

			if _row_max_f <= 0 or _col_max_f <= 0 or \
					_row_min_f >= _rows_o or _col_min_f >= _cols_o:
				continue

			_row_min = int(math.floor(_row_min_f))
			_row_min = max(0, _row_min)

			_col_min = int(math.floor(_col_min_f))
			_col_min = max(0, _col_min)

			_row_max = int(math.ceil(_row_max_f))
			_row_max = min(_rows_o, _row_max)

			_col_max = int(math.ceil(_col_max_f))
			_col_max = min(_cols_o, _col_max)

			_vs = {}
			_ns = 0
			_tp = _nodata

			for _row_o from _row_min<=_row_o<_row_max:
				for _col_o from _col_min<=_col_o<_col_max:
					_a = (min(_row_o + 1, _row_max_f) - \
							max(_row_o, _row_min_f)) * \
							(min(_col_o + 1, _col_max_f) - \
							(max(_col_o, _col_min_f)))

					if _a < 0.5:
						continue

					_v = dat[_row_o, _col_o]
					if _v == _nodata:
						continue

					# if _v not in _valid_vals:
					# 	continue

					_ns += _a

					if _v in _chang_vals:
						_tp = _v
						break

					_vs[_v] = _vs.get(_v, 0) + 1

				if _tp in _chang_vals:
					break

			if _tp == _nodata:
				if _ns <= 2:
					continue

				_mx = 0
				for _kk in _vs:
					if _vs[_kk] > _mx:
						_mx = _vs[_kk]
						_tp = _kk

			_dat[_row_n, _col_n] = _tp

	return _dat

def average_percent(bnd_in, bnd_ot, nodata=255):
	# only support float32 data type
	assert(bnd_in.data.dtype == np.uint8)

	_geo_in = list(bnd_in.geo_transform)
	_geo_ot = list(bnd_ot.geo_transform)

	cdef float _cell_in = _geo_in[1]
	cdef float _cell_ot = _geo_ot[1]

	cdef float _dive = _cell_ot / _cell_in
	_size = [bnd_ot.height, bnd_ot.width]
	_offs = [(_geo_ot[3] - _geo_in[3]) / _geo_in[5],
				(_geo_ot[0] - _geo_in[0]) / _geo_in[1]]

	cdef np.ndarray[np.uint8_t, ndim=2] _dat = average_fcc_percent_pixels(bnd_in.data,
			_offs[0], _offs[1], _dive,
			nodata, _size[0], _size[1])

	return ge.geo_band_cache(_dat, _geo_ot, bnd_ot.proj,
				nodata, ge.pixel_type())

cdef np.ndarray[np.uint8_t, ndim=2] average_fcc_percent_pixels(np.ndarray[np.uint8_t, ndim=2] dat,
		float off_y, float off_x, float scale,
		int nodata, unsigned int rows, unsigned int cols):

	cdef unsigned int _rows_o, _cols_o
	cdef unsigned int _rows_n, _cols_n

	_rows_o = dat.shape[0]
	_cols_o = dat.shape[1]

	_rows_n = rows
	_cols_n = cols

	cdef unsigned int _row_o, _col_o
	cdef unsigned int _row_n, _col_n

	cdef int _row_min, _row_max
	cdef int _col_min, _col_max

	cdef float _row_min_f, _row_max_f
	cdef float _col_min_f, _col_max_f

	cdef double _ns
	cdef float _a
	cdef int _tp

	cdef int _nodata
	_nodata = nodata

	_dat = np.empty([_rows_n, _cols_n], np.uint8)
	_dat.fill(_nodata)

	_valid_vals = (11, 99, 19, 91)

	_row_min_f = off_y - scale
	for _row_n from 0<=_row_n<_rows_n:
		_row_min_f = _row_min_f + scale
		_row_max_f = _row_min_f + scale

		_col_min_f = off_x - scale
		for _col_n from 0<=_col_n<_cols_n:
			_col_min_f = _col_min_f + scale
			_col_max_f = _col_min_f + scale

			if _row_max_f <= 0 or _col_max_f <= 0 or \
					_row_min_f >= _rows_o or _col_min_f >= _cols_o:
				continue

			_row_min = int(math.floor(_row_min_f))
			_row_min = max(0, _row_min)

			_col_min = int(math.floor(_col_min_f))
			_col_min = max(0, _col_min)

			_row_max = int(math.ceil(_row_max_f))
			_row_max = min(_rows_o, _row_max)

			_col_max = int(math.ceil(_col_max_f))
			_col_max = min(_cols_o, _col_max)

			_vs = {}
			_ns = 0

			for _row_o from _row_min<=_row_o<_row_max:
				for _col_o from _col_min<=_col_o<_col_max:
					_a = (min(_row_o + 1, _row_max_f) - \
							max(_row_o, _row_min_f)) * \
							(min(_col_o + 1, _col_max_f) - \
							(max(_col_o, _col_min_f)))

					if _a < 0.5:
						continue

					_v = dat[_row_o, _col_o]
					if _v == _nodata:
						continue

					if _v not in _valid_vals:
						continue

					_ns += _a
					_vs[_v] = _vs.get(_v, 0.0) + _a

			if _ns <= 2:
				continue

			# calculate the precent forest for each aggregated pixel coverage
			_dat[_row_n, _col_n] = int(((_vs.get(11, 0.0) + _vs.get(19, 0.0)) * 100.0) / _ns)

	return _dat


def update_band(np.ndarray[np.uint8_t, ndim=2] dat, np.ndarray[np.uint8_t, ndim=2] ref, int off_x, int off_y, int nodata):
	cdef int _row, _col
	cdef int _v1, _v2

	_chang_vals = (19, 91)

	_row_s = max(off_y, 0)
	_row_e = min(off_y + ref.shape[0], dat.shape[0])

	_col_s = max(off_x, 0)
	_col_e = min(off_x + ref.shape[1], dat.shape[1])

	for _row in xrange(_row_s, _row_e):
		for _col in xrange(_col_s, _col_e):
			_v1 = dat[_row, _col]
			_v2 = ref[_row - off_y, _col - off_x]

			if _v1 == nodata or _v2 in _chang_vals:
				dat[_row, _col] = _v2

def update_band_fill(np.ndarray[np.uint8_t, ndim=2] dat, np.ndarray[np.uint8_t, ndim=2] ref, int off_x, int off_y, int nodata):
	cdef int _row, _col
	cdef int _v1, _v2

	_row_s = max(off_y, 0)
	_row_e = min(off_y + ref.shape[0], dat.shape[0])

	_col_s = max(off_x, 0)
	_col_e = min(off_x + ref.shape[1], dat.shape[1])

	for _row in xrange(_row_s, _row_e):
		for _col in xrange(_col_s, _col_e):
			_v1 = dat[_row, _col]
			_v2 = ref[_row - off_y, _col - off_x]

			if _v1 == nodata:
				dat[_row, _col] = _v2

def update_band_fill_int16(np.ndarray[np.int16_t, ndim=2] dat, np.ndarray[np.int16_t, ndim=2] ref, int off_x, int off_y, int nodata):
	cdef int _row, _col
	cdef int _v1, _v2

	_row_s = max(off_y, 0)
	_row_e = min(off_y + ref.shape[0], dat.shape[0])

	_col_s = max(off_x, 0)
	_col_e = min(off_x + ref.shape[1], dat.shape[1])

	for _row in xrange(_row_s, _row_e):
		for _col in xrange(_col_s, _col_e):
			_v1 = dat[_row, _col]
			_v2 = ref[_row - off_y, _col - off_x]

			if _v1 == nodata:
				dat[_row, _col] = _v2

