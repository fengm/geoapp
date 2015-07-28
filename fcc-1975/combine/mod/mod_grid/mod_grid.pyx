
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def update( np.ndarray[np.int16_t, ndim=2] dat1, np.ndarray[np.uint8_t, cast=True, ndim=2] dati, np.ndarray[np.int16_t, ndim=2] dat2):
	cdef int _rows = dat1.shape[0], _cols = dat1.shape[1]
	cdef int _row, _col, _v1, _v2, _vi

	# cdef float _n
	# _n = 0
	for _row in xrange(_rows):
		for _col in xrange(_cols):
			_vi = dati[_row, _col]
			# if _row == 300 and _col == 300:
			# 	print _vi, dati[0, 0]
			if _vi != 1:
				continue

			_v2 = dat2[_row, _col]
			dat1[_row, _col] = _v2
			# _n += 1

	# print 'found', _n

cdef int cal_ndvi(int v1, int v2, int v3):
	if v1 > 3000:
		return -9999

	if v2 > 2500:
		return -9999

	if v3 > 1600:
		return -9999

	return cal_index(v1, v2)

cdef int cal_ndwi(int v1, int v2):
	if v1 > 2500 or v2 > 2000:
		return -9999

	return cal_index(v1, v2)

cdef int cal_index(int v1, int v2):
	if v1 < 0 or v2 < 0:
		return -9999

	if v1 + v2 == 0:
		return -9999

	cdef float _v = float(v1 - v2) / (v2 + v1)
	if _v < -1.0:
		_v = -1.0
	if _v > 1.0:
		_v = 1.0

	return int(_v * 1000)

def read_pixel(bnds, int row, int col):
	cdef np.ndarray[np.int16_t, ndim=2] _dat
	cdef int _v, _b

	_vs = []
	_nu = 0
	for _b in xrange(len(bnds)):
		assert(row < bnds[_b].height and row >= 0)
		assert(col < bnds[_b].width and col >= 0)

		_dat = bnds[_b].data
		_v = _dat[row, col]

		if _v >= 0:
			_nu += 1
			_vs.append(_v)
		else:
			_vs.append(-9999)

	return _vs, _nu

def save_pixel(bnds, int row, int col, vs):
	cdef np.ndarray[np.int16_t, ndim=2] _dat
	cdef int _v, _b

	_vs = []
	for _b in xrange(len(bnds)):
		assert(row < bnds[_b].height and row >= 0)
		assert(col < bnds[_b].width and col >= 0)

		_dat = bnds[_b].data
		_dat[row, col] = vs[_b]

def pick_values(vs, th, min_val):
	_vs = []
	for _v in vs:
		if _v > 0:
			_vs.append(_v)
	
	if len(_vs) == 0: return -9999, -1
	if len(_vs) == 1: return _vs[0], vs.index(_vs[0])

	# _v = _vs[int(th * (len(_vs) - 1))]
	_v = max(_vs)
	# if len(_vs) > 4:
	# 	_v = _vs[int(th * (len(_vs) - 1))]
	# 	# _v = sorted(vs)[len(vs) - 2]
	# else:
	# 	_v = max(_vs)

	_idx = -1
	if _v > min_val:
		_idx = vs.index(_v)

	return _v, _idx

def combine(imgs):
	assert(len(imgs) > 0)

	_bnds = []
	_bnd_out = None
	for i in xrange(len(imgs)):
		_bbb = [imgs[i].get_band(_b+1) for _b in xrange(imgs[i].band_num)]
		if _bnd_out == None:
			_bnds.append([_bnd.cache() for _bnd in _bbb])
			_bnd_out = _bnds[-1]
		else:
			_bnds.append([_bnd.read_block(_bnd_out[0]) for _bnd in _bbb])

	if len(imgs) == 1:
		return _bnd_out

	cdef int _rows = _bnd_out[0].height, _cols = _bnd_out[0].width
	cdef int _row, _col, _v1, _v2

	cdef np.ndarray[np.int16_t, ndim=2] _dat_ndvi = np.empty((_rows, _cols), np.int16)
	cdef np.ndarray[np.int16_t, ndim=2] _dat_ndwi = np.empty((_rows, _cols), np.int16)
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_indx = np.empty((_rows, _cols), np.uint8)
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_lyrs = np.empty((_rows, _cols), np.uint8)

	_dat_ndvi.fill(-9999)
	_dat_ndwi.fill(-9999)
	_dat_indx.fill(0)
	_dat_lyrs.fill(0)

	# cdef float _n
	# _n = 0

	import progress_percentage
	_ppp = progress_percentage.progress_percentage(_rows)

	for _row in xrange(_rows):
		_ppp.next()

		for _col in xrange(_cols):
			# if not(_row == 3404 and _col == 2849):
			# 	continue

			_vi, _nu = read_pixel(_bnd_out, _row, _col)
			if _nu == 0:
				_dat_indx[_row, _col] = 255
				_dat_lyrs[_row, _col] = 255
				continue

			# inavailable pixels
			_vr = []
			_is = []

			_vr.append(_vi)
			_is.append(0)

			for _i in xrange(1, len(_bnds)):
				_vi, _nu = read_pixel(_bnds[_i], _row, _col)
				if _nu == 4:
					_vr.append(_vi)
					_is.append(_i)

			# check NDVI
			_b_ndvi = [cal_ndvi(_vr[_b][2], _vr[_b][1], _vr[_b][0]) for _b in xrange(len(_vr))]
			_v_ndvi, _i_ndvi = pick_values(_b_ndvi, 0.8, 200)
			_dat_ndvi[_row, _col] = _v_ndvi

			# check NDWI
			_b_ndwi = [cal_ndwi(_vr[_b][0], _vr[_b][2]) for _b in xrange(len(_vr))]
			_v_ndwi, _i_ndwi = pick_values(_b_ndwi, 0.7, 200)
			_dat_ndwi[_row, _col] = _v_ndwi

			# if _row == 3404 and _col == 2849:
			# if _col == 2821-1 and _row == 2813-1:
			# 	for _rr in _vr:
			# 		print 'SR', _rr
			# 	print 'NDWI', _b_ndwi, _v_ndwi
			# 	print 'NDVI', _b_ndvi, _v_ndvi

			if max(_v_ndvi, _v_ndwi) > 100:
				_i_idx = _i_ndvi if _v_ndvi > _v_ndwi else _i_ndwi
				if _i_idx >= 0:
					if _i_idx > 0:
						save_pixel(_bnd_out, _row, _col, _vr[_i_idx])
						_dat_lyrs[_row, _col] = _is[_i_idx]
					_dat_indx[_row, _col] = 1
					continue


			# pick the mediean value
			_b_swir = [_vr[_b][3] for _b in xrange(len(_vr))]
			_v_swir, _i_swir = pick_values(_b_ndvi, 0.8, 0)
			if _i_swir >= 0:
				if _i_swir > 0:
					save_pixel(_bnd_out, _row, _col, _vr[_i_swir])
					_dat_lyrs[_row, _col] = _is[_i_swir]
				_dat_indx[_row, _col] = 3
				continue

			_dat_indx[_row, _col] = 4

	_ppp.done()

	return _bnd_out, \
				_bnd_out[0].from_grid(_dat_ndvi, nodata=-9999), \
				_bnd_out[0].from_grid(_dat_ndwi, nodata=-9999), \
				_bnd_out[0].from_grid(_dat_indx, nodata=255), \
				_bnd_out[0].from_grid(_dat_lyrs, nodata=255)

def combine_bak(bs1, bs2):
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


