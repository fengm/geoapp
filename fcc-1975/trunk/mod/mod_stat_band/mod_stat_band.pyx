import collections
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def stat(bnd, bnd_qa, int qa):
	cdef np.ndarray[np.uint8_t, ndim=2] _dat_qa = bnd_qa.data
	cdef np.ndarray[np.int16_t, ndim=2] _dat_va = bnd.data
	cdef int _rows = bnd_qa.height, _cols = bnd_qa.width
	cdef int _row, _col, _v_va, _v_qa

	_vs = collections.defaultdict(lambda: 0.0)
	for _row in xrange(_rows):
		for _col in xrange(_cols):
			_v_qa = _dat_qa[_row, _col]
			if _v_qa != qa:
				continue

			_v_va = _dat_va[_row, _col]
			_vs[_v_va] += 1.0
	
	return _vs

