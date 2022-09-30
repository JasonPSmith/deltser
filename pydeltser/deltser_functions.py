import numpy as np
from pydelt import run_deltser

def deltser(face_list, approx=False, approx_val=100000):

    return run_deltser(face_list, approx, approx_val, False, 'null')


def deltser_file(in_file, approx=False, approx_val=100000):

    return run_deltser(np.array([[[]]]), approx, approx_val, True, in_file)
