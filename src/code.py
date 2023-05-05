from constants import *
from reed_muller import reed_muller_encode, reed_muller_decode
from reed_salomon import reed_solomon_encode, reed_solomon_decode

import numpy as np


def code_encode(em, m):
    tmp = np.zeros(VEC_N1_SIZE_64, dtype=np.uint64)

    encoded_message = reed_solomon_encode(tmp, m)
    encoded_tensor = reed_muller_encode(em, tmp)
    return encoded_tensor, encoded_message


def code_decode(m, em):
    tmp = np.zeros(VEC_N1_SIZE_64, dtype=np.uint64)

    decoded_message = reed_muller_decode(tmp, em)
    decoded_tensor = reed_solomon_decode(m, tmp)

    return decoded_message, decoded_tensor

