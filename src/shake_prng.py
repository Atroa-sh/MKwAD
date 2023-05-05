from constants import *
import hashlib



def seedexpander(state, output, outlen):
    bsize = 8  # size of uint64_t in bytes
    remainder = outlen % bsize
    tmp = bytearray(bsize)

    hashlib.
    shake256_inc_squeeze(output, outlen - remainder, state)

    if remainder != 0:
        shake256_inc_squeeze(tmp, bsize, state)
        output[-remainder:] = tmp[:remainder]

    return output