from constants import *
import numpy as np

TABLE = 16
WORD = 64


def swap(tab, elt1, elt2):
    tmp = tab[elt1]
    tab[elt1] = tab[elt2]
    tab[elt2] = tmp
    return tab


def reduce(o, a):
    for i in range(VEC_N_SIZE_64):
        r = a[i + VEC_N_SIZE_64 - 1] >> (PARAM_N & 0x3F)
        carry = (a[i + VEC_N_SIZE_64] << (64 - (PARAM_N & 0x3F))) & 0xFFFFFFFFFFFFFFFF
        o[i] = a[i] ^ r ^ carry
    o[VEC_N_SIZE_64 - 1] &= RED_MASK
    return o, a


def fast_convolution_mult(o, a1, a2, weight, ctx):
    carry = 0
    dec, s = 0, 0
    table = [0] * (TABLE * (VEC_N_SIZE_64 + 1))
    permuted_table = [0] * TABLE
    permutation_table = [i for i in range(TABLE)]
    permuted_sparse_vect = [0] * PARAM_OMEGA_E
    permutation_sparse_vect = [0] * PARAM_OMEGA_E

    seedexpander(ctx,  permutation_table, TABLE << 1) #tocheck

    for i in range(TABLE - 1):
        swap(permuted_table + i, 0, permutation_table[i] % (TABLE - i))

    pt = table[permuted_table[0] * (VEC_N_SIZE_64 + 1):]
    pt[:VEC_N_SIZE_64] = a2[:]
    pt[VEC_N_SIZE_64] = 0

    for i in range(1, TABLE):
        carry = 0
        idx = permuted_table[i] * (VEC_N_SIZE_64 + 1)
        pt = table[idx:idx + VEC_N_SIZE_64 + 1]

        for j in range(VEC_N_SIZE_64):
            pt[j] = (a2[j] << i) ^ carry
            carry = (a2[j] >> (WORD - i))

        pt[VEC_N_SIZE_64] = carry

    for i in range(weight):
        permuted_sparse_vect[i] = i

    seedexpander(ctx, permutation_sparse_vect, weight << 1) #tocheck

    for i in range(weight - 1):
        swap(permuted_sparse_vect + i, 0, permutation_sparse_vect[i] % (weight - i))

    res_16 = bytearray(o)
    # for i in range(weight):
    #     carry = 0
    #     dec, s = divmod(a1[permuted_sparse_vect[i]], 16)
    #
    #     pt = table[permuted_table[dec] * (VEC_N_SIZE_64 + 1):]
    #
    #     for j in range(VEC_N_SIZE_64 + 1):
    #         tmp = int.from_bytes(res_16[s * 16:s * 16 + 8], byteorder='little')
    #         tmp ^= pt[j]
    #         res_16[s * 16:s * 16 + 8] = tmp.to_bytes(8, byteorder='little')
    #         s += 4
    for i in range(weight):
        carry = 0x0
        dec = a1[permuted_sparse_vect[i]] & 0xf
        s = a1[permuted_sparse_vect[i]] >> 4

        res_16 = o.view(dtype=np.uint16)[s * 8:] #tocheck
        pt = table[permuted_table[dec] * (VEC_N_SIZE_64 + 1):]

        for j in range(VEC_N_SIZE_64 + 1):
            tmp = int.from_bytes(res_16[:8], byteorder='little')
            tmp ^= pt[j]
            res_16[:8] = tmp.to_bytes(8, byteorder='little')
            res_16 = res_16[4:]
    return o


def vect_mul(o, a1, a2, weight, ctx):
    tmp = [0] * ((VEC_N_SIZE_64 << 1) + 1)

    tmp = fast_convolution_mult(tmp, a1, a2, weight, ctx)
    o, tmp = reduce(o, tmp)
    return o
