from constants import *
from shake_prng import seekexpander
from typing import List
import numpy as np
import random


def vect_set_random_fixed_weight_by_coordinates(ctx, v, weight):
    random_bytes_size = 3 * weight
    rand_bytes = bytearray([0] * (3 * PARAM_OMEGA_R)) # weight is expected to be <= PARAM_OMEGA_R
    inc = 0
    i, j = 0, random_bytes_size

    while i < weight:
        while True:
            if j == random_bytes_size:
                seedexpander(ctx, rand_bytes, random_bytes_size)
                j = 0

            v[i] = ((rand_bytes[j]) << 16) | ((rand_bytes[j+1]) << 8) | (rand_bytes[j+2])
            j += 3
            #tocheck
            if v[i] < UTILS_REJECTION_THRESHOLD:
                break

        v[i] = v[i] % PARAM_N

        inc = 1
        for k in range(i):
            if v[k] == v[i]:
                inc = 0
                break

        i += inc

    return ctx, v

def vect_set_random_fixed_weight(ctx, v, weight):
    # PARAM_OMEGA_R is assumed to be a constant defined elsewhere
    tmp = np.zeros((PARAM_OMEGA_R,), dtype=np.uint32)

    ctx, tmp = vect_set_random_fixed_weight_by_coordinates(ctx, tmp, weight)

    for i in range(weight):
        index = tmp[i] // 64
        pos = tmp[i] % 64
        v[index] |= (1 << pos)

    return ctx, v


def vect_set_random(ctx, v):
    rand_bytes = np.zeros(VEC_N_SIZE_BYTES, dtype=np.uint8)

    seedexpander(ctx, rand_bytes, VEC_N_SIZE_BYTES)

    v = np.frombuffer(rand_bytes, dtype=np.uint64)
    v[VEC_N_SIZE_64 - 1] &= BITMASK(PARAM_N, 64)

    return ctx, v


def vect_set_random_from_prng(v: bytearray) -> None:#tocheck
    rand_bytes = bytearray(VEC_K_SIZE_BYTES)
    shake_prng(rand_bytes, VEC_K_SIZE_BYTES)
    v[:VEC_K_SIZE_BYTES] = rand_bytes


def vect_add(v1, v2, size):
    o = [0]*size
    for i in range(size):
        o[i] = v1[i] ^ v2[i]
    return o


def vect_compare(v1, v2, size):
    r = 0

    for i in range(size):
        r |= v1[i] ^ v2[i]

    r = (~r + 1) >> 63
    return int(r & 0x01)#tocheck


def vect_resize(o, size_o, v, size_v):
    mask = 0x7FFFFFFFFFFFFFFF
    val = 0
    if size_o < size_v:
        if size_o % 64:
            val = 64 - (size_o % 64)

        o[:VEC_N1N2_SIZE_64] = v[:VEC_N1N2_SIZE_64]

        for i in range(val):
            o[VEC_N1N2_SIZE_64 - 1] &= (mask >> i)
    else:
        o[:CEIL_DIVIDE(size_v, 8)] = v[:CEIL_DIVIDE(size_v, 8)]

    return o


def vect_print(v, size):
    if size == VEC_K_SIZE_BYTES:
        tmp = bytearray(v)
        print(''.join(format(x, '02x') for x in tmp))#tocheck
    elif size == VEC_N_SIZE_BYTES:
        tmp = bytearray(v)
        print(''.join(format(x, '02x') for x in tmp))
    elif size == VEC_N1N2_SIZE_BYTES:
        tmp = bytearray(v)
        print(''.join(format(x, '02x') for x in tmp))
    elif size == VEC_N1_SIZE_BYTES:
        tmp = bytearray(v)
        print(''.join(format(x, '02x') for x in tmp))


def vect_print_sparse(v, weight):
    for i in range(weight-1):
        print(v[i], end=", ")
    print(v[weight-1])
