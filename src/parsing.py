from constants import *
from vector import vect_set_random_fixed_weight, vect_set_random_fixed_weight_by_coordinates, vect_set_random


def hqc_secret_key_to_string(sk, sk_seed, pk):
    sk[:SEED_BYTES] = sk_seed
    sk[SEED_BYTES:] = pk
    return sk


def hqc_secret_key_from_string(x, y, pk, sk):
    sk_seedexpander = 0 #tocheck nie wiem o co chodzi z tym typem
    sk_seed = sk[:SEED_BYTES]

    seedexpander_init(sk_seedexpander, sk_seed, SEED_BYTES)#tocheck

    sk_seedexpander, x = vect_set_random_fixed_weight(sk_seedexpander, x, PARAM_OMEGA)
    sk_seedexpander, y = vect_set_random_fixed_weight_by_coordinates(sk_seedexpander, y, PARAM_OMEGA)
    pk[:] = sk[SEED_BYTES:SEED_BYTES + PUBLIC_KEY_BYTES]#tocheck czy nie wychodzi poza, za malo
    return x, y, pk


def hqc_public_key_to_string(pk, pk_seed, s):
    pk[:SEED_BYTES] = pk_seed[:SEED_BYTES]
    pk[SEED_BYTES:SEED_BYTES + VEC_N_SIZE_BYTES] = s[VEC_N_SIZE_BYTES:]
    return pk


def hqc_public_key_from_string(h, s, pk):
    pk_seedexpander = 0 #tocheck nie wiem o co chodzi z tym typem
    pk_seed = [0] * SEED_BYTES

    pk_seed[:SEED_BYTES] = pk[:SEED_BYTES]
    seedexpander_init(pk_seedexpander, pk_seed, SEED_BYTES)

    pk_seedexpander, h = vect_set_random(pk_seedexpander, h)
    s[:VEC_N_SIZE_BYTES] = pk[SEED_BYTES:SEED_BYTES + VEC_N_SIZE_BYTES]
    return h, s


def hqc_ciphertext_to_string(ct, u, v, d):
    ct[:VEC_N_SIZE_BYTES] = u
    ct[VEC_N_SIZE_BYTES:VEC_N_SIZE_BYTES+VEC_N1N2_SIZE_BYTES] = v
    ct[VEC_N_SIZE_BYTES+VEC_N1N2_SIZE_BYTES:] = d
    return ct


def hqc_ciphertext_from_string(u, v, d, ct):
    u[:VEC_N_SIZE_BYTES] = ct[:VEC_N_SIZE_BYTES]
    v[:] = ct[VEC_N_SIZE_BYTES:VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES]
    d[:] = ct[VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES:VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES + SHAKE256_512_BYTES]
    return u, v, d