from constants import *



def gf_generate(exp, log, m):
    elt = 1
    alpha = 2  # primitive element of GF(2^PARAM_M)
    gf_poly = PARAM_GF_POLY

    for i in range((1 << m) - 1):#tocheck czy na pewno dobrze przetlumaczylo tego bit
        exp[i] = elt
        log[elt] = i

        elt *= alpha
        if elt >= 1 << m:
            elt ^= gf_poly

    exp[(1 << m) - 1] = 1
    exp[1 << m] = 2
    exp[(1 << m) + 1] = 4
    log[0] = 0  # by convention
    return exp, log


def gf_mul(a: int, b: int):
    mask = ((a != 0) & (b != 0)) * 0xFFFF
    return mask & gf_exp[(gf_log[a] + gf_log[b]) % (PARAM_GF_MUL_ORDER - 1)]


def gf_square(a):
    mask = (a != 0) * (-1)  # equivalent to (uint16_t) (-((int32_t) a) >> 31)
    return mask & gf_exp[gf_mod(2 * gf_log[a])]


def gf_inverse(a):
    mask = (a != 0) * 0xffff
    return mask & gf_exp[PARAM_GF_MUL_ORDER - gf_log[a]]


def gf_mod(i):
    tmp = i - PARAM_GF_MUL_ORDER
    mask = -((tmp >> 15) & 1)  # mask = 0xffff if i < PARAM_GF_MUL_ORDER
    return tmp + (mask & PARAM_GF_MUL_ORDER)