from constants import *
from gf import gf_mul, gf_inverse
from fft import fft, fft_retrieve_error_poly


def mod(i, modulus):
    tmp = i - modulus
    mask = -(tmp >> 15) # mask is -1 if tmp is negative, 0 otherwise
    return tmp + (mask & modulus)


def compute_generator_poly(poly):
    poly[0] = 1
    tmp_degree = 0

    for i in range(1, 2 * PARAM_DELTA + 1):
        for j in range(tmp_degree, 0, -1):
            poly[j] = gf_exp[mod(gf_log[poly[j]] + i, PARAM_GF_MUL_ORDER)] ^ poly[j - 1]

        poly[0] = gf_exp[mod(gf_log[poly[0]] + i, PARAM_GF_MUL_ORDER)]
        tmp_degree += 1
        poly[tmp_degree] = 1

    print("")
    for i in range(PARAM_G):
        print(poly[i], end=", ")
    print("")
    return poly


def reed_solomon_encode(cdw, msg):
    gate_value = 0
    tmp = [0] * PARAM_G
    PARAM_RS_POLY = RS_POLY_COEFS

    msg_bytes = bytearray(msg.to_bytes(PARAM_K, 'little'))#tocheck dunno if little
    cdw_bytes = bytearray(PARAM_N1)

    for i in range(PARAM_K):
        gate_value = msg_bytes[PARAM_K - 1 - i] ^ cdw_bytes[PARAM_N1 - PARAM_K - 1]

        for j in range(PARAM_G):
            tmp[j] = gf_mul(gate_value, PARAM_RS_POLY[j])

        for k in range(PARAM_N1 - PARAM_K - 1, 0, -1):
            cdw_bytes[k] = cdw_bytes[k - 1] ^ tmp[k]

        cdw_bytes[0] = tmp[0]

    cdw_bytes[PARAM_N1 - PARAM_K:PARAM_N1] = msg_bytes[:PARAM_K]
    cdw[:] = int.from_bytes(cdw_bytes, 'little')
    return cdw


def compute_syndromes(syndromes, cdw):
    for i in range(2 * PARAM_DELTA):
        for j in range(1, PARAM_N1):
            syndromes[i] ^= gf_mul(cdw[j], alpha_ij_pow[i][j-1])
        syndromes[i] ^= cdw[0]
    return syndromes


def compute_elp(sigma, syndromes):
    deg_sigma = 0
    deg_sigma_p = 0
    deg_sigma_copy = 0
    sigma_copy = [0] * (PARAM_DELTA + 1)
    X_sigma_p = [0, 1] + [0] * PARAM_DELTA
    pp = (1 << 16) - 1  # 2*rho
    d_p = 1
    d = syndromes[0]

    for mu in range(2 * PARAM_DELTA):
        # Save sigma in case we need it to update X_sigma_p
        sigma_copy[:2 * (PARAM_DELTA)] = sigma[:2 * (PARAM_DELTA)]
        deg_sigma_copy = deg_sigma

        dd = gf_mul(d, gf_inverse(d_p))

        for i in range(1, min(mu+2, PARAM_DELTA+1)):
            sigma[i] ^= gf_mul(dd, X_sigma_p[i])

        deg_X = mu - pp
        deg_X_sigma_p = deg_X + deg_sigma_p

        # mask1 = 0xffff if(d != 0) and 0 otherwise
        mask1 = -(-d << 15)

        # mask2 = 0xffff if(deg_X_sigma_p > deg_sigma) and 0 otherwise
        mask2 = -((deg_sigma - deg_X_sigma_p) >> 15)

        # mask12 = 0xffff if the deg_sigma increased and 0 otherwise
        mask12 = mask1 & mask2
        deg_sigma ^= mask12 & (deg_X_sigma_p ^ deg_sigma)

        if mu == (2 * PARAM_DELTA - 1):
            break

        pp ^= mask12 & (mu ^ pp)
        d_p ^= mask12 & (d ^ d_p)

        for j in range(PARAM_DELTA, 0, -1):
            X_sigma_p[j] = (mask12 & sigma_copy[j-1]) ^ (~mask12 & X_sigma_p[j-1])

        deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p)
        d = syndromes[mu + 1]

        for j in range(1, min(mu+2, PARAM_DELTA+1)):#tocheck nie wiem czy nie +1 i +0 powinno byc
            d ^= gf_mul(sigma[j], syndromes[mu + 1 - j])

    return deg_sigma, sigma


def compute_roots(error, sigma):
    w = [0] * (1 << PARAM_M)

    w = fft(w, sigma, PARAM_DELTA + 1)
    error, w = fft_retrieve_error_poly(error, w)
    return error


def compute_z_poly(z, sigma, degree, syndromes):
    z[0] = 1
    for i in range(PARAM_DELTA + 1):
        # mask = int((i - degree - 1) >> 15) & 0xffff
        mask = -((i - degree - 1) >> 15)
        z[i] = mask & sigma[i]

    z[1] ^= syndromes[0]

    for i in range(2, PARAM_DELTA + 1):#tocheck nie wiem czy nie +0
        # mask = int((i - degree - 1) >> 15) & 0xffff
        mask = -((i - degree - 1) >> 15)
        z[i] ^= mask & syndromes[i - 1]

        for j in range(1, i):
            z[i] ^= mask & gf_mul(sigma[j], syndromes[i - j - 1])
    return z


def compute_error_values(error_values, z, error):
    beta_j = [0] * PARAM_DELTA
    e_j = [0] * PARAM_DELTA

    delta_counter = 0
    for i in range(PARAM_N1):
        found = 0
        mask1 = (-(error[i]) >> 31)
        for j in range(PARAM_DELTA):
            mask2 = ~(-(j ^ delta_counter) >> 31)
            beta_j[j] += mask1 & mask2 & gf_exp[i]
            found += mask1 & mask2 & 1
        delta_counter += found
    delta_real_value = delta_counter

    for i in range(PARAM_DELTA):
        tmp1 = 1
        tmp2 = 1
        inverse = gf_inverse(beta_j[i])
        inverse_power_j = 1

        for j in range(1, PARAM_DELTA + 1):
            inverse_power_j = gf_mul(inverse_power_j, inverse)
            tmp1 ^= gf_mul(inverse_power_j, z[j])
        for k in range(1, PARAM_DELTA):
            tmp2 = gf_mul(tmp2, (1 ^ gf_mul(inverse, beta_j[(i + k) % PARAM_DELTA])))
        mask1 = ((i - delta_real_value) >> 15)
        e_j[i] = mask1 & gf_mul(tmp1, gf_inverse(tmp2))

    delta_counter = 0
    for i in range(PARAM_N1):
        found = 0
        mask1 = (-(error[i]) >> 31)
        for j in range(PARAM_DELTA):
            mask2 = ~(-(j ^ delta_counter) >> 31)
            error_values[i] += mask1 & mask2 & e_j[j]
            found += mask1 & mask2 & 1
        delta_counter += found
    return error_values, z, error


def correct_errors(cdw, error_values):
    for i in range(PARAM_N1):
        cdw[i] ^= error_values[i]
    return cdw


def reed_solomon_decode(msg, cdw):
    cdw_bytes = bytearray(PARAM_N1)
    syndromes = [0] * (2 * PARAM_DELTA)
    sigma = [0] * (1 << PARAM_FFT)
    error = bytearray(1 << PARAM_M)
    z = [0] * PARAM_N1
    error_values = [0] * PARAM_N1

    # Copy the vector in an array of bytes
    cdw_bytes[:] = cdw[:PARAM_N1]

    # Calculate the 2*PARAM_DELTA syndromes
    syndromes = compute_syndromes(syndromes, cdw_bytes)

    # Compute the error locator polynomial sigma
    # Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room
    deg, sigma = compute_elp(sigma, syndromes)

    # Compute the error polynomial error
    error = compute_roots(error, sigma)

    # Compute the polynomial z(x)
    z = compute_z_poly(z, sigma, deg, syndromes)

    # Compute the error values
    error_values, z, error = compute_error_values(error_values, z, error)

    # Correct the errors
    cdw_bytes = correct_errors(cdw_bytes, error_values)

    # Retrieve the message from the decoded codeword
    msg[:] = cdw_bytes[PARAM_G - 1 : PARAM_G - 1 + PARAM_K]
    return msg