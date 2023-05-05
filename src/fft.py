from gf import *


def compute_fft_betas(betas):
    for i in range(PARAM_M - 1):#tocheck czy ++i zle ty wplywa
        betas[i] = 1 << (PARAM_M - 1 - i)
    return betas


def compute_subset_sums(subset_sums, set, set_size: int):
    subset_sums[0] = 0

    for i in range(set_size):
        for j in range(1 << i):
            subset_sums[(1 << i) + j] = set[i] ^ subset_sums[j]
    return subset_sums


def radix(f0, f1, f, m_f: int):
    if m_f == 4:
        f0[4] = f[8] ^ f[12]
        f0[6] = f[12] ^ f[14]
        f0[7] = f[14] ^ f[15]
        f1[5] = f[11] ^ f[13]
        f1[6] = f[13] ^ f[14]
        f1[7] = f[15]
        f0[5] = f[10] ^ f[12] ^ f1[5]
        f1[4] = f[9] ^ f[13] ^ f0[5]

        f0[0] = f[0]
        f1[3] = f[7] ^ f[11] ^ f[15]
        f0[3] = f[6] ^ f[10] ^ f[14] ^ f1[3]
        f0[2] = f[4] ^ f0[4] ^ f0[3] ^ f1[3]
        f1[1] = f[3] ^ f[5] ^ f[9] ^ f[13] ^ f1[3]
        f1[2] = f[3] ^ f1[1] ^ f0[3]
        f0[1] = f[2] ^ f0[2] ^ f1[1]
        f1[0] = f[1] ^ f0[1]

    elif m_f == 3:
        f0[0] = f[0]
        f0[2] = f[4] ^ f[6]
        f0[3] = f[6] ^ f[7]
        f1[1] = f[3] ^ f[5] ^ f[7]
        f1[2] = f[5] ^ f[6]
        f1[3] = f[7]
        f0[1] = f[2] ^ f0[2] ^ f1[1]
        f1[0] = f[1] ^ f0[1]

    elif m_f == 2:
        f0[0] = f[0]
        f0[1] = f[2] ^ f[3]
        f1[0] = f[1] ^ f0[1]
        f1[1] = f[3]

    elif m_f == 1:
        f0[0] = f[0]
        f1[0] = f[1]

    else:
        radix_big(f0, f1, f, m_f)

    return f0, f1


def radix_big(f0, f1, f, m_f):
    Q = [0] * (2 * (1 << (PARAM_FFT - 2)))
    R = [0] * (2 * (1 << (PARAM_FFT - 2)))
    Q0 = [0] * (1 << (PARAM_FFT - 2))
    Q1 = [0] * (1 << (PARAM_FFT - 2))
    R0 = [0] * (1 << (PARAM_FFT - 2))
    R1 = [0] * (1 << (PARAM_FFT - 2))
    n = 1
    n <<= (m_f - 2)
    Q[:2 * n] = f[3 * n:5 * n]
    Q[n:2 * n] = f[3 * n:4 * n]
    R[:4 * n] = f[:4 * n]

    for i in range(n):
        Q[i] ^= f[2 * n + i]
        R[n + i] ^= Q[i]

    Q0, Q1 = radix(Q0, Q1, Q, m_f - 1)
    R0, R1 = radix(R0, R1, R, m_f - 1)

    f0[:2 * n] = R0[:2 * n]
    f0[n:2 * n] = Q0[:2 * n]
    f1[:2 * n] = R1[:2 * n]
    f1[n:2 * n] = Q1[:2 * n]

    return f0, f1


def fft_rec(w, f, f_coeffs, m, m_f, betas):
    f0 = [0] * (1 << (PARAM_FFT - 2))
    f1 = [0] * (1 << (PARAM_FFT - 2))
    gammas = [0] * (PARAM_M - 2)
    deltas = [0] * (PARAM_M - 2)
    gammas_sums = [0] * (1 << (PARAM_M - 2))
    u = [0] * (1 << (PARAM_M - 2))
    v = [0] * (1 << (PARAM_M - 2))
    tmp = [0] * (PARAM_M - (PARAM_FFT - 1))
    beta_m_pow = 0
    x = 0
    # Step 1
    if m_f == 1:
        for i in range(m):
            tmp[i] = gf_mul(betas[i], f[1])

        w[0] = f[0]
        x = 1
        for j in range(m):
            for k in range(x):
                w[x + k] = w[k] ^ tmp[j]
            x <<= 1

        return

    # Step 2: compute g
    if betas[m - 1] != 1:
        beta_m_pow = 1
        x = 1
        x <<= m_f
        for i in range(1, x):
            beta_m_pow = gf_mul(beta_m_pow, betas[m - 1])
            f[i] = gf_mul(beta_m_pow, f[i])

    # Step 3
    f0, f1 = radix(f0, f1, f, m_f)

    # Step 4: compute gammas and deltas
    for i in range(m - 1):
        gammas[i] = gf_mul(betas[i], gf_inverse(betas[m - 1]))
        deltas[i] = gf_square(gammas[i]) ^ gammas[i]

    # Compute gammas sums
    gammas_sums = compute_subset_sums(gammas_sums, gammas, m - 1)

    # Step 5
    u, f0 = fft_rec(u, f0, (f_coeffs + 1) // 2, m - 1, m_f - 1, deltas)

    k = 1
    k <<= (m - 1) & 0xf
    if f_coeffs <= 3:  # 3-coefficient polynomial f case: f1 is constant
        w[0] = u[0]
        w[k] = u[0] ^ f1[0]
        for i in range(1, k):
            w[i] = u[i] ^ gf_mul(gammas_sums[i], f1[0])
            w[k + i] = w[i] ^ f1[0]
    else:
        v, f1 = fft_rec(v, f1, f_coeffs // 2, m - 1, m_f - 1, deltas)

        # Step 6
        w[k:] = v
        w[0] = u[0]
        w[k] ^= u[0]
        for i in range(1, k):
            w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i])
            w[k + i] ^= w[i]
    return w, f


def fft(w, f, f_coeffs):
    PARAM_M = 8  # replace with actual value
    PARAM_FFT = 16  # replace with actual value

    betas = [0] * (PARAM_M - 1)
    betas_sums = [0] * (1 << (PARAM_M - 1))
    f0 = [0] * (1 << (PARAM_FFT - 1))
    f1 = [0] * (1 << (PARAM_FFT - 1))
    deltas = [0] * (PARAM_M - 1)
    u = [0] * (1 << (PARAM_M - 1))
    v = [0] * (1 << (PARAM_M - 1))

    # Follows Gao and Mateer algorithm
    betas = compute_fft_betas(betas)  # replace with actual function call

    # Step 1: PARAM_FFT > 1, nothing to do

    # Compute gammas sums
    betas_sums = compute_subset_sums(betas_sums, betas, PARAM_M - 1)  # replace with actual function call

    # Step 2: beta_m = 1, nothing to do

    # Step 3
    f0, f1 = radix(f0, f1, f, PARAM_FFT)  # replace with actual function call

    # Step 4: Compute deltas
    for i in range(PARAM_M - 1):
        deltas[i] = gf_square(betas[i]) ^ betas[i]  # replace with actual function call

    # Step 5
    u, f0 = fft_rec(u, f0, (f_coeffs + 1) // 2, PARAM_M - 1, PARAM_FFT - 1, deltas)  # replace with actual function call
    v, f1 = fft_rec(v, f1, f_coeffs // 2, PARAM_M - 1, PARAM_FFT - 1, deltas)  # replace with actual function call

    k = 1 << (PARAM_M - 1)
    # Step 6, 7 and error polynomial computation
    w[k:2*k] = v[:k]

    # Check if 0 is root
    w[0] = u[0]

    # Check if 1 is root
    w[k] ^= u[0]

    # Find other roots
    for i in range(1, k):
        w[i] = u[i] ^ gf_mul(betas_sums[i], v[i])  # replace with actual function call
        w[k + i] ^= w[i]
    return w


def fft_retrieve_error_poly(error, w):
    gammas = [0] * (PARAM_M - 1)
    gammas_sums = [0] * (1 << (PARAM_M - 1))
    k = 1 << (PARAM_M - 1)

    gammas = compute_fft_betas(gammas)
    gammas_sums = compute_subset_sums(gammas_sums, gammas, PARAM_M - 1)

    error[0] ^= 1 ^ (-(w[0] >> 15) & 1)
    error[0] ^= 1 ^ (-(w[k] >> 15) & 1)

    for i in range(1, k):
        index = PARAM_GF_MUL_ORDER - gf_log[gammas_sums[i]]
        error[index] ^= 1 ^ (-(w[i] >> 15) & 1)

        index = PARAM_GF_MUL_ORDER - gf_log[gammas_sums[i] ^ 1]
        error[index] ^= 1 ^ (-(w[k + i] >> 15) & 1)
    return error, w
