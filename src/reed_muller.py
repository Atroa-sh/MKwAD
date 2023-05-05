from constants import *


MULTIPLICITY = CEIL_DIVIDE(PARAM_N2, 128)

def BIT0MASK(x):
    return -((x) & 1)

class codeword:
    def __init__(self):
        self.u8 = [0] * 16
        self.u32 = [0] * 4


def encode(word: codeword, message):
    first_word = BIT0MASK(message >> 7)

    first_word ^= BIT0MASK(message >> 0) & 0xaaaaaaaa
    first_word ^= BIT0MASK(message >> 1) & 0xcccccccc
    first_word ^= BIT0MASK(message >> 2) & 0xf0f0f0f0
    first_word ^= BIT0MASK(message >> 3) & 0xff00ff00
    first_word ^= BIT0MASK(message >> 4) & 0xffff0000

    word.u32[0] = first_word

    first_word ^= BIT0MASK(message >> 5)
    word.u32[1] = first_word
    first_word ^= BIT0MASK(message >> 6)
    word.u32[3] = first_word
    first_word ^= BIT0MASK(message >> 5)
    word.u32[2] = first_word
    return word


def hadamard(src, dst):#tocheck nie wiem jak mamy modyfikowac src i dst w ten sposob
    p1 = src
    p2 = dst
    for _ in range(7):
        for i in range(64):
            p2[i] = p1[2*i] + p1[2*i+1]
            p2[i+64] = p1[2*i] - p1[2*i+1] #tocheck dziwnie w oryginale tam zagladaja do tych tablic
        # swap p1, p2 for next round
        p3 = p1
        p1 = p2
        p2 = p3
    return src, dst


def expand_and_sum(dest, src):
    # start with the first copy
    for part in range(4):
        for bit in range(32):
            dest[part * 32 + bit] = src[0].u32[part] >> bit & 1#tocheck przetlumaczylo z nawiasami od src do bit

    # sum the rest of the copies
    for copy in range(1, MULTIPLICITY):
        for part in range(4):
            for bit in range(32):
                dest[part * 32 + bit] += src[copy].u32[part] >> bit & 1
    return dest


def find_peaks(transform):
    peak_abs_value = 0
    peak_value = 0
    peak_pos = 0
    for i in range(128):
        # get absolute value
        t = transform[i]
        pos_mask = -(t > 0)#tocheck looks sus
        absolute = (pos_mask & t) | (~pos_mask & -t)
        # all Python versions have a conditional expression
        peak_value = t if absolute > peak_abs_value else peak_value
        peak_pos = i if absolute > peak_abs_value else peak_pos
        peak_abs_value = absolute if absolute > peak_abs_value else peak_abs_value
    # set bit 7
    peak_pos |= 128 * (peak_value > 0)
    return peak_pos


def reed_muller_encode(cdw, msg):
    message_array = msg #tocheck dziwne te copy
    code_array = cdw
    for i in range(VEC_N1_SIZE_BYTES):
        # fill entries i * MULTIPLICITY to (i+1) * MULTIPLICITY
        pos = i * MULTIPLICITY
        # encode first word
        code_array[pos] = encode(code_array[pos], message_array[i])
        # copy to other identical codewords
        for copy in range(1, MULTIPLICITY):
            code_array[pos + copy] = code_array[pos]#tocheck chyba zadziala ale not sure
    return cdw, msg


def reed_muller_decode(msg, cdw):
    message_array = msg
    codeArray = cdw
    expanded = 0
    for i in range(VEC_N1_SIZE_BYTES):
        # collect the codewords

        expanded = expand_and_sum(expanded, codeArray[i*MULTIPLICITY])
        # apply hadamard transform
        transform = bytearray()
        expanded, transform = hadamard(expanded, transform)
        # fix the first entry to get the half Hadamard transform
        transform[0] -= 64 * MULTIPLICITY
        # finish the decoding
        message_array[i] = find_peaks(transform)
    return msg, cdw