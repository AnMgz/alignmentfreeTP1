
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def encode_kmer(text, k):
    kmer = 0
    for letter in text[0:k]:
        kmer<<=2
        kmer += encode_nucl(letter)
    return kmer

def stream_kmers(text, k):
    """enumerate_kmer"""
    mask = (1<<(2*(k-1)))-1
    kmer = encode_kmer(text,k)
    for i in range(len(text)-(k)):
        yield kmer
        kmer &= mask
        kmer <<= 2
        kmer += encode_nucl(text[i+k])
    yield kmer

def encode_nucl(letter):
    encoding = {'A':0, 'C':1, 'T':2, 'G':3}
    return encoding[letter]
