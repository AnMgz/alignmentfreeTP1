from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, encode_kmer


def compare_kmers(km1,km2):
    for i in range(len(km1)):
        if km1[i]!=km2[i]:
            return False
    return True


def create_index(file, k):
    index = {}
    for kmer in stream_kmers(file,k):
        kmer = kmer2str(kmer,k)
        if kmer in index:
            index[kmer] = index[kmer] + 1
        else :
            index[kmer] = 1
    return index

def jaccard(fileA, fileB, k):
    j = 0
    index = create_index(fileA,k)
    intersect = 0
    union = sum(index.values)
    for idx,kmer in enumerate(stream_kmers(fileB, k)):
        if kmer in index and index[kmer]>0:
            intersect += 1
            index[kmer] -= 1
        else:
            union += 1
    j = intersect/union
    return j



if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---

            j = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], j)
