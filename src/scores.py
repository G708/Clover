import numpy as np

def glint(de_prior, gini):
    """
    By this score, the user can rank the genes which is:
    - Genes with is Low probability of being DEG
    - Reraness (tissue-specific expression)

    Lower score -> Rera to be DEG and tissue-specific expressed gene
    = 0.2 x np.log2(1/0.9) = 0.0649
    Larger score -> used to be DEG and Housekeeping gene
    = 0.9 x np.log2(1/0.2) = 0.7003
    """
    prob = de_prior # close to 1: More DEG
    ts = np.log2(2/(gini +1)) # 1: More housekeeping, 0: More TS
    return prob * ts
    

def dowsing(de_prior, gini, fdr):
    rare_deg = glint(de_prior, gini)
    rareness = np.log2(2/(rare_deg+1))
    return rareness * -np.log10(fdr)

def treasure_hunt(g2p_rank, dowsing):
    a = 1 - g2p_rank
    b = dowsing
    return a * b

def ropeway(g2p_rank, dowsing):
    a = g2p_rank
    b = dowsing
    return a * b