import numpy as np

def clover(de_prior, gini):
    """
    By this score, user can rank the genes which is:
    - Genes which is Low probability to be DEG
    - Reraness (tissue selective expression)

    Lower score -> Rera to be DEG and Tissue selective expression gene
    = 0.2 x np.log1p(1/0.9) = 0.0649
    Biger score -> used to be DEG and House keeping gene
    = 0.9 x np.log1p(1/0.2) = 0.7003
    """
    prob = de_prior # close to 1: More DEG
    ts = np.log2(2/(gini +1)) # 1: More housekeeping, 0: More TS
    return prob * ts
    

def dowsing(de_prior, gini, fdr):
    rare_deg = clover(de_prior, gini)
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