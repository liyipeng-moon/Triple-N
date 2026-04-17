function vec = gen_rdm_vec(INRDM)
% Extracts the lower-triangular (excluding diagonal) elements of a
% 1000×1000 representational dissimilarity matrix (INRDM) and returns them
% as a vector for further analysis.

LOC = find(tril(ones(size(INRDM)),-1));
vec = INRDM(LOC);