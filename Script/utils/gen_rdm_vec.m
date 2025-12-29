function vec = gen_rdm_vec(INRDM)

LOC = find(tril(ones(1000),-1));
vec = INRDM(LOC);