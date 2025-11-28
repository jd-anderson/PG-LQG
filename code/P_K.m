function p_K = P_K(A,B,M_inv,Q,R,K,V,W)
% calculate p_K 
K*M_inv
Q
p_K=dlyap((A-B*K*M_inv)',Q+M_inv'*K'*R*K*M_inv);
end

