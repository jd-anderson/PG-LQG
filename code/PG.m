function [K,costs] = PG(A,B,C,M_inv,Q,R,p,K_ini,stepsize,num_iter,V,W)
% policy gradient without model

K=K_ini;

a =1
sigma_0=eye(4);
costs=[];
r = 0.1;
for i=1:num_iter
    delta_K=gradient_modelfree(A,B,C,Q,R,K,p,r,100,V,W);%r=0.2
    K=K-stepsize*delta_K;
    Acl = A-B*K*M_inv;
    eig_vals = abs(eig(Acl));

    if mod(i,100)==1
        p_K_PG_output = P_K(A,B,M_inv,C'*Q*C,R,K,V,W);
        costs=[costs,cost(p_K_PG_output,sigma_0)]
    end
end
end

