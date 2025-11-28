function [K,z] = warmup(gamma,A,B,C,Q,R,p,M_inv)

A_gamma = sqrt(gamma)*A;
B_gamma = sqrt(gamma)*B;


[~,input_dim]=size(B);
[output_dim,state_dim]=size(C);
K = zeros(input_dim,p*(input_dim+output_dim));
u=zeros(input_dim,p);
x=zeros(state_dim,p);
y=zeros(output_dim,p);


%Initial conditons
x0 = mvnrnd(zeros(4,1),20*eye(4))';

x(:,1)=x0;
u(:,1) = zeros(input_dim,1);
y(:,1)=C*x(:,1);
%initial p steps with stabilizing input
for t=2:p
    u(:,t)=zeros(input_dim,1);
    x(:,t)=sqrt(gamma)*A*x(:,t-1)+sqrt(gamma)*B*u(:,t-1);
    y(:,t)=C*x(:,t);
end
z = [u;y];

spectral_radius = [];
while gamma < 1

A_gamma = sqrt(gamma)*A;
B_gamma = sqrt(gamma)*B;
M_gamma=calculate_M(A_gamma,B_gamma,C,p);
M_inv_gamma=M_gamma'*inv(M_gamma*M_gamma');

[~,S,~]=dlqr(A_gamma,B_gamma,C'*Q*C,R);
cost_optimal=trace(S);    

eta = 10e-5;
costs=[];
K0 = K;
for i=1:100
    delta_K=gradient_modelfree(A_gamma,B_gamma,C,Q,R,K,p,0.01,10);
    K=K-eta*delta_K;
    sr = max(abs(eig(A_gamma - B_gamma*K*M_inv_gamma)));
    if sr >1
        K = K0;
        break;
    end
    p_K_PG_output = P_K(A_gamma,B_gamma,M_inv_gamma,C'*Q*C,R,K);
    costs=[costs,trace(p_K_PG_output)];
end


spectral_radius = [spectral_radius max(abs(eig((A-B*K*M_inv))))];

%update z
u=zeros(input_dim,p+20);
x=zeros(state_dim,p+20);
x(:,1)=x0;
y=zeros(output_dim,p+20);
u(:,1:p) = z(1:input_dim,:);
y(:,1:p) = z(input_dim+1:end,:);

for t=2:p
    x(:,t)=A_gamma*x(:,t-1)+B_gamma*u(:,t-1);
end

for t=(p+1):(p+20)
        u_tp=vec(u(:,(t-1):-1:(t-p)));
        y_tp=vec(y(:,(t-1):-1:(t-p)));
        z_t=[u_tp;y_tp];
        u(:,t)=-K*z_t;
        x(:,t)=A_gamma*x(:,t-1)+B_gamma*u(:,t-1);
        y(:,t)=C*x(:,t);
end


u_temp = u(:,end-p-1:end-p);
y_temp = y(:,end-p-1:end-p);

z = [u_temp;y_temp];

%update gamma
gamma = update_gamma_rule(gamma,A,B,K,M_inv); %need to change it to a update rule that uses the cost and remove M_inv

end



end