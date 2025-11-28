clear
close all

num=1; % experiment iter
state_dim=4;
input_dim=2;
output_dim=2;
maxiter=200000;


A = [-0.2639,    0.5924,   -0.6445,   -0.8047;
    0.5288,    0.4654,    0.6087,    0.0537;
     -0.2803,   -0.4883,    0.1135,   -0.6962;
     -1.0480,    0.2543,   -0.1278,   -0.1279];

B =[-0.9313,   -1.0678;
    2.0774,    0.3084;
   -1.4758,   -0.7451;
   -0.2621,   -1.5536];

C = [2.2795,   -0.6637,   -1.1390,   -0.8495;
    0.4608,    1.2424,    1.4244,  -1.3973];
W = 0.01 * eye(state_dim);
V = 0.01 * eye(output_dim);

%controllability
nc = rank(ctrb(A,B));
no = rank(obsv(A,C));

Q=eye(output_dim);
R=eye(input_dim);
sigma_0 = Sigma_K(A, C, W, V)

cost_PG_output=zeros(num,1);
cost_PG_output_rec=[];
time_PG_output=zeros(num,1);
Qc=C'*Q*C;
sigma_1=eye(state_dim);
[K_state_feedback,S,e]=dlqr(A,B,Qc,R);


%% pg_output
for iter=1:num
tic

p=calculate_p(A,B,C,state_dim);
z_dim=p*(input_dim+output_dim);

M=calculate_S(A, B, C, Q, R, W, V, p);
M_inv=M'*inv(M*M')





%gamma_0 = 0.1;
%[K0,z] = warmup(gamma_0,A,B,C,Q,R,p,M_inv);
%K_ini = K0;
%learned from the warmup - uncoment the warmup to generate another initial
K_ini = [0.2409   -0.5372    0.0746    0.5889   -0.3675   -0.5425    0.4393    0.3226
         0.5353   -0.1826    0.1608    0.3159   -0.0954   -0.1154    0.2875    0.2316];
Acl = A-B*K_ini*M_inv;
eig_vals = abs(eig(Acl));
p_K_PG_output = P_K(A,B,M_inv,C'*Q*C,R,K_ini);
cost_init = trace(p_K_PG_output*sigma_0)






stepsize=5e-9;
[K_PG_output,costs]=PG(A,B,C,M_inv,Q,R,p,K_ini,stepsize,maxiter,V,W);
p_K_PG_output = P_K(A,B,M_inv,Qc,R,K_PG_output,V,W);
cost_PG_output(iter) = cost(p_K_PG_output,sigma_0);

time_PG_output(iter)=toc;
cost_PG_output_rec=[cost_PG_output_rec;costs];
end


%% plot
figure;

mean_cost=mean(cost_PG_output_rec,1);
max_cost=max(cost_PG_output_rec);
min_cost=min(cost_PG_output_rec);

set(gca,'FontSize',14,'YScale','log');
xlabel('Iteration','FontSize',14)
ylabel('Cost','FontSize',14)
set(gcf,'unit','centimeters','position',[1,2,14,8])
hold on


temp=0:100:maxiter*-100;
plot(temp,mean_cost,'LineWidth',1.5,'color',[255 153 18]/255)
hold on
