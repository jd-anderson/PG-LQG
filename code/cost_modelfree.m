function J = cost_modelfree(A,B,C,K,Q,R,p,V,W)
% Calculate cost with a simulator
% Monte Carlo
%tic
num=1; %num of Monte Carlo
step=100; 
[~,input_dim]=size(B);
[output_dim,state_dim]=size(C);


u=zeros(input_dim,p+step);
x=zeros(state_dim,p+step);
y=zeros(output_dim,p+step);


J=zeros(num,1);
for i=1:num
    %Initial conditons
    x(:,1)=0;
    u(:,1)=0*mvnrnd(zeros(input_dim,1),20*eye(input_dim))';
    y(:,1)=C*x(:,1);
    %initial p steps with stabilizing input
    for t=2:p
        u(:,t)=0*mvnrnd(zeros(input_dim,1),eye(input_dim))';
        w = mvnrnd(zeros(1, state_dim), W)';
        v = mvnrnd(zeros(1, output_dim), V)';
        x(:,t)=A*x(:,t-1)+B*u(:,t-1)+w;
        y(:,t)=C*x(:,t)+v;
    end
    % last steps
    for t=(p+1):(p+step)
        w = mvnrnd(zeros(1, state_dim), W)';
        v = mvnrnd(zeros(1, output_dim), V)';
        u_tp = reshape(u(:, (t-1):-1:(t-p)), input_dim * p, 1);
        y_tp = reshape(y(:, (t-1):-1:(t-p)), output_dim * p, 1);

        z_t=[u_tp;y_tp];
        u(:,t)=-K*z_t;
        x(:,t)=A*x(:,t-1)+B*u(:,t-1)+w;
        y(:,t)=C*x(:,t)+v;
        J(i)=J(i)+y(:,t)'*Q*y(:,t)+u(:,t)'*R*u(:,t);
    end


end

J=mean(J)/(step+p);
end

