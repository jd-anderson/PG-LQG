function [gamma_lb]= update_gamma_rule(gamma,A,B,K,M_inv)


gamma_lb = gamma;

gammas = linspace(gamma, 1, 1000);

for i=1:length(gammas)

spectral_radius = max(abs(eig(sqrt(gammas(i))*(A-B*K*M_inv))));
if spectral_radius<1 && gammas(i)>gamma_lb
    gamma_lb = gammas(i);
elseif spectral_radius>1
    break;
end

 


end 

end