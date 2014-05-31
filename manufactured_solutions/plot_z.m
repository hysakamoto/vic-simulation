function [ output_args ] = plot_z( X1,X2,X3, t, tau, stratio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% x,y,z
phix = @(X) X*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
phiy = @(Y) Y*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
phiz = @(Z) Z - Z*((t - tau)^2/tau^2 - 1)*(stratio - 1);

%% pressure
k = 0.1;
p = @(x,y,z) (z^2*(2*t - 2*tau)*(stratio - 1))/(2*k*tau^2*(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1)^2*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(3/2)) - (z^2*(2*t - 2*tau)*(stratio - 1))/(2*k*tau^2*(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1)^2);

%% Cauchy stress
E = 100;
nu = 0.45;
mu = E/(2*(1+nu));
lm = E*nu/((1+nu)*(1-2*nu));
Sigma_zz = @(x,y,z) mu*((t^2 - stratio*t^2 - 2*t*tau + tau^2 + 2*stratio*t*tau)^2/tau^4 - 1) + (z^2*(2*t - 2*tau)*(stratio - 1))/(2*k*tau^2*(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1)^2) - (z^2*(2*t - 2*tau)*(stratio - 1))/(2*k*tau^2*(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1)^2*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(3/2));


x1 = zeros(size(X1));
x2 = zeros(size(X2));
x3 = zeros(size(X3));
szz = zeros(size(X3));
pres = zeros(size(X3));

for i=1:length(X1)
    for j=1:length(X2)
        for k=1:length(X3)
            x1(i,j,k) = phix(X1(i,j,k));
            x2(i,j,k) = phiy(X2(i,j,k));
            x3(i,j,k) = phiz(X3(i,j,k));
            
            szz(i,j,k) = Sigma_zz(x1(i,j,k),x2(i,j,k),x3(i,j,k));
            pres(i,j,k) = p(x1(i,j,k),x2(i,j,k),x3(i,j,k));
        end
    end
end

scatter3(reshape(x1,1,[]),reshape(x2,1,[]),reshape(x3,1,[]),100,reshape(pres,1,[]),'fill');
title(sprintf('t=%.1f',t));
xlim([0,1]);
ylim([0,1]);
zlim([0,2]);
colorbar;



end

