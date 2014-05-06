function [ output_args ] = plot_z( X1,X2,X3, t, tau, stratio )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


phix = @(X) X*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
phiy = @(Y) Y*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
phiz = @(Z) Z - Z*((t - tau)^2/tau^2 - 1)*(stratio - 1);

x1 = zeros(size(X1));
x2 = zeros(size(X2));
x3 = zeros(size(X3));

for i=1:length(X1)
    for j=1:length(X2)
        for k=1:length(X3)
            x1(i,j,k) = phix(X1(i,j,k));
            x2(i,j,k) = phiy(X2(i,j,k));
            x3(i,j,k) = phiz(X3(i,j,k));
        end
    end
end

scatter3(reshape(x1,1,[]),reshape(x2,1,[]),reshape(x3,1,[]),100,reshape(x3,1,[]),'fill');
title(sprintf('t=%.1f',t));
xlim([0,1]);
ylim([0,1]);
zlim([0,2]);




end

