clear all

syms x y z t
syms tau stratio


X = x/(1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) + 1))^(1/2);
Y = y/(1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) + 1))^(1/2);
Z = z/(1 + ((t - tau)^2/tau^2 - 1)*(stratio - 1));

ux = X*(1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) + 1))^(1/2) - X;
uy = Y*(1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) + 1))^(1/2) - Y;
uz = Z*((t - tau)^2/tau^2 - 1)*(stratio - 1);

u = [ux;uy;uz];
v = diff(u,t);

invF = [gradient(X,[x,y,z]).'; gradient(Y,[x,y,z]).'; gradient(Z,[x,y,z]).'];
F = inv(invF);


% permeability
syms k
K = eye(3)*k;

% pressure
pzz = divergence(v,[x,y,z])/k;
pxx = 0; % px=const
pyy = 0; % py=const

pz = int(pzz,z); % pz(t=0) = 0
p = int(pz,z);

% check (==0)
if simplify(divergence(v-K*gradient(p,[x,y,z]),[x,y,z]))==0
    disp ('conservation of mass good!');
end

%% Elastic stress (neo-Hookean)

b = F*F.'; % left Cauchy-Green
J = det(F); % volume change (=1 for incompressible)
I = sym(eye(3)); % identity tensor

syms lm mu

Sigma_E = mu/J*(b-I) + lm/J*(log(J))*I;

Sigma = Sigma_E-p*I;

bf(1) = diff(Sigma(1,1),x);
bf(2) = diff(Sigma(2,2),y);
bf(3) = diff(Sigma(3,3),z);
bf = bf.';



