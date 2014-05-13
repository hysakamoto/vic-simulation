clear all

syms x y z t
syms tau stratio

tau = 20;
stratio = 2.0;
% z = Z*stratio at time t=tau

X = x/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2); 
Y = y/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
Z = z/(1 - ((t - tau)^2/tau^2 - 1)*(stratio - 1));

ux = x-X;
uy = y-Y;
uz = z-Z;

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
gradp = gradient(p,[x,y,z]);

% check (==0)
if simplify(divergence(v-K*gradient(p,[x,y,z]),[x,y,z]))==0
    disp ('conservation of mass good!');
end

%% Elastic stress (compressible neo-Hookean)

b = F*F.'; % left Cauchy-Green
J = det(F); % volume change (=1 for incompressible)
I = sym(eye(3)); % identity tensor

% elastic parameter
syms lm mu

% effective (elastic) Cauchy stress
Sigma_E = mu/J*(b-I) + lm/J*(log(J))*I;

% total Cauchy stress
Sigma = Sigma_E-p*I;

% body force
bf(1) = diff(Sigma(1,1),x);
bf(2) = diff(Sigma(2,2),y);
bf(3) = diff(Sigma(3,3),z);
bf = bf.';

%% Boundary conditions: Neumann
n_top = [0;0;1];
n_bottom = n_top*-1;
n_right = [1;0;0];
n_left = n_right*-1;
n_back = [0;1;0];
n_front = n_back*-1;

gbar_top = (-K*gradient(p,[x,y,z])).'*n_top; % inner product
gbar_bottom = (-K*gradient(p,[x,y,z])).'*n_bottom;
gbar_right = (-K*gradient(p,[x,y,z])).'*n_right; % inner product
gbar_left = (-K*gradient(p,[x,y,z])).'*n_left;
gbar_front = (-K*gradient(p,[x,y,z])).'*n_front; % inner product
gbar_back = (-K*gradient(p,[x,y,z])).'*n_back;

tbar_top = -p*n_top+Sigma_E*n_top;
tbar_bottom = -p*n_top+Sigma_E*n_bottom;
tbar_right = -p*n_right+Sigma_E*n_right;
tbar_left = -p*n_top+Sigma_E*n_left;
tbar_front = -p*n_front+Sigma_E*n_front;
tbar_back = -p*n_top+Sigma_E*n_back;


%% test
% 
% [X1,X2,X3]=ndgrid([0:0.1:1],[0:0.1:1],[0:0.1:1]);
% stratio = 2;
% tau = 1;
% 
% for t=0:0.1:1.2
%     plot_z( X1,X2,X3, t, tau, stratio );
%     pause
% end
%  
% 






