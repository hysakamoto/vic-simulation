clear all

syms x y z t
syms tau stratio

tau = 20;
stratio = 2.0;
% z = Z*stratio at time t=tau

X = x;
Y = y;
Z = z/(x*y*(exp(t)-1)+1);

% X = x/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2); 
% Y = y/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
% Z = z/(1 - ((t - tau)^2/tau^2 - 1)*(stratio - 1));

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
p = x*y*z*(exp(t)-1);

%% source term
source = simplify(divergence(v-K*gradient(p,[x,y,z]),[x,y,z]));


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

%% body force
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

gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_front, gbar_back];
tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_front, tbar_back];

%% initial condition

u_initial = subs(u,t,0);
p_initial = subs(p,t,0);

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


%% Convert to initial representations


syms X_ Y_ Z_

source_ = subs(source, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);
bf_ = subs(bf, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);

u_initial_ = subs(u_initial, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);
p_initial_ = subs(p_initial, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);

gbars_ = subs(gbars, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);
tbars_ = subs(tbars, [x,y,z], [X_,Y_,X_*Y_*Z_*(exp(t)-1)+Z_]);

disp(source_)
disp(bf_)
disp(u_initial_)
disp(p_initial_)
disp(gbars_.')
disp(tbars_.')


