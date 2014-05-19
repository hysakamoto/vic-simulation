clear all

syms x y z t
syms tau stratio

tau = 20;
stratio = 2.0;
% z = Z*stratio at time t=tau

% permeability
syms k
K = eye(3)*k;

% X = x;
% Y = y;
% Z = z/(x*y*(exp(t)-1)+1);

X = x/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2); 
Y = y/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
Z = z/(1 - ((t - tau)^2/tau^2 - 1)*(stratio - 1));
p = -(z^2*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/(2*k);



ux = x-X;
uy = y-Y;
uz = z-Z;

u = [ux;uy;uz];
v = diff(u,t);

invF = [gradient(X,[x,y,z]).'; gradient(Y,[x,y,z]).'; gradient(Z,[x,y,z]).'];
F = inv(invF);



% pressure
% p = 0;
% p = x*y*z*(exp(t)-1);

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
bf = bf.'*J;

%% Boundary conditions: Neumann
n_top = [0;0;1];
n_bottom = n_top*-1;
n_right = [1;0;0];
n_left = n_right*-1;
n_back = [0;1;0];
n_front = n_back*-1;

dsdS = @(n) (J*sqrt((n.'*invF)*(invF.'*n)));
gbar = @(n) (-K*gradient(p,[x,y,z])).'*n * dsdS(n);
tbar = @(n) (-p*n+Sigma_E*n) * dsdS(n);

gbar_top = gbar(n_top);
gbar_bottom = gbar(n_bottom);
gbar_right = gbar(n_right);
gbar_left = gbar(n_left);
gbar_back = gbar(n_back);
gbar_front = gbar(n_front);

tbar_top = tbar(n_top);
tbar_bottom = tbar(n_bottom);
tbar_right = tbar(n_right);
tbar_left = tbar(n_left);
tbar_back = tbar(n_back);
tbar_front = tbar(n_front);

gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front];
tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front];

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

% solve('x/(-1/((t - 20)^2/400 - 2))^(1/2)-X_=0', 'x')
% solve('y/(-1/((t - 20)^2/400 - 2))^(1/2)-Y_=0', 'y')
% solve('-z/((t - 20)^2/400 - 2)-Z_=0', 'z')

x_ = X_ *(-1/((t - 20)^2/400 - 2))^(1/2);
y_ = Y_ *(-1/((t - 20)^2/400 - 2))^(1/2);
z_ = -Z_ * ((t - 20)^2/400 - 2);

u_ = subs(u, [x,y,z], [x_,y_,z_]);
p_ = subs(p, [x,y,z], [x_,y_,z_]);
v_ = subs(v, [x,y,z], [x_,y_,z_]);

source_ = subs(source, [x,y,z], [x_,y_,z_]);
bf_ = subs(bf, [x,y,z], [x_,y_,z_]);

u_initial_ = subs(u_initial, [x,y,z], [x_,y_,z_]);
p_initial_ = subs(p_initial, [x,y,z], [x_,y_,z_]);

gbars_ = subs(gbars, [x,y,z], [x_,y_,z_]);
tbars_ = subs(tbars, [x,y,z], [x_,y_,z_]);

disp('u = ')
disp(u_);
disp('p = ')
disp(p_);
disp('v = ')
disp(v_);
disp('source = ')
disp(source_);
disp('bf = ')
disp(bf_);
disp('u_initial = ')
disp(u_initial_);
disp('p_initial = ')
disp(p_initial_);
disp('gbars = ')
disp(gbars_.');
disp('tbars = ')
disp(tbars_.');


