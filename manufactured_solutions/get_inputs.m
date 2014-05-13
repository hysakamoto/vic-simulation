%% compute inputs from manufatrued_poro_current and manufactured_poro 
%% results

clear all

%% material parameters
% elastic parameters
E = 100;
nu = 0.45;
lm = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
% permeability
k = 1.0; 

% simulation parameters
tau = 20.0; % simulation time
stratio = 2.0; % final stretch ratio

%% initial-to-current
syms X Y Z t

x = X*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
y = Y*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
z = Z - Z*((t - tau)^2/tau^2 - 1)*(stratio - 1);

%% solutions!!

% displacements
ux = x-X;
uy = y-Y;
uz = z-Z;
u = [ux;uy;uz];
% pressure
p  = -(z^2*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/(2*k);
gradp = [0; 0; ...
    -(z*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/k];

%% kinematics
F = [gradient(x,[X,Y,Z]) gradient(y,[X,Y,Z]) gradient(z,[X,Y,Z])].';
invF = inv(F);

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
body_force = [0.0, 0.0,...
(z*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 ...
- (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/k];

body_force_0 = body_force * J; % wood 129

%% boundary conditions - Neumann

% normals
n_top = [0;0;1];
n_bottom = n_top*-1;
n_right = [1;0;0];
n_left = n_right*-1;
n_back = [0;1;0];
n_front = n_back*-1;
ns = [n_top n_bottom n_right n_left n_back n_front];

% area change
dsdS = @(n0) J*sqrt(n0.'*invF*invF.'*n0);

% normal flux
K = eye(3)*k;
gbar = @(n) (-K*gradp).'*n;
gbar_0 = @(n) gbar(n)*dsdS(n); % assuming that n does not change.

gbars = sym(zeros(1,6));
gbars_0 = sym(zeros(1,6));

for i=1:6
   gbars(i) = gbar(ns(:,i));
   gbars_0(i) = gbar_0(ns(:,i));
end

% traction
% dadA = J/(sqrt(n*bn))
% t0 = t(dadA) wood 129
tbar = @(n) Sigma*n;
tbar_0 = @(n0) tbar(n0)*dsdS(n0); % assuming that normal direction does not change (second n0 should be n)

tbars_0 = sym(zeros(3,6));
tbars = sym(zeros(3,6));
for i=1:6
    tbars_0(:,i) = tbar_0(ns(:,i));
    tbars(:,i) = tbar(ns(:,i));
end


%% boundary conditions - Dirichlet

% displacement bcs
u_left = subs(u,X,0);
u_right = subs(u,X,1);
u_front = subs(u,Y,0);
u_back=subs(u,Y,1);
u_bottom = subs(u,Z,0);
u_top = subs(u,Z,1);

% pressure bc
p = -(z^2*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/(2*k);
p_left = subs(p,X,0);
p_right = subs(p,X,1);
p_front = subs(p,Y,0);
p_back = subs(p,Y,1);
p_top = subs(p,Z,0);
p_bottom = subs(p,Y,1);

%% initial conditions
ux_init = subs(ux, t,0);
uy_init = subs(ux, t,0);
uz_init = subs(ux, t,0);
p_init = subs(p, t, 0);



