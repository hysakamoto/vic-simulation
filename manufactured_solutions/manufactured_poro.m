clear all

syms ux uy uz
syms tau
syms uz_top_max
syms t
syms W L H
syms stratio

uz_top = H*(stratio-1) *(1- 1/tau^2*(t-tau)^2);

syms X Y Z
uz = uz_top/H*Z;
z = uz+Z;

z_Z = diff(z,Z);
x_X = sqrt(1/z_Z);
y_Y = x_X;

% zero initial displacement
x = int(x_X,X);
y = int(y_Y,Y); 
% solve('Y*(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2)=y', 'Y') 
% to get Y in terms of y

ux = x-X;
uy = y-Y;

w = subs(x,X,W/2)-subs(x,X,-W/2);
l = subs(y,Y,L/2)-subs(y,Y,-L/2);
h = subs(z,Z,H)-0;

% check
if simplify(w*l*h-W*L*H)==0
    disp( 'incompressibility good!');
end

% displacement
u = [ux;uy;uz];

% velocity
v = diff(u,t);
gradv = [gradient(v(1),[X,Y,Z]).'; gradient(v(2),[X,Y,Z]).';gradient(v(3),[X,Y,Z]).'];

% deformation gradient
F = [gradient(x,[X,Y,Z]).'; gradient(y,[X,Y,Z]).'; gradient(z,[X,Y,Z]).'];
J = det(F);

% inverse deformation gradient
invF = inv(F);

% permeability
syms k
K = eye(3)*k;


A=J*trace(gradv*invF) % ==0 ---> grad_0p constant 

% pressure
pzz = divergence(v,[X,Y,Z])/k;
pxx = 0; % px=const
pyy = 0; % py=const

pz = int(pzz,Z); % pz(t=0) = 0
p = int(pz,Z);
gradp = gradient(p, [X,Y,Z]);

% check (==0)
if simplify(divergence(v-K*gradient(p,[X,Y,Z]),[X,Y,Z]))==0
    disp( 'conservation of mass good!');
end


%% Elastic stress (neo-Hookean)

b = F*F.'; % left Cauchy-Green
J = det(F); % volume change (=1 for incompressible)
I = sym(eye(3)); % identity tensor


% J*trace(gradv*invF) - divergence(invF*J*K*invF.'*gradp,[X,Y,Z])

syms lm mu

Sigma = mu/J*(b-I) + lm/J*(log(J))*I;

bf = divergence(-p*I+Sigma,[X,Y,Z]);

