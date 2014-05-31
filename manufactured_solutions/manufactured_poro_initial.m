clear all

syms X0 Y0 Z0 t

% permeability
syms k
K = eye(3)*k;

% kinematics
u=[ 0,0,Z0 - (Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))/(X0*Y0*(exp(t/20) - 1) + 1) - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)];
p = X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20));

% u = [-t/20*X0, -t/20*Y0, t/4*Z0];
% p = X0*Y0*Z0*sin(t/10*pi);

v = diff(u, t);
x = X0+u(1);
y = Y0+u(2);
z = Z0+u(3);

phi = [x;y;z];

F = [gradient(x,[X0,Y0,Z0]) gradient(y,[X0,Y0,Z0]) gradient(z,[X0,Y0,Z0])].';
invF = inv(F);
J = det(F); % volume change (=1 for incompressible)

gradv = [gradient(v(1),[X0,Y0,Z0]) gradient(v(2),[X0,Y0,Z0]) gradient(v(3),[X0,Y0,Z0])].';
gradp = gradient(p,[X0,Y0,Z0]);

%% source term
curvel = invF*J*K*invF.'*gradp; % current velocity-like parameter
source = J* trace(gradv.'*invF.') - divergence(curvel, [X0,Y0,Z0]);

%% Elastic stress (compressible neo-Hookean)
% elastic parameter
syms lm mu

b = F*F.'; % left Cauchy-Green
I = sym(eye(3)); % identity tensor
C = F.'*F;
invC = inv(C);

% second PK stress tensor (wood 148)
Se = mu*(I-invC)+lm*(log(J))*invC;
Stot = F*Se - p*J*invF.';

% body force
bf(1) = divergence(Stot(1,:), [X0,Y0,Z0]);
bf(2) = divergence(Stot(2,:), [X0,Y0,Z0]);
bf(3) = divergence(Stot(3,:), [X0,Y0,Z0]);
bf = bf.';


%% Boundary conditions: Neumann
% initial boundary normals
n_top = [0;0;1];
n_bottom = n_top*-1;
n_right = [1;0;0];
n_left = n_right*-1;
n_back = [0;1;0];
n_front = n_back*-1;

%% using area change: wood 99
% left Cauchy?Green def tensor: wood 85
% using are change by current n

gbar_0 = @(n0) -curvel.' * n0;
tbar_0 = @(n0) Stot * n0;


gbar_top = gbar_0(n_top);
gbar_bottom = gbar_0(n_bottom);
gbar_right = gbar_0(n_right);
gbar_left = gbar_0(n_left);
gbar_back = gbar_0(n_back);
gbar_front = gbar_0(n_front);

tbar_top = tbar_0(n_top);
tbar_bottom = tbar_0(n_bottom);
tbar_right = tbar_0(n_right);
tbar_left = tbar_0(n_left);
tbar_back = tbar_0(n_back);
tbar_front = tbar_0(n_front);

gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front];
tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front];

%% initial condition

u_initial = limit(u,t,0);
p_initial = subs(p,t,0);
v_initial = subs(v,t,0);

%% Convert to initial representations


%% Output to file
fileID = fopen('../exact_solutions.py','w');


%% u,p
fprintf(fileID, 'U = [''%s'',\n ''%s'',\n ''%s'',\n ''%s'']\n\n', ...
    char(u(1)), char(u(2)), char(u(3)), char(p));

%% v
fprintf(fileID, 'V = [''%s'',\n ''%s'',\n ''%s'']\n\n', ...
    char(v(1)), char(v(2)), char(v(3)));

%% source
fprintf(fileID, 'source = ''%s''\n\n', ...
    char(source));

%% body force
fprintf(fileID, 'bf = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(bf(1)), char(bf(2)), char(bf(3)));

%% u_initial
fprintf(fileID, 'u_initial = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(u_initial(1)), char(u_initial(2)), char(u_initial(3)) );

%% p_initial
fprintf(fileID, 'p_initial = ''%s''\n\n', char(p_initial));

%% v_initial
fprintf(fileID, 'v_initial = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(v_initial(1)), char(v_initial(2)), char(v_initial(3)) );

%% gbars
fprintf(fileID,'gbars = \\\n[');
for i=1:length(gbars)
    fprintf(fileID,'''%s''',char(gbars(i)));
    if i<length(gbars)
        fprintf(fileID,',\n');
    end
end
fprintf(fileID,']\n\n');

%% tbars
[n,m] = size(tbars);
fprintf(fileID,'tbars = \\\n[');
for i=1:m
    fprintf(fileID,'[');
    for j=1:n
        fprintf(fileID,'''%s''',char(tbars(j,i)));
        if j<n
            fprintf(fileID,',\n');
        end
    end
    fprintf(fileID,']');
    if i<m
        fprintf(fileID,', \n');
    end
end
fprintf(fileID,']\n\n');


fclose(fileID);



