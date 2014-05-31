clear all

syms x y z t
syms tau stratio

tau = 20;
stratio = 2.0;
% z = Z*stratio at time t=tau

% permeability
syms k
K = eye(3)*k;

X = x;
Y = y;
% Z = z;
Z = ((4*z*exp(t/40) - 4*z + 1)^(1/2) - 1)/(2*(exp(t/40) - 1));
p = x*y*z*(sin(t/10*pi));

Z = z/(x*y*(exp(t/20)-1)+1);
% Z = z - X*Y*(exp(t/20) - 1); %% working fine
% p = x*y*z*(exp(t/20)-1);



% 
% X = x/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2); 
% Y = y/(-1/(((t - tau)^2/tau^2 - 1)*(stratio - 1) - 1))^(1/2);
% Z = z/(1 - ((t - tau)^2/tau^2 - 1)*(stratio - 1));
% p = -(z^2*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2)))/(2*k);

ux = x-X;
uy = y-Y;
uz = z-Z;

u = [ux;uy;uz];
v = diff(u,t);

invF = [gradient(X,[x,y,z]).'; gradient(Y,[x,y,z]).'; gradient(Z,[x,y,z]).'];
F = inv(invF);


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
bf(1) = divergence(Sigma(1,:), [x,y,z]);
bf(2) = divergence(Sigma(2,:), [x,y,z]);
bf(3) = divergence(Sigma(3,:), [x,y,z]);
bf = bf.';

%% Boundary conditions: Neumann
% assuming that normal direction does not change (second n0 should be n)
n_top = [0;0;1];
n_bottom = n_top*-1;
n_right = [1;0;0];
n_left = n_right*-1;
n_back = [0;1;0];
n_front = n_back*-1;

%% using area change: wood 99
% left Cauchy?Green def tensor: wood 85
% % using are change by current n
% b = F*F.';
% dadA = @(n) J/sqrt(n.'*b*n); % wood 129, 82
% tbar = @(n) (-p*I+Sigma_E)*n *dadA(n); % use current normal 
% gbar = @(n) (-K*gradient(p,[x,y,z])).'*n * dadA(n);
% 
% % using area change by initial n0
% dsdS = @(n0) (J*sqrt((n0.'*inv(F.'*F)*n0))); % Un penetration: (11)
% ncur = @(n0) invF.'*n0 / sqrt(n0.'*inv(F.'*F)*n0);
% tbar_0 = @(n0) (-p*I+Sigma_E)*ncur(n0) *dsdS(n0);
% gbar_0 = @(n0) (-K*gradient(p,[x,y,z])).'*ncur(n0) * dsdS(n0);
%%

%% using wood
%% http://en.wikipedia.org/wiki/Stress_measures, wood 133
% S = J*invF*Sigma*invF.'
% P = J*Signa*invF.'
% t0 = F*S*N = J*Sigma*invF.'*N = P*N
gbar = @(n0) (-K*gradient(p,[x,y,z])).'* (J*invF.'*n0);
tbar = @(n0) (J*(-p*I+Sigma_E)*invF.') *n0;




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

u_initial = limit(u,t,0);
p_initial = subs(p,t,0);
v_initial = subs(v,t,0);

%% Convert to initial representations

syms X0 Y0 Z0

A=solve(strcat(char(X), ' == X'), ...
    strcat(char(Y), ' == Y'),...
    strcat(char(Z), ' == Z'),...
    'x','y','z');
x_char = char(A.x(1));
y_char = char(A.y(1));
z_char = char(A.z(1));
x_char = strrep(strrep(strrep(x_char, 'X', 'X0'), 'Y','Y0'), 'Z', 'Z0');
y_char = strrep(strrep(strrep(y_char, 'X', 'X0'), 'Y','Y0'), 'Z', 'Z0');
z_char = strrep(strrep(strrep(z_char, 'X', 'X0'), 'Y','Y0'), 'Z', 'Z0');

x_ = sym(x_char);
y_ = sym(y_char);
z_ = sym(z_char);

u_ = subs(u, [x,y,z], [x_,y_,z_]);
p_ = subs(p, [x,y,z], [x_,y_,z_]);
v_ = subs(v, [x,y,z], [x_,y_,z_]);

source_ = subs(source, [x,y,z], [x_,y_,z_]);
bf_ = subs(bf, [x,y,z], [x_,y_,z_]);

u_initial_ = subs(subs(u_initial, [x,y,z], [x_,y_,z_]),t,0);
p_initial_ = subs(subs(p_initial, [x,y,z], [x_,y_,z_]),t,0);
v_initial_ = subs(subs(v_initial, [x,y,z], [x_,y_,z_]),t,0);

gbars_ = subs(gbars, [x,y,z], [x_,y_,z_]);
tbars_ = subs(tbars, [x,y,z], [x_,y_,z_]);

%% Output to file
fileID = fopen('../exact_solutions.py','w');


%% u,p
fprintf(fileID, 'U = [''%s'',\n ''%s'',\n ''%s'',\n ''%s'']\n\n', ...
    char(u_(1)), char(u_(2)), char(u_(3)), char(p_));

%% v
fprintf(fileID, 'V = [''%s'',\n ''%s'',\n ''%s'']\n\n', ...
    char(v_(1)), char(v_(2)), char(v_(3)));

%% source
fprintf(fileID, 'source = ''%s''\n\n', ...
    char(source_));

%% body force
fprintf(fileID, 'bf = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(bf_(1)), char(bf_(2)), char(bf_(3)));

%% u_initial
fprintf(fileID, 'u_initial = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(u_initial_(1)), char(u_initial_(2)), char(u_initial_(3)) );

%% p_initial
fprintf(fileID, 'p_initial = ''%s''\n\n', char(p_initial_));

%% v_initial
fprintf(fileID, 'v_initial = [''%s'', ''%s'', ''%s'']\n\n', ...
    char(v_initial_(1)), char(v_initial_(2)), char(v_initial_(3)) );

%% gbars
fprintf(fileID,'gbars = \\\n[');
for i=1:length(gbars_)
    fprintf(fileID,'''%s''',char(gbars_(i)));
    if i<length(gbars_)
        fprintf(fileID,',\n');
    end
end
fprintf(fileID,']\n\n');

%% tbars
[n,m] = size(tbars_);
fprintf(fileID,'tbars = \\\n[');
for i=1:m
    fprintf(fileID,'[');
    for j=1:n
        fprintf(fileID,'''%s''',char(tbars_(j,i)));
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



