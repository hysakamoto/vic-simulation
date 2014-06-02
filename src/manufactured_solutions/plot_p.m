% plot the minimum pressure value over time based on the manufactured solution
% the minimum pressure occurs anywhere Z=1.0

total_T = 10;

k=1;
tt = zeros(100*11,1);
zz = zeros(100*11,1);
pp = zeros(100*11,1);

for i=1:100
    for zv=0:0.1:1
        tt(k) = i/100*total_T;
        zz(k) = zv;
        pp(k) = subs(p, [Z0,t],[zv,tt(k)]);
    end
end

scatter(tt,zz,1000,pp, 'filled');
colorbar;
