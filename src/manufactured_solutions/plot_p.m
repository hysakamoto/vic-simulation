% plot the minimum pressure value over time based on the manufactured solution
% the minimum pressure occurs anywhere Z=1.0

for i=1:100
    tt(i) = i/100*20;
    pp(i) = subs(p, [Z,t],[1,tt(i)]);
end
plot(tt,pp);
