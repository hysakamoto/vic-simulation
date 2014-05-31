U = ['0',
 '0',
 'Z0 - (Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))/(X0*Y0*(exp(t/20) - 1) + 1) - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)',
 'X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))']

V = ['0',
 '0',
 '(X0*Y0*exp(t/20)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(20*(X0*Y0*(exp(t/20) - 1) + 1)^2)']

source = '(X0*Y0*exp(t/20))/(20*(X0*Y0*(exp(t/20) - 1) + 1)^2)'

bf = ['(lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Y0 - Y0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - (lm*(Y0 - Y0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)) + (Y0*mu*(exp(t/20) - 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2', '(lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1)*(X0 - X0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - (lm*(X0 - X0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - X0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)) + (X0*mu*(exp(t/20) - 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2', '(mu*((2*X0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 + (2*Y0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10) + (2*X0*mu*(exp(t/20) - 1)*(X0 - X0*exp(t/20))*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^3 + (2*Y0*mu*(exp(t/20) - 1)*(Y0 - Y0*exp(t/20))*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)^3']

u_initial = ['0', '0', '0']

p_initial = '0'

v_initial = ['0', '0', '(X0*Y0*Z0)/20']

gbars = \
['(X0^2*k*sin((pi*t)/10)*(exp(t/20) - 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*(exp(t/20) - 1) + 1)^2 - (X0*Y0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*(exp(t/20) - 1) + 1) + (Y0^2*k*sin((pi*t)/10)*(exp(t/20) - 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*(exp(t/20) - 1) + 1)^2',
'(X0*Y0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*(exp(t/20) - 1) + 1) - (X0^2*k*sin((pi*t)/10)*(exp(t/20) - 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*(exp(t/20) - 1) + 1)^2 - (Y0^2*k*sin((pi*t)/10)*(exp(t/20) - 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*(exp(t/20) - 1) + 1)^2',
'-Y0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))',
'Y0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))',
'-X0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))',
'X0*k*sin((pi*t)/10)*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))']

tbars = \
[['(Y0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/((X0*Y0*(exp(t/20) - 1) + 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)) - (Y0*(exp(t/20) - 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*(exp(t/20) - 1) + 1)^2',
'(X0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/((X0*Y0*(exp(t/20) - 1) + 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1)) - (X0*(exp(t/20) - 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*(exp(t/20) - 1) + 1)^2',
'((X0*Y0*exp(t/20) - X0*Y0 + 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) + (mu*((X0*Y0*exp(t/20) - X0*Y0 + 1)^2 + (X0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 + (Y0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))))/(X0*Y0*(exp(t/20) - 1) + 1) - (X0^2*mu*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/((X0*Y0*(exp(t/20) - 1) + 1)^2*(X0*Y0*exp(t/20) - X0*Y0 + 1)) - (Y0^2*mu*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/((X0*Y0*(exp(t/20) - 1) + 1)^2*(X0*Y0*exp(t/20) - X0*Y0 + 1))'], 
['(Y0*(exp(t/20) - 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*(exp(t/20) - 1) + 1)^2 - (Y0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/((X0*Y0*(exp(t/20) - 1) + 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1))',
'(X0*(exp(t/20) - 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*(exp(t/20) - 1) + 1)^2 - (X0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/((X0*Y0*(exp(t/20) - 1) + 1)*(X0*Y0*exp(t/20) - X0*Y0 + 1))',
'(X0^2*mu*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/((X0*Y0*(exp(t/20) - 1) + 1)^2*(X0*Y0*exp(t/20) - X0*Y0 + 1)) - ((X0*Y0*exp(t/20) - X0*Y0 + 1)*((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) + (mu*((X0*Y0*exp(t/20) - X0*Y0 + 1)^2 + (X0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 + (Y0^2*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/(X0*Y0*exp(t/20) - X0*Y0 + 1)^2 - 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))))/(X0*Y0*(exp(t/20) - 1) + 1) + (Y0^2*mu*(exp(t/20) - 1)^2*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20))^2)/((X0*Y0*(exp(t/20) - 1) + 1)^2*(X0*Y0*exp(t/20) - X0*Y0 + 1))'], 
['((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)',
'0',
'(Y0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)'], 
['-((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)',
'0',
'-(Y0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)'], 
['0',
'((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)',
'(X0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)'], 
['0',
'-((lm*log(X0*Y0*exp(t/20) - X0*Y0 + 1))/(X0*Y0*exp(t/20) - X0*Y0 + 1) - X0*Y0*sin((pi*t)/10)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))*(X0*Y0*exp(t/20) - X0*Y0 + 1)',
'-(X0*mu*(exp(t/20) - 1)*(Z0 - X0*Y0*Z0 + X0*Y0*Z0*exp(t/20)))/(X0*Y0*exp(t/20) - X0*Y0 + 1)']]

