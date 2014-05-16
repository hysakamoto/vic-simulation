import ast, _ast
import re


#### Parsing Routine ####
ops = {_ast.Mult: "*", _ast.Add: "+", _ast.Sub: "-", 
       _ast.Div: "/", _ast.USub: "-", _ast.Pow: "pow", _ast.Not: "log"}

def rec(n):
    if isinstance(n, _ast.Expr):
        return rec(n.value)
    # function call (compatible with only one argument)
    elif isinstance(n, _ast.Call):
        return "{}({})".format(n.func.id, rec(n.args[0]))
    elif isinstance(n, _ast.UnaryOp):
        return "({} {})".format(ops[type(n.op)], rec(n.operand))
    elif isinstance(n, _ast.BinOp):
        if isinstance(n.op, _ast.Pow):
            return "{}({}, {})".format(ops[type(n.op)], rec(n.left),rec(n.right))
        else: 
            return "({} {} {})".format(rec(n.left), ops[type(n.op)],rec(n.right))
    elif isinstance(n, _ast.Name):
        return n.id
    elif isinstance(n, _ast.Num):
        # make the number float (int may not work)
        return float(n.n)


directions = ["top", "bottom", "right", "left", "back", "front"]

source = """(X_*Y_*exp(t))/(X_*Y_*(exp(t) - 1) + 1)^2"""

bf = """ (lm*log(X_*Y_*exp(t) - X_*Y_ + 1)*(Y_ - Y_*exp(t)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - (lm*(Y_ - Y_*exp(t)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))
                                         (lm*log(X_*Y_*exp(t) - X_*Y_ + 1)*(X_ - X_*exp(t)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - (lm*(X_ - X_*exp(t)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - X_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))
 (mu*((2*X_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 + (2*Y_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2))/(X_*Y_*exp(t) - X_*Y_ + 1) - X_*Y_*(exp(t) - 1)"""

gbars_0 = """                           -X_*Y_*k*(exp(t) - 1)
                            X_*Y_*k*(exp(t) - 1)
 -Y_*k*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))
  Y_*k*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))
  X_*k*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))
 -X_*k*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))"""

tbars_0 = """[                                  (Y_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2,                                  (X_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2,   (lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1) + (mu*((X_*Y_*exp(t) - X_*Y_ + 1)^2 + (X_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1))^2)/(X_*Y_*exp(t) - X_*Y_ + 1)^2 + (Y_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1))^2)/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - 1))/(X_*Y_*exp(t) - X_*Y_ + 1) - X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))]
[                                 -(Y_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2,                                 -(X_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2, - (lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1) - (mu*((X_*Y_*exp(t) - X_*Y_ + 1)^2 + (X_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1))^2)/(X_*Y_*exp(t) - X_*Y_ + 1)^2 + (Y_^2*(exp(t) - 1)^2*(Z_ + X_*Y_*Z_*(exp(t) - 1))^2)/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - 1))/(X_*Y_*exp(t) - X_*Y_ + 1) - X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))]
[ (lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1) - X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)),                                                                                                               0,                                                                                                                                                                                                                                                                                 (Y_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2]
[                                                  -(lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1),                                                                                                               0,                                                                                                                                                                                                                             - X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)) - (Y_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2]
[                                                                                                               0, X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)) - (lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1),                                                                                                                                                                                                                                                                                -(X_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2]
[                                                                                                               0,                                                   (lm*log(X_*Y_*exp(t) - X_*Y_ + 1))/(X_*Y_*exp(t) - X_*Y_ + 1),                                                                                                                                                                                                                               (X_*mu*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1)))/(X_*Y_*exp(t) - X_*Y_ + 1)^2 - X_*Y_*(exp(t) - 1)*(Z_ + X_*Y_*Z_*(exp(t) - 1))]"""


source = source.replace(' ','')
bf = bf.replace(' ','').split('\n')
gbars_0 = gbars_0.replace(' ','').split('\n')
tbars_0 = tbars_0.replace(' ','').split('\n')
for i in range(len(tbars_0)):
    tbars_0[i] = tbars_0[i].replace('[','').replace(']','').split(',')



#### source ####
source = source.translate(None, " ")
source = re.sub('\^', '**', source)
source = str(rec(ast.parse(source).body[0]) )
source = re.sub('X\_', 'x[0]', source)
source = re.sub('Y\_', 'x[1]', source)
source = re.sub('Z\_', 'x[2]', source)
source = source.translate(None, " ")

print "source = Expression(\"" + source + ", t=dt)"


#### body force ####
# bf = "0, 0, ((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))".split(",")

for i in range(3):
    bf[i] = bf[i].translate(None, " ")
    bf[i] = re.sub('\^', '**', bf[i])
    bf[i] = str(rec(ast.parse(bf[i]).body[0]) )
    bf[i] = re.sub('X\_', 'x[0]', bf[i])
    bf[i] = re.sub('Y\_', 'x[1]', bf[i])
    bf[i] = re.sub('Z\_', 'x[2]', bf[i])
    bf[i] = bf[i].translate(None, " ")

print "body_force = Expression(\"" + bf[0] + "\"," + "\"" + bf[1] + "\"," +\
    "\"" + bf[2] + "\"" + ", t=dt)"

#### pressure bc ####
# p_top = "0.0"
# p_bottom   = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
# p_right = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
# p_left = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
# p_back = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
# p_front = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"

# pp = [p_top, p_bottom, p_right, p_left, p_back, p_front]
# pp_c = []
# for p in pp:
#     k = p
#     k = re.sub('\^', '**', k)
#     k = str(rec(ast.parse(k).body[0]) )
#     k = re.sub('X', 'x[0]', k)
#     k = re.sub('Y', 'x[1]', k)
#     k = re.sub('Z', 'x[2]', k)
#     k = k.translate(None, " ")
#     pp_c.append(p)
#     print "p_ = Expression(\"" + k + "\", t=dt)"

# u_top = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "1 - (t - 20)^2/400"]
# u_bottom = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "0"]
# u_right = ["(-1/((t - 20)^2/400 - 2))^(1/2) - 1", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "-Z*((t - 20)^2/400 - 1)"]
# u_left = ["0", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "-Z*((t - 20)^2/400 - 1)"]
# u_back = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "(-1/((t - 20)^2/400 - 2))^(1/2) - 1", "-Z*((t - 20)^2/400 - 1)"]
# u_front = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "0", "-Z*((t - 20)^2/400 - 1)"]

# uu = [u_top, u_bottom, u_right, u_left, u_back, u_front]
# uu_c = []
# for u in uu:
#     k = u
#     for i in range(3):
#         k[i] = re.sub('\^', '**', k[i])
#         k[i] = str(rec(ast.parse(k[i]).body[0]) )
#         k[i] = re.sub('X', 'x[0]', k[i])
#         k[i] = re.sub('Y', 'x[1]', k[i])
#         k[i] = re.sub('Z', 'x[2]', k[i])
#         k[i] = k[i].translate(None, " ")
#     uu_c.append(k)
#     print "u_ = Expression((\"" + k[0] + "\",\n\"" + k[1] + "\",\n\"" + k[2] + "\"), t=dt)"

#### normal flux bc ####
# gbars_0 = "400*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))*(1/(- t^2 + 40*t + 400)^2)^(1/2), -400*((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))*(1/(- t^2 + 40*t + 400)^2)^(1/2), 0, 0, 0, 0".split(",")

gbars_c = []
for i in range(len(gbars_0)):
    k = gbars_0[i]
    k = k.replace(" ", "")
    k = re.sub('\^', '**', k)
    k = str(rec(ast.parse(k).body[0]) )
    k = re.sub('X\_', 'x[0]', k)
    k = re.sub('Y\_', 'x[1]', k)
    k = re.sub('Z\_', 'x[2]', k)
    k = k.translate(None, " ")
    gbars_c.append(k)
    print "gbar_" + directions[i] +"= Expression(\"" + k + "\", t=dt)"


#### traction bc ####
# tbars_0 = ["0,0,400*((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2+mu*(((t-20)^2/400-2)^2-1))*(1/(-t^2+40*t+400)^2)^(1/2)".split(","),
#            "0,0,-400*((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2+mu*(((t-20)^2/400-2)^2-1))*(1/(-t^2+40*t+400)^2)^(1/2)".split(","),
#            "((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2-mu*(1/((t-20)^2/400-2)+1))*(((400/(-t^2+40*t+400))*(-t^2+40*t+400)^2)/160000)^(1/2),0,0".split(","),
#            "-((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2-mu*(1/((t-20)^2/400-2)+1))*(((400/(-t^2+40*t+400))*(-t^2+40*t+400)^2)/160000)^(1/2),0,0".split(","),
#            "0,((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2-mu*(1/((t-20)^2/400-2)+1))*(((400/(-t^2+40*t+400))*(-t^2+40*t+400)^2)/160000)^(1/2),0".split(","),
#            "0,-((((t/200-1/10)/((t-20)^2/400-2)^2-(t/200-1/10)/((-1/((t-20)^2/400-2))^(3/2)*((t-20)^2/400-2)^2))*(Z-Z*((t-20)^2/400-1))^2)/2-mu*(1/((t-20)^2/400-2)+1))*(((400/(-t^2+40*t+400))*(-t^2+40*t+400)^2)/160000)^(1/2),0".split(",")]

tbars_c = []
j=0
for tb in tbars_0:
    k = tb
    for i in range(3):
        k[i] = k[i].replace(" ", "")
        k[i] = re.sub('\^', '**', k[i])
        k[i] = str(rec(ast.parse(k[i]).body[0]) )
        k[i] = re.sub('X\_', 'x[0]', k[i])
        k[i] = re.sub('Y\_', 'x[1]', k[i])
        k[i] = re.sub('Z\_', 'x[2]', k[i])
        k[i] = k[i].translate(None, " ")
    tbars_c.append(k)
    print "tbar_" + directions[j] + "= Expression((\"" + k[0] + "\",\n\"" + k[1] + "\",\n\"" + k[2] + "\"), t=dt)"
    j+=1

#### Exact Solutions ####

## displacement and pressure
ux = ' X*(-1/((t - 20)^2/400 - 2))^(1/2) - X'
uy = ' Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y'
uz =  '-Z*((t - 20)^2/400 - 1)'
p = '-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2'

U = [ux,uy,uz,p]

for i in range(4):
    U[i] = U[i].translate(None, " ")
    U[i] = re.sub('\^', '**', U[i])
    U[i] = str(rec(ast.parse(U[i]).body[0]) )
    U[i] = re.sub('X\_', 'x[0]', U[i])
    U[i] = re.sub('Y\_', 'x[1]', U[i])
    U[i] = re.sub('Z\_', 'x[2]', U[i])
    U[i] = U[i].translate(None, " ")

print "u_e = Expression(\"" + U[0] + "\"," + "\"" + U[1] + "\"," +\
    "\"" + U[2] + "\"" + ", t=dt)"

print "p_e = Expression(\"" + U[3] + "\", d=dt)"

## velocity 

V = ['(X*(t/200 - 1/10))/(2*(-1/((t - 20)^2/400 - 2))^(1/2)*((t - 20)^2/400 - 2)^2)',
     '(Y*(t/200 - 1/10))/(2*(-1/((t - 20)^2/400 - 2))^(1/2)*((t - 20)^2/400 - 2)^2)',
     '-Z*(t/200 - 1/10)']

for i in range(3):
    V[i] = V[i].translate(None, " ")
    V[i] = re.sub('\^', '**', V[i])
    V[i] = str(rec(ast.parse(V[i]).body[0]) )
    V[i] = re.sub('X\_', 'x[0]', V[i])
    V[i] = re.sub('Y\_', 'x[1]', V[i])
    V[i] = re.sub('Z\_', 'x[2]', V[i])
    V[i] = V[i].translate(None, " ")

print "v_e = Expression((\"" + V[0] + "\"," + "\"" + V[1] + "\"," +\
    "\"" + V[2] + "\"" + "), t=dt)"
