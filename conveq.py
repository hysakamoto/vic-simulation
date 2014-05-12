import ast, _ast
import re

ops = {_ast.Mult: "*", _ast.Add: "+", _ast.Sub: "-", 
       _ast.Div: "/", _ast.USub: "-", _ast.Pow: "pow"}

def rec(n):
    if isinstance(n, _ast.Expr):
        return rec(n.value)
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

p_top = "0.0"
p_bottom   = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
p_right = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
p_left = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
p_back = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"
p_front = "-(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2"

pp = [p_top, p_bottom, p_right, p_left, p_back, p_front]
pp_c = []
for p in pp:
    k = p
    k = re.sub('\^', '**', k)
    k = str(rec(ast.parse(k).body[0]) )
    k = re.sub('X', 'x[0]', k)
    k = re.sub('Y', 'x[1]', k)
    k = re.sub('Z', 'x[2]', k)
    k = k.translate(None, " ")
    pp_c.append(p)
    print "p_ = Expression(\"" + k + "\", t=0.0)"

u_top = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "1 - (t - 20)^2/400"]
u_bottom = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "0"]
u_right = ["(-1/((t - 20)^2/400 - 2))^(1/2) - 1", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "-Z*((t - 20)^2/400 - 1)"]
u_left = ["0", "Y*(-1/((t - 20)^2/400 - 2))^(1/2) - Y", "-Z*((t - 20)^2/400 - 1)"]
u_back = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "(-1/((t - 20)^2/400 - 2))^(1/2) - 1", "-Z*((t - 20)^2/400 - 1)"]
u_front = ["X*(-1/((t - 20)^2/400 - 2))^(1/2) - X", "0", "-Z*((t - 20)^2/400 - 1)"]

uu = [u_top, u_bottom, u_right, u_left, u_back, u_front]
uu_c = []
for u in uu:
    k = u
    for i in range(3):
        k[i] = re.sub('\^', '**', k[i])
        k[i] = str(rec(ast.parse(k[i]).body[0]) )
        k[i] = re.sub('X', 'x[0]', k[i])
        k[i] = re.sub('Y', 'x[1]', k[i])
        k[i] = re.sub('Z', 'x[2]', k[i])
        k[i] = k[i].translate(None, " ")
    uu_c.append(k)
    print "u_ = Expression((\"" + k[0] + "\",\n\"" + k[1] + "\",\n\"" + k[2] + "\"), t=0.0)"


gbar_top = "((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))"
gbar_bottom = "-((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))"
gbar_right = "0.0"
gbar_left = "0.0"
gbar_back = "0.0"
gbar_front = "0.0"

gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front]
gbars_c = []
for gb in gbars:
    k = gb
    k = re.sub('\^', '**', k)
    k = str(rec(ast.parse(k).body[0]) )
    k = re.sub('X', 'x[0]', k)
    k = re.sub('Y', 'x[1]', k)
    k = re.sub('Z', 'x[2]', k)
    k = k.translate(None, " ")
    gbars_c.append(k)
    print "gbar_ = Expression(\"" + k + "\", t=0.0)"

#### traction bc ####
tbar_top ="0; 0; (((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2 + (1000*(- t^2/400 + t/10 + 1)^2)/29 - 1000/29";
tbar_bottom ="0; 0; (((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2 - (1000*(- t^2/400 + t/10 + 1)^2)/29 + 1000/29";
tbar_right = "(((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2 + (1000*(- 400*t^2 + 16000*t + 160000))/(29*(- t^2 + 40*t + 400)^2) - 1000/29; 0; 0";
tbar_left = "1000/29 - (1000*(- 400*t^2 + 16000*t + 160000))/(29*(- t^2 + 40*t + 400)^2); 0; (((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2";
tbar_back = "0; (((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2 + (1000*(- 400*t^2 + 16000*t + 160000))/(29*(- t^2 + 40*t + 400)^2) - 1000/29; 0";
tbar_front = "0; 1000/29 - (1000*(- 400*t^2 + 16000*t + 160000))/(29*(- t^2 + 40*t + 400)^2); (((t/200 - 1/10)/((t - 20)^2/400 - 2)^2 - (t/200 - 1/10)/((-1/((t - 20)^2/400 - 2))^(3/2)*((t - 20)^2/400 - 2)^2))*(Z - Z*((t - 20)^2/400 - 1))^2)/2";

tbar_top    = tbar_top.translate(None, " ").split(";")
tbar_bottom = tbar_bottom.translate(None, " ").split(";")
tbar_right  = tbar_right.translate(None, " ").split(";")
tbar_left   = tbar_left.translate(None, " ").split(";")
tbar_back   = tbar_back.translate(None, " ").split(";")
tbar_front  = tbar_front.translate(None, " ").split(";")

tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front]
tbars_c = []
for tb in tbars:
    k = tb
    for i in range(3):
        k[i] = re.sub('\^', '**', k[i])
        k[i] = str(rec(ast.parse(k[i]).body[0]) )
        k[i] = re.sub('X', 'x[0]', k[i])
        k[i] = re.sub('Y', 'x[1]', k[i])
        k[i] = re.sub('Z', 'x[2]', k[i])
        k[i] = k[i].translate(None, " ")
    tbars_c.append(k)
    print "tbar_ = Expression((\"" + k[0] + "\",\n\"" + k[1] + "\",\n\"" + k[2] + "\"), t=0.0)"

