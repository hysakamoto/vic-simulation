import ast, _ast
import re


#### Parsing Routine ####
ops = {_ast.Mult: "*", _ast.Add: "+", _ast.Sub: "-", 
       _ast.Div: "/", _ast.USub: "-", _ast.Pow: "pow", _ast.Not: "log"}

def rec(n):
    # pdb.set_trace()
    if isinstance(n, _ast.Expr):
        return rec(n.value)
    # function call (compatible with only one argument)
    elif isinstance(n, _ast.Call):
        if(len(n.args)==1):
            return "{}({})".format(n.func.id, rec(n.args[0]))
        else:
            print "not compatible with function with more than one argument."
            return None
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


#### Exact Solutions ####
execfile("exact_solutions.py")

for i in range(4):
    U[i] = U[i].translate(None, " ")
    U[i] = re.sub('\^', '**', U[i])
    U[i] = str(rec(ast.parse(U[i]).body[0]) )
    U[i] = re.sub('X0', 'x[0]', U[i])
    U[i] = re.sub('Y0', 'x[1]', U[i])
    U[i] = re.sub('Z0', 'x[2]', U[i])
    U[i] = U[i].translate(None, " ")

print "u_e = Expression((\"" + U[0] + "\"," + "\"" + U[1] + "\"," +\
    "\"" + U[2] + "\"" + "), t=dt)"

print "p_e = Expression(\"" + U[3] + "\", t=dt, k=k)"

## velocity 
for i in range(3):
    V[i] = V[i].translate(None, " ")
    V[i] = re.sub('\^', '**', V[i])
    V[i] = str(rec(ast.parse(V[i]).body[0]) )
    V[i] = re.sub('X0', 'x[0]', V[i])
    V[i] = re.sub('Y0', 'x[1]', V[i])
    V[i] = re.sub('Z0', 'x[2]', V[i])
    V[i] = V[i].translate(None, " ")

print "v_e = Expression((\"" + V[0] + "\"," + "\"" + V[1] + "\"," +\
    "\"" + V[2] + "\"" + "), t=dt)"

#### source ####
source = source.translate(None, " ")
source = re.sub('\^', '**', source)
source = str(rec(ast.parse(source).body[0]) )
source = re.sub('X0', 'x[0]', source)
source = re.sub('Y0', 'x[1]', source)
source = re.sub('Z0', 'x[2]', source)
source = source.translate(None, " ")

print "source = Expression(\"" + source + "\", t=dt)"

#### body force ####
for i in range(3):
    bf[i] = bf[i].translate(None, " ")
    bf[i] = re.sub('\^', '**', bf[i])
    bf[i] = str(rec(ast.parse(bf[i]).body[0]) )
    bf[i] = re.sub('X0', 'x[0]', bf[i])
    bf[i] = re.sub('Y0', 'x[1]', bf[i])
    bf[i] = re.sub('Z0', 'x[2]', bf[i])
    bf[i] = bf[i].translate(None, " ")

print "body_force = Expression((\"" + bf[0] + "\"," + "\"" + bf[1] + "\"," +\
    "\"" + bf[2] + "\"" + "), t=dt, lm=lm, mu=mu, k=k)"

#### normal flux bc ####
gbars_c = []
for i in range(len(gbars)):
    k = gbars[i]
    k = k.replace(" ", "")
    k = re.sub('\^', '**', k)
    k = str(rec(ast.parse(k).body[0]) )
    k = re.sub('X0', 'x[0]', k)
    k = re.sub('Y0', 'x[1]', k)
    k = re.sub('Z0', 'x[2]', k)
    k = k.translate(None, " ")
    gbars_c.append(k)
    print "gbar_" + directions[i] +"= Expression(\"" + k + "\", t=dt, k=k)"


#### traction bc ####
tbars_c = []
j=0
for tb in tbars:
    k = tb
    for i in range(3):
        k[i] = k[i].replace(" ", "")
        k[i] = k[i].replace('^', '**')
        k[i] = str(rec(ast.parse(k[i]).body[0]) )
        k[i] = re.sub('X0', 'x[0]', k[i])
        k[i] = re.sub('Y0', 'x[1]', k[i])
        k[i] = re.sub('Z0', 'x[2]', k[i])
        k[i] = k[i].translate(None, " ")
    tbars_c.append(k)
    print "tbar_" + directions[j] + "= Expression((\"" + k[0] + "\",\n\"" + k[1] + "\",\n\"" + k[2] + "\"), t=dt, mu=mu, lm=lm, k=k)"
    j+=1

#### initial conditions #### 

for i in range(3):
    u_initial[i] = u_initial[i].translate(None, " ")
    u_initial[i] = re.sub('\^', '**', u_initial[i])
    u_initial[i] = str(rec(ast.parse(u_initial[i]).body[0]) )
    u_initial[i] = re.sub('X0', 'x[0]', u_initial[i])
    u_initial[i] = re.sub('Y0', 'x[1]', u_initial[i])
    u_initial[i] = re.sub('Z0', 'x[2]', u_initial[i])
    u_initial[i] = u_initial[i].translate(None, " ")

print "u_initial = Expression((\"" + u_initial[0] + "\"," + "\"" + u_initial[1] + "\"," +\
    "\"" + u_initial[2] + "\"" + "), lm=lm, mu=mu, k=k)"

p_initial = p_initial.translate(None, " ")
p_initial = re.sub('\^', '**', p_initial)
p_initial = str(rec(ast.parse(p_initial).body[0]) )
p_initial = re.sub('X0', 'x[0]', p_initial)
p_initial = re.sub('Y0', 'x[1]', p_initial)
p_initial = re.sub('Z0', 'x[2]', p_initial)
p_initial = p_initial.translate(None, " ")

print "p_initial = Expression(\"" + p_initial + "\", t=dt, k=k)"
