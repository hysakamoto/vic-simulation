""" Simulation of VIC"""

# Simulation of hyperelastic solid in initial configuration

# Begin simulation
import pdb

from dolfin import *
import numpy as np
import newton_solve
from vic_bcs import *

# DEBUG
# set_log_level(DEBUG) #PROGRESS

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["num_threads"] = 2
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}


def manufactured_solutiond(dt, k, mu, lm):



    u_e = Expression(("0.0","0.0","(((x[2]-(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))"), t=dt)
    p_e = Expression("0.0", t=dt, k=k)
    v_e = Expression(("0.0","0.0","((((x[0]*x[1])*exp((t/20.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/(20.0*pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))"), t=dt)
    source = Expression("((((x[0]*x[1])*exp((t/20.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/(20.0*pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))", t=dt)
    body_force = Expression(("((-(((lm*(x[1]-(x[1]*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[1]-(x[1]*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))","((-(((lm*(x[0]-(x[0]*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[0]-(x[0]*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))","(mu*(((((2.0*pow(x[0],2.0))*pow((exp((t/20.0))-1.0),2.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))+((((2.0*pow(x[1],2.0))*pow((exp((t/20.0))-1.0),2.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))))"), t=dt, lm=lm, mu=mu, k=k)
    gbar_top= Expression("0.0", t=dt, k=k)
    gbar_bottom= Expression("0.0", t=dt, k=k)
    gbar_right= Expression("0.0", t=dt, k=k)
    gbar_left= Expression("0.0", t=dt, k=k)
    gbar_back= Expression("0.0", t=dt, k=k)
    gbar_front= Expression("0.0", t=dt, k=k)
    tbar_top= Expression(("(((((x[1]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[1]*lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))",
    "(((((x[0]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[0]*lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))",
    "(((((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((mu*(((pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)+(((pow(x[0],2.0)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))+(((pow(x[1],2.0)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))-1.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-((((pow(x[0],2.0)*mu)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/(pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))-((((pow(x[1],2.0)*mu)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/(pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_bottom= Expression(("((((((x[1]*lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))-((((x[1]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))",
    "((((((x[0]*lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))-((((x[0]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))",
    "((((((pow(x[0],2.0)*mu)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/(pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((mu*(((pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)+(((pow(x[0],2.0)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))+(((pow(x[1],2.0)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))-1.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))+((((pow(x[1],2.0)*mu)*pow((exp((t/20.0))-1.0),2.0))*pow(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))),2.0))/(pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_right= Expression(("(lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))",
    "0.0",
    "((((x[1]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_left= Expression(("((-lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))",
    "0.0",
    "((-(((x[1]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_back= Expression(("0.0",
    "(lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))",
    "((((x[0]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_front= Expression(("0.0",
    "((-lm)*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))",
    "((-(((x[0]*mu)*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))"), t=dt, mu=mu, lm=lm, k=k)
    u_initial = Expression(("0.0","0.0","0.0"), lm=lm, mu=mu, k=k)
    p_initial = Expression("0.0", t=dt, k=k)




    # u_e = Expression(("((x[0]*pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(1.0/2.0)))-x[0])","((x[1]*pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(1.0/2.0)))-x[1])","((-x[2])-(x[2]*((pow((t-20.0),2.0)/400.0)-2.0)))"), t=dt)
    # p_e = Expression("((-((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))/(2.0*k))", t=dt, k=k)
    # v_e = Expression(("((-(x[0]*((t/200.0)-(1.0/10.0))))/(2.0*((pow((t-20.0),2.0)/400.0)-2.0)))","((-(x[1]*((t/200.0)-(1.0/10.0))))/(2.0*((pow((t-20.0),2.0)/400.0)-2.0)))","((x[2]*((t/200.0)-(1.0/10.0)))/((pow((t-20.0),2.0)/400.0)-2.0))"), t=dt)
    # source = Expression("0.0", t=dt)
    # body_force = Expression(("0.0","0.0","((-((x[2]*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*((pow((t-20.0),2.0)/400.0)-2.0)))/k)"), t=dt, lm=lm, mu=mu, k=k)
    # gbar_top= Expression("((((-x[2])*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow((1.0/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)),(1.0/2.0)))*((pow((t-20.0),2.0)/400.0)-2.0))", t=dt, k=k)
    # gbar_bottom= Expression("(((x[2]*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow((1.0/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)),(1.0/2.0)))*((pow((t-20.0),2.0)/400.0)-2.0))", t=dt, k=k)
    # gbar_right= Expression("0.0", t=dt, k=k)
    # gbar_left= Expression("0.0", t=dt, k=k)
    # gbar_back= Expression("0.0", t=dt, k=k)
    # gbar_front= Expression("0.0", t=dt, k=k)
    # tbar_top= Expression(("0.0",
    # "0.0",
    # "(pow((1.0/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)),(1.0/2.0))*((mu*(pow((((t/10.0)-(pow(t,2.0)/400.0))+1.0),2.0)-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k))))"), t=dt, mu=mu, lm=lm, k=k)
    # tbar_bottom= Expression(("0.0",
    # "0.0",
    # "((-pow((1.0/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)),(1.0/2.0)))*((mu*(pow((((t/10.0)-(pow(t,2.0)/400.0))+1.0),2.0)-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k))))"), t=dt, mu=mu, lm=lm, k=k)
    # tbar_right= Expression(("(((mu*(((160000.0*(((t/10.0)-(pow(t,2.0)/400.0))+1.0))/pow((((40.0*t)-pow(t,2.0))+400.0),2.0))-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k)))*pow((2.0-(pow((t-20.0),2.0)/400.0)),(1.0/2.0)))",
    # "0.0",
    # "0.0"), t=dt, mu=mu, lm=lm, k=k)
    # tbar_left= Expression(("((-((mu*(((160000.0*(((t/10.0)-(pow(t,2.0)/400.0))+1.0))/pow((((40.0*t)-pow(t,2.0))+400.0),2.0))-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k))))*pow((2.0-(pow((t-20.0),2.0)/400.0)),(1.0/2.0)))",
    # "0.0",
    # "0.0"), t=dt, mu=mu, lm=lm, k=k)
    # tbar_back= Expression(("0.0",
    # "(((mu*(((160000.0*(((t/10.0)-(pow(t,2.0)/400.0))+1.0))/pow((((40.0*t)-pow(t,2.0))+400.0),2.0))-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k)))*pow((2.0-(pow((t-20.0),2.0)/400.0)),(1.0/2.0)))",
    # "0.0"), t=dt, mu=mu, lm=lm, k=k)
    # tbar_front= Expression(("0.0",
    # "((-((mu*(((160000.0*(((t/10.0)-(pow(t,2.0)/400.0))+1.0))/pow((((40.0*t)-pow(t,2.0))+400.0),2.0))-1.0))+(((pow(x[2],2.0)*((((t/200.0)-(1.0/10.0))/pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))-(((t/200.0)-(1.0/10.0))/(pow((-1.0/((pow((t-20.0),2.0)/400.0)-2.0)),(3.0/2.0))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0)))))*pow(((pow((t-20.0),2.0)/400.0)-2.0),2.0))/(2.0*k))))*pow((2.0-(pow((t-20.0),2.0)/400.0)),(1.0/2.0)))",
    # "0.0"), t=dt, mu=mu, lm=lm, k=k)


    tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front]
    gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front]

    return [u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial]


def vic_sim( m_num, p_order, dt, T_total, max_it, \
             omega, Ee, nu, gamma, tau, perm, \
             top_trac, body_force ):
    """ 
    run the VIC simulation
    """

    ## avoid recompilation
    dt_const = Constant(dt)
    
    # Elasticity parameters
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

    # Exact solutions
    [u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial] \
        = manufactured_solutiond(dt, perm, mu, lmbda)

    
    ## Create mesh and define function space
    mesh = UnitCubeMesh(m_num, m_num, m_num)
    dim  = mesh.topology().dim() 

    ##  Function Spaces
    Pu = VectorFunctionSpace(mesh, "Lagrange", p_order)  # space for displacements
    Pp = FunctionSpace(mesh, "Lagrange", p_order)        # space for pressure
    V  = MixedFunctionSpace([Pu,Pp])                    # mixed space

    ## Functions
    dup  = TrialFunction(V)    # Incremental displacement-pressure
    wq   = TestFunction(V)     # Test function
    up   = Function(V)         # Displacement-pressure from previous iteration
    u, p = split(up)           # Function in each subspace to write the functional
    w, q = split(wq)           # Test Function split

    ## Boundary Conditions

    # Create mesh function over cell facets
    exterior_facet_domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

    # Define Neumann boundaries
    bd_tol = 1E-14   # tolerance for coordinate comparisons
    ds_neumann = neumann_boundaries(bd_tol, exterior_facet_domains)

    # Get the neumann boundary conditions
    gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front \
        = gbars
    tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front \
        = tbars

    # Define Dirichlet boundaries
    bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
        = dirichlet_boundaries(bd_tol, V, dt, u_e, p_e)

    # bcs = [bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
    #        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]

    bcs = [bc_ubottom,\
           bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]


    # bcs = [bc_ubottom, bc_uleft, bc_ufront, bc_pbottom, bc_pleft, bc_pfront]

    ## Initial conditions
    up_1   = Function(V)         # Displacement-pressure from previous iteration
    u_1, p_1 = split(up_1)       # Function in each subspace to write the functional
    assign (up_1.sub(0), interpolate(u_initial,Pu))
    assign (up_1.sub(1), interpolate(p_initial,Pp))

    ## Kinematics
    I    = Identity(dim)           # Identity tensor
    F    = I + grad(u)             # Deformation gradient
    F    = variable(F)             # Make F a variable for tensor differentiations
    C    = F.T*F                   # Right Cauchy-Green tensor
    invF = inv(F)

    # virtual Kinematics
    ddotF = grad(w)
    ddotE = 1.0/2.0 * (ddotF.T*F + F.T*ddotF)

    # Invariants of deformation tensors
    J    = det(F)
    Ic   = tr(C)
    IIIc = det(C)

    # Elasticity parameters
    mu, lmbda = Constant(Ee/(2*(1 + nu))), Constant(Ee*nu/((1 + nu)*(1 - 2*nu)))

    # Permeability
    K_perm = Constant(np.eye(3)*perm) # avoid recompilation

    ## Potential Energy
    # Strain energy density (compressible neo-Hookean model)
    psi = (mu/2.0)*(Ic - 3) - mu*ln(J) + (lmbda/2.0)*(ln(J))**2

    # Strain energy density (nearly incompressible neo-hookean model)
    kappa = mu*10
    C_hat = (IIIc)**(-1.0/3.0)*C
    # psi   = mu/2.0*(tr(C_hat)-3) + 1.0/2.0*kappa*ln(J)**2.0

    # PK1 stress tensor
    P = diff(psi,F)
    # PK2 stress tensor
    S = inv(F)*P

    ## time steps
    t      = 0.0
    tn     = 0         # time step number

    ## Viscoelasticity!!!!

    ## H-value
    HS = TensorFunctionSpace(mesh, "Lagrange", p_order )
    H_1 = Function(HS)
    H_1.interpolate(Constant(np.zeros([3,3])))
    
    ## viscoelastic parameters
    gamma_const = Constant(gamma)
    tau_const = Constant(tau)

    ## Kinematics
    F_1    = I + grad(u_1)             # Deformation gradient
    F_1    = variable(F_1)             # Make F a variable for tensor differentiations
    C_1    = F_1.T*F_1                   # Right Cauchy-Green tensor
    invF_1 = inv(F_1)
    # Invariants of deformation tensors
    J_1    = det(F_1)
    Ic_1   = tr(C_1)
    IIIc_1 = det(C_1)
    # Strain energy density (nearly incompressible neo-hookean model)
    C_hat_1 = (IIIc_1)**(-1.0/3.0)*C_1
    # psi_1   = mu/2.0*(tr(C_hat_1)-3) + 1.0/2.0*kappa*ln(J_1)**2.0
    # compressible neo-Hookean
    psi_1 = (mu/2.0)*(Ic_1 - 3.0) - mu*ln(J_1) + (lmbda/2.0)*(ln(J_1))**2.0
    # PK1 stress tensor
    P_1 = diff(psi_1,F_1)
    # PK2 stress tensor
    S_1 = invF_1*P_1

    H = exp(-dt_const/tau_const)*H_1 + (1-exp(-dt_const/tau_const))*(S-S_1)/(dt_const/tau_const)
    Sc = S+gamma_const*H

    ## Poroelasticity!!!!
    omega_const = Constant(omega)
    v_1 = Function(Pu)
    v_1.interpolate(Constant((0.0, 0.0, 0.0)))
    v = ((u-u_1)/dt_const - (1-omega_const)*v_1)/omega_const

    # Compute residual
    R = (inner(Sc, ddotE) + dot(J*(K_perm*invF.T*grad(p)), invF.T*grad(q)))*dx \
        - p*J*inner(ddotF, invF.T)*dx                                         \
        + (q*J*inner(grad(v),invF.T))*dx                                      \
        - (source*q)*dx \
        + (inner(body_force,w))*dx \
        - (inner(tbar_top,w))*ds_neumann(0) \
        - (inner(tbar_right,w))*ds_neumann(2) \
        - (inner(tbar_left,w))*ds_neumann(3) \
        - (inner(tbar_back,w))*ds_neumann(4) \
        - (inner(tbar_front,w))*ds_neumann(5) \
        # + (inner(gbar_top,q))*ds_neumann(0)\
        # + (inner(gbar_right,q))*ds_neumann(2) \
        # + (inner(gbar_back,q))*ds_neumann(4) \

    # Compute Jacobian of R
    Jac = derivative(R, up, dup)

    # Set up the problem
    problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
    solver = NonlinearVariationalSolver(problem)

    solver.parameters["nonlinear_solver"] = "newton"
    solver.parameters["newton_solver"]["linear_solver"] = "bicgstab"
    solver.parameters["newton_solver"]["preconditioner"] = "ilu"

    # solver.parameters["nonlinear_solver"] = "snes"
    # solver.parameters["snes_solver"]["line_search"] = "bt"
    # solver.parameters["snes_solver"]["linear_solver"] = "gmres"
    # solver.parameters["snes_solver"]["preconditioner"] = "ilu"        
    # solver.parameters["snes_solver"]["method"] = "tr"

    ## Save initial conditions in VTK format
    assign(up, up_1)
    dfile = File("results/displacement.pvd");
    pfile = File("results/pressure.pvd");
    dfile << (up.sub(0),0.0);
    pfile << (up.sub(1),0.0);

    # dup   = Function(V)

    ### Run Simulation
    u_max = []
    Eus = []
    Eps = []
    while tn<max_it:
        print 'time = ', (t+dt)

        # solve
        solver.solve()

        ### Error against exact solutions
        error_u = (u-u_e)**2*dx
        as_tmp = assemble(error_u)
        if as_tmp < 0.0:
            as_tmp = 0.0
        Eu = sqrt(as_tmp)

        error_p = (p-p_e)**2*dx
        as_tmp = assemble(error_p)
        if as_tmp < 0.0:
            as_tmp = 0.0
        Ep = sqrt(as_tmp)

        # Explicit interpolation of u_e onto the same space as u:
        # u_e_V = interpolate(u_e, V)
        # error_u = (u - u_e_V)**2*dx
        # E2 = sqrt(assemble(error_u))
        # Eu_4 = errornorm(u_e, u, normtype='l2', degree=3)

        Eus.append(Eu)
        Eps.append(Ep)

        # pdb.set_trace()
        u_tent, p_tent = up.split(deepcopy=True) 

        u_max.append( np.max(u_tent.vector().array()))

        # update
        t    += dt
        tn   += 1

        # update body force
        body_force.t = t+dt
        # update source term
        source.t = t+dt
        # update boundary conditions
        gbar_top.t = t+dt
        gbar_bottom.t = t+dt
        gbar_right.t = t+dt
        gbar_left.t = t+dt
        gbar_back.t = t+dt
        gbar_front.t = t+dt
        tbar_top.t = t+dt
        tbar_bottom.t = t+dt
        tbar_right.t = t+dt
        tbar_left.t = t+dt
        tbar_back.t = t+dt
        tbar_front.t = t+dt 
        # update exact solutions
        u_e.t = t+dt
        p_e.t = t+dt

        # Define Dirichlet boundaries
        bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
            bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront \
            = dirichlet_boundaries(bd_tol, V, t+dt, u_e, p_e)


        # bcs = [bc_utop, bc_ubottom, bc_uright, bc_uleft, bc_uback, bc_ufront, \
        #        bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]

        bcs = [bc_ubottom,\
               bc_ptop, bc_pbottom, bc_pright, bc_pleft, bc_pback, bc_pfront]

        # define problem
        problem = NonlinearVariationalProblem(R, up, bcs=bcs, J=Jac)
        solver = NonlinearVariationalSolver(problem)

        v_1  = project(v,Pu)
        up_1.assign(up)
        assign(up_1.sub(0),up.sub(0))
        H_1 = project(H, HS)
        S_1 = project(S, HS)

        # Save solution in VTK format
        dfile << (up.sub(0), t);
        pfile << (up.sub(1), t);

        # plot(p_1, title = "pressure", axes=True, interactive = True)


        

    return u_max, Eus, Eps

    # Plot and hold solution
    # plot(u, mode = "displacement", title="displacement", axes=True, interactive = True)
    # plot(p, title  = "pressure", axes=True, interactive = True)
