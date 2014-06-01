from dolfin import *

def getManuSolutions(dt, k, mu, lm):

    u_e = Expression(("0.0","0.0","(((x[2]-(((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0))))"), t=dt)

    p_e = Expression("((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)", t=dt, k=k)

    v_e = Expression(("0.0","0.0","((((((x[0]*x[1])*x[2])*exp((t/20.0)))/20.0)-((((x[0]*x[1])*x[2])*exp((t/20.0)))/(20.0*(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))))+((((x[0]*x[1])*exp((t/20.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/(20.0*pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))))"), t=dt)

    source = Expression("((((((((((((x[0]*x[1])*exp((t/20.0)))/20.0)-((k*((((x[1]*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*(x[1]-(x[1]*exp((t/20.0)))))/100.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))))-((k*((((x[0]*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*(x[0]-(x[0]*exp((t/20.0)))))/100.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))))-(((x[0]*x[1])*exp((t/20.0)))/(20.0*(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))))+(((((x[0]*k)*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))+(((((x[1]*k)*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*(((((2.0*k)*(x[0]-(x[0]*exp((t/20.0)))))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((((2.0*k)*(x[1]-(x[1]*exp((t/20.0)))))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))+((((x[0]*x[1])*exp((t/20.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/(20.0*pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))))+((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*(x[0]-(x[0]*exp((t/20.0)))))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/100.0))+((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*(x[1]-(x[1]*exp((t/20.0)))))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/100.0))", t=dt, k=k)

    body_force = Expression(("((((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[1]-(x[1]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((mu*(x[1]-(x[1]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-((lm*(x[1]-(x[1]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-((((x[1]*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))","((((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[0]-(x[0]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((mu*(x[0]-(x[0]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-((lm*(x[0]-(x[0]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-((((x[0]*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))","(((((((((((((((((lm*(x[1]-(x[1]*exp((t/20.0)))))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))+(((mu*(x[1]-(x[1]*exp((t/20.0)))))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))-((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[1]-(x[1]*exp((t/20.0)))))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))-((2.0*(((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*(((((x[0]-(x[0]*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-x[0])+(x[0]*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-((2.0*(((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*(((((x[1]-(x[1]*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-x[1])+(x[1]*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))+((((mu*(((((((2.0*pow(x[0],2.0))*x[2])+((2.0*pow(x[1],2.0))*x[2]))+(((2.0*pow(x[0],2.0))*x[2])*exp((t/10.0))))-(((4.0*pow(x[0],2.0))*x[2])*exp((t/20.0))))+(((2.0*pow(x[1],2.0))*x[2])*exp((t/10.0))))-(((4.0*pow(x[1],2.0))*x[2])*exp((t/20.0)))))/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(((((((2.0*pow(x[0],2.0))*x[2])+((2.0*pow(x[1],2.0))*x[2]))+(((2.0*pow(x[0],2.0))*x[2])*exp((t/10.0))))-(((4.0*pow(x[0],2.0))*x[2])*exp((t/20.0))))+(((2.0*pow(x[1],2.0))*x[2])*exp((t/10.0))))-(((4.0*pow(x[1],2.0))*x[2])*exp((t/20.0)))))/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0)))-((((mu*(x[0]-(x[0]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[0]-(x[0]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-((((mu*(x[1]-(x[1]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[1]-(x[1]*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))+((((((lm*(x[0]-(x[0]*exp((t/20.0)))))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0))+(((mu*(x[0]-(x[0]*exp((t/20.0)))))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))-((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(x[0]-(x[0]*exp((t/20.0)))))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/pow(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0),2.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0)))-((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(((((2.0*pow(x[0],2.0))*pow((exp((t/20.0))-1.0),2.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),3.0))+((((2.0*x[0])*(exp((t/20.0))-1.0))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(((((2.0*pow(x[1],2.0))*pow((exp((t/20.0))-1.0),2.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),3.0))+((((2.0*x[1])*(exp((t/20.0))-1.0))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-(((lm*(x[0]-(x[0]*exp((t/20.0)))))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((lm*(x[1]-(x[1]*exp((t/20.0)))))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))"), t=dt, lm=lm, mu=mu, k=k)

    gbar_top= Expression("(((((-k)*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((((x[0]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/100.0)))-((k*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((((x[1]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/100.0))))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*(((k/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((k*pow(((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))),2.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))+((k*pow(((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))),2.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)
    gbar_bottom= Expression("((((k*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((((x[0]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/100.0)))+((k*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((((x[1]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/100.0))))+(((((x[0]*x[1])*sin(((pi*t)/10.0)))*(((k/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((k*pow(((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))),2.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))+((k*pow(((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))),2.0))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)
    gbar_right= Expression("((((-k)*((((x[1]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/100.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)
    gbar_left= Expression("(((k*((((x[1]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/100.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)
    gbar_back= Expression("((((-k)*((((x[0]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/100.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)
    gbar_front= Expression("(((k*((((x[0]*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/100.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))+((((((x[0]*x[1])*k)*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))/100.0))", t=dt, k=k)

    tbar_top= Expression(("(((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"(((((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"((((((mu*(((((((((pow(x[0],2.0)*pow(x[2],2.0))+(pow(x[1],2.0)*pow(x[2],2.0)))+((pow(x[0],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[2],2.0))*exp((t/20.0))))+((pow(x[1],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[1],2.0))*pow(x[2],2.0))*exp((t/20.0))))+1.0)/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0))-1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(((((((pow(x[0],2.0)*pow(x[2],2.0))+(pow(x[1],2.0)*pow(x[2],2.0)))+((pow(x[0],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[2],2.0))*exp((t/20.0))))+((pow(x[1],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[1],2.0))*pow(x[2],2.0))*exp((t/20.0))))+1.0))/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))-((((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-((((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))-((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_bottom= Expression(("((((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))+(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"((((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))+(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"(((((((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0))))-(((mu*(((((((((pow(x[0],2.0)*pow(x[2],2.0))+(pow(x[1],2.0)*pow(x[2],2.0)))+((pow(x[0],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[2],2.0))*exp((t/20.0))))+((pow(x[1],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[1],2.0))*pow(x[2],2.0))*exp((t/20.0))))+1.0)/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0))-1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*(((((((pow(x[0],2.0)*pow(x[2],2.0))+(pow(x[1],2.0)*pow(x[2],2.0)))+((pow(x[0],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[2],2.0))*exp((t/20.0))))+((pow(x[1],2.0)*pow(x[2],2.0))*exp((t/10.0))))-(((2.0*pow(x[1],2.0))*pow(x[2],2.0))*exp((t/20.0))))+1.0))/((((((pow(x[0],2.0)*pow(x[1],2.0))-((2.0*x[0])*x[1]))+((pow(x[0],2.0)*pow(x[1],2.0))*exp((t/10.0))))-(((2.0*pow(x[0],2.0))*pow(x[1],2.0))*exp((t/20.0))))+(((2.0*x[0])*x[1])*exp((t/20.0))))+1.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0)))+((((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))+((((x[0]*x[1])*sin(((pi*t)/10.0)))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_right= Expression(("((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"0.0",
"(((((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))+((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_left= Expression(("((((((x[0]*x[1])*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-(lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))",
"0.0",
"(((-(((mu*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))-((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[1]*x[2])-((x[1]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[1]*x[2]))+((x[1]*x[2])*exp((t/20.0))))+(((x[1]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_back= Expression(("0.0",
"((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))-(((((x[0]*x[1])*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0))",
"(((((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))+((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))"), t=dt, mu=mu, lm=lm, k=k)
    tbar_front= Expression(("0.0",
"((((((x[0]*x[1])*sin(((pi*t)/10.0)))*((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/100.0)-(lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))",
"(((-(((mu*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))-(((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0)))))/((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0))))*((((x[0]*x[1])+(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0)))-((x[0]*x[1])*exp((t/20.0))))-2.0))-((lm*log(((((x[0]*x[1])*exp((t/20.0)))-(x[0]*x[1]))+1.0)))*((((((x[0]*x[2])-((x[0]*x[2])*exp((t/20.0))))/(((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0))-(x[0]*x[2]))+((x[0]*x[2])*exp((t/20.0))))+(((x[0]*(exp((t/20.0))-1.0))*((x[2]-((x[0]*x[1])*x[2]))+(((x[0]*x[1])*x[2])*exp((t/20.0)))))/pow((((x[0]*x[1])*(exp((t/20.0))-1.0))+1.0),2.0)))))"), t=dt, mu=mu, lm=lm, k=k)

    u_initial = Expression(("0.0","0.0","0.0"), lm=lm, mu=mu, k=k)

    p_initial = Expression("0.0", t=dt, k=k)

    v_initial = Expression(("0.0","0.0","(((x[0]*x[1])*x[2])/20.0)"), lm=lm, mu=mu, k=k)

    tbars = [tbar_top, tbar_bottom, tbar_right, tbar_left, tbar_back, tbar_front]

    gbars = [gbar_top, gbar_bottom, gbar_right, gbar_left, gbar_back, gbar_front]

    return [u_e, p_e, v_e, source, body_force, tbars, gbars, u_initial, p_initial, v_initial]

