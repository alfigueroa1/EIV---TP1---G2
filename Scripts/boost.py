import numpy as np

def CCM_IL_on(Vd, L, D, Ts, Io, t):
    Ix = Io/(1-D)
    IL = (Vd/L)*(t-D*Ts/2) + Ix
    return IL     
    
def CCM_IL_off(Vd, Vo, L, D, Ts, Io, t):
    Ix = Io/(1-D)
    IL = -((Vo-Vd)/L)*(t-(1+D)*Ts/2) + Ix
    return IL

def DCM_IL_on(Vd, L, D, Ts, t):
    Delta_Il = (Vd/L)*D*Ts
    IL = (Vd/L)*(t-D*Ts/2) + Delta_Il/2
    return IL

def DCM_IL_off(Vd, Vo, L, D, Ts, Io, t):
    Delta_1 = 2*L*Io / (Vd*D*Ts)
    Delta_Il = (Vd/L)*D*Ts
    if (t/Ts) % 1 < (D + Delta_1):
        IL = -((Vo-Vd)/L)*(t-D*Ts) + Delta_Il
    else:
        IL = 0
    return IL

def get_IL(Vd, Vo, L, Ts, Io, t):
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 < D:
            IL = CCM_IL_on(Vd, L, D, Ts, Io, (t%Ts))
        else:
            IL = CCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        if (t/Ts) % 1 < D:
            IL = DCM_IL_on(Vd, L, D, Ts, (t%Ts))
        else:
            IL = DCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    return IL
            

def CCM_get_Delta_Q(Vd, Vo, D, Ts, L, Io, Delta_Vo):
    """Calculate the capacitor for a continuous conduction mode boost converter

    Args:
        Vd (float)
        Vo (float)
        D (float)
        Ts (float)
        L (float)
        Io (float)
        Delta_Vo (float)

    Returns:
        float: C
    """
    Ix = Io/(1-D)
    Delta_IL = (Vd/L)*D*Ts
    if Ix - Delta_IL/2 - Io >= 0:
        print("Capacitor caso rectángulo")
        Delta_Q = Io*D*Ts
    else:
        print("Capacitor caso triángulo")
        h = Ix + Delta_IL/2 - Io
        m = -(Vo-Vd)/L
        Delta_Q = -h**2/(2*m)
        
    return Delta_Q

def DCM_get_Delta_Q(Vd, Vo, D, Ts, L, Io, Delta_Vo):
    """Calculate the capacitor for a diccontinuous conduction mode boost converter

    Args:
        Vd (float)
        Vo (float)
        D (float)
        Ts (float)
        L (float)
        Io (float)
        Delta_Vo (float)

    Returns:
        float: C
    """

    Delta_IL = (Vd/L)*D*Ts
    # En DCM, siempre vamos a estar en el caso triángulo
    h = Delta_IL - Io
    m = -(Vo-Vd)/L
    Delta_Q = -h**2/(2*m)
    
    return Delta_Q

def CCM_Pd(Vd, fsw, Io, t_ri, t_fv, t_rv, t_fi):
    Pd = 0.5*Vd*Io*(t_ri + t_fv + t_rv + t_fi)*fsw
    return Pd