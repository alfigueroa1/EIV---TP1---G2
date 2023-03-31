import numpy as np 

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