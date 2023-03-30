import numpy as np 

def CCM_get_C(Vd, Vo, D, Ts, L, Io, Delta_Vo):
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
        C = Io*D*Ts/Delta_Vo
    else:
        print("Capacitor caso triángulo")
        h = Ix + Delta_IL/2 - Io
        m = -(Vo-Vd)/L
        Delta_Q = -h**2/(2*m)
        C = Delta_Q/Delta_Vo
    return C

def DCM_get_C(Vd, Vo, D, Ts, L, Io, Delta_Vo):
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
    C = Delta_Q/Delta_Vo
    return C
