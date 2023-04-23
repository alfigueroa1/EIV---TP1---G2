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
    if (t/Ts) % 1 <= D + Delta_1:
        IL = -((Vo-Vd)/L)*(t-D*Ts) + Delta_Il
    else:
        IL = 0
    return IL

def get_VL(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the inductor voltage in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: VL
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 <= D:
            VL = Vd
        else:
            VL = -(Vo - Vd)
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        Delta_1 = 2*L*Io / (Vd*D*Ts)
        if (t/Ts) % 1 <= D:
            VL = Vd
        elif (t/Ts) % 1 < D + Delta_1:
            VL = -(Vo - Vd)
        else:
            VL = 0
    return VL

def get_IL(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the inductor current in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: IL
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 <= D:
            IL = CCM_IL_on(Vd, L, D, Ts, Io, (t%Ts))
        else:
            IL = CCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        if (t/Ts) % 1 <= D:
            IL = DCM_IL_on(Vd, L, D, Ts, (t%Ts))
        else:
            IL = DCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    return IL
            
def get_Vsw(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the switch voltage in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: Vsw
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 <= D:
            Vsw = 0
        else:
            Vsw = Vo
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        Delta_1 = 2*L*Io / (Vd*D*Ts)
        if (t/Ts) % 1 <= D:
            Vsw = 0
        elif (t/Ts) % 1 <= D + Delta_1:
            Vsw = Vo
        else:
            Vsw = Vd
    return Vsw

def get_Isw(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the switch current in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: Isw
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 < D:
            Isw = CCM_IL_on(Vd, L, D, Ts, Io, (t%Ts))
        else:
            Isw = 0
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        if (t/Ts) % 1 < D:
            Isw = DCM_IL_on(Vd, L, D, Ts, (t%Ts))
        else:
            Isw = 0
    return Isw

def get_VD(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the diode voltage in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: VD
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 < D:
            VD = Vo
        else:
            VD = 0
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        Delta_1 = 2*L*Io / (Vd*D*Ts)
        if (t/Ts) % 1 < D:
            VD = Vo
        elif (t/Ts) % 1 < D + Delta_1:
            VD = 0
        else:
            VD = Vo - Vd
    return VD

def get_ID(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the diode current in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: ID
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 < D:
            ID = 0
        else:
            ID = CCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        if (t/Ts) % 1 < D:
            ID = 0
        else:
            ID = DCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts))
    return ID

def get_VC(Vd, Vo, L, C, Ts, Io, t):
    """Obtain the value of the capacitor voltage in a Boost converter in an instant t.
    MEAN VALUE MUST BE ADJUSTED

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: VC
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 <= D:
            VC = VC_on(Vo, Vd, L, C, D, Ts, Io, t%Ts - D*Ts/2)
        else:
            VC = CCM_VC_off(Vo, Vd, L, C, D, Ts, Io, t%Ts)
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        Delta_1 = 2*L*Io / (Vd*D*Ts)
        if (t/Ts) % 1 <= D:
            VC = VC_on(Vo, Vd, L, C, D, Ts, Io, t%Ts)
        elif (t/Ts) % 1 <= D + Delta_1:
            VC = DCM_VC_off(Vo, Vd, L, C, D, Ts, Io, t%Ts)
        else:
            VC = VC_on(Vo, Vd, L, C, D, Ts, Io, t%Ts - Ts)
    return VC

def CCM_VC_off(Vo, Vd, L, C, D, Ts, Io, t):
    VC = -((1/C)*((t-D*Ts)*((D-1)*(Vo-Vd)*t - (D-1)*Ts*Vo + (D-1)*Ts*Vd + 2*D*Io*L))/(2*(D-1)*L) + (Io/C)*(D/2)*Ts)
    return VC

def VC_on(Vo, Vd, L, C, D, Ts, Io, t):
    VC = -(Io/C)*t
    return VC

def DCM_VC_off(Vo, Vd, L, C, D, Ts, Io, t):
    Delta_Il = (Vd/L)*D*Ts
    VC = (-((Vo-Vd)/L)*(t**2/2-D*Ts*t) + (Delta_Il -Io)*t)/C + (-Io*D*Ts-((Vo-Vd)/L)*(((D*Ts)**2)/2) - (Delta_Il - Io)*D*Ts)/C 
    return VC

def get_IC(Vd, Vo, L, Ts, Io, t):
    """Obtain the value of the capacitor current in a Boost converter in an instant t

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: IC
    """
    D = 1 - Vd/Vo
    Iob = Vo*Ts*((1-D)**2)*D / (2*L)
    if Io >= Iob:       # If CCM
        if (t/Ts) % 1 <= D:
            IC = -Io
        else:
            IC = CCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts)) - Io
    else:
        D = np.sqrt(( 2*L*Io / (Vd*Ts)) * (Vo/Vd - 1))
        if (t/Ts) % 1 <= D:
            IC = -Io
        else:
            IC = DCM_IL_off(Vd, Vo, L, D, Ts, Io, (t%Ts)) - Io
    return IC

def get_VO(Vd, Vo, L, C, ESR, Ts, Io, t):
    """Obtain the value of the output voltage in a Boost converter in an instant t.
    MEAN VALUE MUST BE ADJUSTED

    Args:
        Vd (float)
        Vo (float)
        L (float)
        Ts (float)
        Io (float)
        t (float): time passed since the Boost converter was turned on

    Returns:
        float: VO
    """
    VO = get_VC(Vd, Vo, L, C, Ts, Io, t) + ESR*get_IC(Vd, Vo, L, Ts, Io, t)
    return VO

def CCM_get_Delta_Q(Vd, Vo, D, Ts, L, Io):
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

def DCM_get_Delta_Q(Vd, Vo, D, Ts, L, Io):
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