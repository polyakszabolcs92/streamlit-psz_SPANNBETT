import numpy as np

# 'Lambda' factor for effective height of compressive zone
def cal_lambda(fck):
    if fck <= 50:
        lam = 0.8
    else:
        lam = 0.8 - (fck-50)/400
    return lam

# 'eta' factor for effective compressive strength
def cal_eta(fck):   
    if fck <= 50:
        eta = 1.0
    else:
        eta = 1.0 - (fck-50)/200
    return eta

# Ecm mean value of concrete elastic modulus [MPa]
def cal_Ecm(fcm):
    E = 22*(fcm/10)**(0.3) *1000
    return E

# GET CONCRETE PROPERTIES AS A LIST
def concrete_props(concrete, alpha, pfactor, df):
    fck = df.loc[concrete, 'fck']          
    eta = cal_eta(fck)                               
    fcd = alpha * eta * fck / pfactor
    eps_cu = df.loc[concrete, 'eps_cu']
    lambda_x = cal_lambda(fck)
    fctm = df.loc[concrete, 'fctm']                   
    fcm = df.loc[concrete, 'fcm']
    Ecm = cal_Ecm(fcm)

    return fck, fcd, eps_cu, lambda_x, fctm, Ecm, fcm

# GET REBAR PROPERTIES AS A LIST
def rebar_props(rebar, pfactor, df):
    fyk = df.loc[rebar, 'fyk']
    eps_usk = df.loc[rebar, 'eps.us']
    eps_usd = 0.9*eps_usk
    Es = df.loc[rebar, 'Es']
    fyd = fyk / pfactor
    k = df.loc[rebar, 'k']
    ftd = k*fyd                                     
    eps_sy = fyd / Es

    return fyk, fyd, eps_usd, eps_sy, Es, ftd

# GET STRAND PROPERTIES AS A LIST
def strand_props(strand, pfactor, df):
    dnom_p = df.loc[strand, 'd.nom']
    fp_01k = df.loc[strand, 'fp.0.1k']
    fpk = df.loc[strand, 'fpk']
    Ap_i = df.loc[strand, 'Ap']
    eps_upk = df.loc[strand, 'eps.upk']
    eps_upd = 0.9*eps_upk
    Ep = df.loc[strand, 'Ep']
    fpd = fp_01k / pfactor
    fp_ud = fpd / 0.9
    eps_py = fpd / Ep

    return dnom_p, Ap_i, fp_01k, fpd, fp_ud, eps_upd, eps_py, Ep