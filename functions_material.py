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


# LINEAR CREEP COEFFICIENT
def phi_lin(cement, t0, t, RH, fcm, Ac, u):
    """
    cement: cement type ['N' or 'R']
    t0: time of prestress (days)
    t: time where creep factor is determined (days)
    RH: relative humidity
    fcm: mean value of concrete compressive strength
    Ac: concrete cross-sectional area
    u: preimeter of cross-section
    """
    # Factor of cement type
    if cement == 'N':
        alfa_cem = 0
    if cement == 'R':
        alfa_cem = 1
    
    t0_mod = max(0.5, t0*(9/(2 + t0**1.2) +1)**alfa_cem)
    
    # Factors of compressive strength
    alfa1 = (35/fcm)**0.7
    alfa2 = (35/fcm)**0.2
    alfa3 = (35/fcm)**0.5
    
    # phi0 calculation
    h0 = 2*Ac / u
    
    if fcm <= 35:
        phi_RH = 1 + (1-RH/100) / (0.1 * (h0)**(1/3))
    else:
        phi_RH = (1 + (1-RH/100) / (0.1 * (h0)**(1/3))*alfa1) * alfa2
        
    beta_fcm = 16.8 / np.sqrt(fcm)
    beta_t0 = 1 / (0.1 + t0_mod**0.2)
    
    phi_0 = phi_RH * beta_fcm * beta_t0
    
    # beta_c calculation
    if fcm <= 35:
        beta_H = min(1.5*(1+(0.012*RH)**18)*h0 + 250, 1500)
    else:
        beta_H = min(1.5*(1+(0.012*RH)**18)*h0 + 250*alfa3, 1500*alfa3)
    
    beta_c = ((t-t0) / (t-t0+beta_H))**0.3
    
    # creep factor
    phi = phi_0 * beta_c
    
    return phi


# FACTOR FOR NON-LINEAR CREEP
def cal_k_phiNL(fck, sigma):
    if sigma <= 0.45*fck:
        k = 1
    else:
        k = np.exp(1.5*(sigma/fck - 0.45))
    return k


# alpha.E RATIO OF E-MODULI
def alfa_e(Ecm, Es, Ep, phi):
    Ec_eff = Ecm / (1+phi)
    alfa_es = Es / Ec_eff
    alfa_ep = Ep / Ec_eff
    return alfa_es, alfa_ep
