import streamlit as st
import pandas as pd
import numpy as np

# IMPORT MATERIAL PROPERTIES FROM CSV FILES
@st.cache_data
def csv_import():
    df_concrete = pd.read_csv("static/concrete_data.csv", delimiter=";", header=0, index_col=0)
    df_rebar = pd.read_csv("static/rebar_data.csv", delimiter=";", header=0, index_col=0)
    df_strand = pd.read_csv("static/strand_data.csv", delimiter=";", header=0, index_col=0)
    return df_concrete, df_rebar, df_strand

df_concrete, df_rebar, df_strand = csv_import()

# LIST OF STRENGTH CLASSES OF MATERIALS FOR DROPDOWN LISTS
concrete_classes = df_concrete.index.values.tolist()
rebar_classes = df_rebar.index.values.tolist()
strand_classes = df_strand.index.values.tolist()

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
def concrete_props(concrete, alpha, pfactor):
    fck = df_concrete.loc[concrete, 'fck']          
    eta = cal_eta(fck)                               
    fcd = alpha * eta * fck / pfactor
    eps_cu = df_concrete.loc[concrete, 'eps_cu']
    lambda_x = cal_lambda(fck)
    fctm = df_concrete.loc[concrete, 'fctm']                   
    fcm = df_concrete.loc[concrete, 'fcm']
    Ecm = cal_Ecm(fcm)

    return fck, fcd, eps_cu, lambda_x, fctm, Ecm

# GET REBAR PROPERTIES AS A LIST
def rebar_props(rebar, pfactor):
    fyk = df_rebar.loc[rebar, 'fyk']
    eps_usk = df_rebar.loc[rebar, 'eps.us']
    eps_usd = 0.9*eps_usk
    Es = df_rebar.loc[rebar, 'Es']
    fyd = fyk / pfactor
    k = df_rebar.loc[rebar, 'k']
    ftd = k*fyd                                     
    eps_sy = fyd / Es

    return fyk, fyd, eps_usd, eps_sy, Es, ftd

# GET STRAND PROPERTIES AS A LIST
def strand_props(strand, pfactor):
    dnom_p = df_strand.loc[strand, 'd.nom']
    fp_01k = df_strand.loc[strand, 'fp.0.1k']
    fpk = df_strand.loc[strand, 'fpk']
    Ap_i = df_strand.loc[strand, 'Ap']
    eps_upk = df_strand.loc[strand, 'eps.upk']
    eps_upd = 0.9*eps_upk
    Ep = df_strand.loc[strand, 'Ep']
    fpd = fp_01k / pfactor
    fp_ud = fpd / 0.9
    eps_py = fpd / Ep

    return dnom_p, Ap_i, fp_01k, fpd, fp_ud, eps_upd, eps_py, Ep