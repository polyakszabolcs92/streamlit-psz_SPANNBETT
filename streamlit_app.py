import streamlit as st
import pandas as pd
import numpy as np
import functions_material as rc

# PAGE CONFIG
st.set_page_config(page_title="Spannbett-psz", layout="centered")

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


# PROJECT DATA INPUT-----------------------------------------------------------------
with st.popover("Project data"):
    project = st.text_input("Project name", placeholder="Type project name here")
    element = st.text_input("Element ID", placeholder="Type element ID here")
    designer = st.text_input("Designer", placeholder="Type designer name here")


tab1, tab2, tab3,  tab4, tab5 = st.tabs(["MATERIAL", 
                                         "GEOMETRY", 
                                         "REINFORCEMENT & PRESTRESS", 
                                         "LOADS", 
                                         "ANALYSIS RESULTS"])

# MATERIAL INPUT
with tab1:

    # CONCRETE
    st.subheader("CONCRETE")

    ccol1, ccol2 = st.columns(2)
    with ccol1:
        concrete = st.selectbox("Concrete strength class",
                                options=concrete_classes,
                                index=5)
        alpha_cc = st.number_input(r"Factor for long-term effects ($\alpha_{cc}$) ",
                                   value= 1.00,
                                   format="%0.2f",
                                   step=0.05,
                                   help="Modification factor of compressive strength, considering long-term effects")
        rho_rc = st.number_input(r"Unit weight ($\rho_{rc}$) [kN/m3]", 
                                 value=25.0,
                                 format="%0.1f",
                                 step=0.1,
                                 help="Unit weight of reinforced concrete, default: 25.0 kN/m3")
        pfactor_rc = st.number_input(r"Partial factor - concrete ($\gamma_{c}$)", 
                                     value=1.50,
                                    format="%0.2f",
                                    step=0.05)
        
    with ccol2:
        fck = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[0]
        fcd = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[1]
        eps_cu = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[2]
        lambda_x = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[3]
        fctm = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[4]
        Ecm = rc.concrete_props(concrete, alpha_cc, pfactor_rc, df_concrete)[5]

        st.markdown("$f_{ck}$"+" = {} N/mm$^2$ (comp. strength, char. value)".format(round(fck, 1)))
        st.markdown("$f_{cd}$"+" = {} N/mm$^2$ (comp. strength, design value)".format(round(fcd, 1)))
        st.markdown(r"$\epsilon_{cu}$"+" = {} (limit compressive strain)".format(round(eps_cu, 4)))
        st.markdown(r"$\lambda$"+" = {}  (factor for height of comp. zone)".format(round(lambda_x, 2)))
        st.markdown("$f_{ctm}$"+" = {} N/mm$^2$ (mean tensile strength)".format(round(fctm, 2)))
        st.markdown("$E_{cm}$"+" = {} N/mm$^2$ (Modulus of elasticity, mean value)".format(round(Ecm, 0)))

    # REBAR
    st.divider()
    st.subheader("REBAR")
    scol1, scol2 = st.columns(2)
    
    with scol1:
        rebar = st.selectbox("Steel strength class",
                             options=rebar_classes,
                             index=1)
        pfactor_rebar = st.number_input(r"Partial factor - steel ($\gamma_{s}$)", 
                                        value=1.15,
                                        format="%0.2f",
                                        step=0.05)
    
    with scol2:
        fyk = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[0]
        fyd = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[1]
        eps_usd = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[2]
        eps_sy = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[3]
        Es = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[4]
        ftd = rc.rebar_props(rebar, pfactor_rebar, df_rebar)[5]

        # st.markdown("Rebar properties")
        st.markdown("$f_{yk}$"+" = {} N/mm$^2$ (steel strength, char. value)".format(round(fyk, 1)))
        st.markdown("$f_{yd}$"+" = {} N/mm$^2$ (steel strength, design value)".format(round(fyd, 1)))
        st.markdown("$f_{td}$"+" = {} N/mm$^2$ (steel rupture strength)".format(round(ftd, 1)))
        st.markdown(r"$\epsilon_{us,d}$"+" = {} (limit strain)".format(round(eps_usd, 4)))
        st.markdown(r"$\epsilon_{sy}$"+" = {} (yield strain)".format(round(eps_sy, 5)))
        st.markdown("$E_{s}$"+" = {} N/mm$^2$ (Young modulus)".format(round(Es, 0)))
    


    # PRESTRESSING STRAND
    st.divider()
    st.subheader("PRESTRESSING STRAND")

    pcol1, pcol2 = st.columns(2)
    
    with pcol1:
        strand = st.selectbox("Strand strength class",
                             options=strand_classes,
                             index=1)
        pfactor_strand = st.number_input(r"Partial factor - strand ($\gamma_{p}$)", 
                                        value=1.15,
                                        format="%0.2f",
                                        step=0.05)
    
    with pcol2:
        dnom_p = rc.strand_props(strand, pfactor_strand, df_strand)[0]
        Ap_i = rc.strand_props(strand, pfactor_strand, df_strand)[1]
        fp_01k = rc.strand_props(strand, pfactor_strand, df_strand)[2]
        fpd = rc.strand_props(strand, pfactor_strand, df_strand)[3]
        fp_ud = rc.strand_props(strand, pfactor_strand, df_strand)[4]
        eps_upd = rc.strand_props(strand, pfactor_strand, df_strand)[5]
        eps_py = rc.strand_props(strand, pfactor_strand, df_strand)[6]
        Ep = rc.strand_props(strand, pfactor_strand, df_strand)[7]

        # st.markdown("Rebar properties")
        st.markdown("$d_{nom}$"+" = {} mm (nominal diameter)".format(round(dnom_p, 1)))
        st.markdown("$A_{pi}$"+" = {} mm$^2$ (strand cross-section area)".format(round(Ap_i, 1)))
        st.markdown("$f_{p,0.1,k}$"+" = {} N/mm$^2$ (strand strength, char. value)".format(round(fp_01k, 1)))
        st.markdown("$f_{pd}$"+" = {} N/mm$^2$ (strand strength, design value)".format(round(fpd, 1)))
        st.markdown("$f_{ud}$"+" = {} N/mm$^2$ (strand rupture strength)".format(round(fp_ud, 1)))
        st.markdown(r"$\epsilon_{up,d}$"+" = {} (limit strain)".format(round(eps_upd, 4)))
        st.markdown(r"$\epsilon_{py}$"+" = {} (yield strain)".format(round(eps_py, 5)))
        st.markdown("$E_{p}$"+" = {} N/mm$^2$ (Young modulus)".format(round(Ep, 0)))


# GEOMETRY INPUT---------------------------------------------------------------------------
with tab2:
    gcol1, gcol2 = st.columns(2)
    
    with gcol1:
        CS_shape = st.selectbox("CROSS-SECTION SHAPE",
                                options=['Rect', 'T', 'I', 'Trapz', 'invT'],
                                index=0)
        h = st.number_input("Beam height (h) [cm]", min_value=0., value= 100.0, format="%0.1f")
        bw = st.number_input("Web width (bw) [cm]", min_value=0., value= 16.0, format="%0.1f")
        b_top = 0
        t_t1 = 0
        t_t2 = 0
        b_bot = 0
        t_b1 = 0
        t_b2 = 0
        b_top_trapz = 0
        h_end = h
        bw_end = bw
        h_sub = h

        if CS_shape == "T":
            b_top = st.number_input("Top flange width (b.top) [cm]", min_value=0., value= 50.0, format="%0.1f")
            t_t1 = st.number_input("Top flange height (t.t1) [cm]", min_value=0., value= 15.0, format="%0.1f")
            t_t2 = st.number_input("Top flange transition (t.t2) [cm]", min_value=0., value= 3.0, format="%0.1f")
        
        if CS_shape == "I":
            b_top = st.number_input("Top flange width (b.top) [cm]", min_value=0., value= 50.0, format="%0.1f")
            t_t1 = st.number_input("Top flange height (t.t1) [cm]", min_value=0., value= 15.0, format="%0.1f")
            t_t2 = st.number_input("Top flange transition (t.t2) [cm]", min_value=0., value= 3.0, format="%0.1f")
            b_bot = st.number_input("Bottom flange width (b.bot) [cm]", min_value=0., value= 30.0, format="%0.1f")
            t_b1 = st.number_input("Bottom flange height (t.b1) [cm]", min_value=0., value= 18.0, format="%0.1f")
            t_b2 = st.number_input("Bottom flange transition (t.b2) [cm]", min_value=0., 
                                   value= (b_bot - bw)/2, format="%0.1f")
        
        if CS_shape == "Trapz":
            b_top_trapz = st.number_input("Top flange width (b.top.trapz) [cm]", min_value=0., value= 30.0, format="%0.1f")
        
        if CS_shape == "invT":
            b_bot = st.number_input("Bottom flange width (b.bot) [cm]", min_value=0., value= 50.0, format="%0.1f")
            t_b1 = st.number_input("Bottom flange height (t.b1) [cm]", min_value=0., value= 25.0, format="%0.1f")


    with gcol2:
        L = st.number_input("Beam length [m]", 
                            min_value=0.,
                            value= 12.00,
                            format="%0.2f")
        var_height = st.checkbox("Variable height", value=False)
        
        if var_height:
            slope = st.number_input("Slope [%]", min_value=0.,
                                    value=3.0,
                                    format="%0.1f")
            h_end = h - (L*100/2)*(slope/100)
            k_height = 0.85      # interpolation factor between beam end and center
            h_sub = k_height*h + (1-k_height)*h_end
                
        if CS_shape == "I":
            web_widening = st.checkbox("Widened web at beam end", value=False)
            if web_widening:
                bw_end = st.number_input("Web width at beam end", min_value=0.,
                                        value=b_bot,
                                        format="%0.1f")
        
    # Gathering cross-section dimensions in a list and converting to [mm]
    cs_dims = [x*10 for x in [h, bw, b_top, t_t1, t_t2, b_bot, t_b1, t_b2, b_top_trapz, h_end]]

with tab3:
    rcol1, rcol2 = st.columns(2)
    with rcol1:
        st.markdown("Tensile reinforcement (bottom)")
        df_As_tens = st.data_editor(pd.DataFrame(data= np.array([[20, 2, 45]]),
                                                columns=["d [mm]", "pcs", "zi [mm]"],
                                                index=[0]),
                                    num_rows='dynamic')
        data_As_tens = df_As_tens.to_numpy()

    with rcol2:
        st.markdown("Compressive reinforcement (top)")
        df_As_comp = st.data_editor(pd.DataFrame(data= np.array([[16, 2, 45]]),
                                                columns=["d [mm]", "pcs", "zi [mm]"],
                                                index=[0]),
                                    num_rows='dynamic')
        data_As_comp = df_As_comp.to_numpy()
    

    st.divider()
    st.markdown("Prestressing")
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        df_Ap = st.data_editor(pd.DataFrame(data= np.array([[100, 3],
                                                            [140, 0],
                                                            [180, 0],
                                                            [220, 0],
                                                            [260, 0],
                                                            [300, 0],
                                                            [340, 0],
                                                            [380, 0]]),
                                                columns=["zi [mm]", "pcs"],
                                                index=np.arange(0,8,1)),
                                    num_rows='fixed',
                                    hide_index=True)
        data_Ap = df_Ap.to_numpy()
    
    with pcol2:
        sig_p0 = st.number_input(r"Prestress (N/mm$^2$)", value=1000, min_value=0, step=50)
        kp_red = st.number_input("Reduction factor for estimated prestress losses", value=0.85, min_value=0., step=0.05)
        concrete_t0 = st.selectbox("Concrete strength at prestress (t0)",
                                    options=concrete_classes,
                                    index=2)
        cement = st.selectbox("Cement type", options= ['R', 'N'], index=0)
        
    with pcol3:
        t0 = st.number_input(r"Time of prestressing ($t_{0}$)", value=3, min_value=0, step=1)
        t1 = st.number_input(r"Time of load application ($t_{1}$)", value=28, min_value=7, step=1)
        t_inf = st.number_input(r"Design working life ($t_{inf}$)", value=18250, min_value=100, step=50)
        RH_t0 = st.number_input("Humidity in temporary state [%]", value=80, min_value=0, step=5)
        RH_inf = st.number_input("Humidity in final state [%]", value=50, min_value=0, step=5)

# Concrete properties in temporary conditions
fck_t0 = rc.concrete_props(concrete_t0, alpha_cc, pfactor_rc, df_concrete)[0]
fctm_t0 = rc.concrete_props(concrete_t0, alpha_cc, pfactor_rc, df_concrete)[4]
fcm_t0 = rc.concrete_props(concrete_t0, alpha_cc, pfactor_rc, df_concrete)[6]
Ecm_t0 = rc.concrete_props(concrete_t0, alpha_cc, pfactor_rc, df_concrete)[5]


# LOADS INPUT---------------------------------------------------------------------------
