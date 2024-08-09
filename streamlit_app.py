import streamlit as st

import functions_material as rc

# PROJECT DATA INPUT
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
        concrete = st.selectbox("Concrete strength class ($t_{inf}$)",
                                options=rc.concrete_classes,
                                index=5)
        alpha_cc = st.number_input(r"$\alpha_{cc}$",
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
        fck = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[0]
        fcd = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[1]
        eps_cu = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[2]
        lambda_x = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[3]
        fctm = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[4]
        Ecm = rc.concrete_props(concrete, alpha_cc, pfactor_rc)[5]

        st.markdown("$f_{ck}$"+" = {} N/mm$^2$ (comp. strength, char. value)".format(round(fck, 1)))
        st.markdown("$f_{cd}$"+" = {} N/mm$^2$ (comp. strength, design value)".format(round(fcd, 1)))
        st.markdown(r"$\epsilon_{cu}$"+" = {} (limit compressive strain)".format(round(eps_cu, 4)))
        st.markdown(r"$\lambda$"+" = {}  (factor for height of comp. zone)".format(round(lambda_x, 2)))
        st.markdown("$f_{ctm}$"+" = {} N/mm$^2$ (mean tensile strength)".format(round(fctm, 2)))
        st.markdown("$E_{cm}$"+" = {} N/mm$^2$ (Young modulus, mean value)".format(round(Ecm, 0)))

    # REBAR
    st.divider()
    st.subheader("REBAR")
    scol1, scol2 = st.columns(2)
    
    with scol1:
        rebar = st.selectbox("Steel strength class",
                             options=rc.rebar_classes,
                             index=1)
        pfactor_rebar = st.number_input(r"Partial factor - steel ($\gamma_{s}$)", 
                                        value=1.15,
                                        format="%0.2f",
                                        step=0.05)
    
    with scol2:
        fyk = rc.rebar_props(rebar, pfactor_rebar)[0]
        fyd = rc.rebar_props(rebar, pfactor_rebar)[1]
        eps_usd = rc.rebar_props(rebar, pfactor_rebar)[2]
        eps_sy = rc.rebar_props(rebar, pfactor_rebar)[3]
        Es = rc.rebar_props(rebar, pfactor_rebar)[4]
        ftd = rc.rebar_props(rebar, pfactor_rebar)[5]

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
                             options=rc.strand_classes,
                             index=1)
        pfactor_strand = st.number_input(r"Partial factor - strand ($\gamma_{p}$)", 
                                        value=1.15,
                                        format="%0.2f",
                                        step=0.05)
    
    with pcol2:
        dnom_p = rc.strand_props(strand, pfactor_strand)[0]
        Ap_i = rc.strand_props(strand, pfactor_strand)[1]
        fp_01k = rc.strand_props(strand, pfactor_strand)[2]
        fpd = rc.strand_props(strand, pfactor_strand)[3]
        fp_ud = rc.strand_props(strand, pfactor_strand)[4]
        eps_upd = rc.strand_props(strand, pfactor_strand)[5]
        eps_py = rc.strand_props(strand, pfactor_strand)[6]
        Ep = rc.strand_props(strand, pfactor_strand)[7]

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
        t_t2 = 3
        b_bot = 0
        t_b1 = 0
        t_b2 = 0

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
        var_height = st.radio("Variable height", options=["Yes", "No"])
        
    