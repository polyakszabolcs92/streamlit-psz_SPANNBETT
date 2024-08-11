import numpy as np

# cs_dims = [h, bw, b_top, t_t1, t_t2, b_bot, t_b1, t_b2, b_top_trapz, h_end]

# BEAM NAME
def beam_name(cs_dims_array, CS_shape, bool_var_height):
    h = cs_dims_array[0]
    bw = cs_dims_array[1]
    btop = cs_dims_array[2]
    tt1 = cs_dims_array[3]
    bbot = cs_dims_array[5]
    tb1 = cs_dims_array[6]
    b_top_trapz = cs_dims_array[8]
    h_end = cs_dims_array[9]
    
    if bool_var_height == False:
        if CS_shape == 'Rect':
            name = "Rectangle {} x {}".format(int(bw/10), int(h/10))
        if CS_shape == 'Trapz':
            name = "Trap {}({}) x {}".format(int(b_top_trapz/10), int(bw/10), int(h/10))
        if CS_shape == 'T':
            name = "T {}({}) / {}({})".format(int(btop/10), int(bw/10), int(h/10), int(tt1/10))
        if CS_shape == 'I':
            name = "I {}({}-{}) / {}({}-{})".format(int(btop/10), int(bw/10), int(bbot/10), int(h/10), int(tt1/10), int(tb1/10))
    else:
        if CS_shape == 'Rect':
            name = "Rectangle {} x {}-{}".format(int(bw/10), int(h/10), int(h_end/10))
        if CS_shape == 'Trapz':
            name = "Trap {}({}) x {}-{}".format(int(b_top_trapz/10), int(bw/10), int(h/10), int(h_end/10))
        if CS_shape == 'T':
            name = "T {}({}) / {}-{}({})".format(int(btop/10), int(bw/10), int(h/10), int(h_end/10), int(tt1/10))
        if CS_shape == 'I':
            name = "I {}({}-{}) / {}-{}({}-{})".format(int(btop/10), int(bw/10), int(bbot/10), int(h/10), int(h_end/10), int(tt1/10), int(tb1/10))
    
    return name


# 'd_mean' MEAN EFFECTIVE DEPTH
def cal_dmean(H, zt, zp, Ast, Ap):
    dt = H - zt
    dp = H - zp
    d =(np.sum(Ast * dt) + np.sum(Ap * dp)) / (np.sum(Ast) + np.sum(Ap))
    return d


# CONCRETE CROSS-SECTION AREA
def cs_area_concrete(cs_shape, h, bw, cs_dims_array, bool_widening):
    btop = cs_dims_array[2]
    tt1 = cs_dims_array[3]
    tt2 = cs_dims_array[4]
    bbot = cs_dims_array[5]
    tb1 = cs_dims_array[6]
    tb2 = cs_dims_array[7]
    b_top_trapz = cs_dims_array[8]
    
    if cs_shape == 'Rect':
        Ac = h * bw
        
    if cs_shape == 'Trapz':
        Ac = 0.5*(bbot + b_top_trapz) * h
        
    if cs_shape == 'fordT':
        Ac = (h-tb1)*bw + tb1*bbot
        
    if cs_shape == 'T':
        t_web = h - (tt1 + tt2)
        Ac = btop*tt1 + 0.5*(btop+bw)*tt2 + t_web*bw
        
    if cs_shape == 'I':
        if bool_widening == False:
            t_web = h - (tt1 + tt2 + tb1 + tb2)
            Ac = btop*tt1 + 0.5*(btop+bw)*tt2 + t_web*bw + 0.5*(bbot+bw)*tb2 + bbot*tb1
        else:
            tf_slope = (tt2 - tt1) / ((btop - bw)/2)
            t_t2_int = tf_slope * ((btop - bw)/2)
            t_b2_new = (bbot-bw/2)
            t_web = h - (tt1+t_t2_int+tb1+t_b2_new)
            Ac = btop*tt1 + 0.5*(btop+bw)*t_t2_int + t_web*bw + 0.5*(bbot+bw)*t_b2_new + bbot*tb1
    
    return Ac


# PERIMETER LENGTH
def u_perim(cs_shape, cs_dims_array):
    h = cs_dims_array[0]
    bw = cs_dims_array[1]
    btop = cs_dims_array[2]
    tt1 = cs_dims_array[3]
    tt2 = cs_dims_array[4]
    bbot = cs_dims_array[5]
    tb1 = cs_dims_array[6]
    tb2 = cs_dims_array[7]
    b_top_trapz = cs_dims_array[8]
    
    if cs_shape == 'Rect':
        u = 2 * (h + bw)
        
    if cs_shape == 'Trapz':
        l_side = np.sqrt(h**2 + (0.5*(b_top_trapz-bw))**2)
        u = b_top_trapz + bw + 2*l_side
        
    if cs_shape == 'fordT':
        u = bw + bbot + (bbot-bw) + 2*h
        
    if cs_shape == 'T':
        tw = h - (tt1+tt2)
        l_slope = np.sqrt(tt2**2 + (0.5*(btop - bw))**2)
        u = btop + bw + 2*(tt1 + l_slope + tw)
        
    if cs_shape == 'I':
        tw = h - (tt1+tt2+tb1+tb2)
        l_slope_top = np.sqrt(tt2**2 + (0.5*(btop - bw))**2)
        l_slope_bot = np.sqrt(tb2**2 + (0.5*(bbot - bw))**2)
        u = btop + bbot + 2*(tt1 + l_slope_top + tw + l_slope_bot + tb1)
    
    return u

