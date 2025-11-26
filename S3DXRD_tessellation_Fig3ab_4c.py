#%% Input parameters

# input file structure
# x y u_11 u_12 u_13 u_21 u_22 u_23 u_31 u_32 u_33

input_file = 'example_1.txt'
misori_threshold = 9 #degree unit, equal or greater
ipf_sel = 'ipf_z' #'ipf_x', 'ipf_y', 'ipf_z'
phase = 0     # 0 / 1 / 2 : cubic / hexagonal / tetragonal

#%% 

import numpy as np
import matplotlib.pyplot as plt
import math
from skimage import measure
import networkx as nx

class IPF_cal:
    # direc = 0 :  z-axis
    # direc = 1 :  x-axis
    # direc = 2 :  y-axis
    # phase = 0 :  cubic
    # phase = 1 :  hexagonal
    # phase = 2 :  tetragonal

    def __init__(self, u11, u12, u13, u21, u22, u23, u31, u32, u33, direct, phase):
        self.u11 = u11
        self.u12 = u12
        self.u13 = u13
        self.u21 = u21
        self.u22 = u22
        self.u23 = u23
        self.u31 = u31
        self.u32 = u32
        self.u33 = u33
        self.direct = direct
        self.phase = phase
        self.n_div = 500
    
    def set_axis(self):
        u2z = [self.u31, self.u32, self.u33]
        u2x = [self.u11, self.u12, self.u13]
        u2y = [self.u21, self.u22, self.u23]

        if self.direct == 1:
            axis = u2x 
        elif self.direct == 2 :
            axis = u2y
        else:
            axis = u2z
        return axis

    def tri_sym(self):
        x, y, z = self.set_axis()
        x0 = math.fabs(x)
        y0 = math.fabs(y)
        z0 = math.fabs(z)
        if self.phase == 0:
            if x0 < y0:
                x1 = y0
                y1 = x0
            else :
                x1 = x0
                y1 = y0
            z1 = z0
            if x1 > z1:
                if y1 > z1 :
                    x2 = y1
                    y2 = z1
                    z2 = x1
                else:
                    x2 = z1
                    y2 = y1
                    z2 = x1
            else:
                x2 = x1
                y2 = y1
                z2 = z1
            xi = 2 * x2 / (1 + z2)
            eta = 2 * y2 / (1 + z2)
            r = xi, eta
            
        elif self.phase == 1:
            if z >= 0:
                xi0= x / (1 + z)
                eta0= y / (1 + z)
            else :
                xi0 = -x / (1 - z)
                eta0 = -y / (1 - z)
            xi1 = math.fabs(xi0)
            eta1 = math.fabs(eta0)
            tan3 = math.tan(math.pi / 3)
            tan6 = math.tan(math.pi / 6)
            sin_p = math.sin(math.pi / 3)
            sin_m = math.sin(-math.pi / 3)
            cos_p = math.cos(math.pi / 3)
            cos_m = math.cos(-math.pi / 3)
            if eta1 > xi1 * tan3:
                xi2 = cos_m * xi1 - sin_m * eta1
                eta2 = sin_m * xi1 + cos_m * eta1
            elif eta1 > xi1 * tan6:
                xi2 = cos_p * xi1 + sin_p * eta1
                eta2 = sin_p * xi1 - cos_p * eta1
            else:
                xi2 = xi1
                eta2 = eta1
            r = xi2, eta2
        elif self.phase == 2:
            if x0 < y0 :
                x1 = y0
                y1 = x0
            else :
                x1 = x0
                y1 = y0
            z1 = z0
            xi = x1 / (1.0 + z1)
            eta = y1 / (1.0 + z1)
            r = xi, eta
        return r
    
    def tri2RGB_sym(self):
        x, y = self.tri_sym()
        if self.phase == 0:
            x1 = 0.0
            y1 = 0.0
            x2 = 0.82843
            y2 = 0.0
            x3 = 0.73059
            y3 = 0.73059
            
        elif self.phase == 1:
            x1 = 0.0
            y1 = 0.0
            x2 = 1.0
            y2 = 0.0
            x3 = 0.866
            y3 = 0.50

        elif self.phase == 2:
            x1 = 0.0
            y1 = 0.0
            x2 = 1.0
            y2 = 0.0
            x3 = 0.7071068
            y3 = x3

        R12 = x2
        R31 = math.sqrt(x3 * x3 + y3 * y3)
        r1 = math.sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1))
        r2 = math.sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2))
        r3 = math.sqrt((x - x3) * (x - x3) + (y - y3) * (y - y3))
        R = 1.0 - r1 / R31
        G = 1.0 - r2 / R12
        B = 1.0 - r3 / R31
        return math.fabs(R), math.fabs(G), math.fabs(B)
        

    def tri2RGB_sym_2(self, x, y):
        if self.phase == 0:
            x1 = 0.0
            y1 = 0.0
            x2 = 0.82843
            y2 = 0.0
            x3 = 0.73059
            y3 = 0.73059
            
        elif self.phase == 1:
            x1 = 0.0
            y1 = 0.0
            x2 = 1.0
            y2 = 0.0
            x3 = 0.866
            y3 = 0.50

        elif self.phase == 2:
            x1 = 0.0
            y1 = 0.0
            x2 = 1.0
            y2 = 0.0
            x3 = 0.7071068
            y3 = x3

        R12 = x2
        R31 = math.sqrt(x3 * x3 + y3 * y3)
        r1 = math.sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1))
        r2 = math.sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2))
        r3 = math.sqrt((x - x3) * (x - x3) + (y - y3) * (y - y3))
        R = 1.0 - r1 / R31
        G = 1.0 - r2 / R12
        B = 1.0 - r3 / R31
        
        return math.fabs(R), math.fabs(G), math.fabs(B)

    def show_RGB(self):
        n_div = 3
        rgb = 1.0,1.0,1.0
        c = [[rgb for j in range (n_div)] for i in range (n_div)]
        RGB = self.tri2RGB_sym()
        R = RGB[0]
        G = RGB[1]
        B = RGB[2]
        c[1][1] = R, G, B
        print(R, G, B)
        im = plt.imshow(c, origin='lower',interpolation='nearest')
        plt.show()
        return im
    
    def show_RGB_basic_tri_sym(self):
        n_div = self.n_div
        n_shift=int(n_div*0.1)        
        rgb = 1.0,1.0,1.0
        if self.phase == 0:
            x3=0.73059     
            z = [[rgb for j in range (n_div)] for i in range (n_div)]
            for ix in range(n_div):
                for iy in range(n_div):
                    x=float(ix) / float(n_div)
                    y=float(iy) / float(n_div)
                    if x <= x3:
                        if y <= x:
                            z[iy+n_shift][ix+n_shift]= self.tri2RGB_sym_2(x, y)
                    else :
                        if y * y <= 8 - (x + 2) * (x + 2):
                            z[iy+n_shift][ix+n_shift]= self.tri2RGB_sym_2(x, y)
        elif self.phase == 1:
            x3 = 0.866
            z = [[rgb for j in range (n_div + n_shift * 2)] for i in range (n_div + n_shift * 2)]

            for ix in range(n_div):
                for iy in range(n_div):
                    x=float(ix) / float(n_div)
                    y=float(iy) / float(n_div)
                    if x <= x3:
                        if y <= x * math.tan(math.pi / 6):
                            z[iy + n_shift][ix + n_shift] = self.tri2RGB_sym_2(x, y)
                    else :
                        if y * y <= 1 - x * x :
                            z[iy + n_shift][ix + n_shift]= self.tri2RGB_sym_2(x, y)
            
        elif self.phase == 2:
            x3 = 0.7071068
            z = [[rgb for j in range (n_div + n_shift * 2)] for i in range (n_div + n_shift * 2)]
            for ix in range(n_div):
                for iy in range(n_div):
                    x=float(ix)/float(n_div)
                    y=float(iy)/float(n_div)
                    if x <= x3:
                        if y <= x:
                            z[iy + n_shift][ix + n_shift]= self.tri2RGB_sym_2(x, y)
                    else :
                        if y * y <= 1 - x * x :
                            z[iy + n_shift][ix + n_shift]= self.tri2RGB_sym_2(x, y)
        im = plt.imshow(z,origin='lower', interpolation='nearest')
        plt.show()
        return im


def tri2xyz(xi, eta):
    psi = xi * xi + eta * eta
    z = (4.0 - psi) / (4.0 + psi)
    sin_alpha = math.sin(math.acos(z))
    try:
        beta = math.atan(eta / xi)
    except:
        beta = 0
    sin_beta = math.sin(beta)
    cos_beta = math.cos(beta)
    y = sin_alpha * sin_beta
    x = sin_alpha * cos_beta
    return x, y, z

def tri2misori(xi0, eta0, xi1, eta1):
    x0, y0, z0 = tri2xyz( xi0, eta0 )
    x1, y1, z1 = tri2xyz( xi1, eta1 )
    r = x0 * x1 + y0 * y1 + z0 * z1
    if r > 1:
        r = 1
    misori_rad = math.acos(r)
    misori_deg = misori_rad * 180.0 / math.pi
    return misori_deg

#%%

s = IPF_cal(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, phase)
fid = open(input_file, 'r')
header = fid.readline().split()
x = []
y = []
u_mat = []
rgb_x = []
rgb_y = []
rgb_z = []
tri_cubic_x0 = []
tri_cubic_y0 = []
tri_cubic_z0 = []
while True:
    line = fid.readline()
    if not line: 
        break
    line_s = line.split()
    x.append(float(line_s[0]))
    y.append(float(line_s[1]))
    u = np.array([[float(line_s[2]), float(line_s[3]), float(line_s[4])],
                  [float(line_s[5]), float(line_s[6]), float(line_s[7])],
                  [float(line_s[8]), float(line_s[9]), float(line_s[10])]])
    u_mat.append(u)
    s.u11 = u[0][0]
    s.u12 = u[0][1]
    s.u13 = u[0][2]
    s.u21 = u[1][0]
    s.u22 = u[1][1]
    s.u23 = u[1][2]
    s.u31 = u[2][0]
    s.u32 = u[2][1]
    s.u33 = u[2][2]
    s.direct = 1
    tc_x0 = s.tri_sym()
    rgb_x.append(s.tri2RGB_sym())
    tri_cubic_x0.append([tc_x0[0], tc_x0[1]])
    
    s.direct = 2
    rgb_y.append(s.tri2RGB_sym())
    tc_y0 = s.tri_sym()
    tri_cubic_y0.append([tc_y0[0], tc_y0[1]])
    
    s.direct = 0
    rgb_z.append(s.tri2RGB_sym())
    tc_z0 = s.tri_sym()
    tri_cubic_z0.append([tc_z0[0], tc_z0[1]])

fid.close()

x_a = np.array(x, dtype=np.int32)
y_a = np.array(y, dtype=np.int32)
u_mat_a = np.array(u_mat)
rgb_x_a = np.array(rgb_x)
rgb_y_a = np.array(rgb_y)
rgb_z_a = np.array(rgb_z)
x_min = np.array(x_a).min()
y_min = np.array(y_a).min()
x_max = np.array(x_a).max()
y_max = np.array(y_a).max()
img_lv = np.zeros((y_max + 2, x_max + 2))
img_in = np.zeros((y_max + 2, x_max + 2))*np.nan

for m, n in enumerate(x):
    img_in[y_a[m] + 1, x_a[m] + 1] = m

r = []
s = 1
for row in range(1, img_in.shape[0]):
    for col in range(1, img_in.shape[1]):
        if ~np.isnan(img_in[row, col]):
            addr_c = int(img_in[row, col])
            
            if np.isnan(img_in[row - 1, col]) & np.isnan(img_in[row, col - 1]):
                img_lv[row, col] = s
                s = s + 1
                r.append((img_lv[row, col], img_lv[row,col]))

            elif ~np.isnan(img_in[row - 1, col]) & np.isnan(img_in[row, col - 1]):
                addr_u = int(img_in[row - 1, col])
                tc_x0_c = tri_cubic_x0[addr_c]
                tc_x0_u = tri_cubic_x0[addr_u]
                tc_y0_c = tri_cubic_y0[addr_c]
                tc_y0_u = tri_cubic_y0[addr_u]
                tc_z0_c = tri_cubic_z0[addr_c]
                tc_z0_u = tri_cubic_z0[addr_u]
                tri2misori_x_cu = tri2misori(tc_x0_c[0], tc_x0_c[1], tc_x0_u[0], tc_x0_u[1])
                tri2misori_y_cu = tri2misori(tc_y0_c[0], tc_y0_c[1], tc_y0_u[0], tc_y0_u[1])
                tri2misori_z_cu = tri2misori(tc_z0_c[0], tc_z0_c[1], tc_z0_u[0], tc_z0_u[1])
                tri2misori_cu = max(tri2misori_x_cu, tri2misori_y_cu, tri2misori_z_cu)
                
                if (tri2misori_cu > misori_threshold):
                    img_lv[row, col] = s
                    s = s + 1
                    r.append((img_lv[row, col], img_lv[row,col]))
                elif (tri2misori_cu <= misori_threshold):
                    img_lv[row, col] = img_lv[row - 1, col]
            
            elif np.isnan(img_in[row - 1, col]) & ~np.isnan(img_in[row, col - 1]):
                addr_l = int(img_in[row, col - 1])
                tc_x0_c = tri_cubic_x0[addr_c]
                tc_x0_l = tri_cubic_x0[addr_l]
                tc_y0_c = tri_cubic_y0[addr_c]
                tc_y0_l = tri_cubic_y0[addr_l]
                tc_z0_c = tri_cubic_z0[addr_c]
                tc_z0_l = tri_cubic_z0[addr_l]
                tri2misori_x_cl = tri2misori(tc_x0_c[0], tc_x0_c[1], tc_x0_l[0], tc_x0_l[1])
                tri2misori_y_cl = tri2misori(tc_y0_c[0], tc_y0_c[1], tc_y0_l[0], tc_y0_l[1])
                tri2misori_z_cl = tri2misori(tc_z0_c[0], tc_z0_c[1], tc_z0_l[0], tc_z0_l[1])
                tri2misori_cl = max(tri2misori_x_cl, tri2misori_y_cl, tri2misori_z_cl)
                 
                if (tri2misori_cl > misori_threshold):
                    img_lv[row, col] = s
                    s = s + 1
                    r.append((img_lv[row, col], img_lv[row,col]))
                elif (tri2misori_cl <= misori_threshold):
                    img_lv[row, col] = img_lv[row, col - 1]
            
            elif ~np.isnan(img_in[row - 1, col]) & ~np.isnan(img_in[row, col - 1]):                
                addr_u = int(img_in[row - 1, col])
                addr_l = int(img_in[row, col - 1])
                tc_x0_c = tri_cubic_x0[addr_c]
                tc_x0_u = tri_cubic_x0[addr_u]
                tc_x0_l = tri_cubic_x0[addr_l]
                tc_y0_c = tri_cubic_y0[addr_c]
                tc_y0_u = tri_cubic_y0[addr_u]
                tc_y0_l = tri_cubic_y0[addr_l]
                tc_z0_c = tri_cubic_z0[addr_c]
                tc_z0_u = tri_cubic_z0[addr_u]
                tc_z0_l = tri_cubic_z0[addr_l]
                tri2misori_x_cu = tri2misori(tc_x0_c[0], tc_x0_c[1], tc_x0_u[0], tc_x0_u[1])
                tri2misori_x_cl = tri2misori(tc_x0_c[0], tc_x0_c[1], tc_x0_l[0], tc_x0_l[1])
                tri2misori_y_cu = tri2misori(tc_y0_c[0], tc_y0_c[1], tc_y0_u[0], tc_y0_u[1])
                tri2misori_y_cl = tri2misori(tc_y0_c[0], tc_y0_c[1], tc_y0_l[0], tc_y0_l[1])
                tri2misori_z_cu = tri2misori(tc_z0_c[0], tc_z0_c[1], tc_z0_u[0], tc_z0_u[1])
                tri2misori_z_cl = tri2misori(tc_z0_c[0], tc_z0_c[1], tc_z0_l[0], tc_z0_l[1])
                tri2misori_cu = max(tri2misori_x_cu, tri2misori_y_cu, tri2misori_z_cu)
                tri2misori_cl = max(tri2misori_x_cl, tri2misori_y_cl, tri2misori_z_cl)
                
                if (tri2misori_cu > misori_threshold) & (tri2misori_cl > misori_threshold):
                    img_lv[row, col] = s
                    s = s + 1
                    r.append((img_lv[row, col], img_lv[row,col]))
                elif (tri2misori_cu <= misori_threshold) & (tri2misori_cl > misori_threshold):
                    img_lv[row, col] = img_lv[row - 1, col]
                elif (tri2misori_cu > misori_threshold) & (tri2misori_cl <= misori_threshold):
                    img_lv[row, col] = img_lv[row, col - 1]
                elif (tri2misori_cu <= misori_threshold) & (tri2misori_cl <= misori_threshold):
                    img_lv[row, col] = img_lv[row, col - 1]
                    if img_lv[row - 1, col] != img_lv[row, col - 1]:
                            r.append((img_lv[row - 1, col], img_lv[row, col - 1]))

rs = set(tuple(r))
rs_list = list(rs)
rs_list_2 = []
for n in rs_list:
    rs_list_2.append(n)

graph = nx.Graph(rs_list_2) 
graph_connected = [tuple(c) for c in nx.connected_components(graph)]
img_lv = img_lv[1:, 1:]
img_lv_2 = np.zeros_like(img_lv)
for m_2, n_2 in enumerate(graph_connected):
    for m_1, n_1 in enumerate(n_2):
        label_where = (img_lv == n_1).nonzero()
        img_lv_2[label_where] = m_2 + 1

lv_val = np.unique(img_lv_2)
grain_list = []        

for m, n in enumerate(lv_val):
    try: #ignore nan
        grain_find = img_lv_2 == n
        grain_addr = np.int64(img_in[1:, 1:][np.where(grain_find)])
        gr_x      = x_a[grain_addr]
        gr_y      = y_a[grain_addr]
        gr_u_mat  = u_mat_a[grain_addr]
        gr_rgb_x  = rgb_x_a[grain_addr]
        gr_rgb_y  = rgb_y_a[grain_addr]
        gr_rgb_z  = rgb_z_a[grain_addr]
        gr_com_x  = np.mean(x_a[grain_addr])
        gr_com_y  = np.mean(y_a[grain_addr])
        gr_radius = (len(grain_addr) / np.pi)**0.5
        gr_contours = measure.find_contours(grain_find * 1, 0.5)
        grain_extracted = [gr_x,     gr_y,     gr_u_mat, \
                           gr_rgb_x, gr_rgb_y, gr_rgb_z, \
                           gr_com_x, gr_com_y, gr_radius, \
                           gr_contours]
        grain_list.append(grain_extracted)
    except:
        pass

#%%
com = []
weights = []
rgb = []
for n in grain_list:
    com.append([n[6], n[7]])
    weights.append(n[8])
    rgb.append([n[3], n[4], n[5]])

weights = np.array(weights)
grid_points = np.c_[x, y]

d_AWV = np.linalg.norm((grid_points[:, None, :] - com), axis=2) - weights #Additively Weighted Voronoi
d_L   = np.linalg.norm((grid_points[:, None, :] - com), axis=2)**2 - weights**2 # Laguerre
d_V   = np.linalg.norm((grid_points[:, None, :] - com), axis=2) # Voronoi

indices_AWV = np.argmin(d_AWV, axis=1)
ipfx_AWV = np.ones((len(indices_AWV), 3))
ipfy_AWV = np.ones((len(indices_AWV), 3))
ipfz_AWV = np.ones((len(indices_AWV), 3))
ipfn_AWV = np.zeros(len(indices_AWV))

indices_L = np.argmin(d_L, axis=1)
ipfx_L = np.ones((len(indices_L), 3))
ipfy_L = np.ones((len(indices_L), 3))
ipfz_L = np.ones((len(indices_L), 3))
ipfn_L = np.zeros(len(indices_L))

indices_V = np.argmin(d_V, axis=1)
ipfx_V = np.ones((len(indices_V), 3))
ipfy_V = np.ones((len(indices_V), 3))
ipfz_V = np.ones((len(indices_V), 3))
ipfn_V = np.zeros(len(indices_V))

for m, n in enumerate(com):
    f_AWV = (indices_AWV == m).nonzero()[0]
    f_L = (indices_L == m).nonzero()[0]
    f_V = (indices_V == m).nonzero()[0]

    if (len(f_AWV) != 0):
        for p in f_AWV:
            ipfx_AWV[p][0] = rgb[m][0][0][0]
            ipfx_AWV[p][1] = rgb[m][0][0][1]
            ipfx_AWV[p][2] = rgb[m][0][0][2]
            ipfy_AWV[p][0] = rgb[m][1][0][0]
            ipfy_AWV[p][1] = rgb[m][1][0][1]
            ipfy_AWV[p][2] = rgb[m][1][0][2]
            ipfz_AWV[p][0] = rgb[m][2][0][0]
            ipfz_AWV[p][1] = rgb[m][2][0][1]
            ipfz_AWV[p][2] = rgb[m][2][0][2]
            ipfn_AWV[p] = m

    if (len(f_L) != 0):
        for p in f_L:
            ipfx_L[p][0] = rgb[m][0][0][0]
            ipfx_L[p][1] = rgb[m][0][0][1]
            ipfx_L[p][2] = rgb[m][0][0][2]
            ipfy_L[p][0] = rgb[m][1][0][0]
            ipfy_L[p][1] = rgb[m][1][0][1]
            ipfy_L[p][2] = rgb[m][1][0][2]
            ipfz_L[p][0] = rgb[m][2][0][0]
            ipfz_L[p][1] = rgb[m][2][0][1]
            ipfz_L[p][2] = rgb[m][2][0][2]
            ipfn_L[p] = m

    if (len(f_V) != 0):
        for p in f_V:
            ipfx_V[p][0] = rgb[m][0][0][0]
            ipfx_V[p][1] = rgb[m][0][0][1]
            ipfx_V[p][2] = rgb[m][0][0][2]
            ipfy_V[p][0] = rgb[m][1][0][0]
            ipfy_V[p][1] = rgb[m][1][0][1]
            ipfy_V[p][2] = rgb[m][1][0][2]
            ipfz_V[p][0] = rgb[m][2][0][0]
            ipfz_V[p][1] = rgb[m][2][0][1]
            ipfz_V[p][2] = rgb[m][2][0][2]
            ipfn_V[p] = m

if ipf_sel == 'ipf_x':
    sel = 3
    plt_dir_AWV = ipfx_AWV
    plt_dir_L = ipfx_L
    plt_dir_V = ipfx_V
elif ipf_sel == 'ipf_y':
    sel = 4
    plt_dir_AWV = ipfy_AWV
    plt_dir_L = ipfy_L
    plt_dir_V = ipfy_V
elif ipf_sel == 'ipf_z':
    sel = 5
    plt_dir_AWV = ipfz_AWV
    plt_dir_L = ipfz_L
    plt_dir_V = ipfz_V
else:
    print('Wrong direction')

#%%
grain_boundary_AWV = []
grain_boundary_L = []
grain_boundary_V = []

img_s = np.zeros_like(img_lv).astype(int)
for n in range(int(ipfn_AWV.max())):
    img_a = img_s.copy()
    gbxy = ipfn_AWV == n
    img_a[np.int64(y_a[gbxy]), np.int64(x_a[gbxy])] = 1
    grain_contours_AVW = measure.find_contours(img_a, 0.5)
    grain_boundary_AWV.append(grain_contours_AVW)

for n in range(int(ipfn_L.max())):
    img_a = img_s.copy()
    gbxy = ipfn_L == n
    img_a[np.int64(y_a[gbxy]), np.int64(x_a[gbxy])] = 1
    grain_contours_L = measure.find_contours(img_a, 0.5)
    grain_boundary_L.append(grain_contours_L)

for n in range(int(ipfn_V.max())):
    img_a = img_s.copy()
    gbxy = ipfn_V == n
    img_a[np.int64(y_a[gbxy]), np.int64(x_a[gbxy])] = 1
    grain_contours_V = measure.find_contours(img_a, 0.5)
    grain_boundary_V.append(grain_contours_V)
    
#%%

fig = plt.figure(figsize=(8, 5))
ax_a = plt.subplot(2, 3, 1)
ax_b = plt.subplot(2, 3, 2)
ax_c = plt.subplot(2, 3, 3)
ax_d = plt.subplot(2, 3, 4)
ax_e = plt.subplot(2, 3, 5)
ax_f = plt.subplot(2, 3, 6)

for n in grain_list:
    ax_a.scatter(n[0], n[1], c=n[sel], marker='s',s=1)

ax_b.imshow(img_lv_2)
ax_b.invert_yaxis()

for n in grain_list:
    gb_ct = n[9]
    for m in gb_ct:
        ax_c.plot(m[:,1], m[:,0],'k',linewidth=0.5)

for n in grain_boundary_AWV:
    for m in n:
        ax_d.plot(m[:,1], m[:,0],'k',linewidth=0.5)
        
for n in grain_boundary_L:
    for m in n:
        ax_e.plot(m[:,1], m[:,0],'k',linewidth=0.5)
        
for n in grain_boundary_V:
    for m in n:
        ax_f.plot(m[:,1], m[:,0],'k',linewidth=0.5)
        
ax_d.scatter(grid_points[:, 0], grid_points[:, 1], c=plt_dir_AWV, marker='s', s=1)
ax_e.scatter(grid_points[:, 0], grid_points[:, 1], c=plt_dir_L, marker='s', s=1)
ax_f.scatter(grid_points[:, 0], grid_points[:, 1], c=plt_dir_V, marker='s', s=1)

ax_a.set_title('Measurement')
ax_b.set_title('Labeled grain')
ax_c.set_title('Extrated grain boundary')
ax_d.set_title('AWVT')
ax_e.set_title('Laguerre Tessellation')
ax_f.set_title('Voronoi Tessellation')

ax_b.set_xlim(ax_a.get_xlim())
ax_b.set_ylim(ax_a.get_ylim())

ax_a.axis('equal')
ax_b.axis('equal')
ax_c.axis('equal')
ax_d.axis('equal')
ax_e.axis('equal')
ax_f.axis('equal')

fig.tight_layout()


