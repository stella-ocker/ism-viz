'''
Python script to create a 3D interactive visualization of scattering screens inferred from scintillation arcs 
and various large-scale ISM features of interest, as in Ocker et al. (arXiv:2309.13809). 


Required input data:

pulsar_properties_ism3d.txt
screen_properties_ism3d.txt
L19_map-inner_final.fits
MolCloud_Distances_Zucker2020.dat

To add new sightlines, update pulsar_properties_ism3d.txt and screen_properties_ism3d.txt, then re-run script below.

Version history:

11-20-23 SKO created Github repo; script replicates figure at https://stella-ocker.github.io/scattering_ism3d_ocker2023

'''

from numpy import *
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

### load pulsar & screen properties

psr_jnames = loadtxt('pulsar_properties_ism3d.txt',usecols=[0],dtype=str)
psr_ra = loadtxt('pulsar_properties_ism3d.txt',usecols=[1],dtype=str)
psr_dec = loadtxt('pulsar_properties_ism3d.txt',usecols=[2],dtype=str)
psr_dist = loadtxt('pulsar_properties_ism3d.txt',usecols=[3])

screen_jnames = loadtxt('screen_properties_ism3d.txt',usecols=[0],dtype=str)
screen_dists = loadtxt('screen_properties_ism3d.txt',usecols=[1])
screen_refs = loadtxt('screen_properties_ism3d.txt',usecols=[2],dtype=str)

### grab unique list of literature refs to include in plot legend

unique_screen_psrs,ref_inds = unique(screen_jnames,return_index=True)
reflist = screen_refs[sort(ref_inds)]

### gather screen coordinates into single set of arrays (one per X,Y,Z)
### need them grouped together for plotting purposes

all_screen_X = zeros(shape(screen_dists))
all_screen_Y = zeros(shape(screen_dists))
all_screen_Z = zeros(shape(screen_dists))

for i in range(len(screen_jnames)):

    # grab parameters for a given pulsar
    psr_i = screen_jnames[i]
    screen = screen_dists[i]
    ind = where(psr_jnames==psr_i)[0]
    rai = psr_ra[ind]
    deci = psr_dec[ind]
    disti = psr_dist[ind]

    # define coordinates for pulsar and screens
    screen_dist = disti*(1.-screen) # screen distance in kpc
    screen_coord = SkyCoord(ra=rai,dec=deci,distance=screen_dist*u.kpc)
    screen_x = screen_coord.galactic.cartesian.x.value # heliocentric Galactic Cartesian coordinates
    screen_y = screen_coord.galactic.cartesian.y.value
    screen_z = screen_coord.galactic.cartesian.z.value

    all_screen_X[i] = screen_x[0]
    all_screen_Y[i] = screen_y[0]
    all_screen_Z[i] = screen_z[0]

### ISM Features ###

# load Pelgrims (2020) dust map of local bubble
L19 = hp.read_map('L19_map-inner_final.fits',field=3)
# convert to physical coordinates (Galactic x,y,z in pc with Sun at origin)
vec = asarray(hp.pix2vec(128,arange(hp.nside2npix(128))))
xyz_6pol = vec * L19
xl = xyz_6pol[0]
yl = xyz_6pol[1]
zl = xyz_6pol[2]

# HII regions
Sh2_27 = SkyCoord(l=8*u.degree,b=23.5*u.degree,frame='galactic',distance=0.112*u.kpc)
Sh2_27_diam = 34. # pc
Sh2_7 = SkyCoord(l=-9*u.degree,b=24*u.degree,frame='galactic',distance=0.136*u.kpc)
Sh2_7_diam = 28. # pc
Sh2_205 = SkyCoord(l=148.4*u.degree,b=-0.2*u.degree,frame='galactic',distance=1*u.kpc) # based on Romero & Cappa 2008; treating here two possibly distinct regions as one
Sh2_205_diam = 24. # pc - adding diameters of two regions in Romero & Kappa 2008
Sh147 = SkyCoord(ra='05h40m01s',dec='+27d48m09s',frame='icrs',distance=1.2*u.kpc)
Sh147_diam = 32*2.

uc = linspace(0,2*pi)
vc = linspace(0,2*pi)

r27 = Sh2_27_diam/2.
Sh2_27x = r27*outer(cos(uc),sin(vc)) + Sh2_27.galactic.cartesian.x.value*1000
Sh2_27y = r27*outer(sin(uc),sin(vc)) + Sh2_27.galactic.cartesian.y.value*1000
Sh2_27z = r27*outer(ones(size(uc)),cos(vc)) + Sh2_27.galactic.cartesian.z.value*1000

r7 = Sh2_7_diam/2.
Sh2_7x = r7*outer(cos(uc),sin(vc)) + Sh2_7.galactic.cartesian.x.value*1000
Sh2_7y = r7*outer(sin(uc),sin(vc)) + Sh2_7.galactic.cartesian.y.value*1000
Sh2_7z = r7*outer(ones(size(uc)),cos(vc)) + Sh2_7.galactic.cartesian.z.value*1000

r205 = Sh2_205_diam/2.
Sh2_205x = r205*outer(cos(uc),sin(vc)) + Sh2_205.galactic.cartesian.x.value*1000
Sh2_205y = r205*outer(sin(uc),sin(vc)) + Sh2_205.galactic.cartesian.y.value*1000
Sh2_205z = r205*outer(ones(size(uc)),cos(vc)) + Sh2_205.galactic.cartesian.z.value*1000

r147 = Sh147_diam/2.
Sh147x = r147*outer(cos(uc),sin(vc)) + Sh147.galactic.cartesian.x.value*1000
Sh147y = r147*outer(sin(uc),sin(vc)) + Sh147.galactic.cartesian.y.value*1000
Sh147z = r147*outer(ones(size(uc)),cos(vc)) + Sh147.galactic.cartesian.z.value*1000

# GSH 238+00+09 superbubble
GSH238 = SkyCoord(l=238.0*u.degree,b=0.*u.degree,distance=0.8*u.kpc,frame='galactic') # center of bubble, Heiles 98 
GSH238_r1 = 500. # pc, Puspitarini 2014
GSH238_r2 = 300. # pc, Puspitarini 2014
GSH238x = GSH238_r1*outer(cos(uc),sin(vc)) + GSH238.galactic.cartesian.x.value*1000 
GSH238y = GSH238_r1*outer(sin(uc),sin(vc)) + GSH238.galactic.cartesian.y.value*1000
GSH238z = GSH238_r2*outer(ones(size(uc)),cos(vc)) 

# plasma filament from Wang et al. 2021
Wang_filament_center = SkyCoord(ra='00h57m46s',dec='-24d23m00s',distance=4*u.pc)
Wang_filament_length = deg2rad(1.7)*4 # length in pc (spans 1.7 deg)
Wang_filament_width = deg2rad(1./60)*4 # width in pc (~1 arcmin)
Wang_filament_lons = linspace(145.,160.) # full range of galactic coords
Wang_filament_lats = linspace(-86.2,-87.2)
# draw filament as a line between two points
Wang_filament1 = SkyCoord(l=Wang_filament_lons[0]*u.degree,b=Wang_filament_lats[0]*u.degree,frame='galactic',distance=4*u.pc)
Wang_filament2 = SkyCoord(l=Wang_filament_lons[-1]*u.degree,b=Wang_filament_lats[-1]*u.degree,frame='galactic',distance=4*u.pc)
Wang_filamentx = array([Wang_filament1.galactic.cartesian.x.value ,Wang_filament2.galactic.cartesian.x.value])
Wang_filamenty = array([Wang_filament1.galactic.cartesian.y.value ,Wang_filament2.galactic.cartesian.y.value])
Wang_filamentz = array([Wang_filament1.galactic.cartesian.z.value ,Wang_filament2.galactic.cartesian.z.value])

# Per-Tau shell (Bialy+2021)
pertau_xcenter = -190. # pc
pertau_ycenter = 65 # pc
pertau_zcenter = -84 # pc
pertau_r = 78. # pc
pertau_x = pertau_r*outer(cos(uc),sin(vc)) + pertau_xcenter
pertau_y = pertau_r*outer(sin(uc),sin(vc)) + pertau_ycenter
pertau_z = pertau_r*outer(ones(size(uc)),cos(vc)) + pertau_zcenter

# major local molecular clouds (Zucker+2020)- Galactic heliocentric cartesian 
cloud_name = loadtxt('MolCloud_Distances_Zucker2020.dat',usecols=[0],dtype=str)
cloud_dist = loadtxt('MolCloud_Distances_Zucker2020.dat',usecols=[13]) # 50th percentile
cloud_coordx = loadtxt('MolCloud_Distances_Zucker2020.dat',usecols=[-3])
cloud_coordy = loadtxt('MolCloud_Distances_Zucker2020.dat',usecols=[-2])
cloud_coordz = loadtxt('MolCloud_Distances_Zucker2020.dat',usecols=[-1])

# define colors
prism_colors= px.colors.qualitative.Prism
pastel_colors = px.colors.qualitative.Pastel1
plotly_colors = px.colors.qualitative.Dark24

# create list of colors for pulsar LOSs
#psr_colors = concatenate((plotly_colors[0:5],plotly_colors[6:9],plotly_colors[9:16])) # colors selected for Ocker+2023
psr_colors = copy(plotly_colors) # longer list of colors to accommodate any new lines of sight

### create figure ###

fig = go.Figure()
# Local Bubble 
fig.add_trace(go.Scatter3d(x=xl, y=yl, z=zl,mode='markers',
                           marker=dict(size=1,color=prism_colors[2],opacity=0.03),name='Local Bubble'))
# molecular clouds
fig.add_trace(go.Scatter3d(x=cloud_coordx,y=cloud_coordy,z=cloud_coordz,mode='markers',
                           marker=dict(size=3,color=prism_colors[0],opacity=0.4),name='Molecular Clouds'))
# superbubble
fig.add_trace(go.Surface(z=GSH238z,x=GSH238x,y=GSH238y,
                         colorscale=[[0, prism_colors[3]],[1,prism_colors[3]]],opacity=0.2,showscale=False,
                         name='GSH 238+00+09',showlegend=True))
# Per-Tau shell
fig.add_trace(go.Surface(z=pertau_z,x=pertau_x,y=pertau_y,
                         colorscale=[[0, prism_colors[3]],[1,prism_colors[3]]],opacity=0.2,showscale=False,
                         name='Per-Tau Shell',showlegend=True))
# HII regions
fig.add_trace(go.Surface(z=Sh2_27z,x=Sh2_27x,y=Sh2_27y,
                         colorscale=[[0,prism_colors[5]],[1,prism_colors[5]]],opacity=0.2,showscale=False,
                         name='Sh 2-27',showlegend=True))
fig.add_trace(go.Surface(z=Sh2_7z,x=Sh2_7x,y=Sh2_7y,
                         colorscale=[[0,prism_colors[5]],[1,prism_colors[5]]],opacity=0.2,showscale=False,
                         name='Sh 2-7',showlegend=True))
fig.add_trace(go.Surface(z=Sh2_205z,x=Sh2_205x,y=Sh2_205y,
                         colorscale=[[0,prism_colors[5]],[1,prism_colors[5]]],opacity=0.2,showscale=False,
                         name='Sh 2-205',showlegend=True))
fig.add_trace(go.Surface(z=Sh147z,x=Sh147x,y=Sh147y,
                         colorscale=[[0,prism_colors[5]],[1,prism_colors[5]]],opacity=0.2,showscale=False,
                         name='S147',showlegend=True))

# Wang filament
fig.add_trace(go.Scatter3d(x=Wang_filamentx,y=Wang_filamenty,z=Wang_filamentz,mode='markers',
                           marker=dict(size=1,color=prism_colors[5]),name='Local Filament'))

# screens
fig.add_trace(go.Scatter3d(x=all_screen_X*1000,y=all_screen_Y*1000,z=all_screen_Z*1000,mode='markers',
                           marker_symbol='square',marker=dict(color='black',size=7),name='Scattering Screens',
                           showlegend=True))

# pulsars -- loop through all 
for i in range(len(psr_jnames)):
    psri = psr_jnames[i]
    rai = psr_ra[i]
    deci = psr_dec[i]
    disti = psr_dist[i]
    psr_coord = SkyCoord(ra=rai,dec=deci,distance=disti*u.kpc)
    psr_x = psr_coord.galactic.cartesian.x.value*1000 # pc
    psr_y = psr_coord.galactic.cartesian.y.value*1000 # pc
    psr_z = psr_coord.galactic.cartesian.z.value*1000 # pc
    fig.add_trace(go.Scatter3d(x=[0,psr_x],y=[0,psr_y],z=[0,psr_z],
                           marker=dict(color=psr_colors[i],size=5,opacity=1),
                           line=dict(color=psr_colors[i],width=5),name=psri+' ('+reflist[i]+')'))

# set axis range and titles
fig.update_layout(
    scene = dict(
        xaxis = dict(range=[-1500,1500],),
        yaxis = dict(range=[-1500,1500],),
        zaxis = dict(range=[-1500,1500],),
        xaxis_title='X (pc)',
        yaxis_title='Y (pc)',
        zaxis_title='Z (pc)'),
        width=1000,
        height=700,
        margin=dict(r=20, l=10, b=10, t=10))
fig.update_layout(scene_aspectmode='cube',legend={'itemsizing':'constant'})
fig.update_layout(legend=dict(yanchor="middle",y=0.5)) 

fig.show() # show figure 
#fig.write_html('ism3d_out.html') # writes to html file
