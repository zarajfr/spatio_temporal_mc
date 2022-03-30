"""
The simulation of Markov chain analysis for
anti-refugee incidents in Germany from 01/01/2015 to 30/06/2020

Parameters
----------

FILE_NAME
    Name of the file that contains spatio temporal data

D_bins
    Cut off thresholds to define the possible states (k) of the Markov chains
    These can be either chosen manually or can be obtained using
    methods implemented in mapclassify library; such as mc.Quantiles.

w
    The spatial weight matrix based upon geographical distances,
    calculated using KNN from pysal and K is set to 4.

"""
import numpy as np
import matplotlib.pyplot as plt
import math
from pprint import pprint
import csv
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import geopandas as gpd
import pandas as pd
import pysal
from libpysal.weights.contiguity import Queen
from libpysal.weights import KNN
import libpysal
from splot.libpysal import plot_spatial_weights
from shapely.geometry import Point
import geopy.distance
import giddy
import mapclassify as mc
import pyproj
from pyproj import Proj, CRS,transform, Transformer
from scipy.stats import chi2_contingency

# Choose the name of the file which contains your data
FILE_NAME = "SM.csv"

def read_spatiotemporal_data(file_name):
    incidents = []
    j = 0
    # Change file_name in case the file is in another fold, such as "D:\\myfiles\SM.csv"
    with open(file_name) as f:
        myf = csv.reader(f, delimiter = ',')
        for row in myf:
            x = []
            for i in range(len(row)):
                if(i>5):
                    x.append(float(row[i]))
            incidents.append(x)
            j+=1
    f.close()
    incidents.pop(0)
    return incidents

def pass_spatial_weight(somename):
    transformer = Transformer.from_crs( "epsg:25832", "epsg:4326")

    df = pd.read_csv(somename, delimiter = ',', header = 0)
    d = {"VERX" : df.VERX, "VERY" : df.VERY}
    df2 = pd.DataFrame(data=d)
    for index, row in df2.iterrows():
        x = transformer.transform(row['VERX'],row['VERY'])
        row['VERX'] = x[0]
        row['VERY'] = x[1]
    gdf = gpd.GeoDataFrame(df2, geometry=gpd.points_from_xy(df2.VERX, df2.VERY))
    w = KNN.from_dataframe(gdf, k = 4)
    return w

# Classic Markov chain
def classic_markov_analysis(somename):
    incidents = read_spatiotemporal_data(FILE_NAME)
    incidents1 = np.array(incidents)
    allnumbers = []
    c=0
    qincidents=[]
    for i in range(len(incidents)):
        for j in range(len(incidents[0])):
            allnumbers.append(incidents[i][j])

    # There are numerous ways to compute the classes - such as the ones below:
    # q = mc.Quantiles(allnumbers, k = 3)
    # q = mc.FisherJenks(allnumbers,k=3)
    # q = mc.EqualInterval(allnumbers, k = 3)

    # Here, the classes are customised in a way that the first class refer to the state of no violence
    D_bins = { "SM.csv": [0, 12.0, max(allnumbers)], "SW.csv": [0, 3.0, max(allnumbers)], "DM.csv": [0, 5.0, max(allnumbers)], "DW.csv": [0, 3.0, max(allnumbers)] }
    q = mc.UserDefined(allnumbers, D_bins[somename])

    for i in range(len(incidents)):
        x = []
        for j in range(len(incidents[0])):
            x.append(q.yb[c])
            c+=1
        qincidents.append(x)

    m = giddy.markov.Markov(np.array(qincidents))

    # Transition probability matrix, classes, steady state, and first mean time passages can be obtained
    # print(m.classes)
    print(m.p)
    print(m.steady_state)
    # print(giddy.ergodic.fmpt(m.p))

    # The chi squared test of homogeneity for transition matrix
    stat, p, dof, expected = chi2_contingency(m.transitions)
    print('stat=%.5f, p = %.5f'%(stat,p))
    if p>0.05:
        print("indep")
    else:
        print("dep")

# Spatial Markov chain
def spatial_markov_analysis(somename):
    incidents = read_spatiotemporal_data(somename)
    w = pass_spatial_weight(somename)

    D_bins = { "SM.csv": [0, 12.0], "SW.csv": [0, 3.0], "DM.csv": [0, 5.0], "DW.csv": [0, 3.0] }
    cc = np.array(D_bins[somename])

    sm = giddy.markov.Spatial_Markov(incidents, w, cutoffs = cc, lag_cutoffs = cc, fixed = True , k= 3, m =3,fill_empty_classes=True)
    "**** probability matrices ****"
    print(sm.p)
    "**** summary ****"
    sm.summary()
    # print("**** first mean passage time ****")
    # print(sm.F)
    # print("**** steady state ****")
    # print(sm.S)
    # print("**** Test of spatial dependence ****")
    # print(giddy.markov.kullback(sm.T))

    boxlabel = ["No", "Moderate", "High"]
    sns.set()
    fig, axes = plt.subplots(2,2,figsize = (7,7))
    for i in range(2):
        for j in range(2):
            ax = axes[i,j]
            if i==0 and j==0:
                # ax.axis('off')
                p_temp = sm.p
                im = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0, vmax=1,square=True, cmap="Reds",fmt='.3f')
                ax.set_title("Pooled",fontsize=13)
                # continue
            else:
                p_temp = sm.P[i*2+j-1]
                im = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0, vmax=1,square=True, cmap="Reds",fmt='.3f')
                # ax.set_title("Spatial Lag %d"%(i*2+j-1),fontsize=13)
                ax.set_title("Spatial Lag: "+boxlabel[i*2+j-1],fontsize=13)
    plt.show()


# LISA Markov chain
def lisa_markov_analysis(somename):
    incidents = read_spatiotemporal_data(somename)
    incidents = np.array(incidents)
    w = pass_spatial_weight(somename)
    lm = giddy.markov.LISA_Markov(incidents, w)
    "**** LISA classes ****"
    print(lm.classes)
    "**** LISA transition times ****"
    T_h1 = lm.transitions
    T_h0 = lm.expected_t
    "**** LISA probability matrix ****"
    print(lm.p)
    "**** LISA steady state ****"
    print(lm.steady_state)
    # print(giddy.ergodic.fmpt(lm.p))
    "**** LISA chi squared test of homogeneity ****"
    print(lm.chi_2)
    return(lm.p)


# Plot LISA probability matrices for all four data sets
def plot_LISA():
    P = []
    names = ["SM.csv", "SW.csv", "DM.csv", "DW.csv"]
    for i in names:
        P.append(lisa_markov_analysis(i))
    sns.set()
    fig, axes = plt.subplots(2,2,figsize = (7,7))
    for i in range(2):
        for j in range(2):
            ax = axes[i,j]
            if i==0 and j==0:
                p_temp = P[0]
                im = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0, vmax=1,square=True, cmap="Reds",fmt='.3f')
                ax.set_title("SM",fontsize=10)
            else:
                p_temp = P[i*2+j]
                im = sns.heatmap(p_temp, annot=True, linewidths=.5, ax=ax, cbar=True, vmin=0, vmax=1,square=True, cmap="Reds",fmt='.3f')
                ax.set_title(names[i*2+j][:2],fontsize=10)
    plt.show()



# The functions below produce the results presented in the three analysis of our paper

classic_markov_analysis(FILE_NAME)
spatial_markov_analysis(FILE_NAME)
lisa_markov_analysis(FILE_NAME)
plot_LISA()
