# Author: Geoffroy Chaussonnet

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt
from scipy.signal import savgol_filter
from pathlib import Path
from covid_utils import *

############### Basic use #############################
# example: plot_country("France",dataParam,fitParam,'3/17/20',ax)
# Argument 1: string, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: fitting parameters, AUTOMATICALLY generated
# Argument 4: date of the confinement start. If no confinement, enter a date in the future
# Argument 5: axis (matplotlib object) where to plot the curves


################ Parameters to define manually ######################
# Path to the folder containing the time series:
path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
figures_path = "../FIGURES"
daysInterval = 7   # To set Major x-axis
startDate = datetime.date(2020, 3,1)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

#yscale = 'linear'
yscale = 'log'

#field = "Confirmed"
field = "Deaths"
#field = "Active"
#field = "DeathRate"

#evolutionType = "cumulative"
evolutionType = "daily"
#evolutionType = "curvature"
#evolutionType = "smoothedCurvature"
#evolutionType = "R0"

iExtrapol = 0

vSmoothing = [5,3]  # [window size,order of fitting polynomial]
################ Parameters to define manually ######################


######################## Definition of Functions ############################

def evolution_country(strCountry,dataParam):

    evolution = evolution_country_aux(field, strCountry, dataParam)

    if dataParam['EvolutionType'] == "cumulative":
        evol =  evolution[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "daily":
        dedt = np.zeros(len(evolution))
        dedt[1:] = np.diff(evolution)
        evol = dedt[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "curvature":
        d2edt2 = np.zeros(len(evolution))
        d2edt2[2:] = np.diff(evolution,2)
        evol = d2edt2[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "smoothedCurvature":
        dedt = np.diff(evolution)
        evol = savgol_filter(dedt, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order
        d2edt2 = np.zeros(len(evolution))
        d2edt2[2:] = np.diff(evol)/evol[-1]
        evol = d2edt2[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "R0":
        R0 = np.zeros(len(evolution))
        delta0 = np.diff(evolution)
        delta = savgol_filter(delta0, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order
        R0[1:] = delta/np.roll(delta,5)
        evol = R0[dataParam['FilterDate']]

    return evol


def get_trend(dates,evol1,fitParam,extParam):
    dtFitBeg = fitParam[0]
    dtFitEnd = fitParam[1]
    dtExtBeg = extParam[0]
    dtExtEnd = extParam[1]
    print("Time windows for fitting: ", dateOut(dtFitBeg), " - ", dateOut(dtFitEnd))
    print("Time windows for extrapo: ", dateOut(dtExtBeg), " - ", dateOut(dtExtEnd))
    bfitDate = (dates>=dtFitBeg) * (dates<=dtFitEnd)
    fitDate = dates[bfitDate]
    Ndfit = (dtFitEnd - dtFitBeg).days + 1
    Ndext = (dtExtEnd - dtExtBeg).days + 1
    Ndtot = (dtExtEnd - dtFitBeg).days + 1
    xfit = np.arange(Ndfit)
    xext = np.arange(Ndtot-Ndext,Ndtot)

    yfit = evol1[bfitDate]
    nz = (yfit>0)
    if sum(nz)>0:
        p1=polyfit(xfit[nz],log(yfit[nz]),1)
        yext = exp(polyval(p1,xext))
    else:
        p1=polyfit(xfit,log(-yfit),1)
        yext = exp(-polyval(p1,xext))
    print(p1)
    correl1 = yext

    xcorrel1 = []
    for i in range(Ndext):
        xcorrel1.append(dateOut(dtExtBeg + dt.timedelta(days=i)))

    rate=correl1[-1]/correl1[-2]-1
    if rate>0: strRate='+%.1f%%' %(rate*100)
    else:strRate='%.1f%%' %(rate*100)

    return xcorrel1, correl1, strRate


def plot_country(strCountry,dataParam,displayParam,fitParam,quarParam,ax):
    print("########## Treating country: %12s ###########" %strCountry)
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]
    iExtrapol = fitParam[2]

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,dataParam)

    # find the quarantine date 
    iQuar = np.where(dataParam['Dates']>=dateIn(quarDate))
    iQuar = dataParam['Dates']>=dateIn(quarDate)

    fitParam1 = []
    extParam1 = []
    # Define the period for the trend
    if sum(iQuar) > 3: # Quarantine found
        dtFitEnd = dateIn(quarDate)

        fitParam2 = []
        extParam2 = []
        fitParam2.append(dateIn(quarDate))
        fitParam2.append(dt.date.today() - dt.timedelta(days=1))
        extParam2.append(dtFitEnd)
        extParam2.append(dtFitEnd + dt.timedelta(days=extrapolPeriod+1))
    else:
        dtFitEnd = dt.date.today() - dt.timedelta(days=1)
    dtFitBeg = dtFitEnd - dt.timedelta(days=fittingPeriod+1)
    dtExtBeg = dtFitEnd
    dtExtEnd = dtExtBeg + dt.timedelta(days=extrapolPeriod+1)
    fitParam1.append(dtFitBeg)
    fitParam1.append(dtFitEnd)
    extParam1.append(dtExtBeg)
    extParam1.append(dtExtEnd)

    evol1s = savgol_filter(evol1, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 

    if displayParam['YScale'] == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
        evol1s = np.ma.masked_where(evol1s<=0,evol1s)
    p = ax[0].semilogy(dataParam['DateAxis'],evol1s,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()
    ax[1].semilogy(dataParam['DateAxis'],evol1,ls='-',lw=4.0,label=strCountry)

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax[0].scatter(dataParam['DateAxis'][iQuar][0],evol1[iQuar][0],c=col,s=300,marker="X")
        ax[1].scatter(dataParam['DateAxis'][iQuar][0],evol1[iQuar][0],c=col,s=300,marker="X")

    if (iExtrapol==0): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam1,extParam1)
    ax[0].semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax[1].semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax[0].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')
    ax[1].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam2,extParam2)
        ax[0].semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax[0].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')
        ax[1].semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax[1].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

def setDisplayParam(field,evolutionType,yscale):
    displayParam = {}

    strUnit = "[-]"
    if field=="Confirmed":
        txtField = "confirmed cases"
    elif field=="Deaths":
        txtField = "deaths"
    elif field=="Active":
        txtField = "active cases"
    elif field=="DeathRate":
        txtField = "death rate"
        strUnit = "[%]"

    if evolutionType == 'cumulative':
        txtEvol = "Cumulative"
    elif evolutionType == 'daily':
        txtEvol = 'Daily'
    elif evolutionType == 'curvature':
        txtEvol = 'Derivative of daily'
    elif evolutionType == 'smoothedCurvature':
        txtEvol = 'Derivative of smoothed daily'
    elif evolutionType == 'R0':
        txtEvol = 'R0 from'

    txtTitle = "%s %s in some Western countries\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s %s" %(txtEvol,txtField,strUnit)
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    Path(figures_path).mkdir(parents=True, exist_ok=True)
    fname = figures_path + "/%s_evolCovid19_%s_%s_with_without_smooth.png" %(strDateToday,txtEvol,txtField)
    displayParam['FileName'] = fname.replace(" ","_")
    displayParam['YScale'] = yscale
    return displayParam


######################## Definition of Functions ############################

dataParam = loadData(path,field,evolutionType,vSmoothing,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType,yscale)
fitParam = setFitExtraParam(field,fittingPeriod, extrapolPeriod,dataParam,iExtrapol)

close(1)
fig = figure(num=1,figsize=(10,6))
ax = []
ax.append(fig.add_subplot(121))
ax.append(fig.add_subplot(122))

#plot_country("World",dataParam,displayParam,fitParam,'3/22/21',ax)
#plot_country("EU",dataParam,displayParam,fitParam,'3/22/21',ax)
#plot_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
plot_country("Italy",dataParam,displayParam,fitParam,'3/9/20',ax)
plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
plot_country("Spain",dataParam,displayParam,fitParam,'3/14/20',ax)
plot_country("Germany",dataParam,displayParam,fitParam,'3/19/20',ax)
plot_country("France",dataParam,displayParam,fitParam,'3/17/20',ax)
#plot_country("Iran",dataParam,displayParam,fitParam,'8/17/20',ax)
#plot_country("Korea, South",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Japan",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Switzerland",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("United Kingdom",dataParam,displayParam,fitParam,'3/22/20',ax)
#plot_country("Denmark",dataParam,displayParam,fitParam,'3/13/20',ax)
#plot_country("Norway",dataParam,displayParam,fitParam,'3/12/20',ax)
#plot_country("Sweden",dataParam,displayParam,fitParam,'3/28/20',ax)
#plot_country("Finland",dataParam,displayParam,fitParam,'3/19/20',ax)
#plot_country("Canada",dataParam,displayParam,fitParam,'5/22/20',ax)
plot_country("Belgium",dataParam,displayParam,fitParam,'3/18/20',ax)
#plot_country("Ireland",dataParam,displayParam,fitParam,'3/28/20',ax)

for lax in ax:
    if dataParam['EvolutionType'] == "R0": lax.axhline(1)
    lax.set_title(displayParam['title'])
    lax.set_yscale(displayParam['YScale'])
    lax.set_xlabel("Date")
    lax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
    lax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    lax.set_ylabel(displayParam['YaxisLabel'])
    lax.legend(loc=2)
    lax.grid(which='major',color='grey', linestyle='-', linewidth=1)
    lax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

fig.tight_layout()
savefig(displayParam['FileName'],dpi=600,bbox='tight')
show()
