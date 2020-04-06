# Author: Geoffroy Chaussonnet, Karlsruhe Institute of Technology, Germany
# Script to load the data related to the COVID-19 gathered by the Johns Hopkins University and to plot the evolution per country, versus time. Observables are confirmed cases, deaths, active cases or death rate. The observables can be cumulative, or daily (first time derivative), or even the curvature (very noisy curve)
# Source of the data: https://github.com/CSSEGISandData/COVID-19

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt
from scipy.signal import savgol_filter
from pathlib import Path
from covid_utils import *

############### Basic use #############################
# example: plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
# Argument 1: name of the country, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: display parameters, AUTOMATICALLY generated
# Argument 4: fitting parameters, AUTOMATICALLY generated
# Argument 5: date of the confinement start. If no confinement, enter a date in the future
# Argument 6: axis (matplotlib object) where to plot the curves

################ Parameters to define manually (BEGIN) ######################
# Path to the folder containing the time series:
path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
figures_path = "../FIGURES"
daysInterval = 7   # To set Major x-axis
startDate = datetime.date(2020, 3,1)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

#yscale = 'linear'
yscale = 'log'

# Type of value to analyse:
field = "Confirmed"
#field = "Deaths"
#field = "Active"
#field = "DeathRate"

# Type of evolution:
#evolutionType = "cumulative"
evolutionType = "daily"
#evolutionType = "curvature"
#evolutionType = "smoothedCurvature"
#evolutionType = "R0"  # (Experimental)

# Extrapolate data before and after lockdown (0=no, 1=yes) 
iExtrapol = 1

# Smoothing: (set window size to 0 to deactivate)
vSmoothing = [7,3]  # [window size,order of fitting polynomial]

# Type of zones (see in the execution section)
#zone = "continents"
zone = "countries"
################ Parameters to define manually (END) ######################

######################## Definition of Functions (BEGIN) ############################

def get_trend(dates,evol1,fitParam,extParam):
    dtFitBeg = fitParam[0]
    dtFitEnd = fitParam[1]
    dtExtBeg = extParam[0]
    dtExtEnd = extParam[1]
    print("Time windows for fitting: ", dateOut(dtFitBeg), " - ", dateOut(dtFitEnd))
    print("Time windows for extrapo: ", dateOut(dtExtBeg), " - ", dateOut(dtExtEnd))
    bfitDate = (dates>=dtFitBeg) * (dates<=dtFitEnd)
    fitDate = dates[bfitDate]
    Ndfit = sum(bfitDate)
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

def dateOut(date):
    return date.strftime('%m/%d/%y').lstrip("0").replace("/0", "/")
def dateIn(strDate):
    spl = strDate.split('/')
    month = int(spl[0])
    day = int(spl[1])
    year = int("20%s" %spl[2])
    return datetime.date(year, month,day)

def plot_country(strCountry,dataParam,displayParam,fitParam,quarParam,ax):
    print("########## Treating country: %12s ###########" %strCountry)
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]
    iExtrapol = fitParam[2]

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,dataParam,displayParam)

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

    if dataParam['Smoothing'][0] != 0:
        evol1 = savgol_filter(evol1, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 

    if displayParam['YScale'] == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
    p = ax.semilogy(dataParam['DateAxis'],evol1,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(dataParam['DateAxis'][iQuar][0],evol1[iQuar][0],c=col,s=300,marker="X")

    if (iExtrapol==0): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam1,extParam1)
    ax.semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam2,extParam2)
        ax.semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

def setDisplayParam(field,evolutionType,yscale,zone):
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

    txtTitle = "%s %s\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s %s" %(txtEvol,txtField,strUnit)
    displayParam['Field'] = field
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    Path(figures_path).mkdir(parents=True, exist_ok=True)
    fname = figures_path + "/%s_evolCovid19_%s_%s_for_%s.png" %(strDateToday,txtEvol,txtField,zone)
    displayParam['FileName'] = fname.replace(" ","_")
    displayParam['YScale'] = yscale
    return displayParam

def setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, iExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="Active":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, iExtrapol]

######################## Definition of Functions (END) ############################


############## Execution section ################

# Initialisation
dataParam = loadData(path,field,evolutionType,vSmoothing,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType,yscale,zone)
fitParam = setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol)

# Set graphic objects
close(1)
fig = figure(num=1,figsize=(10,6))
ax = fig.add_subplot(111)

#plot_country("World",dataParam,displayParam,fitParam,'3/22/21',ax)
if zone == "continents":
    plot_country("EU",dataParam,displayParam,fitParam,'3/22/21',ax)
    #plot_country("European continent",dataParam,displayParam,fitParam,'3/22/21',ax)
    plot_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
elif zone == "countries":
    #plot_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
    plot_country("Italy",dataParam,displayParam,fitParam,'3/9/20',ax)
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
    #plot_country("Belgium",dataParam,displayParam,fitParam,'3/18/20',ax)
    #plot_country("Ireland",dataParam,displayParam,fitParam,'3/28/20',ax)

# Add graph decorations
if dataParam['EvolutionType'] == "R0": ax.axhline(1)
ax.set_title(displayParam['title'])
ax.set_yscale(displayParam['YScale'])
ax.set_xlabel("Date")
ax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.set_ylabel(displayParam['YaxisLabel'])
ax.legend(loc=2)
ax.grid(which='major',color='grey', linestyle='-', linewidth=1)
ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

fig.tight_layout()
savefig(displayParam['FileName'],dpi=600,bbox='tight')
show()
