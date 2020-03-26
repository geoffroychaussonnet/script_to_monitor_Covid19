# Author: Geoffroy Chaussonnet

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt
from scipy.signal import savgol_filter

############### Basic use #############################
# example: plot_country("France",dataParam,fitParam,'3/17/20',ax)
# Argument 1: string, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: fitting parameters, AUTOMATICALLY generated
# Argument 4: date of the confinement start. If no confinement, enter a date in the future
# Argument 5: axis (matplotlib object) where to plot the curves


################ Parameters to define manually ######################
# Path to the folder containing the time series:
path="../csse_covid_19_data/csse_covid_19_time_series/"
daysInterval = 7   # To set Major x-axis:
startDate = datetime.date(2020, 2,22)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

#yscale = 'linear'
yscale = 'log'

#field = "Confirmed"
field = "Deaths"
#field = "DeathRate"

#evolutionType = "cumulative"
evolutionType = "daily"

bExtrapol = True
#bExtrapol = False

bSmoothing = True
#bSmoothing = False
################ Parameters to define manually ######################


######################## Definition of Functions ############################

def evolution_single(strCountry,data):
    if strCountry == "World":
        icountry = range(len(data))
    else:
        icountry = data[data["Country/Region"]==strCountry].index.values
    size=len(data.iloc[icountry[0]].values[4:])
    evolution = zeros(size,dtype=int)
    for ic in icountry:
        locRegion = data.iloc[ic].values[4:]
        locRegion[isnan(locRegion.tolist())] = 0
        evolution[:] += locRegion.astype(int)
    return evolution

def evolution_country(strCountry,dataParam):

    if field=="Confirmed":
        evolution = evolution_single(strCountry,dataParam['Confirmed'])
    elif field=="Deaths":
        evolution = evolution_single(strCountry,dataParam['Deaths'])
    elif field=="DeathRate":
        evolC = evolution_single(strCountry,dataParam['Confirmed'])
        evolD = evolution_single(strCountry,dataParam['Deaths'])
        evolution = evolD/evolC*100

    if dataParam['EvolutionType'] == "cumulative":
        evol =  evolution[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "daily":
        dedt = np.zeros(len(evolution))
        dedt[1:] = (np.roll(evolution,-1) - evolution)[:-1]
        evol = dedt[dataParam['FilterDate']]

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

def dateOut(date):
    return date.strftime('%m/%d/%y').lstrip("0").replace("/0", "/")
def dateIn(strDate):
    spl = strDate.split('/')
    month = int(spl[0])
    day = int(spl[1])
    year = int("20%s" %spl[2])
    return datetime.date(year, month,day)

def plot_country(strCountry,dataParam,displayParam,fitParam,quarParam,ax):
    print("########## Treating country: ", strCountry, " #############")
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]
    bExtrapol = fitParam[2]

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,dataParam)

    # find the quarantine date 
    iQuar = np.where(dataParam['Dates']>=dateIn(quarDate))
    iQuar = dataParam['Dates']>=dateIn(quarDate)
    #print("iQuar:", iQuar)
    #print("quarDAte:", quarDate)
    #print("dataParam['DateAxis']:", dataParam['DateAxis'])

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

    if dataParam['Smoothing']:
        evol1 = savgol_filter(evol1, 7, 3) # window size 51, polynomial order 3

    if displayParam['YScale'] == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
    p = ax.semilogy(dataParam['DateAxis'],evol1,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(dataParam['DateAxis'][iQuar][0],evol1[iQuar][0],c=col,s=300,marker="X")

    if not(bExtrapol): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam1,extParam1)
    ax.semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(dataParam['Dates'],evol1,fitParam2,extParam2)
        ax.semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

def setDisplayParam(field,evolutionType,yscale):
    displayParam = {}

    strUnit = "[-]"
    if field=="Confirmed":
        txtField = "confirmed cases"
    elif field=="Deaths":
        txtField = "deaths"
    elif field=="DeathRate":
        txtField = "death rate"
        strUnit = "[%]"

    if evolutionType == 'cumulative':
        txtEvol = "Cumulative"
    elif evolutionType == 'daily':
        txtEvol = 'Daily'

    txtTitle = "%s %s in some Western countries\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s %s" %(txtEvol,txtField,strUnit)
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    fname = "../FIGURES/%s_evolCovid19_%s_%s.png" %(strDateToday,txtEvol,txtField)
    displayParam['FileName'] = fname.replace(" ","_")
    displayParam['YScale'] = yscale
    return displayParam

def loadData(path,field,evolutionType,bSmoothing,startDate=datetime.date(2020, 1,1)):
    dataParam = {}
    dataParam['Confirmed'] = pd.read_csv(path+"time_series_covid19_confirmed_global.csv")
    dataParam['Deaths'] = pd.read_csv(path+"time_series_covid19_deaths_global.csv")
    dataParam['Field'] = field
    dataParam['EvolutionType'] = evolutionType
    dataParam['Smoothing'] = bSmoothing
    #dateax = dataParam['Confirmed'].columns[4:].values.astype(str)
    dateax = dataParam['Deaths'].columns[4:].values.astype(str)

    # Convert date axis to date vector
    dates = np.array([dt.datetime.strptime(plof,'%m/%d/%y').date() for plof in dateax])
    #dates = []
    #for plof in dateax: 
    #    print("YOP: ", plof)
    #    yop = dt.datetime.strptime(plof,'%m/%d/%y').date()
    #    dates.append(yop)

    # Filter axe of dates
    filterDate = (dates>=startDate)
    dateax = dateax[filterDate]
    dates = dates[filterDate]

    dataParam['FilterDate'] = filterDate
    dataParam['DateAxis'] = dateax
    dataParam['Dates'] = dates

    return dataParam

def setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,bExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, bExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, bExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, bExtrapol]

######################## Definition of Functions ############################

dataParam = loadData(path,field,evolutionType,bSmoothing,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType,yscale)
fitParam = setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,bExtrapol)

close(1)
fig = figure(num=1,figsize=(10,6))
ax = fig.add_subplot(111)

plot_country("World",dataParam,displayParam,fitParam,'3/22/21',ax)
#plot_country("China",dataParam,displayParam,fitParam,'1/22/20',ax)
plot_country("Italy",dataParam,displayParam,fitParam,'3/9/20',ax)
plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
plot_country("Spain",dataParam,displayParam,fitParam,'3/14/20',ax)
plot_country("Germany",dataParam,displayParam,fitParam,'3/19/20',ax)
#plot_country("Iran",dataParam,displayParam,fitParam,'8/17/20',ax)
plot_country("France",dataParam,displayParam,fitParam,'3/17/20',ax)
#plot_country("Korea, South",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Switzerland",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("United Kingdom",dataParam,displayParam,fitParam,'3/22/20',ax)
#plot_country("Norway",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Sweden",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Finland",dataParam,displayParam,fitParam,'5/22/20',ax)
#plot_country("Canada",dataParam,displayParam,fitParam,'5/22/20',ax)

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
