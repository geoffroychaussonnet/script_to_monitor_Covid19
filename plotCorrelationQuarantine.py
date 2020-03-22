# Author: Geoffroy Chaussonnet

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt

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
startDate = datetime.date(2020, 1,22)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

#field = "Confirmed"
#field = "Deaths"
field = "Active"

evolutionType = "cumulative"
#evolutionType = "daily"
################ Parameters to define manually ######################



######################## Definition of Functions ############################

def evolution_single(strCountry,data):
    icountry = data[data["Country/Region"]==strCountry].index.values
    size=len(data.iloc[icountry[0]].values[4:])
    evolution = zeros(size,dtype=int)
    for ic in icountry:
        evolution[:] += data.iloc[ic].values[4:].astype(int)
    return evolution

def evolution_country(strCountry,dataParam):

    if field=="Confirmed":
        evolution = evolution_single(strCountry,dataParam['Confirmed'])
    elif field=="Deaths":
        evolution = evolution_single(strCountry,dataParam['Deaths'])
    elif field=="Active":
        evolC = evolution_single(strCountry,dataParam['Confirmed'])
        evolD = evolution_single(strCountry,dataParam['Deaths'])
        evolR = evolution_single(strCountry,dataParam['Recovered'])
        evolution = evolC - evolR - evolD

    if dataParam['EvolutionType'] == "cumulative":
        return evolution[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "daily":
        dedt = np.zeros(len(evolution))
        dedt[1:] = (np.roll(evolution,-1) - evolution)[:-1]
        return dedt[dataParam['FilterDate']]


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
    p1=polyfit(xfit,log(yfit),1)
    print(p1)
    yext = exp(polyval(p1,xext))
    correl1 = yext

    xcorrel1 = []
    for i in range(Ndext):
        xcorrel1.append(dateOut(dtExtBeg + dt.timedelta(days=i)))

    return xcorrel1, correl1

def dateOut(date):
    return date.strftime('%m/%d/%y').lstrip("0").replace("/0", "/")
def dateIn(strDate):
    spl = strDate.split('/')
    month = int(spl[0])
    day = int(spl[1])
    year = int("20%s" %spl[2])
    return datetime.date(year, month,day)

def plot_country(strCountry,dataParam,fitParam,quarParam,ax):
    print("########## Treating country: ", strCountry, " #############")
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,dataParam)

    # find the quarantine date 
    iQuar = np.where(dataParam['DateAxis']==quarDate)

    fitParam1 = []
    extParam1 = []
    # Define the period for the trend
    if sum(iQuar) > 0: # Quarantine found
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

    p = ax.semilogy(dataParam['DateAxis'],evol1,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()

    # Get the trend
    xcorrel1, correl1 = get_trend(dataParam['Dates'],evol1,fitParam1,extParam1)
    ax.semilogy(xcorrel1,correl1,ls='--',lw=2.0,c=col)
    
    if sum(iQuar) > 0: # Quarantine found
        xcorrel2, correl2 = get_trend(dataParam['Dates'],evol1,fitParam2,extParam2)
        ax.semilogy(xcorrel2,correl2,ls='-',lw=2.0,c=col)

    # Plot the quarantine date
    ax.scatter(dataParam['DateAxis'][iQuar[0]],evol1[iQuar[0]],c=col,s=300,marker="X")

def setDisplayParam(field,evolutionType):
    displayParam = {}
    if field=="Confirmed":
        txtField = "confirmed cases"
    elif field=="Deaths":
        txtField = "deaths"
    elif field=="Active":
        txtField = "active cases"

    if evolutionType == 'cumulative':
        txtEvol = "Cumulative"
    elif evolutionType == 'daily':
        txtEvol = 'Daily'

    txtTitle = "%s %s in some Western countries\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s [-]" %(txtEvol,txtField)
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    fname = "../evolCovid19_%s_%s_stand%s.pdf" %(txtEvol,txtField,strDateToday)
    displayParam['FileName'] = fname
    return displayParam

def loadData(path,field,evolutionType,startDate=datetime.date(2020, 1,1)):
    dataParam = {}
    dataParam['Confirmed'] = pd.read_csv(path+"time_series_19-covid-Confirmed.csv")
    dataParam['Deaths'] = pd.read_csv(path+"time_series_19-covid-Deaths.csv")
    dataParam['Recovered'] = pd.read_csv(path+"time_series_19-covid-Recovered.csv")
    dataParam['Field'] = field
    dataParam['EvolutionType'] = evolutionType
    dateax = dataParam['Confirmed'].columns[4:].values.astype(str)

    # Convert date axis to date vector
    dates = np.array([dt.datetime.strptime(plof,'%m/%d/%y').date() for plof in dateax])

    # Filter axe of dates
    filterDate = (dates>=startDate)
    dateax = dateax[filterDate]
    dates = dates[filterDate]

    dataParam['FilterDate'] = filterDate
    dataParam['DateAxis'] = dateax
    dataParam['Dates'] = dates

    return dataParam

######################## Definition of Functions ############################

dataParam = loadData(path,field,evolutionType,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType)

fitParam = [fittingPeriod, extrapolPeriod]

close(1)
fig = figure(1)
ax = fig.add_subplot(111)


plot_country("France",dataParam,fitParam,'3/17/20',ax)
plot_country("Germany",dataParam,fitParam,'3/19/20',ax)
plot_country("Italy",dataParam,fitParam,'3/9/20',ax)
plot_country("Spain",dataParam,fitParam,'3/14/20',ax)
#plot_country("United Kingdom",dataParam,fitParam,'5/22/20',ax)
#plot_country("US",dataParam,fitParam,'5/22/20',ax)
plot_country("China",dataParam,fitParam,'1/22/20',ax)
plot_country("Korea, South",dataParam,fitParam,'5/22/20',ax)

ax.set_title(displayParam['title'])
#ax.set_yscale('linear')
ax.set_xlabel("Date")
ax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.set_ylabel(displayParam['YaxisLabel'])
ax.legend(loc=0)
ax.grid(which='major',color='grey', linestyle='-', linewidth=1)
ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

savefig(displayParam['FileName'])
show()
