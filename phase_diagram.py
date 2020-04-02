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
daysInterval = 7   # To set Major x-axis
startDate = datetime.date(2020, 1,1)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

yscale = 'linear'
#yscale = 'log'

field = "Confirmed"
#field = "Deaths"
#field = "Active"
#field = "DeathRate"

evolutionType = "cumulative"
#evolutionType = "daily"
#evolutionType = "curvature"
#evolutionType = "smoothedCurvature"
#evolutionType = "R0"

iExtrapol = 0

vSmoothing = [9,3]  # [window size,order of fitting polynomial]

zone = "continents"
#zone = "countries"

################ Parameters to define manually ######################


######################## Definition of Functions ############################

def evolution_single(strCountry,data):

    size=len(data.iloc[0].values[4:])
    evolution = zeros(size,dtype=int)

    lstCountry = [strCountry]
    if strCountry == "EUW":
        lstCountry = ["France", "Germany", "Spain", "Italy", "Netherlands", "Portugal", "Belgium", "Sweden", "Finland", "Greece", "Ireland", "United Kingdom", "Norway","Switzerland", "Poland", "Andorra","Luxembourg", "Liechtenstein", "Malta", "San Marino", "Holy See","Monaco"]
    elif strCountry == "EU":
        lstCountry = ["France", "Germany", "Spain", "Italy", "Netherlands", "Portugal", "Belgium", "Sweden", "Finland", "Greece", "Ireland", "Poland", "Luxembourg", "Malta","Slovenia", "Austria", "Croatia", "Hungary", "Czechia", "Slovakia", "Hungary", "Romania", "Bulgaria", "Cyprus", "Lithuania","Latvia","Estonia"]
    elif strCountry == "European continent":
        lstCountry = ["France", "Germany", "Spain", "Italy", "Netherlands", "Portugal", "Belgium", "Sweden", "Finland", "Greece", "Ireland", "United Kingdom", "Norway","Switzerland", "Poland", "Andorra","Luxembourg", "Liechtenstein", "Malta", "San Marino", "Holy See","Monaco","Hungary", "Czechia","Slovakia", "Slovenia", "Croatia","Bosnia and Herzegovina", "Serbia", "Albania", "Romania", "Bulgaria", "Ukraine", "Belarus", "Latvia", "Estonia", "Lithuania","Moldova","North Macedonia", "Kosovo","Montenegro","Iceland","Cyprus"]

    for ic,cntry in enumerate(data['Country/Region']):
        if (cntry in lstCountry) or (strCountry=="World"):
            locRegion = data.iloc[ic].values[4:]
            locRegion[isnan(locRegion.tolist())] = 0
            evolution[:] += locRegion.astype(int)

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
    elif field=="DeathRate":
        evolC = evolution_single(strCountry,dataParam['Confirmed'])
        evolD = evolution_single(strCountry,dataParam['Deaths'])
        evolution = evolD/evolC*100

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
        #dedt = np.zeros(len(evolution))
        #dedt[1:] = np.diff(evolution)
        dedt = np.diff(evolution)
        evol = savgol_filter(dedt, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        d2edt2 = np.zeros(len(evolution))
        d2edt2[2:] = np.diff(evol)#/evol[-1]
        evol = d2edt2[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "R0":
        R0 = np.zeros(len(evolution))
        delta0 = np.diff(evolution)
        delta = savgol_filter(delta0, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        R0[1:] = delta/np.roll(delta,5)
        ## TEST:
        #print((delta/np.roll(delta,5))[-10:])
        #print((delta*5/np.roll(delta,1))[-10:])
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
    #Ndfit = (dtFitEnd - dtFitBeg).days + 1
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

def plot_phase_country(strCountry,dataParam,displayParam,fitParam,quarParam,ax):
    print("########## Treating country: %12s ###########" %strCountry)
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]
    iExtrapol = fitParam[2]

    # Extract evolution for this country
    dataParam['EvolutionType'] = "smoothedCurvature"
    curvature = evolution_country(strCountry,dataParam)
    dataParam['EvolutionType'] = "daily"
    gradient = evolution_country(strCountry,dataParam)

    # Filter data for more than 100 cases
    dataParam['EvolutionType'] = "cumulative"
    cumul = evolution_country(strCountry,dataParam)

    gooddata = (cumul>100)
    cumul = cumul[gooddata]
    gradient = gradient[gooddata]
    curvature = curvature[gooddata]
    locDates = dataParam['Dates'][gooddata]
    

    # find the quarantine date 
    iQuar = locDates>=dateIn(quarDate)

    scurv = curvature
    sgrad = gradient
    if dataParam['Smoothing'][0] != 0:
        scurv = savgol_filter(curvature, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        sgrad = savgol_filter(gradient, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 

    if displayParam['YScale'] == 'log':
        scurv = np.ma.masked_where(scurv<=0,scurv)
    p = ax.semilogy(sgrad,scurv,ls='-',marker='o',lw=4.0,label=strCountry)
    col = p[0].get_color()
    ax.scatter(sgrad[-1],scurv[-1],c=col,s=100,marker="s")

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(sgrad[iQuar][0],scurv[iQuar][0],c=col,s=300,marker="X")

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

    txtEvol = "Phase portrait from"

#    if evolutionType == 'cumulative':
#        txtEvol = "Cumulative"
#    elif evolutionType == 'daily':
#        txtEvol = 'Daily'
#    elif evolutionType == 'curvature':
#        txtEvol = 'Derivative of daily'
#    elif evolutionType == 'smoothedCurvature':
#        txtEvol = 'Derivative of smoothed daily'
#    elif evolutionType == 'R0':
#        txtEvol = 'R0 from'

    txtTitle = "%s %s\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s %s" %(txtEvol,txtField,strUnit)
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    fname = "../FIGURES/%s_evolCovid19_%s_%s_for_%s.png" %(strDateToday,txtEvol,txtField,zone)
    displayParam['FileName'] = fname.replace(" ","_")
    displayParam['YScale'] = yscale
    return displayParam

def loadData(path,field,evolutionType,vSmoothing,startDate=datetime.date(2020, 1,1)):
    dataParam = {}
    dataParam['Confirmed'] = pd.read_csv(path+"time_series_covid19_confirmed_global.csv")
    dataParam['Deaths'] = pd.read_csv(path+"time_series_covid19_deaths_global.csv")
    dataParam['Recovered'] = pd.read_csv(path+"time_series_covid19_recovered_global.csv")
    dataParam['Field'] = field
    dataParam['EvolutionType'] = evolutionType
    dataParam['Smoothing'] = vSmoothing
    #dateax = dataParam['Confirmed'].columns[4:].values.astype(str)
    dateax = dataParam['Deaths'].columns[4:].values.astype(str)

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

def setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, iExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="Active":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, iExtrapol]

######################## Definition of Functions ############################

dataParam = loadData(path,field,evolutionType,vSmoothing,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType,yscale,zone)
fitParam = setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol)

close(1)
fig = figure(num=1,figsize=(10,6))
ax = fig.add_subplot(111)

i=0

#plot_country("World",dataParam,displayParam,fitParam,'3/22/21',ax)
if zone == "continents":
    plot_phase_country("EU",dataParam,displayParam,fitParam,'3/22/21',ax)
    #plot_phase_country("European continent",dataParam,displayParam,fitParam,'3/22/21',ax)
    plot_phase_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    plot_phase_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
elif zone == "countries":
    #plot_phase_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    #plot_phase_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
    plot_phase_country("Italy",dataParam,displayParam,fitParam,'3/9/20',ax)
    plot_phase_country("Spain",dataParam,displayParam,fitParam,'3/14/20',ax)
    plot_phase_country("Germany",dataParam,displayParam,fitParam,'3/19/20',ax)
    plot_phase_country("France",dataParam,displayParam,fitParam,'3/17/20',ax)
    #plot_phase_country("Iran",dataParam,displayParam,fitParam,'8/17/20',ax)
    plot_phase_country("Korea, South",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Japan",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Switzerland",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("United Kingdom",dataParam,displayParam,fitParam,'3/22/20',ax)
    #plot_phase_country("Denmark",dataParam,displayParam,fitParam,'3/13/20',ax)
    #plot_phase_country("Norway",dataParam,displayParam,fitParam,'3/12/20',ax)
    #plot_phase_country("Sweden",dataParam,displayParam,fitParam,'3/28/20',ax)
    #plot_phase_country("Finland",dataParam,displayParam,fitParam,'3/19/20',ax)
    #plot_phase_country("Canada",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Belgium",dataParam,displayParam,fitParam,'3/18/20',ax)
    #plot_phase_country("Ireland",dataParam,displayParam,fitParam,'3/28/20',ax)

if dataParam['EvolutionType'] == "R0": ax.axhline(1)
ax.set_title(displayParam['title'])
ax.set_yscale(displayParam['YScale'])
#ax.set_xscale('log')
ax.set_xlabel(r'Gradient [new case/day]')
ax.set_ylabel(r'Curvature [new case/day$^2$]')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.legend(loc=2)
ax.grid(which='major',color='grey', linestyle='-', linewidth=1)
ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

fig.tight_layout()
savefig(displayParam['FileName'],dpi=600,bbox='tight')
show()
