# Author: Geoffroy Chaussonnet

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt

# Path to the folder containing the time series:
path="../csse_covid_19_data/csse_covid_19_time_series/"
# To set Major x-axis:
daysInterval = 7   
# Start date of the plot:
startDate = datetime.date(2020, 2,22)   
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

# basic use:
# example: plot_country("France",filterDate,dates,'3/17/20',ax)
# Argument 1: string, as it appears in the CSV file
# Argument 2: array of boolean, AUTOMATICALLY generated from the variable "startDate"
# Argument 3: date axis, AUTOMATICALLY generated 
# Argument 4: date of the confinement start. If no confinement, enter a date of the future
# Argument 5: axis (matplotlib object) where to plot the curves

###################### Function ############################

def evolution_country(strCountry,filterDate):
    icountry = data[data["Country/Region"]==strCountry].index.values
    size=len(data.iloc[icountry[0]].values[4:])
    evolution = zeros(size,dtype=int)
    for ic in icountry:
        evolution[:] += data.iloc[ic].values[4:].astype(int)
    return evolution[filterDate]


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

def plot_country(strCountry,filterDate,dates,fitParam,quarParam,ax):
    print("########## Treating country: ", strCountry, " #############")
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,filterDate)

    # find the quarantine date 
    iQuar = np.where(dateax==quarDate)

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

    p = ax.semilogy(dateax,evol1,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()

    # Get the trend
    xcorrel1, correl1 = get_trend(dates,evol1,fitParam1,extParam1)
    ax.semilogy(xcorrel1,correl1,ls='--',lw=2.0,c=col)
    
    if sum(iQuar) > 0: # Quarantine found
        xcorrel2, correl2 = get_trend(dates,evol1,fitParam2,extParam2)
        ax.semilogy(xcorrel2,correl2,ls='-',lw=2.0,c=col)

    # Plot the quarantine date
    ax.scatter(dateax[iQuar[0]],evol1[iQuar[0]],c=col,s=300,marker="X")

###################### Function ############################

data = pd.read_csv(path+"time_series_19-covid-Deaths.csv")
dateax=data.columns[4:].values.astype(str)
Ndays = len(dateax)
# Convert date axis to date vector
dates = np.array([dt.datetime.strptime(plof,'%m/%d/%y').date() for plof in dateax])

fitParam = [fittingPeriod, extrapolPeriod]

# Filter axe of dates
filterDate = (dates>=startDate)
dateax = dateax[filterDate]
dates = dates[filterDate]

close(1)
fig = figure(1)
ax = fig.add_subplot(111)

plot_country("France",filterDate,dates,fitParam,'3/17/20',ax)
plot_country("Germany",filterDate,dates,fitParam,'3/19/20',ax)
plot_country("Italy",filterDate,dates,fitParam,'3/9/20',ax)
plot_country("Spain",filterDate,dates,fitParam,'3/14/20',ax)
plot_country("United Kingdom",filterDate,dates,fitParam,'5/22/20',ax)
plot_country("US",filterDate,dates,fitParam,'5/22/20',ax)
#plot_country("Iran",filterDate,dates,'3/22/20',ax)

ax.set_title("Cumulated deaths in some Western countries\n (Source: Johns Hopkins University)")

#ax.set_yscale('linear')
ax.set_xlabel("Date")
ax.xaxis.set_major_locator(ticker.MultipleLocator(7))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.set_ylabel("Cumulated number of deaths [-]")
ax.legend(loc=0)
ax.grid(which='major',color='grey', linestyle='-', linewidth=1)
ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

strDateToday = dt.date.today().strftime("%Y%m%d")
fname = "../evolCovid19_Europe_Death_Prediction%idays_stand%s.pdf" %(14,strDateToday)
savefig(fname)
show()
