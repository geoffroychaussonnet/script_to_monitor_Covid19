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

def get_trend(dates,evol1,startMonth,startDay,deltaDays):
    startDate = datetime.date(2020, startMonth,startDay)
    endDate = startDate + dt.timedelta(days=9)
    bcorrelDate = (dates>=startDate) * (dates<=endDate)
    correlDate = dates[bcorrelDate]
    x = np.arange(sum(bcorrelDate))
    y = evol1[bcorrelDate]
    p1=polyfit(x,log(y),1)
    print(p1)
    #xpredict = np.arange(sum(bcorrelDate)+deltaDays)
    xpredict = np.arange(sum(bcorrelDate)-1,sum(bcorrelDate)+deltaDays)
    pol = exp(polyval(p1,xpredict))
    correl1 = pol
    xcorrel1 = dateax[bcorrelDate].tolist()
    xcorrel1 = [dateax[bcorrelDate][-1]]
    correlDate = correlDate.tolist()
    for i in range(deltaDays):
        correlDate.append(correlDate[-1] + dt.timedelta(days=1))
        xcorrel1.append(dateOut(correlDate[-1]))

    return xcorrel1, correl1

def dateOut(date):
    return date.strftime('%m/%d/%y').lstrip("0").replace("/0", "/")
def dateIn(strDate):
    spl = strDate.split('/')
    month = int(spl[0])
    day = int(spl[1])
    year = int("20%s" %spl[2])
    return datetime.date(year, month,day)

def plot_country(strCountry,filterDate,dates,quarParam,ax):
    quarDate = quarParam

    # Extract evolution for this country
    evol1 = evolution_country(strCountry,filterDate)

    # find the quarantine date 
    iQuar = np.where(dateax==quarDate)

    # Define the period for the trend
    if sum(iQuar) > 0: # Quarantine found
        quarDateIn = dateIn(quarDate)
    else:
        quarDateIn = dt.date.today()
    correlDateStart = quarDateIn - dt.timedelta(days=9)
    print("Start correl date before quanratine (or today) :", dateOut(correlDateStart))
    startMonth = correlDateStart.month
    startDay = correlDateStart.day
    deltaDays = 14

    # Get the trend
    xcorrel1, correl1 = get_trend(dates,evol1,startMonth,startDay,deltaDays)
    p = ax.semilogy(dateax,evol1,ls='-',lw=4.0,label=strCountry)
    col = p[0].get_color()
    ax.semilogy(xcorrel1,correl1,ls='-',lw=2.0,c=col)
    # Plot the quarantine date
    ax.scatter(dateax[iQuar[0]],evol1[iQuar[0]],c=col,s=300,marker="X")

###################### Function ############################

data = pd.read_csv(path+"time_series_19-covid-Deaths.csv")
dateax=data.columns[4:].values.astype(str)
Ndays = len(dateax)
# Convert date axis to date vector
dates = np.array([dt.datetime.strptime(plof,'%m/%d/%y').date() for plof in dateax])

# Filter axe of dates
filterDate = (dates>=startDate)
dateax = dateax[filterDate]
dates = dates[filterDate]

close(1)
fig = figure(1)
ax = fig.add_subplot(111)

plot_country("France",filterDate,dates,'3/17/20',ax)
plot_country("Germany",filterDate,dates,'3/19/20',ax)
plot_country("Italy",filterDate,dates,'3/9/20',ax)
plot_country("Spain",filterDate,dates,'3/14/20',ax)
plot_country("United Kingdom",filterDate,dates,'3/22/20',ax)
plot_country("US",filterDate,dates,'3/22/20',ax)
#plot_country("Korea, South",filterDate,dates,'3/22/20',ax)

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
fname = "../evolCovid19_Europe_Death_Prediction%idays_stand%s.pdf" %(deltaDays,strDateToday)
savefig(fname)
show()
