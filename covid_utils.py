import numpy as np
import pandas as pd
import datetime as dt
from scipy.signal import savgol_filter


def loadData(path,field,evolutionType,vSmoothing,startDate=dt.date(2020, 1,1)):
    dataParam = {}
    dataParam['Confirmed'] = pd.read_csv(path+"time_series_covid19_confirmed_global.csv")
    dataParam['Deaths'] = pd.read_csv(path+"time_series_covid19_deaths_global.csv")
    dataParam['Recovered'] = pd.read_csv(path+"time_series_covid19_recovered_global.csv")
    dataParam['Field'] = field
    dataParam['EvolutionType'] = evolutionType
    dataParam['Smoothing'] = vSmoothing
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


def evolution_single(strCountry,data):

    size=len(data.iloc[0].values[4:])
    evolution = np.zeros(size,dtype=int)

    lstCountry = [strCountry]
    if strCountry == "EU":
        lstCountry = ["France", "Germany", "Spain", "Italy", "Netherlands", "Portugal", "Belgium", "Sweden", "Finland", "Greece", "Ireland", "Poland", "Luxembourg", "Malta","Slovenia", "Austria", "Croatia", "Hungary", "Czechia", "Slovakia", "Hungary", "Romania", "Bulgaria", "Cyprus", "Lithuania","Latvia","Estonia"]
    elif strCountry == "European continent":
        lstCountry = ["France", "Germany", "Spain", "Italy", "Netherlands", "Portugal", "Belgium", "Sweden", "Finland", "Greece", "Ireland", "United Kingdom", "Norway","Switzerland", "Poland", "Andorra","Luxembourg", "Liechtenstein", "Malta", "San Marino", "Holy See","Monaco","Hungary", "Czechia","Slovakia", "Slovenia", "Croatia","Bosnia and Herzegovina", "Serbia", "Albania", "Romania", "Bulgaria", "Ukraine", "Belarus", "Latvia", "Estonia", "Lithuania","Moldova","North Macedonia", "Kosovo","Montenegro","Iceland","Cyprus"]

    for ic,cntry in enumerate(data['Country/Region']):
        if (cntry in lstCountry) or (strCountry=="World"):
            locRegion = data.iloc[ic].values[4:]
            locRegion[np.isnan(locRegion.tolist())] = 0
            evolution[:] += locRegion.astype(int)

    return evolution


def evolution_country(strCountry,dataParam,displayParam):

    field = displayParam['Field']
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
        dedt = np.diff(evolution)
        evol = savgol_filter(dedt, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        d2edt2 = np.zeros(len(evolution))
        d2edt2[2:] = np.diff(evol)
        evol = d2edt2[dataParam['FilterDate']]
    elif dataParam['EvolutionType'] == "R0":
        R0 = np.zeros(len(evolution))
        delta0 = np.diff(evolution)
        delta = savgol_filter(delta0, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        R0[1:] = delta/np.roll(delta,5)
        evol = R0[dataParam['FilterDate']]

    return evol


def dateOut(date):
    return date.strftime('%m/%d/%y').lstrip("0").replace("/0", "/")


def dateIn(strDate):
    spl = strDate.split('/')
    month = int(spl[0])
    day = int(spl[1])
    year = int("20%s" %spl[2])
    return dt.date(year, month,day)


def setFitExtraParam(field, fittingPeriod, extrapolPeriod,dataParam,iExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, iExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="Active":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, iExtrapol]
