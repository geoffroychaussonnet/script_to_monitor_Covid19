import sys
from pathlib import Path

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
    dataParam['Countries'] = set(dataParam['Confirmed']['Country/Region'])
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
    elif strCountry == "Africa":
        lstCountry = "Morocco, Tunisia, Algeria, Lybia, Egypt, Mali, Niger, Chad, Sudan ,Ethiopia, Mauritania, Senegal, Guinea, Liberia, Ghana, Benin, Togo, Nigeria, Sierra Leone, Cameroon, Central African Republic, Gabon, Congo (Brazzaville), Congo (Kinshasa), Angola, Namibia, Botswana, Lesotho, South Africa, Eswatini, Zimbabwe, Mozambique, Zambia, Madagascar, Burundi, Kenya, Uganda, Somalia, South Sudan, Cote d'Ivoire, Rwanda, Djibouti"
        lstCountry = lstCountry.split(",")
        lstCountry = [text.strip() for text in lstCountry]

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


def dateIn(str_date):
    """
    >>> dateIn("3/7/20")
    datetime.date(2020, 3, 7)
    >>> dateIn("3/7/2020")
    datetime.date(2020, 3, 7)

    :param str_date: date as M/D/Y
    :return: a datetime.date object
    """
    month, day, year = map(int, str_date.split('/'))
    if year < 100:
        year += 2000
    return dt.date(year, month, day)


def setFitExtraParam(field, fittingPeriod, extrapolPeriod,dataParam,iExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, iExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="Active":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, iExtrapol]


def unit_and_field(field):
    strUnit = "[-]"
    if field == "Confirmed":
        txtField = "confirmed cases"
    elif field == "Deaths":
        txtField = "deaths"
    elif field == "Active":
        txtField = "active cases"
    elif field == "DeathRate":
        txtField = "death rate"
        strUnit = "[%]"
    return strUnit, txtField


def txt_evol(evolutionType):
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
    return txtEvol


def title_and_y_axis(displayParam, field, strUnit, txtEvol, txtField,
                     txt_title_format):
    txtTitle = txt_title_format % (txtEvol, txtField)
    txtYaxis = "%s %s %s" % (txtEvol, txtField, strUnit)
    displayParam['Field'] = field
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis


def file_yscale(displayParam, figures_path, png_format, txtEvol, txtField, yscale, zone):
    strDateToday = dt.date.today().strftime("%Y%m%d")
    Path(figures_path).mkdir(parents=True, exist_ok=True)
    if zone:
        fname = figures_path + "/" + png_format % (
        strDateToday, txtEvol, txtField, zone)
    else:
        fname = figures_path + "/" + png_format % (
            strDateToday, txtEvol, txtField)
    displayParam['FileName'] = fname.replace(" ", "_")
    displayParam['YScale'] = yscale

def parse_confinement(file):
    """
    Parse the confinement.dat file.

    :param file: a file like object, lines are:
                 country, type, date[, type, date, ...]
    :return: a mapping country -> type -> list of dates
    """
    quar_dates_by_type_by_country = {}
    for row in file:
        row = row.strip()
        if not row or row[0] == '#':
            continue

        country, *tds = row.split(",")
        for t, d in zip(tds[::2], tds[1::2]):
            t = t.strip()
            try:
                d = dateIn(d)
            except ValueError:
                print("Ignore row: {}".format(row), file=sys.stderr)
            else:
                quar_dates_by_type_by_country.setdefault(
                    country, {}).setdefault(t, []).append(d)

    return quar_dates_by_type_by_country
