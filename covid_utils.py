import sys
from pathlib import Path

import numpy as np
import pandas as pd
import datetime as dt
from scipy.signal import savgol_filter

# the beginning of values in the csv file
start_of_values = 4


def load_data(path, start_date=dt.date(2020, 1, 1)):
    confirmed = pd.read_csv(path + "time_series_covid19_confirmed_global.csv")
    deaths = pd.read_csv(path + "time_series_covid19_deaths_global.csv")
    recovered = pd.read_csv(path + "time_series_covid19_recovered_global.csv")
    data = {'Confirmed': confirmed, 'Deaths': deaths, 'Recovered': recovered,
            'Countries': set(confirmed['Country/Region'])}
    with Path("confinement.dat").open() as f:
        data['Confinement'] = extract_confinement(parse_confinement(f))

    dates_str = deaths.columns[start_of_values:].values.astype(str)
    # Convert date axis to date vector
    dates = np.array([dt.datetime.strptime(date_str, '%m/%d/%y').date()
                      for date_str in dates_str])

    # Filter axe of dates
    filter_date = dates >= start_date
    data['FilterDate'] = filter_date
    data['DateAxis'] = dates_str[filter_date]
    data['Dates'] = dates[filter_date]

    return data


def cumulative_evolution_single(area, data):
    """
    :param area: a country or the name set of countries
    :param data: the data from J. Hopkins University
    :return: the cumulative evolution of the area
    """
    size = len(data.iloc[0].values[start_of_values:])
    evolution = np.zeros(size, dtype=int)
    countries = _get_countries(area)

    for i, country in enumerate(data['Country/Region']):
        if area == "World" or country in countries:
            country_values = data.iloc[i].values[start_of_values:]
            # replace missing values by zeroes
            country_values[np.isnan(country_values.tolist())] = 0
            evolution[:] += country_values.astype(int)

    return evolution


def _get_countries(area):
    """
    :param area: a country or the name set of countries
    :return: countries that belong to the area
    """
    if area == "EU":
        countries = {"France", "Germany", "Spain", "Italy", "Netherlands",
                     "Portugal", "Belgium", "Sweden", "Finland", "Greece",
                     "Ireland", "Poland", "Luxembourg", "Malta", "Slovenia",
                     "Austria", "Croatia", "Hungary", "Czechia", "Slovakia",
                     "Hungary", "Romania", "Bulgaria", "Cyprus", "Lithuania",
                     "Latvia", "Estonia"}
    elif area == "European continent":
        countries = {"France", "Germany", "Spain", "Italy", "Netherlands",
                     "Portugal", "Belgium", "Sweden", "Finland", "Greece",
                     "Ireland", "United Kingdom", "Norway", "Switzerland",
                     "Poland", "Andorra", "Luxembourg", "Liechtenstein",
                     "Malta", "San Marino", "Holy See", "Monaco", "Hungary",
                     "Czechia", "Slovakia", "Slovenia", "Croatia",
                     "Bosnia and Herzegovina", "Serbia", "Albania", "Romania",
                     "Bulgaria", "Ukraine", "Belarus", "Latvia", "Estonia",
                     "Lithuania", "Moldova", "North Macedonia", "Kosovo",
                     "Montenegro", "Iceland", "Cyprus"}
    elif area == "Africa":
        countries = {"Morocco", "Tunisia", "Algeria", "Lybia", "Egypt", "Mali",
                     "Niger", "Chad", "Sudan", "Ethiopia", "Mauritania",
                     "Senegal", "Guinea", "Liberia", "Ghana", "Benin", "Togo",
                     "Nigeria", "Sierra Leone", "Cameroon",
                     "Central African Republic", "Gabon",
                     "Congo (Brazzaville)", "Congo (Kinshasa)", "Angola",
                     "Namibia", "Botswana", "Lesotho", "South Africa",
                     "Eswatini", "Zimbabwe", "Mozambique", "Zambia",
                     "Madagascar", "Burundi", "Kenya", "Uganda", "Somalia",
                     "South Sudan", "Cote d'Ivoire", "Rwanda", "Djibouti"}
    else:
        countries = {area}
    return countries


def evolution_country(area, data, field, evolution_type, filter_date,
                      smoothing):
    if field == "Confirmed":
        cumulative_evolution = cumulative_evolution_single(area, data['Confirmed'])
    elif field == "Deaths":
        cumulative_evolution = cumulative_evolution_single(area, data['Deaths'])
    elif field == "Active":
        evol_c = cumulative_evolution_single(area, data['Confirmed'])
        evol_d = cumulative_evolution_single(area, data['Deaths'])
        evol_r = cumulative_evolution_single(area, data['Recovered'])
        cumulative_evolution = evol_c - evol_r - evol_d
    elif field == "DeathRate":
        evol_c = cumulative_evolution_single(area, data['Confirmed'])
        evol_d = cumulative_evolution_single(area, data['Deaths'])
        cumulative_evolution = evol_d/evol_c*100
    else:
        raise ValueError(field)

    if evolution_type == "cumulative":
        evolution = cumulative_evolution
    elif evolution_type == "daily":
        evolution = np.zeros(len(cumulative_evolution))
        evolution[1:] = np.diff(cumulative_evolution)
    elif evolution_type == "curvature":
        evolution = np.zeros(len(cumulative_evolution))
        evolution[2:] = np.diff(cumulative_evolution, 2)
    elif evolution_type == "smoothedCurvature":
        # np.diff #1
        window_length, polyorder = smoothing
        evolution = np.diff(cumulative_evolution)
        smoothed_evolution = savgol_filter(evolution, window_length, polyorder)
        # np.diff #2
        evolution = np.zeros(len(cumulative_evolution))
        evolution[2:] = np.diff(smoothed_evolution)
    elif evolution_type == "R0":
        window_length, polyorder = smoothing
        evolution = np.zeros(len(cumulative_evolution))
        delta = np.diff(cumulative_evolution)
        smoothed_delta = savgol_filter(delta, window_length, polyorder)
        evolution[1:] = smoothed_delta/np.roll(smoothed_delta, 5)
    else:
        raise ValueError(evolution_type)

    return evolution[filter_date]


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


extrapol_period_by_field = {
    "Confirmed": 14,
    "Deaths": 21,
    "Active": 21,
    "DeathRate": 21
}


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


def title_and_y_axis(displayParam, strUnit, txtEvol, txtField,
                     txt_title_format):
    txtTitle = txt_title_format.format(txtEvol, txtField)
    displayParam['title'] = txtTitle
    txtYaxis = "{} {} {}".format(txtEvol, txtField, strUnit)
    displayParam['YaxisLabel'] = txtYaxis


def file_name(figures_path, png_format, txt_evol, txt_field, zone=None,
              str_date_today=dt.date.today().strftime("%Y%m%d")):
    """
    Create the directoy if necessary.

    :param figures_path:
    :param png_format:
    :param txt_evol:
    :param txt_field:
    :param zone:
    :param str_date_today: date of today YYYYMMDD
    :return:
    """
    if zone:
        fname = png_format.format(str_date_today, txt_evol, txt_field, zone)
    else:
        fname = png_format.format(str_date_today, txt_evol, txt_field)
    fullpath = str(Path(figures_path, fname))
    return fullpath.replace(" ", "_")


def ensure_figures_directory_exists(figures_path):
    path = Path(figures_path)
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)


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


def extract_confinement(quar_dates_by_type_by_country,
                        future=dt.date(2099, 1, 1)):
    """
    Extract the last partial confinement, or the first total confinement,
    or a date in the future

    :param quar_dates_by_type_by_country: a mapping `country -> type -> [date]`
    :return: a mapping `country -> date`.
    """
    def first_confinement(q_by_t):
        if 'T' in q_by_t:
            d = min(q_by_t['T'], default=future)
        elif 'P' in q_by_t:
            d = max(q_by_t['P'], default=future)
        else:
            d = future
        return dateOut(d)

    return {c: first_confinement(q_by_t)
            for c, q_by_t in quar_dates_by_type_by_country.items()}


# temporary data
quar_date_by_area = {"World": '3/22/21', "EU": '3/22/21', "China": '1/22/22',
                     "US": '3/22/20', "European continent": '3/22/21',
                     "Italy": '3/9/20', "Spain": '3/14/20',
                     "Germany": '3/19/20', "France": '3/17/20',
                     "Iran": '8/17/20', "Korea, South": '5/22/20',
                     "Japan": '5/22/20', "Switzerland": '5/22/20',
                     "United Kingdom": '3/22/20', "Denmark": '3/13/20',
                     "Norway": '3/12/20', "Sweden": '3/28/20',
                     "Finland": '3/19/20', "Canada": '5/22/20',
                     "Belgium": '3/18/20', "Ireland": '3/28/20',
                     }
