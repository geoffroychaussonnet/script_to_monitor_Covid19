import unittest
import io
import sys
from unittest.mock import patch

from covid_utils import *


class TestParseConfinement(unittest.TestCase):
    def test_basic_file(self):
        f = io.StringIO("""### source: ...

France, T, 3/17/2020
Spain, T, 3/14/2020
Ireland, P, 3/15/2020, P, 3/24/2020, T, 3/28/2020""")

        self.assertEqual(parse_confinement(f),
                         {'France': {'T': [dt.date(2020, 3, 17)]},
                          'Spain': {'T': [dt.date(2020, 3, 14)]},
                          'Ireland': {'P': [dt.date(2020, 3, 15),
                                            dt.date(2020, 3, 24)],
                                      'T': [dt.date(2020, 3, 28)]}})

    def test_err(self):
        f = io.StringIO("""### source: ...

Germany, P,
""")

        with patch('sys.stderr', new=io.StringIO()) as fake_out:
            parse_confinement(f)
            self.assertEqual(fake_out.getvalue(), "Ignore row: Germany, P,\n")

    def test_real_file(self):
        with Path("confinement.dat").open() as f:
            parse_confinement(f)


class TestExtractConfinement(unittest.TestCase):
    def test_partial(self):
        self.assertEqual(
            {'c': '3/15/20'},
            extract_confinement({"c": {'P': [dt.date(2020, 3, 10),
                                             dt.date(2020, 3, 15)]}})
        )

    def test_partial_and_total(self):
        self.assertEqual(
            {'c': '3/20/20'},
            extract_confinement(
                {"c": {'P': [dt.date(2020, 3, 10), dt.date(2020, 3, 15)],
                       'T': [dt.date(2020, 3, 20), dt.date(2020, 3, 25)]}})
        )

    def test_empty_total(self):
        self.assertEqual(
            {'c': '1/1/99'},
            extract_confinement({"c": {'P': [dt.date(2020, 3, 10)], 'T': []}})
        )

    def test_empty(self):
        self.assertEqual(
            {'c': '1/1/99'},
            extract_confinement({"c": {}})
        )
