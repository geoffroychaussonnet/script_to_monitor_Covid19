import unittest
import io
import sys

from covid_utils import *


class TestParseConfinement(unittest.TestCase):
    def test(self):
        f = io.StringIO("""### source: ...

France, T, 3/17/2020
Spain, T, 3/14/2020
Ireland, P, 3/15/2020, P, 3/24/2020, T, 3/28/2020""")

        self.assertEqual(parse_confinement(f),
                         {'France': {'T': [datetime.date(2020, 3, 17)]},
                          'Spain': {'T': [datetime.date(2020, 3, 14)]},
                          'Ireland': {'P': [datetime.date(2020, 3, 15),
                                            datetime.date(2020, 3, 24)],
                                      'T': [datetime.date(2020, 3, 28)]}})
