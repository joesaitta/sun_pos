"""
We probably already have this, but need days from epoch to use Vallado's stuff
"""
from dateutil import parser


def jd(ut1):
    """
    Convert a datetime to Julian Date (UT1, days from epoch 1/1, 4713 B.C. 1200)
    Algorithm 14 from Vallado, 3rd edition
    :param ut1: datetime
    :return: days from epoch
    """
    a = 367 * ut1.year
    b = int(7 * (ut1.year + int((ut1.month + 9)/12)) / 4)
    c = int(275 * ut1.month / 9)
    d = (((ut1.second / 60 + ut1.minute) / 60) + ut1.hour) / 24

    return a - b + c + ut1.day + 1721013.5 + d


if __name__ == "__main__":
    t = parser.parse('April 2, 2006, 00:00 UTC')
    print(jd(t))

    t = parser.parse('October 26, 1996, 14:20 UTC')
    print(jd(t))

