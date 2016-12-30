import numpy as np

from ..statistics import stats

def test_simple_stats():
    a = np.ones(10)
    truth = {'mean': 1,
            'median': 1,
            'stddev': 0,
            'total': 9,
            'npoints': 10}

    calc = stats(a)

    assert calc == truth, "Crazy wrong stats"
