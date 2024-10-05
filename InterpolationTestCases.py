import unittest
from InitialData.interpolate import *

# from InitialData.interpolate_grid import *
import numpy as np
import pandas as pd


class TestInterpolation(unittest.TestCase):
    def test_InterpolationTest(self):
        """
        Test the interpolation functions
        """
        # testing if the locating point functionality works
        testr = np.linspace(0, 10, 11)
        testtheta = np.linspace(0, np.pi, 11)
        testphi = np.linspace(0, 2 * np.pi, 11)
        testcombs = [
            [np.float64(0), np.float64(0), np.float64(0)],
            [np.float64(1), np.float64(0.6283185307179586), np.float64(0.6283185307179586)],
            [np.float64(1), np.float64(0.6283185307179586), np.float64(1.2566370614359172)],
            [np.float64(1), np.float64(0.9424777960769379), np.float64(0.6283185307179586)],
            [np.float64(1), np.float64(0.9424777960769379), np.float64(1.2566370614359172)],
        ]

        self.assertEqual(
            locate_point(
                testr,
                testtheta,
                testphi,
                (0.5196152422706632, 0.5196152422706632, 0.5196152422706632),
            ),
            testcombs,
        )

        tmpcombs_2 = [
            [np.float64(3), np.float64(0.3141592653589793), np.float64(0.6283185307179586)],
            [np.float64(3), np.float64(0.3141592653589793), np.float64(1.2566370614359172)],
            [np.float64(3), np.float64(0.6283185307179586), np.float64(0.6283185307179586)],
            [np.float64(3), np.float64(0.6283185307179586), np.float64(1.2566370614359172)],
            [np.float64(4), np.float64(0.3141592653589793), np.float64(0.6283185307179586)],
            [np.float64(4), np.float64(0.3141592653589793), np.float64(1.2566370614359172)],
            [np.float64(4), np.float64(0.6283185307179586), np.float64(0.6283185307179586)],
            [np.float64(4), np.float64(0.6283185307179586), np.float64(1.2566370614359172)],
        ]
        self.assertEqual(
            locate_point(
                testr,
                testtheta,
                testphi,
                (np.float64(3.0), np.float64(1.0), np.float64(2.2)),
            ),
            tmpcombs_2,
        )


if __name__ == "__main__":
    unittest.main()  # this will run all the tests in the file when the file is run directly
