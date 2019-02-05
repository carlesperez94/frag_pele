from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.clustering import thresholdcalculator
import unittest
import numpy as np


class thresholdCalculatorTest(unittest.TestCase):
    def testConstantDefaultParams(self):
        thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()

        clusteringBlock = {
            "type": "irrelevant",
            "thresholdCalculator": {
                "type": "constant"
            }
        }

        thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
        value = thresholdCalculator.value

        goldenValue = 2.
        self.assertAlmostEqual(value, goldenValue, 10)

    def testConstantOtherParams(self):
        thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()

        clusteringBlock = {
            "type": "irrelevant",
            "thresholdCalculator": {
                "type": "constant",
                "params": {
                    "value": 3
                }
            }
        }

        thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
        value = thresholdCalculator.value

        goldenValue = 3.
        self.assertAlmostEqual(value, goldenValue, 10)

    def testDefaultThresholdCalculator(self):
        thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()

        clusteringBlock = {
            "type": "irrelevant"
        }

        thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
        values = thresholdCalculator.values
        conditions = thresholdCalculator.conditions

        goldenValues = [2., 3., 4, 5]
        goldenConditions = [1.0, 0.75, 0.5]
        np.testing.assert_almost_equal(np.array(values), np.array(goldenValues), 10)
        np.testing.assert_almost_equal(np.array(conditions), np.array(goldenConditions), 10)

    def testHeavisideDefaultParams(self):
        thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()

        clusteringBlock = {
            "type": "irrelevant",
            "thresholdCalculator": {
                "type": "heaviside"
            }
        }

        thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
        values = thresholdCalculator.values
        conditions = thresholdCalculator.conditions

        goldenValues = [2., 3, 4., 5.]
        goldenConditions = [1.0, 0.75, 0.5]
        np.testing.assert_almost_equal(np.array(values), np.array(goldenValues), 10)
        np.testing.assert_almost_equal(np.array(conditions), np.array(goldenConditions), 10)

    def testHeavisideParams(self):
        thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()

        clusteringBlock = {
            "type": "irrelevant",
            "thresholdCalculator": {
                "type": "heaviside",
                "params": {
                    "values": [2, 5, 10],
                    "conditions": [3, 10]
                }
            }
        }

        thresholdCalculator = thresholdCalculatorBuilder.build(clusteringBlock)
        values = thresholdCalculator.values
        conditions = thresholdCalculator.conditions

        goldenValues = [2., 5., 10.]
        goldenConditions = [3., 10.]
        self.assertAlmostEqual(values, goldenValues, 10)
        self.assertAlmostEqual(conditions, goldenConditions, 10)
