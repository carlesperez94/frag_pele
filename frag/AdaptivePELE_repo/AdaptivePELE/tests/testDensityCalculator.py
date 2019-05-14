from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.spawning import densitycalculator
from AdaptivePELE.spawning import densitycalculatortypes
import unittest


class densityCalculatorTest(unittest.TestCase):
    def testDensityCalculatorUndefinedBlock(self):
        spawningBlock = {}
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)
        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null)

    def testDensityCalculatorNullCalculator(self):
        spawningBlock = {
            type: "irrelevant",
            "density": {
                "type": "null"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)
        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null)

    def testDensityCalculatorHeavisideNoParams(self):
        spawningBlock = {
            type: "irrelevant",
            "density": {
                "type": "heaviside"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)

        goldenValues = [1.]
        goldenConditions = []

        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside)
        self.assertEqual(densityCalculator.values, goldenValues)
        self.assertEqual(densityCalculator.conditions, goldenConditions)

    def testDensityCalculatorHeavisideParams(self):
        spawningBlock = {
            type: "irrelevant",
            "density": {
                "type": "heaviside",
                "params": {
                    "conditions": [1, 2],
                    "values": [1, 2]
                }
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)

        goldenValues = [1., 2.]
        goldenConditions = [1., 2.]

        self.assertEqual(densityCalculator.type, densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside)
        self.assertEqual(densityCalculator.values, goldenValues)
        self.assertEqual(densityCalculator.conditions, goldenConditions)

    def testDensityCalculatorContinuousParams(self):
        spawningBlock = {
            type: "irrelevant",
            "density": {
                "type": "continuous"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)

        self.assertAlmostEqual(densityCalculator.calculate(0.5, 8), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 8), 8)
        self.assertAlmostEqual(densityCalculator.calculate(0.5, 6), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 6), 8)
        self.assertAlmostEqual(densityCalculator.calculate(0.5, 4), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 4), 8)
        # Changed behaviour, now densities are always set with contact threshold
        # 8
        # self.assertAlmostEqual(densityCalculator.calculate(0.2, 6), 1)
        # self.assertAlmostEqual(densityCalculator.calculate(0.5, 6), 8)
        # self.assertAlmostEqual(densityCalculator.calculate(0.05, 4), 1)
        # self.assertAlmostEqual(densityCalculator.calculate(0.2, 4), 8)

    def testDensityCalculatorInverseContinuousParams(self):
        spawningBlock = {
            type: "irrelevant",
            "density": {
                "type": "exitContinuous"
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(spawningBlock)

        self.assertAlmostEqual(densityCalculator.calculate(0.5, 8), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 8), 0.125)
        self.assertAlmostEqual(densityCalculator.calculate(0.0, 8), 3.375)
        self.assertAlmostEqual(densityCalculator.calculate(0.5, 6), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 6), 0.125)
        self.assertAlmostEqual(densityCalculator.calculate(0.5, 4), 1)
        self.assertAlmostEqual(densityCalculator.calculate(1.5, 4), 0.125)
