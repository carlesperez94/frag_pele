# Python Imports
import sys
import argparse
import unittest
from unittest.mock import Mock, patch

# Third-Party Imports

# Project Imports
from frag_pele.Helpers import argument_parser as ap
import frag_pele.constants as const


class TestArgumentParser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.class_parser = argparse.ArgumentParser()

    @patch.object(ap, '_add_all_arguments')
    def test_create_parser(self, mock_argparser_lib):
        parser = ap._create_parser()

        mock_argparser_lib.assert_called_once_with(parser)
        self.assertEqual(type(parser), argparse.ArgumentParser)

    @patch.object(ap, "_add_frag_required_named_arguments")
    @patch.object(ap, "_add_frag_standard_arguments")
    @patch.object(ap, "_add_plop_arguments")
    @patch.object(ap, "_add_pele_conf_arguments")
    @patch.object(ap, "_add_clustering_arguments")
    @patch.object(ap, "_add_protocol_arguments")
    @patch.object(ap, "_add_output_format_arguments")
    @patch.object(ap, "_add_others_arguments")
    def test_add_all_arguments(self, mock_a, mock_b, mock_c, mock_d, mock_e, mock_f, mock_g, mock_h):
        parser = argparse.ArgumentParser()

        ap._add_all_arguments(parser)
        mock_a.assert_called_once_with(parser)
        mock_b.assert_called_once_with(parser)
        mock_c.assert_called_once_with(parser)
        mock_d.assert_called_once_with(parser)
        mock_e.assert_called_once_with(parser)
        mock_f.assert_called_once_with(parser)
        mock_g.assert_called_once_with(parser)
        mock_h.assert_called_once_with(parser)

    def test_add_required_named_arguments(self):
        flags_to_check = ['-cp', '--complex_pdb', '-sef', '--serie_file']

        parser = argparse.ArgumentParser()
        ap._add_frag_required_named_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(6, len(actions))  # 6 and not 4 as in flag_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_standard_arguments(self):
        flags_to_check = ['--core', '-nc', '--no_check', '-x', '--growing_steps', '-cr', '--criteria', '-rst',
                          '--restart', '-cc', '--c_chain', '-fc', '--f_chain', '-tc', '--clash_thr', '-sc',
                          '--sampling_control', '-op', '--only_prepare', '-og', '--only_grow', '-EX', '--explorative']

        parser = argparse.ArgumentParser()
        ap._add_frag_standard_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(25, len(actions))  # 25 and not 23 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_plop_arguments(self):
        flags_to_check = ['-pl', '--plop_path', '-sp', '--sch_python', '-rot', '--rotamers']

        parser = argparse.ArgumentParser()
        ap._add_plop_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(8, len(actions))  # 8 and not 6 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_pele_conf_arguments(self):
        flags_to_check = ['-d', '--pele_dir', '-c', '--contrl', '-l', '--license', '-r', '--resfold', '-rp', '--report',
                          '-tj', '--traject', '-cs', '--cpus', '-stp', '--steps', '-es', '--pele_eq_steps', '-miov',
                          '--min_overlap', '-maov', '--max_overlap', '-tmp', '--temperature', '-sd', '--seed', '-st',
                          '--steering', '-trh', '--translation_high', '-roth', '--rotation_high', '-trl',
                          '--translation_low', '-rotl', '--rotation_low', '-rad', '--radius_box', '-dat', '--data',
                          '-doc', '--documents']

        parser = argparse.ArgumentParser()
        ap._add_pele_conf_arguments(parser)
        actions = parser._option_string_actions
        self.assertEqual(44, len(actions))  # 44 and not 42 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_clustering_arguments(self):
        flags_to_check = ['-dis', '--distcont', '-ct', '--threshold', '-e', '--epsilon', '-cn', '--condition', '-mw',
                          '--metricweights', '-ncl', '--nclusters', '-pdbf', '--pdbout', '-ban', '--banned', '-lim',
                          '--limit']

        parser = argparse.ArgumentParser()
        ap._add_clustering_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(20, len(actions))  # 20 and not 18 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_protocol_arguments(self):
        flags_to_check = ['-HT', '--highthroughput', '--test']

        parser = argparse.ArgumentParser()
        ap._add_protocol_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(5, len(actions))  # 5 and not 3 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_output_format_arguments(self):
        flags_to_check = ['--mae']

        parser = argparse.ArgumentParser()
        ap._add_output_format_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(3, len(actions))  # 3 and not 1 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_add_other_arguments(self):
        flags_to_check = ['--rename']

        parser = argparse.ArgumentParser()
        ap._add_others_arguments(parser)
        actions = parser._option_string_actions

        self.assertEqual(3, len(actions))  # 3 and not 1 as in flags_to_check, because automatically adds -h --help
        self._loop_assert_in(flags_to_check, actions)

    def test_check_highthroughput_in_args_true(self):
        mock = Mock()
        mock.highthroughput = True
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0

        ap._check_highthroughput_in_args(mock)

        self.assertEqual(mock.growing_steps, 1)
        self.assertEqual(mock.steps, 3)
        self.assertEqual(mock.pele_eq_steps, 10)

    def test_check_highthroughput_in_args_false(self):
        mock = Mock()
        mock.highthroughput = False
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0

        ap._check_highthroughput_in_args(mock)

        self.assertEqual(mock.growing_steps, 0)
        self.assertEqual(mock.steps, 0)
        self.assertEqual(mock.pele_eq_steps, 0)

    def test_check_test_in_args_true(self):
        mock = Mock()
        mock.test = True
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0
        mock.temp = 0

        ap._check_test_in_args(mock)

        self.assertEqual(mock.growing_steps, 1)
        self.assertEqual(mock.steps, 1)
        self.assertEqual(mock.pele_eq_steps, 1)
        self.assertEqual(mock.temperature, 1000000)

    def test_check_test_in_args_false(self):
        mock = Mock()
        mock.test = False
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0
        mock.temp = 0

        ap._check_test_in_args(mock)

        self.assertEqual(mock.growing_steps, 0)
        self.assertEqual(mock.steps, 0)
        self.assertEqual(mock.pele_eq_steps, 0)
        self.assertEqual(mock.temp, 0)

    def test_parse_arguments(self):
        test_args = ['TestFragArgParser', '-cp', 'TEST_CP', '-sef', 'TEST_SEF', '-nc', '-rst', '-op', '-og',
                     '-EX', '--mae', '--rename']

        groups = ['positional arguments', 'optional arguments', 'required named arguments', 'Frag Standard Arguments',
                  'Plop Related Arguments', 'PELE Related Arguments', 'Clustering Related Arguments']

        dict_expected = {'positional arguments': argparse.Namespace(),
                         'optional arguments': argparse.Namespace(help=None, highthroughput=False, mae=True,
                                                                  rename=True, test=False),
                         'required named arguments': argparse.Namespace(complex_pdb='TEST_CP', serie_file='TEST_SEF'),
                         'Frag Standard Arguments': argparse.Namespace(c_chain=const.C_CHAIN,
                                                                       clash_thr=const.CLASH_THRESHOLD, core=None,
                                                                       criteria=const.SELECTION_CRITERIA,
                                                                       explorative=True,
                                                                       f_chain=const.F_CHAIN,
                                                                       growing_steps=const.GROWING_STEPS, no_check=True,
                                                                       only_grow=True, only_prepare=True, restart=True,
                                                                       sampling_control=None),
                         'Plop Related Arguments': argparse.Namespace(plop_path=const.PLOP_PATH,
                                                                      rotamers=const.ROTRES,
                                                                      sch_python=const.SCHRODINGER_PY_PATH),
                         'PELE Related Arguments': argparse.Namespace(contrl=const.CONTROL_TEMPLATE,
                                                                      cpus=const.CPUS, data=const.PATH_TO_PELE_DATA
                                                                      , documents=const.PATH_TO_PELE_DOCUMENTS,
                                                                      license=const.PATH_TO_LICENSE,
                                                                      max_overlap=const.MAX_OVERLAP,
                                                                      min_overlap=const.MIN_OVERLAP,
                                                                      pele_dir=const.PATH_TO_PELE,
                                                                      pele_eq_steps=const.PELE_EQ_STEPS,
                                                                      radius_box=const.RADIUS_BOX,
                                                                      report=const.REPORT_NAME,
                                                                      resfold=const.RESULTS_FOLDER,
                                                                      rotation_high=const.ROTATION_HIGH,
                                                                      rotation_low=const.ROTATION_LOW,
                                                                      seed=const.SEED, steering=const.STEERING,
                                                                      steps=const.STEPS,
                                                                      temperature=const.TEMPERATURE,
                                                                      traject=const.TRAJECTORY_NAME,
                                                                      translation_high=const.TRANSLATION_HIGH,
                                                                      translation_low=const.TRANSLATION_LOW),
                         'Clustering Related Arguments': argparse.Namespace(banned=const.BANNED_ANGLE_THRESHOLD,
                                                                            condition=const.CONDITION,
                                                                            distcont=const.DISTANCE_COUNTER,
                                                                            epsilon=const.EPSILON,
                                                                            limit=const.BANNED_ANGLE_THRESHOLD,
                                                                            metricweights=const.METRICS_WEIGHTS,
                                                                            nclusters=const.NUM_CLUSTERS,
                                                                            pdbout=const.PDBS_OUTPUT_FOLDER,
                                                                            threshold=const.CONTACT_THRESHOLD)}

        with patch.object(sys, 'argv', test_args):
            result = ap.parse_arguments()
            self.assertEqual(len(result.keys()), 7)
            for group in groups:
                self.assertEqual(dict_expected[group], result[group])

    def test_parse_arguments_highthroughput_flag(self):
        test_args = ['TestFragArgParser', '-cp', 'TEST_CP', '-sef', 'TEST_SEF', '-nc', '-rst', '-op', '-og',
                     '-EX', '--mae', '--rename', '-HT']

        groups = ['positional arguments', 'optional arguments', 'required named arguments', 'Frag Standard Arguments',
                  'Plop Related Arguments', 'PELE Related Arguments', 'Clustering Related Arguments']

        dict_expected = {'positional arguments': argparse.Namespace(),
                         'optional arguments': argparse.Namespace(help=None, highthroughput=True, mae=True,
                                                                  rename=True, test=False),
                         'required named arguments': argparse.Namespace(complex_pdb='TEST_CP', serie_file='TEST_SEF'),
                         'Frag Standard Arguments': argparse.Namespace(c_chain=const.C_CHAIN,
                                                                       clash_thr=const.CLASH_THRESHOLD, core=None,
                                                                       criteria=const.SELECTION_CRITERIA,
                                                                       explorative=True,
                                                                       f_chain=const.F_CHAIN,
                                                                       growing_steps=1, no_check=True,
                                                                       only_grow=True, only_prepare=True, restart=True,
                                                                       sampling_control=None),
                         'Plop Related Arguments': argparse.Namespace(plop_path=const.PLOP_PATH,
                                                                      rotamers=const.ROTRES,
                                                                      sch_python=const.SCHRODINGER_PY_PATH),
                         'PELE Related Arguments': argparse.Namespace(contrl=const.CONTROL_TEMPLATE,
                                                                      cpus=const.CPUS, data=const.PATH_TO_PELE_DATA
                                                                      , documents=const.PATH_TO_PELE_DOCUMENTS,
                                                                      license=const.PATH_TO_LICENSE,
                                                                      max_overlap=const.MAX_OVERLAP,
                                                                      min_overlap=const.MIN_OVERLAP,
                                                                      pele_dir=const.PATH_TO_PELE,
                                                                      pele_eq_steps=10,
                                                                      radius_box=const.RADIUS_BOX,
                                                                      report=const.REPORT_NAME,
                                                                      resfold=const.RESULTS_FOLDER,
                                                                      rotation_high=const.ROTATION_HIGH,
                                                                      rotation_low=const.ROTATION_LOW,
                                                                      seed=const.SEED, steering=const.STEERING,
                                                                      steps=3,
                                                                      temperature=const.TEMPERATURE,
                                                                      traject=const.TRAJECTORY_NAME,
                                                                      translation_high=const.TRANSLATION_HIGH,
                                                                      translation_low=const.TRANSLATION_LOW),
                         'Clustering Related Arguments': argparse.Namespace(banned=const.BANNED_ANGLE_THRESHOLD,
                                                                            condition=const.CONDITION,
                                                                            distcont=const.DISTANCE_COUNTER,
                                                                            epsilon=const.EPSILON,
                                                                            limit=const.BANNED_ANGLE_THRESHOLD,
                                                                            metricweights=const.METRICS_WEIGHTS,
                                                                            nclusters=const.NUM_CLUSTERS,
                                                                            pdbout=const.PDBS_OUTPUT_FOLDER,
                                                                            threshold=const.CONTACT_THRESHOLD)}

        with patch.object(sys, 'argv', test_args):
            result = ap.parse_arguments()
            self.assertEqual(len(result.keys()), 7)
            for group in groups:
                self.assertEqual(dict_expected[group], result[group])

    def test_parse_arguments_test_flag(self):
        test_args = ['TestFragArgParser', '-cp', 'TEST_CP', '-sef', 'TEST_SEF', '-nc', '-rst', '-op', '-og',
                     '-EX', '--mae', '--rename', '--test']

        groups = ['positional arguments', 'optional arguments', 'required named arguments', 'Frag Standard Arguments',
                  'Plop Related Arguments', 'PELE Related Arguments', 'Clustering Related Arguments']

        dict_expected = {'positional arguments': argparse.Namespace(),
                         'optional arguments': argparse.Namespace(help=None, highthroughput=False, mae=True,
                                                                  rename=True, test=True),
                         'required named arguments': argparse.Namespace(complex_pdb='TEST_CP', serie_file='TEST_SEF'),
                         'Frag Standard Arguments': argparse.Namespace(c_chain=const.C_CHAIN,
                                                                       clash_thr=const.CLASH_THRESHOLD, core=None,
                                                                       criteria=const.SELECTION_CRITERIA,
                                                                       explorative=True,
                                                                       f_chain=const.F_CHAIN,
                                                                       growing_steps=1, no_check=True,
                                                                       only_grow=True, only_prepare=True, restart=True,
                                                                       sampling_control=None),
                         'Plop Related Arguments': argparse.Namespace(plop_path=const.PLOP_PATH,
                                                                      rotamers=const.ROTRES,
                                                                      sch_python=const.SCHRODINGER_PY_PATH),
                         'PELE Related Arguments': argparse.Namespace(contrl=const.CONTROL_TEMPLATE,
                                                                      cpus=const.CPUS, data=const.PATH_TO_PELE_DATA
                                                                      , documents=const.PATH_TO_PELE_DOCUMENTS,
                                                                      license=const.PATH_TO_LICENSE,
                                                                      max_overlap=const.MAX_OVERLAP,
                                                                      min_overlap=const.MIN_OVERLAP,
                                                                      pele_dir=const.PATH_TO_PELE,
                                                                      pele_eq_steps=1,
                                                                      radius_box=const.RADIUS_BOX,
                                                                      report=const.REPORT_NAME,
                                                                      resfold=const.RESULTS_FOLDER,
                                                                      rotation_high=const.ROTATION_HIGH,
                                                                      rotation_low=const.ROTATION_LOW,
                                                                      seed=const.SEED, steering=const.STEERING,
                                                                      steps=1, temperature=1000000,
                                                                      traject=const.TRAJECTORY_NAME,
                                                                      translation_high=const.TRANSLATION_HIGH,
                                                                      translation_low=const.TRANSLATION_LOW),
                         'Clustering Related Arguments': argparse.Namespace(banned=const.BANNED_ANGLE_THRESHOLD,
                                                                            condition=const.CONDITION,
                                                                            distcont=const.DISTANCE_COUNTER,
                                                                            epsilon=const.EPSILON,
                                                                            limit=const.BANNED_ANGLE_THRESHOLD,
                                                                            metricweights=const.METRICS_WEIGHTS,
                                                                            nclusters=const.NUM_CLUSTERS,
                                                                            pdbout=const.PDBS_OUTPUT_FOLDER,
                                                                            threshold=const.CONTACT_THRESHOLD)}

        with patch.object(sys, 'argv', test_args):
            result = ap.parse_arguments()
            self.assertEqual(len(result.keys()), 7)
            for group in groups:
                self.assertEqual(dict_expected[group], result[group])

    def _loop_assert_in(self, flag_list, actions):
        for flag in flag_list:
            self.assertIn(flag, actions)
