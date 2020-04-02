# Python Imports
import argparse
import unittest
from unittest.mock import Mock, patch

# Third-Party Imports

# Project Imports
from frag_pele.Helpers import argument_parser as ap


class TestArgumentParser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.class_parser = argparse.ArgumentParser()

    def test_create_parser(self):
        pass

    @patch.object(ap, "_add_required_named_arguments")
    @patch.object(ap, "_add_standard_arguments")
    @patch.object(ap, "_add_plop_arguments")
    @patch.object(ap, "_add_pele_conf_arguments")
    @patch.object(ap, "_add_clustering_arguments")
    @patch.object(ap, "_add_protocol_arguments")
    @patch.object(ap, "_add_output_format_arguments")
    @patch.object(ap, "_add_others_arguments")
    def test_add_all_arguments(self, mock_a, mock_b, mock_c, mock_d, mock_e, mock_f, mock_g, mock_h):
        parser = argparse.ArgumentParser()

        ap._add_all_arguments(parser)
        mock_a.assert_called_with(parser)
        mock_b.assert_called_with(parser)
        mock_c.assert_called_with(parser)
        mock_d.assert_called_with(parser)
        mock_e.assert_called_with(parser)
        mock_f.assert_called_with(parser)
        mock_g.assert_called_with(parser)
        mock_h.assert_called_with(parser)

    def test_add_required_named_arguments(self):
        pass

    def test_add_arguments_to_required_named(self):
        mock = Mock()

        ap._add_arguments_to_required_named(mock)

        self.assertEqual(len(mock.add_argument.call_args_list), 2)
        self.assertEqual(mock.add_argument.call_args_list[0][0], ('-cp', '--complex_pdb'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[1][0], ('-sef', '--serie_file'))  # todo: check a better way

    def test_standard_arguments(self):
        mock = Mock()

        ap._add_standard_arguments(mock)

        self.assertEqual(len(mock.add_argument.call_args_list), 11)
        self.assertEqual(mock.add_argument.call_args_list[0][0], ('--core',))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[1][0], ('-nc', '--no_check'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[2][0], ('-x', '--growing_steps'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[3][0], ('-cr', '--criteria'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[4][0], ('-rst', '--restart'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[5][0], ('-cc', '--c_chain'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[6][0], ('-fc', '--f_chain'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[7][0], ('-tc', '--clash_thr'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[8][0], ('-sc', '--sampling_control'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[9][0], ('-op', '--only_prepare'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[10][0], ('-og', '--only_grow'))  # todo: check a better way

    def test_add_plop_arguments(self):
        mock = Mock()

        ap._add_plop_arguments(mock)

        self.assertEqual(len(mock.add_argument.call_args_list), 3)
        self.assertEqual(mock.add_argument.call_args_list[0][0], ('-pl', '--plop_path'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[1][0], ('-sp', '--sch_python'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[2][0], ('-rot', '--rotamers'))  # todo: check a better way

    def test_add_pele_conf_arguments(self):
        mock = Mock()

        ap._add_pele_conf_arguments(mock)

        self.assertEqual(len(mock.add_argument.call_args_list), 21)
        self.assertEqual(mock.add_argument.call_args_list[0][0], ('-d', '--pele_dir'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[1][0], ('-c', '--contrl'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[2][0], ('-l', '--license'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[3][0], ('-r', '--resfold'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[4][0], ('-rp', '--report'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[5][0], ('-tj', '--traject'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[6][0], ('-cs', '--cpus'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[7][0], ('-stp', '--steps'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[8][0], ('-es', '--pele_eq_steps'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[9][0], ('-miov', '--min_overlap'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[10][0], ('-maov', '--max_overlap'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[11][0], ('-tmp', '--temperature'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[12][0], ('-sd', '--seed'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[13][0], ('-st', '--steering'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[14][0], ('-trh', '--translation_high'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[15][0], ('-roth', '--rotation_high'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[16][0], ('-trl', '--translation_low'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[17][0], ('-rotl', '--rotation_low'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[18][0], ('-rad', '--radius_box'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[19][0], ('-dat', '--data'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[20][0], ('-doc', '--documents'))  # todo: check a better way

    def test_add_clustering_arguments(self):
        parser = argparse.ArgumentParser()
        ap._add_clustering_arguments(parser)

        actions = parser._option_string_actions
        self.assertIn("-dis", actions)
        self.assertIn("--distcont", actions)
        self.assertEqual(20, len(actions))  # 20 cause of auto --help flag

        # self.assertEqual(len(mock.add_argument.call_args_list), 9)
        # self.assertEqual(mock.add_argument.call_args_list[0][0], ('-dis', '--distcont'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[1][0], ('-ct', '--threshold'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[2][0], ('-e', '--epsilon'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[3][0], ('-cn', '--condition'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[4][0], ('-mw', '--metricweights'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[5][0], ('-ncl', '--nclusters'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[6][0], ('-pdbf', '--pdbout'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[7][0], ('-ban', '--banned'))  # todo: check a better way
        # self.assertEqual(mock.add_argument.call_args_list[8][0], ('-lim', '--limit'))  # todo: check a better way

    def test_add_protocol_arguments(self):
        mock = Mock()

        ap._add_protocol_arguments(mock)

        self.assertEqual(mock.add_argument.call_args_list[0][0], ('-HT', '--highthroughput'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[1][0], ('-EX', '--explorative'))  # todo: check a better way
        self.assertEqual(mock.add_argument.call_args_list[2][0], ('--test',))  # todo: check a better way

    def test_add_output_format_arguments(self):
        mock = Mock()

        ap._add_output_format_arguments(mock)

        self.assertEqual(mock.add_argument.call_args[0][0], '--mae')  # todo: check a better way

    def test_add_other_arguments(self):
        mock = Mock()

        ap._add_others_arguments(mock)

        self.assertEqual(mock.add_argument.call_args[0][0], '--rename')  # todo: check a better way

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
        self.assertEqual(mock.temp, 1000000)

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
        pass
