# Python Imports
import argparse
from unittest.mock import Mock

# Third-Party Imports
from nose.tools import assert_equal

# Project Imports
from frag_pele.Helpers import argument_parser as ap


class TestArgumentParser:

    @classmethod
    def setup_class(cls):
        cls.class_parser = argparse.ArgumentParser()

    def test_create_parser(self):
        pass

    def test_add_all_arguments(self):
        pass

    def test_add_required_named_arguments(self):
        pass

    def test_standard_arguments(self):
        mock = Mock()

        ap._add_standard_arguments(mock)

        assert_equal(len(mock.add_argument.call_args_list), 11)
        assert_equal(mock.add_argument.call_args_list[0][0], ('--core',))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[1][0], ('-nc', '--no_check'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[2][0], ('-x', '--growing_steps'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[3][0], ('-cr', '--criteria'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[4][0], ('-rst', '--restart'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[5][0], ('-cc', '--c_chain'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[6][0], ('-fc', '--f_chain'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[7][0], ('-tc', '--clash_thr'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[8][0], ('-sc', '--sampling_control'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[9][0], ('-op', '--only_prepare'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[10][0], ('-og', '--only_grow'))  # todo: check a better way

    def test_add_plop_arguments(self):
        mock = Mock()

        ap._add_plop_arguments(mock)

        assert_equal(len(mock.add_argument.call_args_list), 3)
        assert_equal(mock.add_argument.call_args_list[0][0], ('-pl', '--plop_path'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[1][0], ('-sp', '--sch_python'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[2][0], ('-rot', '--rotamers'))  # todo: check a better way

    def test_add_pele_conf_arguments(self):
        mock = Mock()

        ap._add_pele_conf_arguments(mock)

        assert_equal(len(mock.add_argument.call_args_list), 21)
        assert_equal(mock.add_argument.call_args_list[0][0], ('-d', '--pele_dir'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[1][0], ('-c', '--contrl'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[2][0], ('-l', '--license'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[3][0], ('-r', '--resfold'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[4][0], ('-rp', '--report'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[5][0], ('-tj', '--traject'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[6][0], ('-cs', '--cpus'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[7][0], ('-stp', '--steps'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[8][0], ('-es', '--pele_eq_steps'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[9][0], ('-miov', '--min_overlap'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[10][0], ('-maov', '--max_overlap'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[11][0], ('-tmp', '--temperature'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[12][0], ('-sd', '--seed'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[13][0], ('-st', '--steering'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[14][0], ('-trh', '--translation_high'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[15][0], ('-roth', '--rotation_high'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[16][0], ('-trl', '--translation_low'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[17][0], ('-rotl', '--rotation_low'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[18][0], ('-rad', '--radius_box'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[19][0], ('-dat', '--data'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[20][0], ('-doc', '--documents'))  # todo: check a better way

    def test_add_clustering_arguments(self):
        mock = Mock()

        ap._add_clustering_arguments(mock)

        assert_equal(len(mock.add_argument.call_args_list), 9)
        assert_equal(mock.add_argument.call_args_list[0][0], ('-dis', '--distcont'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[1][0], ('-ct', '--threshold'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[2][0], ('-e', '--epsilon'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[3][0], ('-cn', '--condition'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[4][0], ('-mw', '--metricweights'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[5][0], ('-ncl', '--nclusters'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[6][0], ('-pdbf', '--pdbout'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[7][0], ('-ban', '--banned'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[8][0], ('-lim', '--limit'))  # todo: check a better way

    def test_add_protocol_arguments(self):
        mock = Mock()

        ap._add_protocol_arguments(mock)

        assert_equal(mock.add_argument.call_args_list[0][0], ('-HT', '--highthroughput'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[1][0], ('-EX', '--explorative'))  # todo: check a better way
        assert_equal(mock.add_argument.call_args_list[2][0], ('--test',))  # todo: check a better way

    def test_add_output_format_arguments(self):
        mock = Mock()

        ap._add_output_format_arguments(mock)

        assert_equal(mock.add_argument.call_args[0][0], '--mae')  # todo: check a better way

    def test_add_other_arguments(self):
        mock = Mock()

        ap._add_others_arguments(mock)

        assert_equal(mock.add_argument.call_args[0][0], '--rename')  # todo: check a better way

    def test_check_highthroughput_in_args_true(self):
        mock = Mock()
        mock.highthroughput = True
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0

        ap._check_highthroughput_in_args(mock)

        assert_equal(mock.growing_steps, 1)
        assert_equal(mock.steps, 3)
        assert_equal(mock.pele_eq_steps, 10)

    def test_check_highthroughput_in_args_false(self):
        mock = Mock()
        mock.highthroughput = False
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0

        ap._check_highthroughput_in_args(mock)

        assert_equal(mock.growing_steps, 0)
        assert_equal(mock.steps, 0)
        assert_equal(mock.pele_eq_steps, 0)

    def test_check_test_in_args_true(self):
        mock = Mock()
        mock.test = True
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0
        mock.temp = 0

        ap._check_test_in_args(mock)

        assert_equal(mock.growing_steps, 1)
        assert_equal(mock.steps, 1)
        assert_equal(mock.pele_eq_steps, 1)
        assert_equal(mock.temp, 1000000)

    def test_check_test_in_args_false(self):
        mock = Mock()
        mock.test = False
        mock.growing_steps = 0
        mock.steps = 0
        mock.pele_eq_steps = 0
        mock.temp = 0

        ap._check_test_in_args(mock)

        assert_equal(mock.growing_steps, 0)
        assert_equal(mock.steps, 0)
        assert_equal(mock.pele_eq_steps, 0)
        assert_equal(mock.temp, 0)

    def test_parse_arguments(self):
        pass
