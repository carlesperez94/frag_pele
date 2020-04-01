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
        pass

    def _test_add_plop_arguments(self):
        pass

    def test_add_pele_conf_arguments(self):
        pass

    def test_add_clustering_arguments(self):
        pass

    def test_add_protocol_arguments(self):
        pass

    def test_add_output_format_arguments(self):
        pass

    def test_add_other_arguments(self):
        pass

    def test_check_highthroughput_in_args(self):
        pass

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
