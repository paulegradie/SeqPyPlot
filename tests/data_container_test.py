import tempfile
from collections import OrderedDict
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import os

from pandas.testing import assert_frame_equal

from SeqPyPlot.seqpyplot.container.data_container import DataContainer
from SeqPyPlot.seqpyplot.parsers.config_parser import config_parser


def foo(a, b):
    return a + b

def test_foo():

    res = foo(1, 1)
    assert res == 2

@pytest.fixture(scope='module')
def load_container():

    conf = os.path.join('tests', 'test_config.ini')
    config_obj = config_parser(conf)
    dc = DataContainer(config_obj)
    return dc

def test_config_load():

    conf = os.path.join('tests', 'test_config.ini')
    config_obj = config_parser(conf)
    dc = DataContainer(config_obj)
    assert True


def test_split(load_container):

    fake_file_pairs = [('a1', 'a2'), ('b1', 'b2')]
    fake_df = pd.DataFrame({
        'a1': [1, 2, 3, 4],
        'a2': [1, 2 ,3, 4],
        'b1': [1, 2, 3, 4],
        'b2': [1, 2, 3, 4]
    })

    dc = load_container
    dc.file_pairs = fake_file_pairs
    res = dc.split(fake_df)

    expected = [
        pd.DataFrame({
            'a1': [1, 2, 3, 4],
            'a2': [1, 2 ,3, 4]}),
        pd.DataFrame({
            'b1': [1, 2, 3, 4],
            'b2': [1, 2, 3, 4]})
        ]

    assert len(res) == 2
    for framea, frameb in zip(res, expected):
        assert_frame_equal(framea, frameb)

def test_make_col_pairs(load_container):
    
    dc = load_container
    dc.num_file_pairs = 2

    expected_file_pairs = [('a1', 'b1'), ('a2', 'b2')]
    fake_df = pd.DataFrame({
        'a1': [1, 2, 3, 4],
        'a2': [1, 2 ,3, 4],
        'b1': [1, 2, 3, 4],
        'b2': [1, 2, 3, 4]
    })

    res = dc.make_col_pairs(fake_df)

    assert list(res) == expected_file_pairs


def test_normalize_file_pairs(load_container):

    dc = load_container

    fake_df = pd.DataFrame({
        'a1': [10, 0, 1, 4],
        'a2': [1,  0, 3, 4],
        'b1': [1,  0, 1, 4],
        'b2': [5,  0, 3, 0]
    })
    
    exp_dict = OrderedDict({
        'a1': [10.0, 0.0, 1.0, 4.0],
        'b1': [0.650138,  0.000000, 0.650138, 2.600550],
        'a2': [1.0,  0.0, 3.0, 4.0],
        'b2': [6.568001, 0.000000, 3.940801, 0.000000]
    })

    expected = pd.DataFrame(exp_dict, columns=exp_dict.keys())
    res = dc.normalize_file_pairs(fake_df)

    assert_frame_equal(res, expected)


@patch('SeqPyPlot.seqpyplot.container.data_container.TMM')
def test_execute_normalization(mock_tmm, load_container):

    dc = load_container
    fake_df = pd.DataFrame({
        'a1': [10, 0, 1, 4],
        'a2': [1,  0, 3, 4],
        'b1': [1,  0, 1, 4],
        'b2': [5,  0, 3, 0]
    })
    dc.execute_normalization(fake_df)
    mock_tmm.assert_called_once()

def test_merge_dfs(load_container):
    dc = load_container

    fake_df1 = pd.DataFrame({
        'a1': [10, 1, 4],
        'a2': [1, 3, 4],
    }, index=[0, 2, 3])

    fake_df2 = pd.DataFrame({
        'b1': [1,  0, 4],
        'b2': [5,  0, 0]
    }, index=[0, 1, 3])

    res = dc.merge_dfs([fake_df1, fake_df2])
    expected = pd.DataFrame({
        'a1': [10.0, np.NaN, 1.0, 4.0],
        'a2': [1.0,  np.NaN, 3.0, 4.0],
        'b1': [1.0,  0.0, np.NaN, 4.0],
        'b2': [5.0,  0.0, np.NaN, 0.0]
    })

    assert_frame_equal(res, expected)

def test_reorder_cols(load_container):

    dc = load_container
    fake_df = pd.DataFrame({
        'a1': [10, 0, 1, 4],
        'b1': [1,  0, 3, 4],
        'a2': [1,  0, 1, 4],
        'b2': [5,  0, 3, 0]
    })
    res = dc.reorder_cols(fake_df)
    assert res.columns.tolist() == ['a1', 'b1', 'a2', 'b2']