import os

import pytest
from SeqPyPlot.seqpyplot.analyzer.paired_sample_filter import \
    PairedSampleFilter
from SeqPyPlot.seqpyplot.container.data_container import DataContainer
from SeqPyPlot.seqpyplot.parsers.config_parser import config_parser
import pandas as pd
from pandas.testing import assert_frame_equal

# @pytest.fixture(scope='module')
# def load_container():

#     conf = os.path.join('tests', 'test_config.ini')
#     config_obj = config_parser(conf)
#     dc = DataContainer(config_obj)
#     return dc, config_obj

@pytest.fixture(scope='module')
def load_config():

    conf = os.path.join('tests', 'test_config.ini')
    return config_parser(conf)


def test_load_config(load_config):

    conf_obj = load_config
    psf = PairedSampleFilter(conf_obj)

    assert True

def test_count_by_gene(load_config):

    conf_obj= load_config
    psf = PairedSampleFilter(conf_obj)

    
    fake_df = [pd.DataFrame({
        'a': [1, 2, 3, 4],
        'b': [1, 2, 3, 4]
    }, index=['a', 'b', 'c', 'd'])]

    res = psf.count_by_gene(fake_df)
    assert sorted(res.items()) == [('a', 1),
                                   ('b', 1),
                                   ('c', 1),
                                   ('d', 1)]

def test_apply_fold_change(load_config):

    conf_obj= load_config

    psf = PairedSampleFilter(conf_obj, log2fold=1.0)
    psf.file_pairs = [('a1', 'b1'), ('a2', 'b2')]
    
    fake_df_list = [
        pd.DataFrame({
            'a1': [1, 1, 10, 10],
            'b1': [1, 4 ,1, 10]}, index=[0, 1, 2, 3]),
        pd.DataFrame({
            'a2': [1, 2, 3, 4],
            'b2': [1, 5, 3, 4], }, index=[0, 1, 2, 3])
        ]
    
    expected = [
        pd.DataFrame({
            'a1': [1, 10],
            'b1': [4, 1]
        }, index=[1, 2]),

        pd.DataFrame({
            'a2': [2],
            'b2': [5]
        }, index=[1])
    ]

    res = psf.apply_fold_change(fake_df_list)
    for df, exp in zip(res, expected):
        assert_frame_equal(df, exp)


def test_apply_diff(load_config):
    conf_obj= load_config

    psf = PairedSampleFilter(conf_obj, diff=[2.0, 6.0])
    psf.file_pairs = [('a1', 'b1'), ('a2', 'b2')]
    
    fake_df_list = [
        pd.DataFrame({
            'a1': [1, 1, 10, 10],
            'b1': [1, 4 ,1, 8]}, index=[0, 1, 2, 3]),
        pd.DataFrame({
            'a2': [1, 2, 3, 4],
            'b2': [1, 50, 5, 7], }, index=[0, 1, 2, 3])
        ]

    res = psf.apply_diff(fake_df_list)

    expected = [
        pd.DataFrame({
            'a1': [1, 10],
            'b1': [4, 8]
        }, index=[1, 3]),

        pd.DataFrame({
            'a2': [3, 4],
            'b2': [5, 7]
        }, index=[2, 3])
    ]

    for df, exp in zip(res, expected):
        assert_frame_equal(df, exp)


def test_apply_low(load_config):
    conf_obj= load_config

    psf = PairedSampleFilter(conf_obj, low=10.)
    psf.file_pairs = [('a1', 'b1'), ('a2', 'b2')]
    
    fake_df_list = [
        pd.DataFrame({
            'a1': [1, 1, 10, 0 ],
            'b1': [1, 4, 1,  20]}, index=[0, 1, 2, 3]),
        pd.DataFrame({
            'a2': [1, 2, 13, 4],
            'b2': [1, 50, 5, 7], }, index=[0, 1, 2, 3])
        ]

    expected = [
        pd.DataFrame({
            'a1': [10, 0],
            'b1': [1, 20]
        }, index=[2, 3]),

        pd.DataFrame({
            'a2': [2, 13],
            'b2': [50, 5]
        }, index=[1, 2])
    ]

    res = psf.apply_low(fake_df_list)
    for df, exp in zip(res, expected):
        assert_frame_equal(df, exp)

def test_apply_hi(load_config):
    conf_obj= load_config

    psf = PairedSampleFilter(conf_obj, hi=10.)
    psf.file_pairs = [('a1', 'b1'), ('a2', 'b2')]
    
    fake_df_list = [
        pd.DataFrame({
            'a1': [1, 11, 11, 10 ],
            'b1': [1, 4, 1,  12]}, index=[0, 1, 2, 3]),
        pd.DataFrame({
            'a2': [100, 2, 11, 14],
            'b2': [311, 5, 15, 7], }, index=[0, 1, 2, 3])
        ]

    expected = [
        pd.DataFrame({
            'a1': [1, 11, 11],
            'b1': [1, 4, 1]
        }, index=[0, 1, 2]),

        pd.DataFrame({
            'a2': [2, 14],
            'b2': [5, 7]
        }, index=[1, 3])
    ]

    res = psf.apply_hi(fake_df_list)
    for df, exp in zip(res, expected):
        assert_frame_equal(df, exp)
