import os

import pytest
from SeqPyPlot.seqpyplot.analyzer.paired_sample_filter import \
    PairedSampleFilter
from SeqPyPlot.seqpyplot.container.data_container import DataContainer
from SeqPyPlot.seqpyplot.parsers.config_parser import config_parser
import pandas as pd

@pytest.fixture(scope='module')
def load_container():

    conf = os.path.join('tests', 'test_config.ini')
    config_obj = config_parser(conf)
    dc = DataContainer(config_obj)
    return dc, config_obj

def test_load_config(load_container):

    dc, conf_obj= load_container
    psf = PairedSampleFilter(conf_obj, dc)

    assert True

def test_count_by_gene(load_container):

    dc, conf_obj= load_container
    psf = PairedSampleFilter(conf_obj, dc)

    
    fake_df = [pd.DataFrame({
        'a': [1, 2, 3, 4],
        'b': [1, 2, 3, 4]
    }, index=['a', 'b', 'c', 'd'])]

    res = psf.count_by_gene(fake_df)
    assert sorted(res.items()) == [('a', 1),
                                   ('b', 1),
                                   ('c', 1),
                                   ('d', 1),]