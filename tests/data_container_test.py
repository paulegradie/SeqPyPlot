import pytest
import tempfile
from ConfigParser import ConfigParser

from .SeqPyPlotLib.container.data_container import DataContainer

# This is a compromise in testing ediquette
from SeqPyPlotLib.parsers.config_parser import config_parser
CONFIG = config_parser('./test_config.ini')



def test_split():

    dc = DataContainer(CONFIG)
    result = dc.split()

    