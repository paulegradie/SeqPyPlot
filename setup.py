from setuptools import setup, find_packages
# from distutils.core import setup

setup(
    name='SeqPyPlot',
    version='0.4.0',
    author='Paul E Gradie',
    author_email='paul.e.gradie@gmail.com',
    packages=find_packages(exclude=['tests']),
    scripts=['bin/SeqPyPlot.py'],
    #url='http://pypi.python.org/pypi/seqpyplot',
    license='LICENSE.txt',
    description='Package for descriptive analysis.',
    long_description=open('docs/README.txt').read(),
    install_requires=[
        "pandas == 0.22.0",
        "numpy == 1.14.3",
        "pathos",
        "seaborn",
        "sklearn",
        "tqdm"
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

)

