
from setuptools import setup, find_packages

setup(
    name='hoppy',
    version='0.0.2',
    packages=find_packages(exclude=['tests*']),
    license='MIT',
    description=' Libararies and scripts for X-ray analyses of High-energy Observatory Pipelines of PYthon (HOPPY).',
    long_description=open('README.md').read(),
    install_requires=['numpy','pyyaml'],
    url='https://github.com/tenoto/hoppy',
    author='Teru Enoto',
    author_email='teruaki.enoto@gmail.com'
)
