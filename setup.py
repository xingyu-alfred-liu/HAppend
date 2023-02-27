from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
        name='pah101tools',
        version='0.1.0',
        description='Tools to reproduce plots in PAH101 dataset',
        long_description=long_description,
        long_description_content_type="text/markdown",
        author='Xingyu (Alfred) Liu',
        author_email='xingyu.alfred.liu@gmail.com',
        url="https://github.com/BLABABA/PAH101Plot.git",
        packages=find_packages(),
        install_requires=['pymatgen', 'ase', 'matplotlib'],
)