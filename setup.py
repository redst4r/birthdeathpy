from setuptools import setup, find_packages

setup(name='birthdeathpy',
    version=0.1,
    description='Simulating a simple discrete time birth death process',
    url='http://github.com/redst4r/birthdeathpy/',
    author='redst4r',
    maintainer='redst4r',
    maintainer_email='redst4r@web.de',
    license='GNU GPL 3',
    keywords='birth death',
    packages=find_packages(),
    install_requires=[
        'toolz',
        'numpy',
        'tqdm',
        'scipy',
        'pandas',
        'fire',
        'matplotlib',
        'seaborn'
        ],
    entry_points = {
        'console_scripts': ['birthdeath=birthdeathpy.cli:main'],
    },
    zip_safe=False)
