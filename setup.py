from setuptools import setup, find_packages
# include=['spacerplacer', 'run_experiments', 'model', 'model.*']
setup(
    name='spacerplacer',
    version='1.0.0',
    packages=find_packages(),
    package_data={'model': ['mafft_scripts/mafft-linux64/*', 'mafft_scripts/mafft-mac/*', 'mafft_scripts/mafft-win/*']},
    py_modules=['spacerplacer', 'run_experiments', 'input_parser'],
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'spacerplacer = spacerplacer:main'
        ]
    }
)