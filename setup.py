from setuptools import setup, find_packages
# include=['spacerplacer', 'run_experiments', 'model', 'model.*']
setup(
    name='spacerplacer',
    version='1.0.0',
    packages=find_packages(),
    package_data={'model': ['mafft_scripts/mafft-linux64/*', 'mafft_scripts/mafft-mac/*', 'mafft_scripts/mafft-win/*',
                            'bdm_likelihood_computation/sympy_bdm_lh_fct/*'],
                  },
    py_modules=['spacerplacer'],
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'spacerplacer = spacerplacer:main'
        ]
    }
)