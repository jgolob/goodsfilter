from setuptools import setup

setup(
    name='goodsfilter',
    version='0.1.0',
    description="""Filter DADA2-style sequence tables using
    Good's coverage saturation to remove likely contaminant and PCR error sequence variants
    """,
    url='https://github.com/jgolob/goodsfilter/',
    author='Jonathan Golob',
    author_email='j-dev@golob.org',
    license='MIT',
    packages=['goodsfilter'],
    zip_safe=False,
    install_requires=[
    ],
    entry_points={
        'console_scripts': [
            'goodsfilter=goodsfilter.goodsfilter:main',
        ],
    }
)
