try:
    from skbuild import setup
except ImportError:
    raise ImportError('scikit-build is required for installing')

setup(
    name="traji",
    version="0.2.0",
    description="A C++/Python library designed for trajectory calculation in Cartesian / Frenet coordinate",
    long_description='(see project homepage)',
    author='Jacob Zhong',
    author_email='cmpute@gmail.com',
    url='https://github.com/cmpute/traji',
    download_url='https://github.com/cmpute/traji/archive/master.zip',
    license='MIT',
    packages=['traji'],
    install_requires=['numpy>=1.11'],
    setup_requires=['scikit-build'],
    # extras_require={'test': ['pytest']},
    # package_data={'pcl': [p.lstrip("pcl/") for p in (glob.glob("pcl/*.pyi") + glob.glob("pcl/**/*.pyi"))]},
    classifiers=[
        'Programming Language :: C++',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering'
    ],
    keywords=['trajectory', 'frenet', 'numpy', 'binding'],
)