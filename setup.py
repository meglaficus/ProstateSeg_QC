from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ProstateSeg_QC',
      packages=['psqc_tools'],
      py_modules=['psqc'],
      package_dir={'': 'psqc'},
      version='0.0.2',
      description='Automated quality control algorithm for prostate segmentation',
      long_description=long_description,
      url='https: // github.com/meglaficus/ProstateSeg_QC',
      author='Jakob Meglič, under MR Cancer group, NTNU',
      author_email='jakobmeglic123@gmail.com',
      license='Apache License 2.0',
      classifiers=['Programming Language :: Python :: 3',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Image Processing',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: OS Independent'],
      install_requires=[
              "tqdm",
              "scipy",
              "numpy",
              "SimpleITK",
              "pandas",
              "connected-components-3d",
              "regex"],
      python_requires=">=3.6",
      keywords=['prostate', 'segmentation', 'quality control', 'automated']
      )
