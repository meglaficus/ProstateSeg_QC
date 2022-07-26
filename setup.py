from setuptools import setup, find_namespace_packages

setup(name='ProstateSeg_QC',
      packages=find_namespace_packages(
          include=["pqc", "pqc.*"]),
      version='0.2',
      description='Automated quality control algorithm for prostate segmentation',
      url='https: // github.com/meglaficus/ProstateSeg_QC',
      author='Jakob Meglič, under MR Cancer group, NTNU',
      author_email='jakobmeglic123@gmail.com',
      license='Apache',
      install_requires=[
              "tqdm",
              "scipy",
              "numpy",
              "SimpleITK",
              "pandas",
              "connected-components-3d",
              "regex"
      ],
      keywords=['prostate', 'segmentation', 'quality control', 'automated']
      )
