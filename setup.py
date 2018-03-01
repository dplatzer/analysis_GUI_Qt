from setuptools import setup, find_packages

setup(name='analysis_GUI_Qt',
      version='1.0',
	  description='The analysis_GUI_Qt Python library',
	  long_description=open('README.md').read(),
	  author='Dominique Platzer',
	  url="https://github.com/dplatzer/analysis_GUI_Qt",
	  install_requires=["pyqt>=5.6"],
	  packages=find_packages(),
	  include_package_data=True)