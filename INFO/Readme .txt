To use the program, you need to install:
	C++ https://sourceforge.net/projects/orwelldevcpp/
	Python 3 https://www.python.org/downloads/
	Anaconda https://repo.anaconda.com/archive/Anaconda3-2020.11-Windows-x86_64.exe (Add checkbox with PATH)

Install libraries using these commands:
	pip install plotly
	pip install plotly-express
	pip install tkintertable
	pip install numpy
	conda create -n geo_env
	conda activate geo_env
	conda config --env --add channels conda-forge
	conda config --env --set channel_priority strict
	conda install python=3 geopandas
	

To execute the program open OptiPoint.py using a Python. 
Then choose the file with the statistics (default statistic is located in program's folder “Statistic.csv”).
If you want to get ellipse parameters you can open file "FinalOutput.txt" which will appear after program's execution.
