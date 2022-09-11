# OptiPoint
Program for finding optimal point based on popularity density.
## User Interface
![image](https://user-images.githubusercontent.com/49571325/189542371-4ec01f5f-a206-4348-a7f9-8c2f03564f70.png)
![image](https://user-images.githubusercontent.com/49571325/189542378-9a52e15a-fe3e-44f4-b4b2-94ec6bc70f76.png)

## To use the program, you need to install:
	Python 3 https://www.python.org/downloads/
	Anaconda https://repo.anaconda.com/archive/Anaconda3-2020.11-Windows-x86_64.exe (Add checkbox with PATH)
## Install libraries using these commands:
	pip install plotly
	pip install plotly-express
	pip install tkintertable
	pip install numpy
	conda create -n geo_env
	conda activate geo_env
	conda config --env --add channels conda-forge
	conda config --env --set channel_priority strict
	conda install python=3 geopandas
	
To execute the program open OptiPoint.py using a Python. Then choose the file with the statistics (default statistic is located in program's folder “Statistic.csv”). If you want to get ellipse parameters you can open file "FinalOutput.txt" which will appear after program's execution.
