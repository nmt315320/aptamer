# aptamer

#License
Copyright (C) 2022 Mengting Niu(yunzeer@gmail.com)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/.


#Type: Package

Files: 1.data

data.xlsx-DNA VS aptamer datasets
ASCpeck.fasta-ASC-speck sequences


2.code

2-1. python code----aptamer identification

Aptagen Spider-aptamer crawler, extract data

Aptamer Data Cleanup.ipynb-python code-aptamer data clean

getData.py 

AnalyseFASTA.py

basic_units.py

Deal_Kmer.py

DProcess.py

attention.py

RicEnnTool.py

train.py(main code)

2-2. matlab code -----new aptamer design

BFOA_my.m(main code)

The tool is developed for aptamer identification using deep hierarchical network
![image](https://github.com/nmt315320/aptamer/blob/main/Architecture.png)
# Requirements
- Python 3.7 (64-bit)
- Keras 2.2.0 in Python
- TensorFlow-GPU 1.14.0 in Python
- Numpy 1.18.0 in Python
- Gensim 3.8.3
- Ubuntu 18.04 (64-bit)
# Usage

command: python train.py 
MATLAB BFOA_my.m


You can train the model of 5-fold cross-validation with a very simple way by the command blow:  
*Python train.py* . The script of if **name == "main"** calls training process which trains several models of each model type and finds the best set of hyperparameters. The main function then trains the models several times (num_final_runs) and saves the best model.


The prediction results will be displayed automatically. If you need to save the results, please specify the path yourself. Thank you and enjoy the tool!

 If you have any suggestions or questions, please email me at *yunzeer@gmail.com*.
