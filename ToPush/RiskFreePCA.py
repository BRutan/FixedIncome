##############################################
# RiskFreePCA.py
##############################################
# Description:
# * Perform PCA using daily treasury data.

from bs4 import BeautifulSoup as Soup
import csv
from datetime import datetime
import os
from pandas import DataFrame
import re
import math 
import numpy as np

class RiskFreePCA(object):
    """
    * Perform principle component analysis on risk-free yield curve.
    """
    def __init__(self, filepath):
        """
        * Pull in html data at path.
        """
        if not isinstance(filepath, str):
            raise Exception('filepath must be a string.')
        elif '.' not in filepath:
            raise Exception('filepath must point to file.')
        elif not os.path.exists(filepath):
            raise Exception('File at filepath does not exist.')

        if not os.path.exists("TreasuryData.csv"):
            self.__PullData(filepath)
            self.__WriteDataToFile()
        else:
            self.__PullDataFromFile()

        self.__CleanData()

    @property
    def Data(self):
        return self.__data
    ##################
    # Interface Methods:
    ##################
    def SelectData(self, startDate, endDate):
        """
        * Pull all rows between two dates.
        """
        mask = (self.__data['Date'] >= startDate) & (self.__data['Date'] <= endDate)
        return self.__data.loc[mask]

    def PerformPCA(self, data, outputFile):
        """
        * Perform PCA on data.
        """
        targetCols = [key for key in self.__data if key != 'Date']
        rsquareMatrix = { key : [] for key in targetCols }
        rsquareMatrix['DepCols'] = []
        maxLen = 0
        for indepCol in targetCols:
            rsquareMatrix['DepCols'].append(indepCol)
            for depCol in targetCols:
                x_data = np.array(data[indepCol]).reshape(1, -1)
                y_data = np.array(data[depCol]).reshape(1, -1)
                y_data = y_data.squeeze()
                x_data = x_data.squeeze()
                mean_y = sum(y_data) / len(y_data)
                mean_x = sum(x_data) / len(x_data)
                cov = sum([(y - mean_y) * x_data[num] for num,y in enumerate(y_data)])
                var_x = sum([(x - mean_x) ** 2 for x in x_data])
                beta_1 = cov / var_x
                beta_0 = mean_y - mean_x * beta_1
                preds = [beta_0 + beta_1 * x for x in x_data]
                TSS = sum([(y - mean_y) ** 2 for y in y_data])
                r_2 = sum([(y_pred - mean_y) ** 2 for y_pred in preds]) / TSS
                rsquareMatrix[indepCol].append(r_2)
        results = DataFrame.from_dict(rsquareMatrix)
        results.set_index('DepCols', inplace=True)
        results.to_csv(outputFile)
        
    ##################
    # Private Helpers:
    ##################
    def __CleanData(self):
        """
        * Clean all rows containing non-numeric values in dataset.
        """
        targetCols = [key for key in self.__data if key != 'Date']
        numMatch = re.compile('[0-9]+(\.[0-9]+)')
        indices = set()
        for col in targetCols:
            for num, val in enumerate(self.__data[col]):
                if not numMatch.match(str(val)):
                    indices.add(num)

        self.__data = self.__data.drop(self.__data.index[list(indices)])

    def __PullData(self, filepath):
        soup = Soup(open(filepath, 'r'), 'lxml')
        colMatch = re.compile('d:bc_.+')
        dateMatch = re.compile('d:new_date')
        dateConv = re.compile('\d{4}-\d{2}-\d{2}')
        # Set all tenors using header columns:
        rows = soup.find_all('m:properties')
        tags = rows[0].find_all(colMatch)
        self.__headers = {}
        skipMaturities = { '30yeardisplay' : True, "1month" : True, "2month" : True }
        data = {"Date" : []}
        for tag in tags:
            match = re.search('d:bc_[0-9a-z]+', str(tag))[0]
            cleaned = re.sub('d:bc_', '', match) 
            if not cleaned in skipMaturities:
                self.__headers[cleaned] = match
                data[cleaned] = []
        for row in rows:
            date = row.find(dateMatch)
            date = dateConv.search(str(date))[0]
            date = datetime.strptime(date, '%Y-%m-%d')
            data["Date"].append(date)
            for header in self.__headers:
                col = row.find(self.__headers[header])
                if col.text:
                    data[header].append(float(col.text.strip()))
                else:
                    data[header].append('na')
        # Header is rate name, row index is for given date:
        self.__data = DataFrame.from_dict(data)
        #self.__data.set_index("Date", inplace=True)
    def __WriteDataToFile(self):
        """
        * Write all data to local file for easy access.
        """
        path = "TreasuryData.csv"
        self.__data.to_csv(path)

    def __PullDataFromFile(self):   
        data = {}
        numType = re.compile('[0-9]+(\.[0-9]+){0,1}')
        dateType = re.compile('[0-9]{1,2}//[0-9]{1,2}//[0-9]{4}')
        with open("TreasuryData.csv", 'r', newline='') as f:
            reader = csv.reader(f)
            numToCol = {}
            for num, row in enumerate(reader):
                for colNum, col in enumerate(row):
                    if num != 0 and colNum in numToCol:
                        if numToCol[colNum] == "Date":
                            val = datetime.strptime(col, '%Y-%m-%d')
                        elif numType.match(col):
                            val = float(col)
                        else:
                            val = col
                        data[numToCol[colNum]].append(val)    
                    elif num == 0:
                        if col.strip() and not col.isspace() and not re.match('^[0-9]$',col.strip()):
                            numToCol[colNum] = col.strip()
                            data[col] = []
            self.__data = DataFrame.from_dict(data)
            #self.__data.set_index("Date", inplace=True)