##############################################
# PerformPCA.py
##############################################
# Description:
# * Perform principal component analysis
# on daily treasury data to determine
# relationships between treasury yields to maturity.

import datetime
from RiskFreePCA import RiskFreePCA
from argparse import ArgumentParser
import os

def GeneratePrincipleComponents():
    parser = ArgumentParser()
    parser.add_argument("datapath", type=str, help="Path to .html file containing yield curve data.")
    args = parser.parse_args()

    if not os.path.exists(args.datapath):
        raise Exception("File does not exist at datapath.")

    # Pull in data from local html file:
    pcaGen = RiskFreePCA(args.datapath)
    # Perform PCA for 1998-2008 period:
    formatStr = '%Y-%m-%d'
    startDate = datetime.datetime.strptime('1998-01-01', formatStr)
    endDate = datetime.datetime.strptime('2008-12-31', formatStr)
    data = pcaGen.SelectData(startDate, endDate)
    pcaGen.PerformPCA(data, "PCA_1998_2008.csv")
    # Perform PCA for 2009-2019 period:
    startDate = datetime.datetime.strptime('2009-01-01', formatStr)
    endDate = datetime.datetime.strptime('2019-12-31', formatStr)
    data = pcaGen.SelectData(startDate, endDate)
    pcaGen.PerformPCA(data, "PCA_2009_2019.csv")

if __name__ == '__main__':
    GeneratePrincipleComponents()




