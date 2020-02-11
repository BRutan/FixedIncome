##############################################
# Bootstrapping.py
##############################################
# Description:
# * Bootstrap yield curve for given credit quality
# using bond prices and attributes.

from FixedIncome import Bond, ZeroCurve, InterpolationMethods, Interpolation
import csv
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from enum import Enum
import pandas
import os
from sortedcontainers import SortedList

class BootStrapper(object):
    """
    * Produces a bootstrapped risk-free yield curve using the Monotone Convex algorithm and
    risk-free bonds.
    """
    def __init__(self, path, method):
        """
        * Generate new Bootstrapper.
        Inputs:
        * path: Path to bond containing bonds trading in market, containing 
        prerequisite columns (string).
        * method: Interpolation method (InterpolationMethods enum). 
        """
        self.__settlement = None
        self.__interpmethod = method
        group = InputBonds(path)
        self.__settlement = group.Settlement
        self.__sorted = SortedList(group.Bonds.values())

    @property
    def Settlement(self):
        return self.__settlement

    def GenerateSpotCurve(self):
        """
        * Generate zero curve incorporating interpolation method. 
        """
        if not self.__sorted:
            raise Exception('Need bonds to bootstrap.')
        bonds = self.__sorted
        # Map tenor to zero discount factor curve:
        yieldCurve = {}
        numBonds = len(bonds)
        method = self.__interpmethod
        num = 0
        # Find and remove all zero coupon bonds, swaps:
        # We are assuming that all bonds have same settlement date (bootstrapping wouldn't make sense otherwise).
        settle = bonds[0].Settlement
        zeroCount = 0
        swaps = []
        while num < numBonds:
            if bonds[num].Type == 'swap':
                swaps.append(bonds[num])
                num += 1
            elif len(bonds[num].CashFlows.keys()) == 1 and bonds[num].YTM < .2:
                zeroCount += 1
                bond = bonds[num]
                yieldCurve[bond.Tenor] = bond.YTM
                del bonds[num]
                numBonds -= 1
            else:
                num += 1
        if zeroCount == 0:
            raise Exception('Need at least one zero coupon bond.')
        for swap in swaps:
            self.__BackoutSwap(swap, yieldCurve, method)
        # Perform the bootstrapping using coupon bonds/swaps:
        for bond in bonds:
            if bond.Type == 'swap':
                # Perform bootstrapping using swaps:
                self.__BackoutSwap(bond, yieldCurve, method)
            elif bond.TenorToCashFlow:
                # Generate present value using known zero rates:
                pv = 0
                cashFlows = bond.TenorToCashFlow
                price = bond.Price
                finalTenor = max(cashFlows.keys())
                for num, tenor in enumerate(cashFlows.keys()):
                    # Stop bootstrapping if cash flow date occurs after observed zero curves:
                    if num >= len(yieldCurve.keys()) or tenor > max(yieldCurve.keys()):
                        break
                    elif tenor not in yieldCurve:
                        # Interpolate using current spot curve:
                        interpolator = Interpolation(method, yieldCurve)
                        _yield = interpolator.Interpolate(tenor)
                    else:
                        _yield = yieldCurve[tenor]
                    pv += cashFlows[tenor] * (1 + _yield) ** (-1 * tenor)
                num += 1
                # Skip bond if cannot bootstrap:
                if len(cashFlows.keys()) - num == 1:
                    # Calculate present value of remaining cash flows using price:
                    terminalPmt = cashFlows[finalTenor]
                    _yield = (((price - pv) / terminalPmt) ** (-1 / tenor)) - 1
                    yieldCurve[finalTenor] = _yield
        
        yields = yieldCurve
        discountFactors = ZeroCurve.ConvertZeroYieldsToDFs(yieldCurve)
            
        return ZeroCurve(discountFactors, yields, self.__interpmethod, self.__settlement)

    def __BackoutSwap(self, swap, yieldCurve, method):
        """
        * Generate zero rates using swap rate.
        """        
        # Skip if bond's tenor is greater than max tenor in yield curve:
        if swap.Tenor > max(yieldCurve.keys()):
            return
        # Fixed swap rate goes into Coupon column:
        fixedRate = swap.Coupon
        # We assume that first fixing occurs at NextPaymentDate:
        dates = list(swap.CashFlows.keys())
        pmtfreq = swap.PmtFreq
        targetTenor = Bond.CalculateTenor(dates[0], dates[1])
        if swap.Tenor not in yieldCurve:
            # Interpolate to get appropriate zero coupon yield factor:
            interpolator = Interpolation(method, yieldCurve)
            _yield = interpolator.Interpolate(swap.Tenor)
        else:
            _yield = yieldCurve[swap.Tenor]
        if '1.26' in str(_yield) or _yield > .2:
            _yield == _yield
        df = (1 + _yield) ** -swap.Tenor
        # Backout target zero from swap rate:
        df_target = (1 - pmtfreq * fixedRate * df) / (1 + pmtfreq * fixedRate)
        yieldCurve[targetTenor] = (df_target ** -targetTenor) - 1

class InputBonds(object):
    """
    * Stores zero/coupon bonds for use in bootstrapping.
    """
    __ReqColumns = { 'name' : False, 'price' : False, 'maturity' : False, 'settlement' : False, 'coupon' : False, 
                    'freq' : False, 'facevalue' : False, 'type' : False, 'daycount' : False }
    __OptColumns = { 'nextpmtdate' : False, 'ytm' : False }
    def __init__(self, path):
        """
        * Load all bonds into the object from file at path.
        """
        errs = []
        if not isinstance(path, str):
            errs.append('path must be a string.')
        elif not os.path.exists(path):
            errs.append('File at path does not exist.')
        elif '.csv' not in path:
            errs.append('path must point to csv file.')
        if errs:
            raise Exception('\n'.join(errs))
        self.__maxorder = 0
        # Map { Name -> Bond }:
        self.__settlement = None
        self.Bonds = {}
        self.__LoadBonds(path)
        if self.Bonds:
            # Assumes that settlement dates are uniform in file:
            first = list(self.Bonds.keys())[0]
            self.__settlement = self.Bonds[first].Settlement

    @property
    def Settlement(self):
        return self.__settlement
        
    ###################
    # Interface Methods:
    ###################
    def CalculatePriceRisks(self, orders):
        """
        * Calculate all price risks.
        """
        if not isinstance(orders, list):
            raise Exception('orders must be a list of positive numeric values.')
        elif [order for order in orders if not isinstance(order, (float, int))]:
            raise Exception('all orders must be positive numeric values.')
        self.__maxorder = max(orders)
        for name in self.Bonds:
            self.Bonds[name].CalculatePriceRisks(orders)

    def PrintAllAttributes(self, path = None):
        """
        * Print all attributes to local file.
        """
        if not path:
            path = 'BondAttributes.csv'
        elif not isinstance(path, str):
            raise Exception('path must be a string if provided.')
        elif not os.path.exists(path[0:path.rfind('\\') + 1]):
            raise Exception('folder containg target file at path does not exist.')
        elif '.' in path and os.path.exists(path):
            raise Exception('File at path already exists.')
        elif not '.' in path:
            path = ''.join([path[0:path.rfind('\\') + 1], 'BondAttributes.csv'])
        reqCols = InputBonds.RequiredColumns()
        columns = list(reqCols)
        columns.append('ytm')
        columns = [col.capitalize() for col in columns]
        allorders = range(1, self.__maxorder + 1) if self.__maxorder > 0 else []
        if allorders:
            endCol = len(columns)
            for order in orders: 
                if order == 1:
                    postfix = 'st'
                elif order == 2:
                    postfix = 'nd'
                elif order == 3:
                    postfix = 'rd'
                else:
                    postfix = 'th'
                columns.append(''.join([str(order),postfix,' Order Term']))
        with open(path, 'w', newline = '') as f:
            writer = csv.writer(f)
            # Write column headers:
            writer.writerow(columns)
            for name in self.Bonds:
                bond = self.Bonds[name]
                row = []
                for colNum in range(0, endCol):
                    col = columns[colNum].lower()
                    val = None
                    if col == 'name':
                        val = name
                    elif col == 'daycount':
                        val = bond.GetAttr(col).ToString()
                    elif col in reqCols:
                        val = bond.GetAttr(col)
                        if isinstance(val, (date, datetime)):
                            val = val.strftime('%m/%d/%Y')
                    elif col == 'ytm':
                        val = bond.GetAttr('ytm')
                    row.append(str(val))
                for order in allorders:
                    val = ''
                    if order in bond.PriceRisks:
                        val = bond.PriceRisks[order]
                    row.append(str(val))
                writer.writerow(row)

    ###################
    # Private Helpers:
    ###################
    def __LoadBonds(self, path):
        req = InputBonds.RequiredColumns()
        opt = InputBonds.OptionalColumns()
        errs = set()
        columns = {}
        # Load all bonds from file:
        with open(path, 'r') as f:
            reader = csv.reader(f)
            data = { val : 0 for val in InputBonds.__ReqColumns }
            for rowNum, row in enumerate(reader):
                if rowNum != 0:
                    # Pull in bond:
                    for colNum, cellvalue in enumerate(row):
                        if not cellvalue.strip():
                            break
                        data[columns[colNum]] = cellvalue.strip().lower()
                    try:
                        bond = Bond(data)
                        self.Bonds[bond.Name] = bond 
                    except BaseException as ex:
                        errs.add(str(ex))
                else:
                    # Determine if all correct headers were entered:
                    for colNum, col in enumerate(row):
                        lowered = col.lower()
                        if lowered in req:
                            req[col] = True
                            columns[colNum] = lowered
                        elif lowered in opt:
                            columns[colNum] = lowered
                    missingCols = []
                    for key in req:
                        if not req[key]:
                            missingCols.append(pair[0])
                    if missingCols:
                        raise Exception(''.join(['The following columns are missing from the file:', ','.join(missingCols)]))
                    numCols = max(columns.keys()) + 1
    ################
    # Static Methods:
    ################
    @classmethod
    def GenerateSampleFile(cls, path = None):
        """
        * Generate sample file to path. If path not provided then
        generate to local file with predetermined name.
        """
        if path and not isinstance(path, str):
            raise Exception('path must be a string or None.')
        elif not path:
            path = 'SampleBondFile.csv'
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(InputBonds.__ReqColumns.keys())
    @staticmethod
    def RequiredColumns():
        return InputBonds.__ReqColumns.copy()
    @staticmethod
    def OptionalColumns():
        return InputBonds.__OptColumns.copy()

def test():
    """
    * Test the above classes.
    """
    InputBonds.GenerateSampleFile("BootstrapperInputSampleFile.csv")


if __name__ == '__main__':
    test()