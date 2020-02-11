##############################################
# FixedIncome.py
##############################################
# Description:
# * Objects for modeling fixed income.

from calendar import monthrange
import csv
from datetime import datetime, date, timedelta
from dateutil.relativedelta import relativedelta
import dateutil.rrule as dategen
from enum import Enum
import pandas
from pandas import DataFrame, concat
import numpy as np
import os
from scipy.interpolate import interp1d, CubicSpline

    
class InterpolationMethods(Enum):
    LINEAR = 1
    CUBIC = 2

class BondGroup(object):
    """
    * Stores zero/coupon/swaps from file using predetermined column headers, for fair value pricing .
    """
    __ReqColumns = { 'name' : False, 'maturity' : False, 'settlement' : False,  
                    'maturity' : False, 'coupon' : False, 'freq' : False, 'facevalue' : False, 'type' : False, 'daycount' : False }
    __OptColumns = { 'price' : False, 'nextpmtdate' : False, 'ytm' : False }
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
        self.Bonds = {}
        # Map { Name -> Bond }, price using zero curve:
        self.FairValues = {}
        self.__LoadBonds(path)
        
    ###################
    # Interface Methods:
    ###################
    def CalculateFairValues(self, zeroCurve):
        """
        * Calculate fair values of all bonds using passed zero curve.
        """
        for name in self.Bonds:
            self.FairValues[name] = self.Bonds[name].CalculateFairValue(zeroCurve)

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

    def PrintPriceComparison(self, outFile = None):
        """
        * Print document detailing if bond is overvalued or not.
        """
        if not outFile:
            outFile = 'FairValues.csv'
        with open(outFile, 'w', newline = '') as f:
            writer = csv.writer(f)

            writer.writerow(['Name', 'Coupon', 'Settlement', 'Maturity', 'QuotedPrice', 'FairValue'])
            for name in self.Bonds:
                bond = self.Bonds[name]
                row = [name, bond.Coupon, bond.Settlement, bond.Maturity, bond.Price, self.FairValues[name]]
                writer.writerow(row)

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
        reqCols = BondGroup.RequiredColumns()
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
        req = BondGroup.RequiredColumns()
        optional = BondGroup.OptionalColumns()
        errs = []
        columns = {}
        # Load all bonds from file:
        with open(path, 'r') as f:
            reader = csv.reader(f)
            data = { val : 0 for val in BondGroup.__ReqColumns if val != 'name' }
            for rowNum, row in enumerate(reader):
                if rowNum != 0 and [col for col in row if col]:
                    # Pull in bond:
                    for colNum, cellvalue in enumerate(row):
                        if colNum in columns:
                            data[columns[colNum]] = cellvalue.strip().lower()
                    try:
                        bond = Bond(data)
                        self.Bonds[bond.Name] = bond 
                    except BaseException as ex:
                        errs.append(str(ex))
                else:
                    # Determine if all correct headers were entered:
                    for colNum, col in enumerate(row):
                        lowered = col.lower()
                        if lowered in req:
                            req[col] = True
                            columns[colNum] = lowered
                        elif lowered in optional:
                            columns[colNum] = lowered
                    missingCols = []
                    for key in req:
                        if not req[key]:
                            missingCols.append(key)
                    if missingCols:
                        raise Exception(''.join(['The following columns are missing from the file:', ','.join(missingCols)]))
                    numCols = max(columns.keys()) + 1
            if errs:
                raise Exception('\n'.join(errs))
    ################
    # Static Methods:
    ################
    @staticmethod
    def GenerateSampleFile(path = None):
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
            writer.writerow(BondGroup.__ReqColumns.keys())
    @staticmethod
    def RequiredColumns():
        return BondGroup.__ReqColumns.copy()
    @staticmethod
    def OptionalColumns():
        return BondGroup.__OptColumns.copy()

class Bond(object):
    """
    * Store all attributes for bond.
    """
    __reqArgs = { 'name' : False, 'maturity' : False, 'settlement' : False, 'coupon' : False, 'freq' : False, 'facevalue' : False, 
                 'daycount' : False, 'type' : False }
    __argToSetterMap = { 'name' : '_Bond__Name', 'facevalue' : '_Bond__FaceValue', 'coupon' : '_Bond__Coupon', 'settlement' : '_Bond__Settlement', 
                'maturity' : '_Bond__Maturity', 'freq' : '_Bond__PmtFreq', 'price' : '_Bond__Price', 
                'nextpmtdate' : '_Bond__NextPaymentDate', 'daycount' : '_Bond__DayCount', 'type' : '_Bond__Type' }
    __argToAttrMap = { 'name' : '_Bond__name', 'facevalue' : '_Bond__facevalue', 'coupon' : '_Bond__coupon', 'settlement' : '_Bond__settlement', 
                'maturity' : '_Bond__maturity', 'freq' : '_Bond__pmtfreq', 'price' : '_Bond__price', 
                'nextpmtdate' : '_Bond__nextpmtdate', 'cashflows' : '_Bond__cashflows', 
                'ytm' : '_Bond__ytm', 'accruedint' : '_Bond__accruedinterest', 'daycount' : '_Bond__daycount' }
    __optionalArgs = { 'cashflows' : '_Bond__TenorToCashFlow', 'ytm' : '_Bond__YTM', 'nextpmtdate' : '_Bond__NextPaymentDate', 
                      'accruedint' : '_Bond__AccruedInterest', 'price' : '_Bond__Price' }
    __validTypes = { 'tbill' : False, 'swap' : False, 'bond' : False, 'rate' : False }
    def __init__(self, args):
        """
        * Instantiate new Bond object using passed arguments.
        """
        self.__tenor = None
        self.__ytm = None
        self.__pricerisks = {}
        self.__nextpmtdate = None 
        self.__prevpmtdate = None
        self.__tenorToCashFlow = None
        self.__accruedinterest = None
        self.__daycount = None
        self.__price = None
        self.__type = None
        if 'type' in args:
            self.__Type(args['type'])
        errs = []
        invalidArgs = []
        reqArgs = Bond.__reqArgs.copy()
        for arg in args:
            lowered = arg.lower()
            if lowered in reqArgs:
                reqArgs[lowered] = True
                targetArg = Bond.__argToSetterMap[lowered]
            elif lowered in Bond.__optionalArgs:
                targetArg = Bond.__optionalArgs[lowered]
            else:
                targetArg = None
                invalidArgs.append(arg)
            if targetArg and targetArg != 'nextpmtdate':
                try:
                    setter = getattr(self, targetArg)
                    setter(args[arg])
                except BaseException as ex:
                    errs.append(str(ex))
            
        missingArgs = [key for key in reqArgs if not reqArgs[key]]
        if missingArgs:
            errs.append(''.join(['The following required arguments are missing: ', ','.join(missingArgs)]))
        if invalidArgs:
            errs.append(''.join(['The following arguments are invalid: ', ','.join(invalidArgs)]))
        # Ensure that dates make sense:
        elif self.__settlement > self.__maturity:
            errs.append('Settlement must be before Maturity.')
        if errs:
            raise Exception('\n'.join(errs))

        if not self.__nextpmtdate:
            self.__CalcNextPaymentDate()

        self.__CalcPrevPaymentDate()
        self.__tenor = Bond.CalculateTenor(self.__settlement, self.__maturity)
        if not self.__tenorToCashFlow:
            self.__MapTenorToCashFlow()
        if not self.__accruedinterest:
            self.__CalculateAccruedInterest()
        if self.__type == 'tbill':
            self.__ConvertQuotation()
        elif not self.__price is None and not self.__ytm:
            self.__YTM()
        
    def __lt__(self, other):
        """
        * For use in sorting Bonds in Bootstrapper.
        """
        return self.Tenor < other.Tenor
        
    ################
    # Properties:
    ################
    @property
    def AccruedInterest(self):
        return self.__accruedinterest
    @property
    def CashFlows(self):
        """
        * Return mapping { PmtDate -> CashFlow }.
        """
        return self.__cashflows.copy()
    @property
    def Coupon(self):
        """
        * Coupon rate in percent form.
        """
        return self.__coupon
    @property
    def CouponPmt(self):
        """
        * Coupon rate as cash payment.
        """
        return self.__coupon * self.FaceValue
    @property
    def DirtyPrice(self):
        return self.__price + self.__accruedinterest
    @property
    def FaceValue(self):
        return self.__facevalue
    @property
    def Maturity(self):
        return self.__maturity
    @property
    def Name(self):
        return self.__name
    @property
    def NextPaymentDate(self):
        return self.__nextpmtdate
    @property
    def PmtFreq(self):
        return self.__pmtfreq
    @property
    def PreviousPaymentDate(self):
        return self.__prevpmtdate
    @property
    def Price(self):
        return self.__price
    @property
    def PriceRisks(self):
        return self.__pricerisks.copy()
    @property
    def Settlement(self):
        return self.__settlement
    @property
    def Tenor(self):
        """
        * Tenor of bond in years.
        """
        return self.__tenor
    @property
    def TenorToCashFlow(self):
        return self.__tenorToCashFlow.copy()
    @property
    def Type(self):
        return self.__type
    @property
    def YTM(self):
        return self.__ytm
    #############
    # Interface methods:
    #############
    def CalculateFairValue(self, zeroCurve):
        """
        * Calculate fair value of bond using existing zero curve.
        """
        if not isinstance(zeroCurve, ZeroCurve):
            raise Exception('zeroCurve must be a ZeroCurve object')
        yields = zeroCurve.Yields
        dfs = zeroCurve.DiscountFactors
        pv = 0
        for tenor in self.__tenorToCashFlow:
            cf = self.__tenorToCashFlow[tenor]
            df = zeroCurve.DiscountFactor(tenor)
            if df is None:
                return 'NA'
            pv += cf * df
        return pv

    #############
    # Private Setters:
    #############
    def __AccruedInterest(self, ai):
        if not isinstance(ai, (int, float, str)):
            raise Exception('AccruedInterest must be numeric or numeric string.')
        elif isinstance(ai, str):
            self.__accruedinterest = float(ai)
        self.__accruedinterest = float(ai)
        if self.__accruedinterest < 0:
            raise Exception('AccruedInterest must be non-negative.')
    def __Coupon(self, coupon):
        if not isinstance(coupon, (int, float, str)):
            raise Exception('Coupon must be numeric or numeric string.')
        self.__coupon = float(coupon)
        if self.__coupon < 0:
            raise Exception('Coupon must be non-negative.')
    def __DayCount(self, dc):
        self.__daycount = DayCount(dc)
    def __FaceValue(self, fv):
        if not isinstance(fv, (float, int, str)):
            raise Exception('fv must be numeric or numeric string.')
        self.__facevalue = float(fv)
        if self.__facevalue < 0:
            raise Exception('fv must be non-negative.')
    def __Maturity(self, mat):
        if not isinstance(mat, (date, datetime, str)):
            raise Exception('maturity must be date/datetime/ MM/DD/YYYY string.')
        elif isinstance(mat, str):
            self.__maturity = datetime.strptime(mat, '%m/%d/%Y')
        else:
            self.__maturity = mat
    def __Name(self, name):
        if not isinstance(name, str):
            raise Exception('name must be a string.')
        self.__name = name.strip()
    def __NextPaymentDate(self, pmtDate):
        if not isinstance(pmtDate, (date, datetime, str)):
            raise Exception('nextpmtdate must be date/datetime/ "MM/DD/YYYY" or blank string.')
        elif isinstance(pmtDate, str) and pmtDate.strip() == '':
            self.__nextpmtdate = None
        elif isinstance(pmtDate, str):
            self.__nextpmtdate = datetime.strptime(pmtDate.strip(), '%m/%d/%Y')
        else:
            self.__nextpmtdate = pmtDate 
    def __PmtFreq(self, freq):
        if isinstance(freq, str):
            self.__pmtfreq= PaymentFreq.FileStringToFreq.value[freq.lower()]
        elif isinstance(freq, PaymentFreq):
            self.__pmtfreq= freq
        elif isinstance(freq, (float, int)):
            self.__pmtfreq= PaymentFreq.IntToFreq.value[int(freq)]
        else:
            raise Exception('PmtFreq must be a PaymetFreq enumeration/valid integer or string.')
    def __Price(self, price):
        if not isinstance(price, (float, int, str)):
            raise Exception('price must be numeric or numeric string.')
        self.__price = float(price)
        if self.__type != 'swap' and self.__price < 0:
            raise Exception('price must be non-negative.')
    def __Type(self, _type):
        if self.__type:
            return
        if not isinstance(_type, str):
            raise Exception('type must be a string.')
        elif not _type in self.__validTypes:
            raise Exception(''.join(['Type must be one of ', ','.join([key for key in self.__validTypes])]))
        self.__type = _type
    def __Settlement(self, settle):
        if not isinstance(settle, (date, datetime, str)):
            raise Exception('Settlement must be date/datetime/ MM/DD/YYYY string.')
        elif isinstance(settle, str):
            self.__settlement = datetime.strptime(settle, '%m/%d/%Y')
        else:
            self.__settlement = settle
    def __CalcPrevPaymentDate(self):
        """
        * Infer the previous payment date based upon frequency and next payment date.
        """
        if not self.__prevpmtdate:
            attr = PaymentFreq.FreqToRelativeDeltaAttr.value[self.PmtFreq]
            adj = PaymentFreq.StepAdjust.value[self.PmtFreq]
            delta = relativedelta()
            setattr(delta, attr, 1 * adj)
            self.__prevpmtdate = self.NextPaymentDate - delta
    def __CalcNextPaymentDate(self):
        """
        * Calculate next payment date.
        """
        if not self.__nextpmtdate:
            attr = PaymentFreq.FreqToRelativeDeltaAttr.value[self.PmtFreq]
            adj = PaymentFreq.StepAdjust.value[self.PmtFreq]
            delta = relativedelta()
            setattr(delta, attr, 1 * adj)
            settle = self.Settlement
            if self.Maturity - delta < settle:
                self.__nextpmtdate = self.Maturity
            else:
                curr = self.Maturity
                while curr >= settle:
                    curr -= delta
                self.__nextpmtdate = curr + delta

    @staticmethod
    def CalculateTenor(startDate, endDate):
        """
        * Calculate bond or payment's tenor.
        """
        start = startDate
        tenor = 0
        while start.year < endDate.year:
            tenor += 1
            start = datetime(year = start.year + 1, month = 1, day = 31)

        startOfYear = datetime(year = start.year, month = 1, day = 1)
        endOfYear = datetime(year = endDate.year, month = 12, day = 31)
        if endDate != start:
            tenor += ((endDate - start).days)/ ((endOfYear - startOfYear).days + 1)

        return tenor

    def __ConvertQuotation(self):
        """
        * Convert discount yield quotation to price, bond equivalent yield.
        """
        pct = self.__price / self.__facevalue
        days = (self.Maturity - self.Settlement).days
        self.__price = self.__facevalue * (1 - pct * days/360)
        self.__ytm = ((self.__facevalue - self.__price) / self.__price) * (365 / days)

    def __MapTenorToCashFlow(self, cfs = None):
        """
        * Generate date to cash flow table if not specified at constructor.
        """
        self.__tenorToCashFlow = {}
        if cfs and not isinstance(cfs, dict):
            raise Exception('TenorToCashFlow must be a dictionary.')
        elif cfs:
            self.__cashflows = cfs
        else:
            # Generate cash flow table:
            self.__cashflows = {}
            freq = int(self.PmtFreq)
            stepAdj = PaymentFreq.StepAdjust.value[self.PmtFreq]
            attr = PaymentFreq.FreqToRelativeDeltaAttr.value[self.PmtFreq]
            step = PaymentFreq.FreqToRRuleFreq.value[self.PmtFreq]
            delta = relativedelta()
            coupon = self.Coupon * self.FaceValue / freq
            start = self.__nextpmtdate if self.__nextpmtdate and self.__nextpmtdate < self.Maturity else self.Settlement
            dates = list(dategen.rrule(freq = step, dtstart = start, interval = stepAdj, until = self.Maturity))
            for count, date in enumerate(dates):
                if count != len(dates) - 1:
                    # Convert the tenor to relevant period units:
                    self.__cashflows[date] = coupon 
                else:
                    # Add the face value payment:
                    self.__cashflows[date] = coupon + self.FaceValue
                
        # Map Tenor (years) to cashflows:
        for count, date in enumerate(self.__cashflows):
            tenor = Bond.CalculateTenor(self.Settlement, date)
            if count != len(dates) - 1:
                # Convert the tenor to relevant period units:
                self.__tenorToCashFlow[tenor] = coupon
            else:
                # Add the face value payment:
                self.__tenorToCashFlow[tenor] = coupon + self.FaceValue
                
    ################
    # Interface Methods
    ################
    def GetAttr(self, attr):
        if attr in Bond.__argToAttrMap:
            return getattr(self, Bond.__argToAttrMap[attr])
        return None
    @staticmethod
    def RequiredArgs():
        """
        * Return all required arguments for constructor.
        """
        return { key : '' for key in Bond.__reqArgs.keys() }
    @staticmethod
    def AllAttributes():
        """
        * Return all possible attributes to attribute name map.
        """
        return Bond.__argToAttrMap.copy()
    @staticmethod
    def PriceSensitivity(ytm, cashflows, freq, order):
        df_raw = (1 + ytm / freq)
        order = int(order)
        direction = (-1 if not order % 2 == 0 else 1) 
        term = 0
        for tenor in cashflows:
            pmt = cashflows[tenor]
            term += direction * pmt * (tenor ** order) * df_raw ** (-tenor * freq - order)
        
        return term

    def __CalculateAccruedInterest(self):
        """
        * Calculate accrued interest for the bond.
        """
        self.__accruedinterest = self.CouponPmt * self.__daycount.DaycountFactor(self.PreviousPaymentDate, self.Settlement)

    def __YTM(self, ytm = None):
        """
        * Calculate the yield to maturity of the security using 
        Newton's Method given all other parameters.
        """
        if ytm and not isinstance(ytm, (int, float, str)):
            raise Exception('YTM must be numeric or numeric string.')
        elif ytm:
            self.__ytm = float(ytm)
            if self.__ytm < 0:
                raise Exception('YTM must be non-negative.')
        elif self.Type == 'rate':
            self.__ytm = self.Price / 100
        elif len(self.__tenorToCashFlow.keys()) == 1:
            tenor = self.__tenor
            finalPmt = self.__tenorToCashFlow[tenor]
            self.__ytm = (finalPmt / self.Price) ** (1/tenor) - 1
        else:
            pass
            # Calculate the yield to maturity:
            tol_consec = .000001
            tol_approx = .000001
            tenorToCashFlow = self.__tenorToCashFlow
            price = lambda x_0, x_1 = tenorToCashFlow, x_2 = self.PmtFreq, x_3 = 0: Bond.PriceSensitivity(x_0, x_1, x_2, x_3)
            duration = lambda x_0, x_1 = tenorToCashFlow, x_2 = self.PmtFreq, x_3 = 1: Bond.PriceSensitivity(x_0, x_1, x_2, x_3)
            self.__ytm = NewtonsMethod(price, duration, self.Coupon, target = self.DirtyPrice)

    def CalculatePriceRisks(self, orders):
        """
        * Calculate Taylor Series terms.
        Inputs:
        * orders: List containing positive integers denoting the 
        order.
        * dateToCashFlow: Map { PaymentDate -> CashFlow }. If not specified then will assume 
        Outputs:
        * Return dictionary mapping order to Taylor Series Expansion term.
        """
        errs = []
        if not isinstance(orders, list):
            errs.append('orders must be list containing positive integers.')
        elif [order for order in orders if not (isinstance(order, (float, int)) and order > 0)]:
            errs.append('orders must be all positive integers.')
        if not self.__price:
            errs.append('Need to calculate price first.')
        if errs:
            raise Exception('\n'.join(errs))

        # Calculate Taylor Series expansion terms:
        freq = int(self.PmtFreq)
        output = { }
        ytm = self.__ytm
        df_raw = (1 + ytm / freq)
        # Calculate price sensitivities using cash flows and yield to maturity:
        for order in orders:
            int_order = int(order)
            output[int_order] = 0
            direction = (-1 if not int_order % 2 == 0 else 1)
            for tenor in self.__tenorToCashFlow:
                pmt = self.__tenorToCashFlow[tenor]
                output[order] +=  direction * pmt * tenor ** int_order * df_raw ** (tenor - int_order)

        # Map order to taylor series expansion term:
        self.__pricerisks = output

    def PriceChangeEst(self, ytm_new):
        """
        * Estimate the change in the price given new yield to maturity.
        """
        errs = []
        if not isinstance(ytm_new, (int, float)):
            errs.append('ytm_new must be numeric.')
        if not self.__price:
            errs.append('Need to calculate Price first.')
        change = 0
        ytm_change = float(ytm_new) - self.YTM 
        terms = self.__pricerisks
        for order in terms:
            change += (1 if order % 2 == 0 else -1) * terms[order] * ytm_change ** order  

        return change / self.Price

class PaymentFreq(Enum):
    ANNUAL = 1
    SEMIANNUAL = 2
    QUARTERLY = 4
    MONTHLY = 12
    WEEKLY = 52
    DAILY = 365
    StepAdjust = { ANNUAL : 1, SEMIANNUAL : 6, QUARTERLY : 4, MONTHLY : 12, WEEKLY : 1, DAILY : 1 }
    FreqToRelativeDeltaAttr = { 1 : 'years', 2 : 'months', 4 : 'months', 12 : 'months', 52 : 'weeks', 365 : 'days'}
    IntToFreq = { 1 : ANNUAL, 2 : SEMIANNUAL, 4 : QUARTERLY, 12 : MONTHLY, 52 : WEEKLY, 365 : DAILY  }
    FileStringToFreq  = { 'annual' : ANNUAL, 'semiannual' : SEMIANNUAL, 'quarterly' : QUARTERLY, 
                             'monthly' : MONTHLY, 'weekly' : WEEKLY, 'daily' : DAILY }
    FreqToRRuleFreq = { ANNUAL : dategen.YEARLY, SEMIANNUAL : dategen.MONTHLY, 
                       QUARTERLY : dategen.MONTHLY, MONTHLY : dategen.MONTHLY, 
                       WEEKLY : dategen.WEEKLY, DAILY : dategen.DAILY }


class DayCount(object):
    """
    * Store and handle functionality of daycounts in calculating fair values of bonds.
    """
    ACT_ACT = 1
    ACT_365 = 2
    ACT_364 = 3
    THIRTY_360 = 4
    ACT_360 = 5
    __validConv = { 1 : ACT_ACT, 2 : ACT_365, 3 : ACT_364, 4 : THIRTY_360, 5 : ACT_360 }
    __validStrings = { "act_act" : 1, "30_360" : 4, "act_365" : 2, "act_364" : 3, "act_360" : 5 }
    __enumToString = { 1 : "act_act", 4 : "30_360", 2 : "act_365", 3 : "act_364", 5 : "act_360" }
    def __init__(self, convention):
        if isinstance(convention, (int, float)) and not int(convention) in DayCount.__validConv:
            raise Exception('convention is must be 2-4.')
        elif isinstance(convention, str) and not convention in DayCount.__validStrings:
            raise Exception(''.join(['convention must be one of ', ','.join([key for key in DayCount.__validStrings.keys()])]))
        elif isinstance(convention, str):
            self.__conv = DayCount.__validStrings[convention]
        elif isinstance(convention, (int, float)):
            self.__conv = int(convention)
    ################
    # Properties:
    ################
    @property
    def Convention(self):
        return self.__conv
    @staticmethod
    def EnumToString():
        return DayCount.__enumToString.copy()
    @staticmethod
    def ValidConventions():
        return DayCount.__validConv.copy()
    @staticmethod
    def ValidStrings():
        return DayCount.__validStrings.copy()
    ################
    # Interface Methods:
    ################
    def DaycountFactor(self, startdate, enddate):
        """
        * Return daycount factor based upon convention.
        Assumes that startdate and enddate are in same year (i.e. not used for calculating tenor).
        """
        if startdate == enddate:
            return 0
        if self.__conv == DayCount.ACT_ACT:
            startOfYear = datetime(year = startdate.year, month = 1, day = 1)
            endOfYear = datetime(year = enddate.year, month = 12, day = 31)
            return ((enddate - startdate).days + 1)/ ((endOfYear - startOfYear).days + 1)
        if self.__conv == DayCount.ACT_365:
            return ((enddate - startdate).days + 1) / 365
        if self.__conv == DayCount.ACT_364:
            return ((enddate - startdate).days + 1) / 364
        if self.__conv == DayCount.THIRTY_360:
            # Start with # of months that passed:
            startMonth = datetime(year = startdate.year, month = startdate.month, day = 1)
            endMonth = datetime(year = enddate.year, month = enddate.month, day = 1)
            days = DayCount.NumberOfMonths(startMonth, endMonth) * 30
            days += ((enddate - endMonth).days + 1)
            days += (startMonth - startdate).days + 1
            return days / 360
        if self.__conv == DayCount.ACT_360:
            return ((enddate - startdate).days + 1) / 360
    def ToInt(self):
        return self.__conv
    def ToString(self):
        return DayCount.__enumToString[self.__conv]

    @staticmethod
    def NumberOfMonths(d1, d2):
        delta = 0
        while True:
            mdays = monthrange(d1.year, d1.month)[1]
            d1 += timedelta(days=mdays)
            if d1 <= d2:
                delta += 1
            else:
                break
        return delta

    @staticmethod
    def IsLeapYear(val):
        if isinstance(val, (datetime, date)):
            year = val.year
        elif isinstance(val, (float, int)) and val >= 0:
            year = int(val)
        else:
            raise Exception('date must be a datetime/date/non-negative numeric')
        if year % 4 != 0:
            return False
        else:
            if year % 100 != 0:
                return True
            elif year % 400 != 0:
                return False
            else:
                return True

def NewtonsMethod(func, firstDeriv, x_0, tol_consec = .01, tol_approx = .01, target = 0):
    """
    * Zerofinding using Newton's Method.
    """
    while firstDeriv(x_0) == 0:
        x_0 += .001
    
    x_new = x_0
    x_old = x_new - .00001
    while not (abs(func(x_new) - target) <= tol_approx and abs(x_new - x_old) <= tol_consec):
        x_old = x_new
        while firstDeriv(x_old) == 0:
            x_old += tol_consec
        x_new -= (func(x_old) - target) / firstDeriv(x_old)

    return x_new 


class ZeroCurve(object):
    """"
    * Risk-free zero coupon curve used for discounting cash flows.
    """
    __reqCols = { "tenor" : -1, "discountfactor" : -1, "yield" : -1 }
    __validInterpTypes = { 'linear' : InterpolationMethods.LINEAR, 'cubic' : InterpolationMethods.CUBIC }
    def __init__(self, interp):
        """
        * Use if wish to load from file.
        """
        self.__dfs = {}
        self.__yields = {}
        self.__interp = None
        
    def __init__(self, dfs, yields, interp, settle):
        """
        * Zero curve used for discounting cash flows.
        Inputs:
        * dfs: Discount factor { Tenor -> DiscountFactor }.
        * yields: { Tenor -> ZeroCouponYield }.
        * interp: Expecting an InterpolationMethods enumeration.
        * settle: Expecting a string/datetime denoting settlement date for discounting.
        """
        # Map { Date -> DiscountFactor }:
        self.__dfs = dfs
        self.__yields = yields
        self.__settlement = settle
        tenors = list(self.__dfs.keys())
        if self.__dfs:
            self.__starttenor = min(tenors)
            self.__endtenor = max(tenors)
        self.__interp = Interpolation(interp, yields)
    
    #############
    # Properties:
    #############
    @property
    def DiscountFactors(self):
        return self.__dfs.copy()
    @property
    def EndTenor(self):
        return self.__endtenor
    @property
    def InterpolationMethod(self):
        return self.__interp
    @property
    def Settlement(self):
        return self.__settlement
    @property
    def StartTenor(self):
        return self.__starttenor
    @property
    def Yields(self):
        return self.__yields.copy()
    @InterpolationMethod.setter
    def InterpolationMethod(self, method):
        """
        * Change the interpolation method using same yields.
        """
        self.__interp = Interpolation(method, self.__yields)
    #############
    # Interface Methods:
    #############
    def DiscountFactor(self, tenor):
        """
        * Return single annualized discount factor 
        corresponding to tenor, using interpolation if 
        necessary.
        """
        if tenor in self.__dfs:
            return self.__dfs[tenor]
        else:
            return (1 + self.__interp.Interpolate(tenor)) ** -tenor

    def CalcDFsForDates(self, path):
        """
        * Calculate discount factors for given dates listed in file, assuming
        settlement is today.
        """
        if not os.path.exists(path):
            raise Exception("File at path does not exist.")
        elif '.csv' not in path:
            raise Exception("File must be a csv.")
        inData = pandas.read_csv(path)
        inData = inData.rename(columns={col : col.lower() for col in inData.columns})
        if 'date' not in inData.columns:
            raise Exception("File requires Date column (case insensitive).")
        try:
            dates = [datetime.strptime(str(date), '%m/%d/%Y') for date in inData['date']]
        except:
            raise Exception('Dates in Date column must follow %m/%d/%Y format.')
        columns = ["Date", "Tenor", "DiscountFactor", "Yield"]
        curr = {col : [] for col in columns}
        settle = self.__settlement
        interpolator = self.__interp
        for date in dates:
            if date > settle:
                tenor = Bond.CalculateTenor(settle,date)
                if tenor not in self.__dfs:
                    _yield = interpolator.Interpolate(tenor)
                    df = self.DiscountFactor(tenor)
                else:
                    _yield = self.__yields[tenor]
                    df = self.__dfs[tenor]
                curr["Date"].append(date)
                curr["Tenor"].append(tenor)
                curr["DiscountFactor"].append(df)
                curr["Yield"].append(_yield)
        outData = DataFrame(curr, columns = columns)
        outData = outData.set_index('Date')
        return outData.sort_values(by="Date")
        
    def LoadFromFile(self, path, interp = None):
        """
        * Load discount factors from file.
        Inputs:
        * path: string path to file.
        * interp: InterpolationMethods enumeration, or None if already loaded or do not
        want to respecify.
        """
        errs = []
        if not isinstance(path, str):
            errs.append('path must be a string.')
        elif '.' not in path:
            errs.append('Path must correspond to file.')
        elif not os.path.exists(path):
            errs.append('File at path does not exist.')
        if not interp and not self.__interp:
            errs.append('interp must be provided.')
        if errs:
            raise Exception(''.join(errs))

        if self.__dfs:
            self.__dfs = {}
            self.__yields = {}
            self.__starttenor = None
            self.__endtenor = None
        reqCols = ZeroCurve.__reqCols.copy()
        numToCol = {}
        keyCols = {'date' : False, 'tenor' : False}
        with open(path, 'r') as f:
            reader = csv.reader(f)
            atHeader = True
            for row in reader:
                for num, col in enumerate(row):
                    if not atHeader:
                        if numToCol[num] == 'tenor':
                            currTenor = float(col.strip())
                        elif numToCol[num] == 'discountfactor':
                            self.__dfs[currTenor] = float(col.strip())
                    elif atHeader and col.lower() in reqCols:
                        reqCols[col.lower()] = num
                        numToCol[num] = col.lower()
                if atHeader:
                    missingCols = [ key for key in reqCols.keys() if reqCols[key] == -1 ]
                    if missingCols:
                        raise Exception(''.join(['The following columns are missing:', ','.join(missingCols)]))
                    atHeader = False
        if self.__dfs:
            self.__yields = ZeroCurve.ConvertZeroYieldsToDFs(self.__dfs)
            self.__starttenor = min(self.__dfs.keys())
            self.__endtenor = max(self.__dfs.keys())
            if interp:
                self.__interp = Interpolation(interp, self.__yields)
            else:
                self.__interp = Interpolation(self.__interp.Method, self.__yields)

    def PrintToFile(self, path = None):
        """
        * Print zero curve to file at path.
        Inputs:
        * path: String to filepath.
        """
        if not path:
            path = ''.join(["ZeroCurve_",self.__interp.String,".csv"])
        sorted = list(self.__dfs.keys())
        sorted.sort()
        reqCols = ZeroCurve.__reqCols.copy()
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([col for col in reqCols.keys()])
            for tenor in sorted:
                row = [tenor, self.__dfs[tenor], self.__yields[tenor]]
                if self.__dfs[tenor] < .5:
                    self.__dfs[tenor] == self.__dfs[tenor]
                writer.writerow(row)

    @classmethod
    def ValidInterpTypes(cls):
        return cls.__validInterpTypes.copy()

    @staticmethod
    def ConvertZeroYieldsToDFs(data):
        """
        * Convert yields to discount factors.
        Inputs:
        * data: Maps {Tenor (years)-> Yield (annualized)}.
        """
        dfs = {}
        for tenor in data.keys():
            dfs[tenor] = (1 + data[tenor]) ** (-tenor)
        return dfs

class Interpolation(object):
    __toString = { InterpolationMethods.LINEAR : "linear", InterpolationMethods.CUBIC : "cubic" }
    def __init__(self, method, data):
        """
        * Interpolation method given dataset.
        """
        self.__method = method
        self.__interpolator = None
        self.__x_data = list(data.keys())
        self.__x_data.sort()
        self.__y_data = []
        for x in self.__x_data:
            self.__y_data.append(data[x])
        self.__x_data = np.asarray(self.__x_data).squeeze()
        self.__y_data = np.asarray(self.__y_data).squeeze()

        if self.__method == InterpolationMethods.CUBIC:
            self.__interpolator = CubicSpline(self.__x_data, self.__y_data)
        elif self.__method == InterpolationMethods.LINEAR:
            self.__interpolator = interp1d(self.__x_data, self.__y_data)
    @property
    def Method(self):
        return self.__method

    @property
    def String(self):
        return Interpolation.__toString[self.__method]

    def Interpolate(self, target):
        """
        * Perform interpolation between start and end inputs.
        Inputs:
        * target: Expecting tenor (floating point).
        """
        # Interpolate zero yield depending on method:
        val = self.__interpolator(target).item()
        return val

if __name__ == '__main__':
    inuts = InputBond

    
    