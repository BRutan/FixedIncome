##############################################
# GenerateSpotCurveCalcFairValues.py
##############################################
# Description:
# * Generate spot curve using securities listed
# in file in predetermined format. Price other securities
# using generated spot curves.

from argparse import ArgumentParser
from Bootstrapping import BootStrapper, InputBonds
from FixedIncome import BondGroup, InterpolationMethods, ZeroCurve
import os
import sys
import warnings
warnings.simplefilter("ignore")

def GenerateSpotCurveCalculateFairValues():
    """
    * Generate spot curve using securities listed
    in file in predetermined format. Price other securities
    using generated spot curves.
    """
    desc = ["Generate spot curve using securities listed in predetermined format in file."]
    desc.append("Optionally calculate fair values of bonds in file at path.")
    parser = ArgumentParser(prog="GenerateSpotCurveCalculateFairValues",description = ''.join(desc))
    # Positional args:
    parser.add_argument("bootstrappath", type=str, help="Path to file containing securities to boostrap with.")
    parser.add_argument("interptype", type=str, help="Interpolation type.")
    parser.add_argument("zerocurveoutput", type=str, help="Output file for zero curve.")
    # Optional args:
    parser.add_argument("--fairvaluepath", type=str, help="Path to file containing securities to price using generated zero curve.")
    parser.add_argument("--outputfile", type=str, help="Output file for price comparison.")
    parser.add_argument("--datesfile", type=str, help="Pass file path containing dates to calculate interpolate discount factors.")
    parser.add_argument("--datesoutfile", type=str, help="Output file for discount factors for supplied dates in datesfile. Must include if using datesfile.")
    args = parser.parse_args()

    # Check arguments:
    errs = []
    if args.datesfile is None != args.datesoutfile is None:
        errs.append("Must include both datesfile and datesoutfile if including one.")
    elif args.datesoutfile and '.csv' not in args.datesoutfile:
        errs.append('datesoutfile must be csv.')
    if not args.interptype in ZeroCurve.ValidInterpTypes():
        errs.append('interptype must be one of ')
        errs.append(','.join(ZeroCurve.ValidInterpTypes.keys()))
    if errs:
        raise Exception('\n'.join(errs))
    interptype = args.interptype.lower() 
    _interptype = ZeroCurve.ValidInterpTypes()[interptype]
    # Load all bonds and bootstrap:
    print("Generating zero curve using input bonds.")
    strapper = BootStrapper(args.bootstrappath, _interptype)
    curve = strapper.GenerateSpotCurve()
    # Print curve to file:
    curve.PrintToFile(args.zerocurveoutput)
    # Calculate fair values, using generated curve:
    if args.fairvaluepath and args.outputfile:
        print("Calculating fair values of provided bonds using curve.")
        targetBonds = BondGroup(args.fairvaluepath)
        outFile = ''.join(['PriceComparison', interptype, '.csv'])
        targetBonds.CalculateFairValues(curve)
        targetBonds.PrintPriceComparison(args.outputfile)
    if args.datesfile:
        data = curve.CalcDFsForDates(args.datesfile)
        data.to_csv(args.datesoutfile)
    
if __name__ == '__main__':
    GenerateSpotCurveCalculateFairValues()
