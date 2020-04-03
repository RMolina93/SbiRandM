import argparse, warnings

try:
	from .modeller_comparison import *
except:
	print ("Modeller is not correctly installed.")

from .mainSuperimp import *

args = parse_arguments()

mainMod(args)
