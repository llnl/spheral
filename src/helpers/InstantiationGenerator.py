#-------------------------------------------------------------------------------
# IntantiationGenerator
# 
# A python script to automatically generate Spheral++ instantion files to be 
# compiled.  Assumed arguments:
#    infile - the file to be read, defining "text"
#   outfile - the file to be written out
#      ndim - an integer value for the dimensionality being generated (1,2,3)
#-------------------------------------------------------------------------------
import re
import sys

assert len(sys.argv) >= 4
infile = sys.argv[1]
outfile = sys.argv[2]
dims = [int(x) for x in sys.argv[3:]]

# Read the input file to get the definition of the string "text",
# which we use to generate the explicit instantiation .cc file
exec(open(infile).read())

# Parse "text" into a header, instantiations, and footer. We need to make sure headers are only included once since many of them do not have include guards.
index = re.search("^[ \t]*template", text, re.MULTILINE).start()
header = text[:index].lstrip()
remainder = text[index:]
index = remainder.rfind("}")
instantiations = remainder[:index].rstrip()
footer = remainder[index:].lstrip()

# Build up the text to write out
outtext = header

for ndim in dims:
    dictionary = {"ndim"      : ndim,
                  "Dim"       : "Dim<%s>" % ndim,
                  "Scalar"    : "Dim<%s>::Scalar" % ndim,
                  "Vector"    : "Dim<%s>::Vector" % ndim,
                  "Tensor"    : "Dim<%s>::Tensor" % ndim,
                  "SymTensor" : "Dim<%s>::SymTensor" % ndim,
    }

    outtext += f"{instantiations % dictionary}\n"

outtext += footer
outtext = outtext.rstrip()

with open(outfile, "w") as f:
    f.write(outtext)
