#-------------------------------------------------------------------------------
# IntantiationGenerator
# 
# A python script to automatically generate Spheral++ instantiation files
# to be compiled.  Assumed arguments:
#   input_file  - the file to be read, defining "text"
#   output_file - the file to be written out
#   dimensions  - a list of dimensions for which to generate instantiations (e.g. 1 2 3)
#-------------------------------------------------------------------------------
import argparse
import re

def generate_instantiations(input_file, output_file, dimensions):
    # Read the input file to get the definition of the string "text",
    # which we use to generate the explicit instantiation .cc file.
    variables = {}

    with open(input_file) as f:
        exec(f.read(), variables)

    # Get the intersection of the given dimensions and the supported dimensions,
    # if specified.
    supported_dimensions = variables.get("dimensions")

    if supported_dimensions:
        dimensions = [dim for dim in dimensions if dim in supported_dimensions]

    # Get specializations, if any.
    specializations = variables.get("specializations", {})

    # Parse "text" into a header, instantiations, and footer.
    text = variables["text"]

    index = re.search("^[ \t]*template", text, re.MULTILINE).start()
    header = text[:index]
    remainder = text[index:]

    index = remainder.rfind("}")
    instantiations = remainder[:index].rstrip()
    footer = remainder[index:].lstrip()

    # Build up the text to write out
    output_text = header

    for ndim in dimensions:
        if ndim in specializations:
            output_text += f"{specializations[ndim].rstrip()}\n\n"

        dictionary = {"ndim"      : ndim,
                      "Dim"       : "Dim<%s>" % ndim,
                      "Scalar"    : "Dim<%s>::Scalar" % ndim,
                      "Vector"    : "Dim<%s>::Vector" % ndim,
                      "Tensor"    : "Dim<%s>::Tensor" % ndim,
                      "SymTensor" : "Dim<%s>::SymTensor" % ndim,
        }

        output_text += f"{instantiations % dictionary}\n"

    output_text += footer
    output_text = output_text.strip()

    with open(output_file, "w") as f:
        f.write(output_text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a Spheral++ instantiation file from a template string."
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Input file containing the template as a Python string variable 'text'."
    )

    parser.add_argument(
        "output_file",
        type=str,
        help="Output file to write the generated instantiations."
    )

    parser.add_argument(
        "dimensions",
        type=int,
        nargs="+",
        help="One or more dimensions for which to generate instantiations (e.g., 1 2 3)."
    )

    args = parser.parse_args()
    generate_instantiations(args.input_file, args.output_file, args.dimensions)
