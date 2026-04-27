# example to call this script: python pubs/pptx_to_pdf.py pubs/ezra_seminar_slides/

from pptxtopdf import convert
import sys

args = sys.argv[1]
convert(args, args)
