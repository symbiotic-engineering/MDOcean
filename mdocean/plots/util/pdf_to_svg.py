import pymupdf
import sys

# example to call this script: python mdocean/plots/util/pdf_to_svg.py pubs/ezra_seminar_slides/figs/N2_diagram

def pdf_to_svg(filename_in,filename_out):
    doc = pymupdf.open(filename+".pdf")
    page = doc[0]
    
    # Convert page to SVG
    svg_content = page.get_svg_image()
    
    # Save to file
    with open(filename+".svg", "w", encoding="utf-8") as f:
        f.write(svg_content)
    
    doc.close()

args = sys.argv[1:]
pdf_to_svg(args[0],args[1])