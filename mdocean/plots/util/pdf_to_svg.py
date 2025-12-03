import pymupdf
import sys
import os
import glob

# example to call this script: python mdocean/plots/util/pdf_to_svg.py pubs/ezra_seminar_slides/figs/pdf/ pubs/ezra_seminar_slides/figs/svg/

def pdf_to_svg(folder_in,folder_out):
    files_in = glob.glob(folder_in + '*.pdf')
    for file in files_in:
        doc = pymupdf.open(file)
        page = doc[0]
        
        # Convert page to SVG
        svg_content = page.get_svg_image()
        
        # Save to file
        base = os.path.basename(file)
        os.makedirs(folder_out, exist_ok=True)
        file_out = folder_out + os.path.splitext(base)[0]
        with open(file_out + ".svg", "w", encoding="utf-8") as f:
            f.write(svg_content)
        
        doc.close()

args = sys.argv[1:]
pdf_to_svg(args[0],args[1])