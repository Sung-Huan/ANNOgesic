import os


def print_track(track_num, svg_out, figure_width):
    id_num = 3067
    x = 2.0744663
    y = 131
    for track in range(track_num):
        if (track % 2) == 0:
            svg_out.write("  <rect\n")
            svg_out.write("     width=\"{0}\"\n".format(figure_width))
            svg_out.write("     height=\"40\"\n")
            svg_out.write("     x=\"{0}\"\n".format(x))
            if track == 0:
                svg_out.write("     y=\"{0}\"\n".format(y))
            else:
                y = y + 40
                svg_out.write("     y=\"{0}\"\n".format(y))
            svg_out.write("     id=\"rect{0}\"\n".format(id_num))
            svg_out.write("     style=\"opacity:0.25;fill:#37c84f;"
                          "fill-opacity:0.25;fill-rule:evenodd;")
            svg_out.write("stroke:#000000;stroke-width:1px;"
                          "stroke-linecap:butt;stroke-linejoin:miter;"
                          "stroke-opacity:0.25\" />\n")
        if (track % 2) == 1:
            svg_out.write("  <rect\n")
            svg_out.write("     width=\"{0}\"\n".format(figure_width))
            svg_out.write("     height=\"40\"\n")
            svg_out.write("     x=\"{0}\"\n".format(x))
            y = y + 40
            svg_out.write("     y=\"{0}\"\n".format(y))
            svg_out.write("     id=\"rect{0}\"\n".format(id_num))
            svg_out.write("     style=\"opacity:0.25;fill:#c8374f;"
                          "fill-opacity:0.25;fill-rule:evenodd;")
            svg_out.write("stroke:#000000;stroke-width:1px;"
                          "stroke-linecap:butt;stroke-linejoin:miter;"
                          "stroke-opacity:0.25\" />\n")
        id_num += 1


def gen_svg(input_png, track_num, figure_height, figure_width):
    svg_out = open(input_png[:-4] + ".svg", "w")
    svg_out.write("""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   version="1.1"
""")
    svg_out.write("   width=\"{0}\"\n".format(figure_width))
    svg_out.write("   height=\"{0}\"\n".format(figure_height))
    svg_out.write("   viewBox=\"0 0 1860 {0}\"\n".format(figure_height))
    svg_out.write("""   id="svg3055">
  <metadata
     id="metadata3061">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <defs
     id="defs3059" />
  <image
""")
    svg_out.write("     xlink:href=\"file://{0}/{1}\"\n".format(
                  os.getcwd(), input_png))
    svg_out.write("""     width="100%"
     height="100%"
     preserveAspectRatio="xMidYMin meet"
     id="image3063" />\n""")
    print_track(track_num, svg_out, figure_width)
    svg_out.write("</svg>")
    svg_out.close()
