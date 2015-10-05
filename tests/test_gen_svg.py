import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
from mock_helper import gen_file, import_data
import annogesiclib.gen_svg as gs


class TestGenSvg(unittest.TestCase):

    def setUp(self):
        self.test_folder = "test_folder"
        self.example = Example()
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_gen_svg(self):
        gs.gen_svg("test_folder/test.png", 4, 1000, 400)
        data = import_data("test_folder/test.svg")
        self.assertEqual("\n".join(data), self.example.svg)


class Example(object):

    svg = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   version="1.1"
   width="400"
   height="1000"
   viewBox="0 0 1860 1000"
   id="svg3055">
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
     xlink:href="file:///home/silas/ANNOgesic/test_folder/test.png"
     width="100%"
     height="100%"
     preserveAspectRatio="xMidYMin meet"
     id="image3063" />
  <rect
     width="400"
     height="40"
     x="2.0744663"
     y="131"
     id="rect3067"
     style="opacity:0.25;fill:#37c84f;fill-opacity:0.25;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:0.25" />
  <rect
     width="400"
     height="40"
     x="2.0744663"
     y="171"
     id="rect3068"
     style="opacity:0.25;fill:#c8374f;fill-opacity:0.25;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:0.25" />
  <rect
     width="400"
     height="40"
     x="2.0744663"
     y="211"
     id="rect3069"
     style="opacity:0.25;fill:#37c84f;fill-opacity:0.25;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:0.25" />
  <rect
     width="400"
     height="40"
     x="2.0744663"
     y="251"
     id="rect3070"
     style="opacity:0.25;fill:#c8374f;fill-opacity:0.25;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:0.25" />
</svg>"""

if __name__ == "__main__":
    unittest.main()
