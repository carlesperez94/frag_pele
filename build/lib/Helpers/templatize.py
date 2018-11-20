import sys
import os
from string import Template


class TemplateBuilder(object):
    """
        Description: Base clase to substitute de keywords
        int the already templetize file.

        e.g.

        If we have a file like that --> dededea $RESIDUENAME $CHAIN
        to replace $RESIDUENAME to ASL and $CHAIN to Z we should do:

        file= /path/to/file
        keywords= { "RESIDUENAME" : "ASL", "CHAIN": "Z" }

        TemplateBuilder(file, keywords)

    """

    def __init__(self, file, keywords):

        self.file = file
        self.keywords = keywords
        self.fill_in()

    def fill_in(self):
        """
        Fill the control file in
        """
        with open(os.path.join(self.file), 'r') as infile:
            confile_data = infile.read()

        confile_template = Template(confile_data) 

        confile_text = confile_template.safe_substitute(self.keywords)

        with open(os.path.join(self.file), 'w') as outfile:
            outfile.write(confile_text)