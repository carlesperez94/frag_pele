""" Wrapper script to call apidoc and sphinx quicstart, usage is:
    python makedocs.py -F -o project_source doc_folder
"""
from sphinx import apidoc
import sys

apidoc.main(sys.argv)
