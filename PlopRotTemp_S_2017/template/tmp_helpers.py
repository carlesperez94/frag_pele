import re


class Helper(object):

	def __init__(self, *args, **kwargs):
		pass
		
	def preproces_file_lines(self, file):
	  with open(file, "r") as f:
	    lines = f.readlines()
	    for i, line in enumerate(lines):
	      line = re.sub(' +',' ',line)
	      line = line.strip('\n').strip()
	      lines[i] = line
	    return lines