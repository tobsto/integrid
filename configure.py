# ****************************************************************************
# This file is part of Multigrid.
#
# Multigrid is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Multigrid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Multigrid.  If not, see <http:#www.gnu.org/licenses/>.
#
# Copyright 2012 Tobias Stollenwerk
# ****************************************************************************

#!/usr/bin/python

import subprocess
from sys import *
import os

install_dir=''
if (len(argv))>1:
	install_dir=argv[1]
else:
	print "Call with install directory. Break"
	exit(1)

# read in makefile
mfile=open("makefile", 'r')
mfileContent=mfile.readlines()
mfile.close()
# write installation directories into the makefile
mfile=open("makefile", 'w')
mfile.write("INSTALL_DIR=%s\n" % install_dir)

for l in mfileContent:
	if not (l.startswith('INSTALL_DIR')):
		mfile.write("%s" % l)
mfile.close()
