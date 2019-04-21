# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (info@kinestat.com)
# *
# * Kinestat Pharma
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'info@kinestat.com'
# *
# **************************************************************************
"""
PKPD functions
"""
import numpy as np
import math
import time
import hashlib
from os.path import (exists, join, splitext, isdir, isfile, islink, expanduser,
                     expandvars, basename, dirname, split, relpath, getmtime)

def parseRange(auxString):
    if auxString=="":
        return None
    elif auxString.startswith('['):
        auxString=auxString.replace('[','')
        auxString=auxString.replace(']','')
        tokens=auxString.split(',')
        auxArray = np.array(tokens, dtype='float')
        return auxArray.astype(np.float)
    elif ':' in auxString:
        tokens=auxString.split(':')
        if len(tokens)!=3:
            raise Exception("The X evaluation string is not well formatted: %s"%auxString)
        fromValue = float(tokens[0])
        step= float(tokens[1])
        toValue = float(tokens[2])
        return np.arange(fromValue,toValue,step)

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def getMD5String(fn):
    fnTime = time.strftime("%Y-%m-%d-%M:%S\n",time.gmtime(getmtime(fn)))
    return hashlib.md5(fnTime+open(fn, 'rb').read()).hexdigest()

def writeMD5(fn):
    if not exists(fn):
        return
    fnMD5=splitext(fn)[0]+".md5"
    fh=open(fnMD5,"w")
    fh.write("%s\n"%getMD5String(fn))
    fh.write("%s\n"%fn)
    fh.close()

def verifyMD5(fn):
    if fn is None or not exists(fn):
        return True
    if fn.endswith(".md5"):
        fnMD5=fn
    else:
        fnMD5=splitext(fn)[0]+".md5"
        if not exists(fnMD5):
            return False
    fh=open(fnMD5,"r")
    md5StringFile=fh.readline().strip()
    fnFile=fh.readline().strip()
    fh.close()
    if fn.endswith(".md5"):
        fn=fnFile
    return getMD5String(fn)==md5StringFile