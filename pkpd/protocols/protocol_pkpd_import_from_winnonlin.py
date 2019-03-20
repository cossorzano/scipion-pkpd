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

from .protocol_pkpd_import_from_csv import ProtPKPDImportFromText


class ProtPKPDImportFromWinnonlin(ProtPKPDImportFromText):
    """ Import experiment from Winnonlin.\n
        Protocol created by http://www.kinestatpharma.com\n"""
    _label = 'import from winnonlin'

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):
        ProtPKPDImportFromText._defineParams(self,form,"Winnonlin")

    #--------------------------- STEPS functions --------------------------------------------
    def readTextFile(self):
        if len(self.experiment.samples)!=1:
            raise Exception("Importing from Winnonlin is designed only for 1 sample")

        sample = self.experiment.samples.values()[0] # First (and only) element of the dictionary
        aux = ['dummy']+self.listOfVariables
        sample.addMeasurementPattern(aux)

        fh=open(self.inputFile.get())
        print("Opening %s"%self.inputFile.get())
        inData = False
        for line in fh.readlines():
            line=line.strip()
            print(line)
            tokens = line.split()
            if not inData and len(tokens)>0:
                if tokens[0].strip().lower()=="data":
                    inData = True
            elif inData:
                if len(tokens)!=len(self.listOfVariables):
                    inData = False
                listOfValues = []
                validList = True
                for token in tokens:
                    token=token.strip().lower()
                    if token=="missing":
                        listOfValues.append("NA")
                    else:
                        try:
                            value = float(token)
                            listOfValues.append(token)
                        except:
                            # Not a float
                            validList = False
                            break
                if validList:
                    sample.addMeasurement(' '.join(listOfValues))
        fh.close()