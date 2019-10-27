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


class PKPDUnit:
    UNIT_TIME_H = 1
    UNIT_TIME_MIN = 2
    UNIT_TIME_SEC = 3
    UNIT_INVTIME_H = 4
    UNIT_INVTIME_MIN = 5
    UNIT_INVTIME_SEC = 6
    UNIT_CONC_g_L= 10
    UNIT_CONC_mg_L= 11
    UNIT_CONC_ug_L= 12
    UNIT_CONC_ng_L= 13
    UNIT_CONC_g_mL= 14
    UNIT_CONC_g_uL= 15
    UNIT_CONC_ug_mL= 16
    UNIT_CONC_mmol_L= 17
    UNIT_CONC_umol_L= 18
    UNIT_CONC_nmol_L= 19
    UNIT_VOLUME_L=30
    UNIT_VOLUME_mL=31
    UNIT_VOLUME_uL=32
    UNIT_VOLUME_nL=33
    UNIT_WEIGHT_kg= 100
    UNIT_WEIGHT_g= 101
    UNIT_WEIGHT_mg= 102
    UNIT_WEIGHT_ug= 103
    UNIT_WEIGHT_ng= 104
    UNIT_WEIGHT_mmol= 105
    UNIT_WEIGHT_umol= 106
    UNIT_WEIGHT_nmol= 107
    UNIT_TIMECONC_H_g_L = 200
    UNIT_TIMECONC_H_mg_L = 201
    UNIT_TIMECONC_H_ug_L = 202
    UNIT_TIMECONC_H_ng_L = 203
    UNIT_TIMECONC_H_g_mL = 204
    UNIT_TIMECONC_H_g_uL = 205
    UNIT_TIMECONC_H_ug_mL = 206
    UNIT_TIMECONC_H_mmol_L = 207
    UNIT_TIMECONC_H_umol_L = 208
    UNIT_TIMECONC_H_nmol_L = 209
    UNIT_TIMECONC_MIN_g_L = 220
    UNIT_TIMECONC_MIN_mg_L = 221
    UNIT_TIMECONC_MIN_ug_L = 222
    UNIT_TIMECONC_MIN_ng_L = 223
    UNIT_TIMECONC_MIN_g_mL = 224
    UNIT_TIMECONC_MIN_g_uL = 225
    UNIT_TIMECONC_MIN_ug_mL = 226
    UNIT_TIMECONC_MIN_mmol_L = 227
    UNIT_TIMECONC_MIN_umol_L = 228
    UNIT_TIMECONC_MIN_nmol_L = 229
    UNIT_TIME2CONC_H2_g_L = 250
    UNIT_TIME2CONC_H2_mg_L = 251
    UNIT_TIME2CONC_H2_ug_L = 252
    UNIT_TIME2CONC_H2_ng_L = 253
    UNIT_TIME2CONC_H2_g_mL = 254
    UNIT_TIME2CONC_H2_g_uL = 255
    UNIT_TIME2CONC_H2_ug_mL = 256
    UNIT_TIME2CONC_H2_mmol_L = 257
    UNIT_TIME2CONC_H2_umol_L = 258
    UNIT_TIME2CONC_H2_nmol_L = 259
    UNIT_TIME2CONC_MIN2_g_L = 270
    UNIT_TIME2CONC_MIN2_mg_L = 271
    UNIT_TIME2CONC_MIN2_ug_L = 272
    UNIT_TIME2CONC_MIN2_ng_L = 273
    UNIT_TIME2CONC_MIN2_g_mL = 274
    UNIT_TIME2CONC_MIN2_g_uL = 275
    UNIT_TIME2CONC_MIN2_ug_mL = 276
    UNIT_TIME2CONC_MIN2_mmol_L = 277
    UNIT_TIME2CONC_MIN2_umol_L = 278
    UNIT_TIME2CONC_MIN2_nmol_L = 279
    UNIT_VOLUMEINVTIME_L_H = 300
    UNIT_VOLUMEINVTIME_mL_H = 301
    UNIT_VOLUMEINVTIME_uL_H = 302
    UNIT_VOLUMEINVTIME_nL_H = 303
    UNIT_VOLUMEINVTIME_L_MIN = 304
    UNIT_VOLUMEINVTIME_mL_MIN = 305
    UNIT_VOLUMEINVTIME_uL_MIN = 306
    UNIT_VOLUMEINVTIME_nL_MIN = 307
    UNIT_VOLUMEINVTIME_L_SEC = 308
    UNIT_VOLUMEINVTIME_mL_SEC = 309
    UNIT_VOLUMEINVTIME_uL_SEC = 310
    UNIT_VOLUMEINVTIME_nL_SEC = 311
    UNIT_VOLUMEINVWEIGHT_L_kg = 312
    UNIT_VOLUMEINVWEIGHT_L_g = 313
    UNIT_WEIGHTINVTIME_kg_H = 350
    UNIT_WEIGHTINVTIME_g_H = 351
    UNIT_WEIGHTINVTIME_mg_H = 352
    UNIT_WEIGHTINVTIME_ug_H = 353
    UNIT_WEIGHTINVTIME_ng_H = 354
    UNIT_WEIGHTINVTIME_kg_MIN = 355
    UNIT_WEIGHTINVTIME_g_MIN = 356
    UNIT_WEIGHTINVTIME_mg_MIN = 357
    UNIT_WEIGHTINVTIME_ug_MIN = 358
    UNIT_WEIGHTINVTIME_ng_MIN = 359
    UNIT_WEIGHTINVTIME_kg_SEC = 360
    UNIT_WEIGHTINVTIME_g_SEC = 361
    UNIT_WEIGHTINVTIME_mg_SEC = 362
    UNIT_WEIGHTINVTIME_ug_SEC = 363
    UNIT_WEIGHTINVTIME_ng_SEC = 364
    UNIT_WEIGHTINVTIME_mmol_H = 365
    UNIT_WEIGHTINVTIME_umol_H = 366
    UNIT_WEIGHTINVTIME_nmol_H = 367
    UNIT_WEIGHTINVTIME_mmol_MIN = 368
    UNIT_WEIGHTINVTIME_umol_MIN = 369
    UNIT_WEIGHTINVTIME_nmol_MIN = 370
    UNIT_WEIGHTINVTIME_mmol_SEC = 371
    UNIT_WEIGHTINVTIME_umol_SEC = 372
    UNIT_WEIGHTINVTIME_nmol_SEC = 373
    UNIT_INVTIME2_MIN2 = 400
    UNIT_PERCENTAGE = 410
    UNIT_NONE = 99999

    unitDictionary = {
        UNIT_TIME_H: "h",
        UNIT_TIME_MIN: "min",
        UNIT_TIME_SEC: "s",
        UNIT_INVTIME_H: "1/h",
        UNIT_INVTIME_MIN: "1/min",
        UNIT_INVTIME_SEC: "1/s",
        UNIT_CONC_g_L: "g/L",
        UNIT_CONC_mg_L: "mg/L",
        UNIT_CONC_ug_L: "ug/L",
        UNIT_CONC_ng_L: "ng/L",
        UNIT_CONC_g_mL: "g/mL",
        UNIT_CONC_g_uL: "g/uL",
        UNIT_CONC_ug_mL: "ug/mL",
        UNIT_CONC_mmol_L: "mmol/L",
        UNIT_CONC_umol_L: "umol/L",
        UNIT_CONC_nmol_L: "nmol/L",
        UNIT_VOLUME_L: "L",
        UNIT_VOLUME_mL: "mL",
        UNIT_VOLUME_uL: "uL",
        UNIT_VOLUME_nL: "nL",
        UNIT_WEIGHT_kg: "kg",
        UNIT_WEIGHT_g: "g",
        UNIT_WEIGHT_mg: "mg",
        UNIT_WEIGHT_ug: "ug",
        UNIT_WEIGHT_ng: "ng",
        UNIT_WEIGHT_mmol: "mmol",
        UNIT_WEIGHT_umol: "umol",
        UNIT_WEIGHT_nmol: "nmol",
        UNIT_NONE: "none",

        UNIT_TIMECONC_H_g_L: "g*h/L",
        UNIT_TIMECONC_H_mg_L: "mg*h/L",
        UNIT_TIMECONC_H_ug_L: "ug*h/L",
        UNIT_TIMECONC_H_ng_L: "ng*h/L",
        UNIT_TIMECONC_H_g_mL: "g*h/mL",
        UNIT_TIMECONC_H_g_uL: "g*h/uL",
        UNIT_TIMECONC_H_ug_mL: "ug*h/mL",
        UNIT_TIMECONC_H_mmol_L: "mmol*h/L",
        UNIT_TIMECONC_H_umol_L:  "umol*h/L",
        UNIT_TIMECONC_H_nmol_L:  "nmol*h/L",
        UNIT_TIMECONC_MIN_g_L: "g*min/L",
        UNIT_TIMECONC_MIN_mg_L: "mg*min/L",
        UNIT_TIMECONC_MIN_ug_L: "ug*min/L",
        UNIT_TIMECONC_MIN_ng_L: "ng*min/L",
        UNIT_TIMECONC_MIN_g_mL: "g*min/mL",
        UNIT_TIMECONC_MIN_g_uL: "g*min/uL",
        UNIT_TIMECONC_MIN_ug_mL: "ug*min/mL",
        UNIT_TIMECONC_MIN_mmol_L: "mmol*min/L",
        UNIT_TIMECONC_MIN_umol_L:  "umol*min/L",
        UNIT_TIMECONC_MIN_nmol_L:  "nmol*min/L",

        UNIT_TIME2CONC_H2_g_L: "g*h^2/L",
        UNIT_TIME2CONC_H2_mg_L: "mg*h^2/L",
        UNIT_TIME2CONC_H2_ug_L: "ug*h^2/L",
        UNIT_TIME2CONC_H2_ng_L: "ng*h^2/L",
        UNIT_TIME2CONC_H2_g_mL: "g*h^2/mL",
        UNIT_TIME2CONC_H2_g_uL: "g*h^2/uL",
        UNIT_TIME2CONC_H2_ug_mL: "ug*h^2/mL",
        UNIT_TIME2CONC_H2_mmol_L: "mmol*h^2/L",
        UNIT_TIME2CONC_H2_umol_L:  "umol*h^2/L",
        UNIT_TIME2CONC_H2_nmol_L:  "nmol*h^2/L",
        UNIT_TIME2CONC_MIN2_g_L: "g*min^2/L",
        UNIT_TIME2CONC_MIN2_mg_L: "mg*min^2/L",
        UNIT_TIME2CONC_MIN2_ug_L: "ug*min^2/L",
        UNIT_TIME2CONC_MIN2_ng_L: "ng*min^2/L",
        UNIT_TIME2CONC_MIN2_g_mL: "g*min^2/mL",
        UNIT_TIME2CONC_MIN2_g_uL: "g*min^2/uL",
        UNIT_TIME2CONC_MIN2_ug_mL: "ug*min^2/mL",
        UNIT_TIME2CONC_MIN2_mmol_L: "mmol*min^2/L",
        UNIT_TIME2CONC_MIN2_umol_L:  "umol*min^2/L",
        UNIT_TIME2CONC_MIN2_nmol_L:  "nmol*min^2/L",

        UNIT_VOLUMEINVTIME_L_H: "L/h",
        UNIT_VOLUMEINVTIME_mL_H: "mL/h",
        UNIT_VOLUMEINVTIME_uL_H: "uL/h",
        UNIT_VOLUMEINVTIME_nL_H: "nL/h",
        UNIT_VOLUMEINVTIME_L_MIN: "L/min",
        UNIT_VOLUMEINVTIME_mL_MIN: "mL/min",
        UNIT_VOLUMEINVTIME_uL_MIN: "uL/min",
        UNIT_VOLUMEINVTIME_nL_MIN: "nL/min",
        UNIT_VOLUMEINVTIME_L_SEC: "L/s",
        UNIT_VOLUMEINVTIME_mL_SEC: "mL/s",
        UNIT_VOLUMEINVTIME_uL_SEC: "uL/s",
        UNIT_VOLUMEINVTIME_nL_SEC: "nL/s",
        UNIT_VOLUMEINVWEIGHT_L_kg: "L/kg",
        UNIT_VOLUMEINVWEIGHT_L_g: "L/g",

        UNIT_WEIGHTINVTIME_kg_H: "kg/h",
        UNIT_WEIGHTINVTIME_g_H: "g/h",
        UNIT_WEIGHTINVTIME_mg_H: "mg/h",
        UNIT_WEIGHTINVTIME_ug_H: "ug/h",
        UNIT_WEIGHTINVTIME_ng_H: "ng/h",
        UNIT_WEIGHTINVTIME_kg_MIN: "kg/min",
        UNIT_WEIGHTINVTIME_g_MIN: "g/min",
        UNIT_WEIGHTINVTIME_mg_MIN: "mg/min",
        UNIT_WEIGHTINVTIME_ug_MIN: "ug/min",
        UNIT_WEIGHTINVTIME_ng_MIN: "ng/min",
        UNIT_WEIGHTINVTIME_kg_SEC: "kg/s",
        UNIT_WEIGHTINVTIME_g_SEC: "g/s",
        UNIT_WEIGHTINVTIME_mg_SEC: "mg/s",
        UNIT_WEIGHTINVTIME_ug_SEC: "ug/s",
        UNIT_WEIGHTINVTIME_ng_SEC: "ng/s",
        UNIT_WEIGHTINVTIME_mmol_H: "mmol/h",
        UNIT_WEIGHTINVTIME_umol_H: "umol/h",
        UNIT_WEIGHTINVTIME_nmol_H: "nmol/h",
        UNIT_WEIGHTINVTIME_mmol_MIN: "mmol/min",
        UNIT_WEIGHTINVTIME_umol_MIN: "umol/min",
        UNIT_WEIGHTINVTIME_nmol_MIN: "nmol/min",
        UNIT_WEIGHTINVTIME_mmol_SEC: "mmol/s",
        UNIT_WEIGHTINVTIME_umol_SEC: "umol/s",
        UNIT_WEIGHTINVTIME_nmol_SEC: "nmol/s",

        UNIT_INVTIME2_MIN2: "1/min^2"
    }

    def __init__(self,unitString=""):
        self.unit = self._fromString(unitString)

    def isTime(self):
        return self.unit>=1 and self.unit<=9

    def isConcentration(self):
        return self.unit>=10 and self.unit<=99

    def isWeight(self):
        return self.unit>=100 and self.unit<=109

    def isWeightInvTime(self):
        return self.unit>=350 and self.unit<=375

    def _fromString(self, unitString):
        return self.stringToCode(unitString)

    def _toString(self):
        return self.codeToString(self.unit)

    @classmethod
    def codeToString(cls, unitCode):
        """ Return the unit string from a given code.
        If code is not recognized, return the empty string.
        """
        return cls.unitDictionary.get(unitCode, "")

    @classmethod
    def stringToCode(cls, unitString):
        """ Return the numeric code from a given string.
        If wrong string value, returns None. """
        if unitString == "":
            return None

        for _unitKey, _unitString in cls.unitDictionary.items():
            if _unitString == unitString:
                return _unitKey
        if unitString == "ug/mL":
            return PKPDUnit.UNIT_CONC_mg_L
        elif unitString == "mg/mL":
            return PKPDUnit.UNIT_CONC_g_L
        elif unitString == "NA" or unitString == "None":
            return PKPDUnit.UNIT_NONE
        return None


def convertUnits(x, unitsIn, unitsOut):
    if unitsIn==unitsOut:
        return x

    if unitsIn == PKPDUnit.UNIT_WEIGHT_g:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e9*x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_mg:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e6*x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_ug:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e3*x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_ng:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-9*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_mmol:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_mmol:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_umol:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_nmol:
            return 1e6*x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_umol:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_mmol:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_umol:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_nmol:
            return 1e3*x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_WEIGHT_nmol:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_mmol:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_umol:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_nmol:
            return x
        else:
            raise Exception("Unknown unit conversion from %s to %s"%(PKPDUnit.unitDictionary[unitsIn],PKPDUnit.unitDictionary[unitsOut]))

    elif unitsIn == PKPDUnit.UNIT_CONC_g_L:
        if unitsOut == PKPDUnit.UNIT_CONC_mg_L or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ug_L:
            return 1e6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ng_L:
            return 1e9*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_mL:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_uL:
            return 1e-6*x

    elif unitsIn == PKPDUnit.UNIT_CONC_mg_L  or unitsIn==PKPDUnit.UNIT_CONC_ug_mL:
        if unitsOut == PKPDUnit.UNIT_CONC_g_L:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_mg_L  or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return 1
        elif unitsOut == PKPDUnit.UNIT_CONC_ug_L:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ng_L:
            return 1e6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_mL:
            return x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_uL:
            return 1e-3*x

    elif unitsIn == PKPDUnit.UNIT_CONC_ug_L:
        if unitsOut == PKPDUnit.UNIT_CONC_g_L:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_mg_L  or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ng_L:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_mL:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_uL:
            return x

    elif unitsIn == PKPDUnit.UNIT_CONC_ng_L:
        if unitsOut == PKPDUnit.UNIT_CONC_g_L:
            return 1e-9*x
        elif unitsOut == PKPDUnit.UNIT_CONC_mg_L or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ug_L:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_mL:
            return x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_uL:
            return 1e3*x

    elif unitsIn == PKPDUnit.UNIT_CONC_g_mL:
        if unitsOut == PKPDUnit.UNIT_CONC_g_L:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_mg_L or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return x
        elif unitsOut == PKPDUnit.UNIT_CONC_ug_L:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ng_L:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_uL:
            return 1e-3*x

    elif unitsIn == PKPDUnit.UNIT_CONC_g_uL:
        if unitsOut == PKPDUnit.UNIT_CONC_g_L:
            return 1e6*x
        elif unitsOut == PKPDUnit.UNIT_CONC_mg_L or unitsOut==PKPDUnit.UNIT_CONC_ug_mL:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_ug_L:
            return x
        elif unitsOut == PKPDUnit.UNIT_CONC_ng_L:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_CONC_g_mL:
            return 1e3*x

    elif unitsIn == PKPDUnit.UNIT_TIME_H:
        if unitsOut == PKPDUnit.UNIT_TIME_MIN:
            return 60*x
        elif unitsOut == PKPDUnit.UNIT_TIME_SEC:
            return 3600*x

    elif unitsIn == PKPDUnit.UNIT_TIME_MIN:
        if unitsOut == PKPDUnit.UNIT_TIME_H:
            return x/60
        elif unitsOut == PKPDUnit.UNIT_TIME_SEC:
            return 60*x

    elif unitsIn == PKPDUnit.UNIT_TIME_SEC:
        if unitsOut == PKPDUnit.UNIT_TIME_H:
            return x/3600
        elif unitsOut == PKPDUnit.UNIT_TIME_MIN:
            return x/60

    elif unitsIn == PKPDUnit.UNIT_TIMECONC_MIN_mg_L:
        if unitsOut == PKPDUnit.UNIT_TIMECONC_H_mg_L:
            return x/60

    elif unitsIn == PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e9*x

    elif unitsIn == PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e6*x

    elif unitsIn == PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return 1e3*x

    elif unitsIn == PKPDUnit.UNIT_WEIGHTINVTIME_ng_MIN:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_g:
            return 1e-9*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_mg:
            return 1e-6*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ug:
            return 1e-3*x
        elif unitsOut == PKPDUnit.UNIT_WEIGHT_ng:
            return x

    elif unitsIn == PKPDUnit.UNIT_WEIGHTINVTIME_umol_MIN:
        if unitsOut == PKPDUnit.UNIT_WEIGHT_umol:
            return x

    return None

def changeRateToMinutes(amount,unit):
    if unit==PKPDUnit.UNIT_WEIGHTINVTIME_g_H:
        return (amount/60,PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN)
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_g_SEC:
        return (amount*60,PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN)

    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_mg_H:
        return (amount/60,PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN)
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_mg_SEC:
        return (amount*60,PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN)

    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ug_H:
        return (amount/60,PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN)
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ug_SEC:
        return (amount*60,PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN)

    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ng_H:
        return (amount/60,PKPDUnit.UNIT_WEIGHTINVTIME_ng_MIN)
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ng_SEC:
        return (amount*60,PKPDUnit.UNIT_WEIGHTINVTIME_ng_MIN)

    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_umol_H:
        return (amount/60,PKPDUnit.UNIT_WEIGHTINVTIME_umol_MIN)
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_umol_SEC:
        return (amount*60,PKPDUnit.UNIT_WEIGHTINVTIME_umol_MIN)

    else:
        return (amount,unit)

def changeRateToWeight(unit):
    if unit==PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN:
        return PKPDUnit.UNIT_WEIGHT_g
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN:
        return PKPDUnit.UNIT_WEIGHT_mg
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN:
        return PKPDUnit.UNIT_WEIGHT_ug
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_ng_MIN:
        return PKPDUnit.UNIT_WEIGHT_ng
    elif unit==PKPDUnit.UNIT_WEIGHTINVTIME_umol_MIN:
        return PKPDUnit.UNIT_WEIGHT_umol
    else:
        return unit

def multiplyUnits(unitX,unitY):
    if unitX==PKPDUnit.UNIT_TIME_H:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_TIMECONC_H_g_L
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_TIMECONC_H_mg_L
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_TIMECONC_H_ug_L
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_TIMECONC_H_ng_L
        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_TIMECONC_H_g_mL
        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_TIMECONC_H_g_uL
        elif unitY==PKPDUnit.UNIT_CONC_ug_mL:
            return PKPDUnit.UNIT_TIMECONC_H_ug_mL
        elif unitY==PKPDUnit.UNIT_CONC_mmol_L:
            return PKPDUnit.UNIT_TIMECONC_H_mmol_L
        elif unitY==PKPDUnit.UNIT_CONC_umol_L:
            return PKPDUnit.UNIT_TIMECONC_H_umol_L
        elif unitY==PKPDUnit.UNIT_CONC_nmol_L:
            return PKPDUnit.UNIT_TIMECONC_H_nmol_L

        elif unitY==PKPDUnit.UNIT_TIMECONC_H_g_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_g_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_mg_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_mg_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_ug_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_ug_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_ng_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_ng_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_g_mL:
            return PKPDUnit.UNIT_TIME2CONC_H2_g_mL
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_g_uL:
            return PKPDUnit.UNIT_TIME2CONC_H2_g_uL
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_ug_mL:
            return PKPDUnit.UNIT_TIME2CONC_H2_ug_mL
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_mmol_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_mmol_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_umol_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_umol_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_H_nmol_L:
            return PKPDUnit.UNIT_TIME2CONC_H2_nmol_L

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_TIME_MIN:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_g_L
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_mg_L
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_ug_L
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_ng_L
        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_TIMECONC_MIN_g_mL
        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_TIMECONC_MIN_g_uL
        elif unitY==PKPDUnit.UNIT_CONC_ug_mL:
            return PKPDUnit.UNIT_TIMECONC_MIN_ug_mL
        elif unitY==PKPDUnit.UNIT_CONC_mmol_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_mmol_L
        elif unitY==PKPDUnit.UNIT_CONC_umol_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_umol_L
        elif unitY==PKPDUnit.UNIT_CONC_nmol_L:
            return PKPDUnit.UNIT_TIMECONC_MIN_nmol_L

        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_g_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_g_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_mg_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_mg_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_ug_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_ug_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_ng_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_ng_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_g_mL:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_g_mL
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_g_uL:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_g_uL
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_ug_mL:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_ug_mL
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_mmol_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_mmol_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_umol_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_umol_L
        elif unitY==PKPDUnit.UNIT_TIMECONC_MIN_nmol_L:
            return PKPDUnit.UNIT_TIME2CONC_MIN2_nmol_L

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUMEINVTIME_mL_MIN:
        if unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN
        elif unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN
        elif unitY==PKPDUnit.UNIT_CONC_ug_mL:
            return PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN

    elif unitX==PKPDUnit.UNIT_VOLUMEINVTIME_L_MIN:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN
        elif unitY==PKPDUnit.UNIT_CONC_umol_L:
            return PKPDUnit.UNIT_WEIGHTINVTIME_umol_MIN

    elif unitX==PKPDUnit.UNIT_INVTIME_MIN:
        if unitY==PKPDUnit.UNIT_INVTIME_MIN:
            return PKPDUnit.UNIT_INVTIME2_MIN2
        else:
            return PKPDUnit.UNIT_NONE


def divideUnits(unitX,unitY):
    if unitX==PKPDUnit.UNIT_WEIGHT_g:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_VOLUME_L
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_NONE
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_NONE
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_NONE

        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_VOLUME_mL

        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_VOLUME_uL

        elif unitY == PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_WEIGHTINVTIME_g_MIN

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_mg:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_VOLUME_mL
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_VOLUME_L
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_NONE
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_NONE

        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_VOLUME_uL

        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_VOLUME_nL

        elif unitY==PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_WEIGHTINVTIME_mg_MIN

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_ug:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_VOLUME_uL
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_VOLUME_mL
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_VOLUME_L
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_NONE

        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_VOLUME_nL
        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_NONE

        elif unitY==PKPDUnit.UNIT_CONC_ug_mL:
            return PKPDUnit.UNIT_VOLUME_mL

        elif unitY == PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_WEIGHTINVTIME_ug_MIN

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_ng:
        if unitY==PKPDUnit.UNIT_CONC_g_L:
            return PKPDUnit.UNIT_VOLUME_nL
        elif unitY==PKPDUnit.UNIT_CONC_mg_L:
            return PKPDUnit.UNIT_VOLUME_uL
        elif unitY==PKPDUnit.UNIT_CONC_ug_L:
            return PKPDUnit.UNIT_VOLUME_mL
        elif unitY==PKPDUnit.UNIT_CONC_ng_L:
            return PKPDUnit.UNIT_VOLUME_L

        elif unitY==PKPDUnit.UNIT_CONC_g_mL:
            return PKPDUnit.UNIT_NONE

        elif unitY==PKPDUnit.UNIT_CONC_g_uL:
            return PKPDUnit.UNIT_NONE

        elif unitY == PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_WEIGHTINVTIME_ng_MIN

        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_mmol:
        if unitY==PKPDUnit.UNIT_CONC_mmol_L:
            return PKPDUnit.UNIT_VOLUME_L
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_umol:
        if unitY==PKPDUnit.UNIT_CONC_umol_L:
            return PKPDUnit.UNIT_VOLUME_L
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_WEIGHT_nmol:
        if unitY==PKPDUnit.UNIT_CONC_nmol_L:
            return PKPDUnit.UNIT_VOLUME_L
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUME_L:
        if unitY==PKPDUnit.UNIT_TIME_H:
            return PKPDUnit.UNIT_VOLUMEINVTIME_L_H
        elif unitY==PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_VOLUMEINVTIME_L_MIN
        elif unitY==PKPDUnit.UNIT_TIME_SEC:
            return PKPDUnit.UNIT_VOLUMEINVTIME_L_SEC
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUME_mL:
        if unitY==PKPDUnit.UNIT_TIME_H:
            return PKPDUnit.UNIT_VOLUMEINVTIME_mL_H
        elif unitY==PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_VOLUMEINVTIME_mL_MIN
        elif unitY==PKPDUnit.UNIT_TIME_SEC:
            return PKPDUnit.UNIT_VOLUMEINVTIME_mL_SEC
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUME_uL:
        if unitY==PKPDUnit.UNIT_TIME_H:
            return PKPDUnit.UNIT_VOLUMEINVTIME_uL_H
        elif unitY==PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_VOLUMEINVTIME_uL_MIN
        elif unitY==PKPDUnit.UNIT_TIME_SEC:
            return PKPDUnit.UNIT_VOLUMEINVTIME_uL_SEC
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUME_nL:
        if unitY==PKPDUnit.UNIT_TIME_H:
            return PKPDUnit.UNIT_VOLUMEINVTIME_nL_H
        elif unitY==PKPDUnit.UNIT_TIME_MIN:
            return PKPDUnit.UNIT_VOLUMEINVTIME_nL_MIN
        elif unitY==PKPDUnit.UNIT_TIME_SEC:
            return PKPDUnit.UNIT_VOLUMEINVTIME_nL_SEC
        else:
            return PKPDUnit.UNIT_NONE

    elif unitX==PKPDUnit.UNIT_VOLUMEINVTIME_L_MIN:
        if unitY==PKPDUnit.UNIT_VOLUME_L:
            return PKPDUnit.UNIT_INVTIME_MIN
        else:
            return PKPDUnit.UNIT_NONE

    else:
        return PKPDUnit.UNIT_NONE


def inverseUnits(unit):
    if unit == PKPDUnit.UNIT_TIME_H:
        return PKPDUnit.UNIT_INVTIME_H
    elif unit == PKPDUnit.UNIT_TIME_MIN:
        return PKPDUnit.UNIT_INVTIME_MIN
    elif unit == PKPDUnit.UNIT_TIME_SEC:
        return PKPDUnit.UNIT_INVTIME_SEC

    else:
        return PKPDUnit.UNIT_NONE

def createUnit(unitNameOrCode):
    unit = PKPDUnit()
    if type(unitNameOrCode)==int:
        unit.unit = unitNameOrCode
    else:
        unit.unit = unit._fromString(unitNameOrCode)
    return unit

# FIXME: Replace the use of this function by PKPDUnit.codeToString
def strUnit(unitCode):
    return PKPDUnit.codeToString(unitCode)

# FIXME: Replace the use of this function by PKPDUnit.stringToCode
def unitFromString(unitString):
    return PKPDUnit.stringToCode(unitString)
