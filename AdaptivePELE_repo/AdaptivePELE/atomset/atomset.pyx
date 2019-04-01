from __future__ import unicode_literals
import numpy as np
import re
from io import StringIO, open
cimport cython
cimport numpy as np
# try:
#     # Check if the basestring type if available, this will fail in python3
#     basestring
# except NameError:
#     basestring = str

# Taken from mdtraj
# https://github.com/mdtraj/mdtraj
_AMINO_ACID_CODES =  {'ACE': None, 'NME':  None, '00C': 'C', '01W':  'X', '02K':
'A', '02L':  'N', '03Y': 'C',  '07O': 'C', '08P':  'C', '0A0': 'D',  '0A1': 'Y',
'0A2': 'K', '0A8':  'C', '0AA': 'V', '0AB': 'V', '0AC':  'G', '0AF': 'W', '0AG':
'L', '0AH':  'S', '0AK': 'D',  '0BN': 'F', '0CS':  'A', '0E5': 'T',  '0EA': 'Y',
'0FL': 'A', '0NC':  'A', '0WZ': 'Y', '0Y8': 'P', '143':  'C', '193': 'X', '1OP':
'Y', '1PA':  'F', '1PI': 'A',  '1TQ': 'W', '1TY':  'Y', '1X6': 'S',  '200': 'F',
'23F': 'F', '23S':  'X', '26B': 'T', '2AD': 'X', '2AG':  'A', '2AO': 'X', '2AS':
'X', '2CO':  'C', '2DO': 'X',  '2FM': 'M', '2HF':  'H', '2KK': 'K',  '2KP': 'K',
'2LU': 'L', '2ML':  'L', '2MR': 'R', '2MT': 'P', '2OR':  'R', '2PI': 'X', '2QZ':
'T', '2R3':  'Y', '2SI': 'X',  '2TL': 'T', '2TY':  'Y', '2VA': 'V',  '2XA': 'C',
'32S': 'X', '32T':  'X', '33X': 'A', '3AH': 'H', '3AR':  'X', '3CF': 'F', '3GA':
'A', '3MD':  'D', '3NF': 'Y',  '3QN': 'K', '3TY':  'X', '3XH': 'G',  '4BF': 'Y',
'4CF': 'F', '4CY':  'M', '4DP': 'W', '4FB': 'P', '4FW':  'W', '4HT': 'W', '4IN':
'W', '4MM':  'X', '4PH': 'F',  '4U7': 'A', '56A':  'H', '5AB': 'A',  '5CS': 'C',
'5CW': 'W', '5HP':  'E', '6CL': 'K', '6CW': 'W', '6GL':  'A', '6HN': 'K', '7JA':
'I', '9NE':  'E', '9NF': 'F',  '9NR': 'R', '9NV':  'V', 'A5N': 'N',  'A66': 'X',
'AA3': 'A', 'AA4':  'A', 'AAR': 'R', 'AB7': 'X', 'ABA':  'A', 'ACB': 'D', 'ACL':
'R', 'ADD':  'X', 'AEA': 'X',  'AEI': 'D', 'AFA':  'N', 'AGM': 'R',  'AGT': 'C',
'AHB': 'N', 'AHH':  'X', 'AHO': 'A', 'AHP': 'A', 'AHS':  'X', 'AHT': 'X', 'AIB':
'A', 'AKL':  'D', 'AKZ': 'D',  'ALA': 'A', 'ALC':  'A', 'ALM': 'A',  'ALN': 'A',
'ALO': 'T', 'ALS':  'A', 'ALT': 'A', 'ALV': 'A', 'ALY':  'K', 'AN8': 'A', 'APE':
'X', 'APH':  'A', 'API': 'K',  'APK': 'K', 'APM':  'X', 'APP': 'X',  'AR2': 'R',
'AR4': 'E', 'AR7':  'R', 'ARG': 'R', 'ARM': 'R', 'ARO':  'R', 'ARV': 'X', 'AS2':
'D', 'AS9':  'X', 'ASA': 'D',  'ASB': 'D', 'ASI':  'D', 'ASK': 'D',  'ASL': 'D',
'ASM': 'X', 'ASN':  'N', 'ASP': 'D', 'ASQ': 'D', 'ASX':  'B', 'AVN': 'X', 'AYA':
'A', 'AZK':  'K', 'AZS': 'S',  'AZY': 'Y', 'B1F':  'F', 'B2A': 'A',  'B2F': 'F',
'B2I': 'I', 'B2V':  'V', 'B3A': 'A', 'B3D': 'D', 'B3E':  'E', 'B3K': 'K', 'B3L':
'X', 'B3M':  'X', 'B3Q': 'X',  'B3S': 'S', 'B3T':  'X', 'B3U': 'H',  'B3X': 'N',
'B3Y': 'Y', 'BB6':  'C', 'BB7': 'C', 'BB8': 'F', 'BB9':  'C', 'BBC': 'C', 'BCS':
'C', 'BE2':  'X', 'BFD': 'D',  'BG1': 'S', 'BH2':  'D', 'BHD': 'D',  'BIF': 'F',
'BIL': 'X', 'BIU':  'I', 'BJH': 'X', 'BL2': 'L', 'BLE':  'L', 'BLY': 'K', 'BMT':
'T', 'BNN':  'F', 'BNO': 'X',  'BOR': 'R', 'BPE':  'C', 'BSE': 'S',  'BTA': 'L',
'BTC': 'C', 'BTR':  'W', 'BUC': 'C', 'BUG': 'V', 'C1X':  'K', 'C22': 'A', 'C3Y':
'C', 'C4R':  'C', 'C5C': 'C',  'C66': 'X', 'C6C':  'C', 'CAF': 'C',  'CAL': 'X',
'CAS': 'C', 'CAV':  'X', 'CAY': 'C', 'CCL': 'K', 'CCS':  'C', 'CDE': 'X', 'CDV':
'X', 'CEA':  'C', 'CGA': 'E',  'CGU': 'E', 'CHF':  'X', 'CHG': 'X',  'CHP': 'G',
'CHS': 'X', 'CIR':  'R', 'CLE': 'L', 'CLG': 'K', 'CLH':  'K', 'CME': 'C', 'CMH':
'C', 'CML':  'C', 'CMT': 'C',  'CPC': 'X', 'CPI':  'X', 'CR5': 'G',  'CS0': 'C',
'CS1': 'C', 'CS3':  'C', 'CS4': 'C', 'CSA': 'C', 'CSB':  'C', 'CSD': 'C', 'CSE':
'C', 'CSJ':  'C', 'CSO': 'C',  'CSP': 'C', 'CSR':  'C', 'CSS': 'C',  'CSU': 'C',
'CSW': 'C', 'CSX':  'C', 'CSZ': 'C', 'CTE': 'W', 'CTH':  'T', 'CUC': 'X', 'CWR':
'S', 'CXM':  'M', 'CY0': 'C',  'CY1': 'C', 'CY3':  'C', 'CY4': 'C',  'CYA': 'C',
'CYD': 'C', 'CYF':  'C', 'CYG': 'C', 'CYJ': 'K', 'CYM':  'C', 'CYQ': 'C', 'CYR':
'C', 'CYS':  'C', 'CZ2': 'C',  'CZZ': 'C', 'D11':  'T', 'D3P': 'G',  'D4P': 'X',
'DA2': 'X', 'DAB':  'A', 'DAH': 'F', 'DAL': 'A', 'DAR':  'R', 'DAS': 'D', 'DBB':
'T', 'DBS':  'S', 'DBU': 'T',  'DBY': 'Y', 'DBZ':  'A', 'DC2': 'C',  'DCL': 'X',
'DCY': 'C', 'DDE':  'H', 'DFI': 'X', 'DFO': 'X', 'DGH':  'G', 'DGL': 'E', 'DGN':
'Q', 'DHA':  'S', 'DHI': 'H',  'DHL': 'X', 'DHN':  'V', 'DHP': 'X',  'DHV': 'V',
'DI7': 'Y', 'DIL':  'I', 'DIR': 'R', 'DIV': 'V', 'DLE':  'L', 'DLS': 'K', 'DLY':
'K', 'DM0':  'K', 'DMH': 'N',  'DMK': 'D', 'DMT':  'X', 'DNE': 'L',  'DNL': 'K',
'DNP': 'A', 'DNS':  'K', 'DOA': 'X', 'DOH': 'D', 'DON':  'L', 'DPL': 'P', 'DPN':
'F', 'DPP':  'A', 'DPQ': 'Y',  'DPR': 'P', 'DSE':  'S', 'DSG': 'N',  'DSN': 'S',
'DSP': 'D', 'DTH':  'T', 'DTR': 'W', 'DTY': 'Y', 'DVA':  'V', 'DYS': 'C', 'ECC':
'Q', 'EFC':  'C', 'EHP': 'F',  'ESB': 'Y', 'ESC':  'M', 'EXY': 'L',  'EYS': 'X',
'F2F': 'F', 'FAK':  'K', 'FB5': 'A', 'FB6': 'A', 'FCL':  'F', 'FGA': 'E', 'FGL':
'G', 'FGP':  'S', 'FH7': 'K',  'FHL': 'K', 'FHO':  'K', 'FLA': 'A',  'FLE': 'L',
'FLT': 'Y', 'FME':  'M', 'FOE': 'C', 'FP9': 'P', 'FRD':  'X', 'FT6': 'W', 'FTR':
'W', 'FTY':  'Y', 'FVA': 'V',  'FZN': 'K', 'GAU':  'E', 'GCM': 'X',  'GFT': 'S',
'GGL': 'E', 'GHG':  'Q', 'GHP': 'G', 'GL3': 'G', 'GLH':  'Q', 'GLJ': 'E', 'GLK':
'E', 'GLM':  'X', 'GLN': 'Q',  'GLQ': 'E', 'GLU':  'E', 'GLX': 'Z',  'GLY': 'G',
'GLZ': 'G', 'GMA':  'E', 'GND': 'X', 'GPL': 'K', 'GSC':  'G', 'GSU': 'E', 'GT9':
'C', 'GVL':  'S', 'H14': 'F',  'H5M': 'P', 'HAC':  'A', 'HAR': 'R',  'HBN': 'H',
'HCS': 'X', 'HFA':  'X', 'HGL': 'X', 'HHI': 'H', 'HIA':  'H', 'HIC': 'H', 'HIP':
'H', 'HIQ':  'H', 'HIE': 'H', 'HIS': 'H',  'HL2': 'L', 'HLU':  'L',  'HMR': 'R',
'HPC': 'F', 'HPE': 'F', 'HPH':  'F', 'HPQ': 'F', 'HQA': 'A', 'HRG':  'R', 'HRP':
    'W', 'HS8': 'H', 'CYT': 'C',
'HS9':  'H', 'HSE': 'S',  'HSL': 'S', 'HSO':  'H', 'HTI': 'C',  'HTN': 'N',
'HTR': 'W', 'HV5':  'A', 'HVA': 'V', 'HY3': 'P', 'HYP':  'P', 'HZP': 'P', 'I2M':
'I', 'I58':  'K', 'IAM': 'A',  'IAR': 'R', 'IAS':  'D', 'IEL': 'K',  'IGL': 'G',
'IIL': 'I', 'ILE':  'I', 'ILG': 'E', 'ILX': 'I', 'IML':  'I', 'IOY': 'F', 'IPG':
'G', 'IT1':  'K', 'IYR': 'Y',  'IYT': 'T', 'IZO':  'M', 'JJJ': 'C',  'JJK': 'C',
'JJL': 'C', 'K1R':  'C', 'KCX': 'K', 'KGC': 'K', 'KNB':  'A', 'KOR': 'M', 'KPI':
'K', 'KST':  'K', 'KYN': 'W',  'KYQ': 'K', 'L2A':  'X', 'LA2': 'K',  'LAA': 'D',
'LAL': 'A', 'LBY':  'K', 'LCK': 'K', 'LCX': 'K', 'LCZ':  'X', 'LDH': 'K', 'LED':
'L', 'LEF':  'L', 'LEH': 'L',  'LEI': 'V', 'LEM':  'L', 'LEN': 'L',  'LET': 'K',
'LEU': 'L', 'LEX':  'L', 'LHC': 'X', 'LLP': 'K', 'LLY':  'K', 'LME': 'E', 'LMF':
'K', 'LMQ':  'Q', 'LP6': 'K',  'LPD': 'P', 'LPG':  'G', 'LPL': 'X',  'LPS': 'S',
'LSO': 'K', 'LTA':  'X', 'LTR': 'W', 'LVG': 'G', 'LVN':  'V', 'LYF': 'K', 'LYK':
'K', 'LYM':  'K', 'LYN': 'K',  'LYR': 'K', 'LYS':  'K', 'LYX': 'K',  'LYZ': 'K',
'M0H': 'C',  'M2L': 'K', 'M2S': 'M',  'M30': 'G', 'M3L': 'K',  'MA': 'A', 'MAA':
'A', 'MAI':  'R', 'MBQ': 'Y',  'MC1': 'S', 'MCG':  'X', 'MCL': 'K',  'MCS': 'C',
'MD3': 'C', 'MD6':  'G', 'MDF': 'Y', 'MDH': 'X', 'MEA':  'F', 'MED': 'M', 'MEG':
'E', 'MEN':  'N', 'MEQ': 'Q',  'MET': 'M', 'MEU':  'G', 'MF3': 'X',  'MGG': 'R',
'MGN': 'Q', 'MGY':  'G', 'MHL': 'L', 'MHO': 'M', 'MHS':  'H', 'MIS': 'S', 'MK8':
'L', 'ML3':  'K', 'MLE': 'L',  'MLL': 'L', 'MLY':  'K', 'MLZ': 'K',  'MME': 'M',
'MMO': 'R', 'MND':  'N', 'MNL': 'L', 'MNV': 'V', 'MOD':  'X', 'MP8': 'P', 'MPH':
'X', 'MPJ':  'X', 'MPQ': 'G',  'MSA': 'G', 'MSE':  'M', 'MSL': 'M',  'MSO': 'M',
'MSP': 'X', 'MT2':  'M', 'MTY': 'Y', 'MVA': 'V', 'N10':  'S', 'N2C': 'X', 'N7P':
'P', 'N80':  'P', 'N8P': 'P',  'NA8': 'A', 'NAL':  'A', 'NAM': 'A',  'NB8': 'N',
'NBQ': 'Y', 'NC1':  'S', 'NCB': 'A', 'NCY': 'X', 'NDF':  'F', 'NEM': 'H', 'NEP':
'H', 'NFA':  'F', 'NHL': 'E',  'NIY': 'Y', 'NLE':  'L', 'NLN': 'L',  'NLO': 'L',
'NLP': 'L', 'NLQ':  'Q', 'NMC': 'G', 'NMM': 'R', 'NNH':  'R', 'NPH': 'C', 'NPI':
'A', 'NSK':  'X', 'NTR': 'Y',  'NTY': 'Y', 'NVA':  'V', 'NYS': 'C',  'NZH': 'H',
'O12': 'X', 'OAR':  'R', 'OAS': 'S', 'OBF': 'X', 'OBS':  'K', 'OCS': 'C', 'OCY':
'C', 'OHI':  'H', 'OHS': 'D',  'OIC': 'X', 'OLE':  'X', 'OLT': 'T',  'OLZ': 'S',
'OMT': 'M', 'ONH':  'A', 'ONL': 'X', 'OPR': 'R', 'ORN':  'A', 'ORQ': 'R', 'OSE':
'S', 'OTB':  'X', 'OTH': 'T',  'OXX': 'D', 'P1L':  'C', 'P2Y': 'P',  'PAQ': 'Y',
'PAS': 'D', 'PAT':  'W', 'PAU': 'A', 'PBB': 'C', 'PBF':  'F', 'PCA': 'E', 'PCC':
'P', 'PCE':  'X', 'PCS': 'F',  'PDL': 'X', 'PEC':  'C', 'PF5': 'F',  'PFF': 'F',
'PFX': 'X', 'PG1':  'S', 'PG9': 'G', 'PGL': 'X', 'PGY':  'G', 'PH6': 'P', 'PHA':
'F', 'PHD':  'D', 'PHE': 'F',  'PHI': 'F', 'PHL':  'F', 'PHM': 'F',  'PIV': 'X',
'PLE': 'L', 'PM3':  'F', 'POM': 'P', 'PPN': 'F', 'PR3':  'C', 'PR9': 'P', 'PRO':
'P', 'PRS':  'P', 'PSA': 'F',  'PSH': 'H', 'PTA':  'X', 'PTH': 'Y',  'PTM': 'Y',
'PTR': 'Y', 'PVH':  'H', 'PVL': 'X', 'PYA': 'A', 'PYL':  'O', 'PYX': 'C', 'QCS':
'C', 'QMM':  'Q', 'QPA': 'C',  'QPH': 'F', 'R1A':  'C', 'R4K': 'W',  'RE0': 'W',
'RE3': 'W', 'RON':  'X', 'RVX': 'S', 'RZ4': 'S', 'S1H':  'S', 'S2C': 'C', 'S2D':
'A', 'S2P':  'A', 'SAC': 'S',  'SAH': 'C', 'SAR':  'G', 'SBL': 'S',  'SCH': 'C',
'SCS': 'C', 'SCY':  'C', 'SD2': 'X', 'SDP': 'S', 'SE7':  'A', 'SEB': 'S', 'SEC':
'U', 'SEG':  'A', 'SEL': 'S',  'SEM': 'S', 'SEN':  'S', 'SEP': 'S',  'SER': 'S',
'SET': 'S', 'SGB':  'S', 'SHC': 'C', 'SHP': 'G', 'SHR':  'K', 'SIB': 'C', 'SLR':
'P', 'SLZ':  'K', 'SMC': 'C',  'SME': 'M', 'SMF':  'F', 'SNC': 'C',  'SNN': 'N',
'SOC': 'C', 'SOY':  'S', 'SRZ': 'S', 'STY': 'Y', 'SUB':  'X', 'SUN': 'S', 'SVA':
'S', 'SVV':  'S', 'SVW': 'S',  'SVX': 'S', 'SVY':  'S', 'SVZ': 'S',  'SYS': 'C',
'T11': 'F', 'T66':  'X', 'TA4': 'X', 'TAV': 'D', 'TBG':  'V', 'TBM': 'T', 'TCQ':
'Y', 'TCR':  'W', 'TDD': 'L',  'TFQ': 'F', 'TH6':  'T', 'THC': 'T',  'THO': 'X',
'THR': 'T', 'THZ':  'R', 'TIH': 'A', 'TMB': 'T', 'TMD':  'T', 'TNB': 'C', 'TNR':
'S', 'TOQ':  'W', 'TPH': 'X',  'TPL': 'W', 'TPO':  'T', 'TPQ': 'Y',  'TQI': 'W',
'TQQ': 'W', 'TRF':  'W', 'TRG': 'K', 'TRN': 'W', 'TRO':  'W', 'TRP': 'W', 'TRQ':
'W', 'TRW':  'W', 'TRX': 'W',  'TRY': 'W', 'TST':  'X', 'TTQ': 'W',  'TTS': 'Y',
'TXY': 'Y', 'TY1':  'Y', 'TY2': 'Y', 'TY3': 'Y', 'TY5':  'Y', 'TYB': 'Y', 'TYI':
'Y', 'TYJ':  'Y', 'TYN': 'Y',  'TYO': 'Y', 'TYQ':  'Y', 'TYR': 'Y',  'TYS': 'Y',
'TYT': 'Y', 'TYW':  'Y', 'TYX': 'X', 'TYY': 'Y', 'TZB':  'X', 'TZO': 'X', 'UMA':
'A', 'UN1':  'X', 'UN2': 'X',  'UNK': 'X', 'VAD':  'V', 'VAF': 'V',  'VAL': 'V',
'VB1': 'K', 'VDL':  'X', 'VLL': 'X', 'VLM': 'X', 'VMS':  'X', 'VOL': 'X', 'WLU':
'L', 'WPA':  'F', 'WRP': 'W',  'WVL': 'V', 'X2W':  'E', 'XCN': 'C',  'XCP': 'X',
'XDT': 'T', 'XPL':  'O', 'XPR': 'P', 'XSN': 'N', 'XX1':  'K', 'YCM': 'C', 'YOF':
'Y', 'YTH':  'T', 'Z01': 'A',  'ZAL': 'A', 'ZCL':  'F', 'ZFB': 'X',  'ZU0': 'T',
'ZZJ': 'A'}

REGEX_PATTERN = re.compile("[0-9]|\+|\-")

cdef class Atom:
    _chargePattern = REGEX_PATTERN
    _ATOM_WEIGHTS = {u"H": 1.00794,
                    u"D": 2.01410178,  # deuterium
                    u"HE": 4.00,
                    u"LI": 6.941,
                    u"BE": 9.01,
                    u"B": 10.811,
                    u"C": 12.0107,
                    u"N": 14.0067,
                    u"O": 15.9994,
                    u"F": 18.998403,
                    u"NE": 20.18,
                    u"NA": 22.989769,
                    u"MG": 24.305,
                    u"AL": 26.98,
                    u"SI": 28.09,
                    u"P": 30.973762,
                    u"S": 32.065,
                    u"CL": 35.453,
                    u"AR": 39.95,
                    u"K": 39.0983,
                    u"CA": 40.078,
                    u"SC": 44.96,
                    u"TI": 47.87,
                    u"V": 50.94,
                    u"CR": 51.9961,
                    u"MN": 54.938045,
                    u"FE": 55.845,
                    u"CO": 58.93,
                    u"NI": 58.6934,
                    u"CU": 63.546,
                    u"ZN": 65.409,
                    u"GA": 69.72,
                    u"GE": 72.64,
                    u"AS": 74.9216,
                    u"SE": 78.96,
                    u"BR": 79.90,
                    u"KR": 83.80,
                    u"RB": 85.47,
                    u"SR": 87.62,
                    u"Y": 88.91,
                    u"ZR": 91.22,
                    u"NB": 92.91,
                    u"W": 95.94,  # Molybdenum?  Not sure why it's not always MO
                    u"MO": 95.94,
                    u"TC": 98.0,
                    u"RU": 101.07,
                    u"RH": 102.91,
                    u"PD": 106.42,
                    u"AG": 107.8682,
                    u"CD": 112.411,
                    u"IN": 114.82,
                    u"SN": 118.71,
                    u"SB": 121.76,
                    u"TE": 127.60,
                    u"I": 126.90447,
                    u"XE": 131.29,
                    u"CS": 132.91,
                    u"BA": 137.33,
                    u"PR": 140.91,
                    u"EU": 151.96,
                    u"GD": 157.25,
                    u"TB": 158.93,
                    u"IR": 192.22,
                    u"PT": 195.084,
                    u"AU": 196.96657,
                    u"HG": 200.59,
                    u"PB": 207.2,
                    u"U": 238.03}
    def __init__(self, basestring atomContent=u""):
        """ Create an atom from a pdb line

            :param atomContent: Line of the pdb from which the atom will be created
            :type atomContent: basestring
        """
        # Force string attributes to be unicode strings
        self.atomSerial = u""
        self.name = u""
        self.resname = u""
        self.resChain = u""
        self.resnum = u""
        self.type = u""
        # atomContent = atomContent.split()
        if len(atomContent) > 6 and (atomContent[:4] == u'ATOM' or atomContent[:6] == u'HETATM'):
            self.atomSerial = atomContent[6:11].strip()
            self.name = atomContent[12:16].strip()
            self.resname = atomContent[17:20].strip()
            self.resChain = atomContent[21]
            self.resnum = atomContent[22:26].strip()
            self.x = float(atomContent[30:38])
            self.y = float(atomContent[38:46])
            self.z = float(atomContent[46:54])

            self.type = re.sub(self._chargePattern, u"", atomContent[76:80]).strip().upper()
            if self.type in self._ATOM_WEIGHTS:
                self.mass = self._ATOM_WEIGHTS[self.type]
            else:
                print("WARNING!: Information about atom type not available, trying to guess from its name, mass properties might be wrong.")
                self.mass = self._ATOM_WEIGHTS[self.name[0]]

            if atomContent.startswith(u'ATOM'):
                self.protein = True
            else:
                self.protein = False

            self.id = self.atomSerial + u":" + self.name + u":" + self.resname
            # self.id = self.atomSerial

    def __getstate__(self):
        # Copy the object's state from
        state = {u"atomSerial": self.atomSerial, u"name": self.name, u"x": self.x,
                 u"y": self.y, u"z": self.z, u"mass": self.mass, u"type": self.type,
                 u"resname": self.resname, u"resChain": self.resChain,
                 u"resnum": self.resnum, u"protein": self.protein, u"id": self.id}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.atomSerial = state[u'atomSerial']
        self.name = state[u'name']
        self.resname = state[u'resname']
        self.resnum = state[u'resnum']
        self.resChain = state[u'resChain']
        self.type = state[u'type']
        self.id = state[u'id']
        self.mass = state[u'mass']
        self.x = state[u'x']
        self.y = state[u'y']
        self.z = state[u'z']
        self.protein = state[u'protein']

    def equivalentResname(self, resname):
        return _AMINO_ACID_CODES[self.resname] == _AMINO_ACID_CODES[resname]

    def set_properties(self, bint isProtein, basestring atomSerial, basestring atomName, basestring resName, basestring resNum, float x, float y, float z, basestring element, basestring resChain):
        self.atomSerial = atomSerial
        self.name = atomName
        self.resname = resName
        self.resChain = resChain
        self.resnum = resNum
        self.x = x
        self.y = y
        self.z = z

        self.type = element
        if self.type in self._ATOM_WEIGHTS:
            self.mass = self._ATOM_WEIGHTS[self.type]
        else:
            print("WARNING!: Information about atom type not available, trying to guess from its name, mass properties might be wrong.")
            self.mass = self._ATOM_WEIGHTS[self.name[0]]

        self.protein = isProtein

        self.id = self.atomSerial + u":" + self.name + u":" + self.resname

    def isHeavyAtom(self):
       """
            Check if Atom is a heavy atom

            :returns: bool -- True if Atom is heavy atom, false otherwise
       """
       return self.type != 'H'

    def isProtein(self):
        """
            Check if Atom is a protein atom

            :returns: bool -- True if Atom is a protein atom, false otherwise
        """
        return self.protein

    def isHeteroAtom(self):
        """
            Check if Atom is an hetero atom

            :returns: bool -- True if Atom is an hetero atom, false otherwise
        """
        return not self.protein

    def printAtom(self):
        """
            Print Atom information
        """
        print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.x, self.y, self.z, self.type, self.mass

    def __richcmp__(self, Atom atom2, int op):
        if op == 2:
            #equality
            if self.name != atom2.name:
                return False
            if self.atomSerial != atom2.atomSerial:
                return False
            if self.protein and atom2.protein:
                if not self.equivalentResname(atom2.resname):
                    return False
            else:
                if self.resname != atom2.resname:
                    return False
            return True
        elif op == 3:
            if self.name == atom2.name:
                return False
            if self.atomSerial == atom2.atomSerial:
                return False
            if self.protein and atom2.protein:
                if self.equivalentResname(atom2.resname):
                    return False
            else:
                if self.resname == atom2.resname:
                    return False
            return True
        elif op == 1:
            if self.id == atom2.id:
                return True
            else:
                return self.serial < atom2.serial
        elif op == 5:
            if self.id == atom2.id:
                return True
            else:
                return self.serial > atom2.serial
        elif op == 0:
            return self.serial < atom2.serial
        elif op == 4:
            return self.serial > atom2.serial

    def __str__(self):
        return u"%s: %s %s %s [%f, %f, %f] %s %f" % (self.id, self.atomSerial,
                                                    self.resChain, self.resnum,
                                                    self.x, self.y, self.z,
                                                    self.type, self.mass)

    def getAtomCoords(self):
        """
            Get the coordinates of the atom

            :returns: numpy.Array -- Array with the coordinate of the atom
        """
        return np.array([self.x, self.y, self.z])

    def squaredDistance(self, Atom atom2):
        """
            Calculate the squared distance between two atoms

            :param atom2: Second Atom to whom the distance will be calculated
            :type atom2: Atom
            :returns: float -- The distance between the atoms
        """
        return (self.x - atom2.x)**2 + (self.y - atom2.y)**2 + (self.z - atom2.z)**2


cdef class PDB:
    _typeProtein = u"PROTEIN"
    _typeHetero = u"HETERO"
    _typeAll = u"ALL"
    _typeCM = u"CM"

    #Atoms to be used in the contact map
    CMAtoms = {u"ALA": u"empty", u"VAL": u"empty", u"LEU": u"empty", u"ILE": u"empty",
               u"MET": u"empty", u"PRO": u"empty", u"PHE": u"CZ", u"TYR": u"OH",
               u"TRP": u"CH2", u"SER": u"empty", u"THR": u"empty", u"CYS": u"empty",
               u"ASN": u"empty", u"GLN": u"empty", u"LYS": u"NZ", u"HIS": u"CE1",
               u"HIE": u"CE1", u"HID": u"CE1", u"HIP": u"CE1", u"ARG": u"NE",
               u"ASP": u"OD1", u"GLU": u"OE1", u"GLY": u"empty"}
    ATOM_LINE_TEMPLATE = u"%s%s %s %s %s%s%s   %.3f%.3f%.3f%.2f%.2f          %s   "

    def __init__(self):
        """
            Object that will contain the information of a PDB file. Has to call
            the initialise method to load the file
        """
        self.atoms = {}
        # {atomId: atom, ...}
        # Where atomId := serial:atomName:resName
        self.totalMass = 0
        # ensure every string is unicode
        self.pdb = None
        self.com = None
        self.centroid = None
        self.ispdb = False

        # Necessary for contactMaps
        self.atomList = []

    def __richcmp__(self, object other, int op):
        """
            Compare two pdb strings, remark lines should be ignored and only the
            atoms and its information should be compared
        """
        cdef list pdb1, pdb2
        if op == 2:
            if self.ispdb and other.ispdb:
                pdb1 = [element.strip() for element in self.pdb.split(u'\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                pdb2 = [element.strip() for element in other.pdb.split(u'\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                return pdb1 == pdb2
            else:
                for atom1, atom2 in zip(self.atomList, other.atomList):
                    if self.atoms[atom1] != other.atoms[atom2]:
                        return False
                return True
        elif op == 3:
            if self.ispdb and other.ispdb:
                pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith(u"ATOM") or element.startswith(u"HETATM")]
                return pdb1 != pdb2
            else:
                for atom1, atom2 in zip(self.atomList, other.atomList):
                    if self.atoms[atom1] == other.atoms[atom2]:
                        return False
                return True
        else:
            print "No boolean operator available for PDB apart from equality"

    def __getstate__(self):
        # Copy the object's state from
        state = {u"atoms": self.atoms, u"atomList": self.atomList,
                 u"com": self.com, u"centroid": self.centroid,
                 u"totalMass": self.totalMass, u"pdb": self.pdb,
                 u"ispdb": self.ispdb}
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.atoms = state[u'atoms']
        self.atomList = state[u'atomList']
        self.com = state.get(u'com')
        self.centroid = state.get(u'centroid')
        self.totalMass = state[u'totalMass']
        self.pdb = state[u'pdb']
        self.ispdb = state.get(u'ispdb', True)

    def isfromPDBFile(self):
        return self.ispdb

    def _initialisePDB(self, basestring PDBstr, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0, basestring element="", dict extra_atoms={}):
        """
            Load the information from a PDB file or a string with the PDB
            contents

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: basestring
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: basestring
            :param atomname: Residue name to select from the pdb (will only select the atoms with that name)
            :type atomname: basestring
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: basestring
            :param chain: Chain name to select from the pdb (will only select the atoms with that name)
            :type chain: basestring
            :param resnum: Residue number to select from the pdb (will only select the atoms with that name)
            :type atomname: int
            :raises: ValueError if the pdb contained no atoms
        """
        cdef object PDBContent
        cdef list stringWithPDBContent
        cdef int atomLineNum
        cdef basestring atomName, resName, atomLine, resnumStr
        cdef Atom atom
        if resnum == 0:
            resnumStr = u""
        else:
            resnumStr = u"%d" % (resnum)
        PDBContent = StringIO(readPDB(PDBstr))  # Using StringIO
        # creates a buffer that can handle a pdb file or a string containing
        # the PDB
        self.pdb = PDBContent.read()  # in case one wants to write it
        if extra_atoms != {}:
            CMAtoms = extra_atoms
        else:
            CMAtoms = self.CMAtoms


        stringWithPDBContent = self.pdb.split(u'\n')
        for atomLine in stringWithPDBContent:
            if not atomLine.startswith(u"ATOM") and not atomLine.startswith(u"HETATM"):
                continue
            if type == self._typeCM:
                atomName = atomLine[12:16].strip()
                resName = atomLine[17:20].strip()
                if resName not in CMAtoms:
                    continue
                if atomName != u"CA" and atomName != CMAtoms[resName]:
                    continue
            else:
                # HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
                if resname != u"" and not atomLine[17:20].strip() == resname:
                    continue
                if atomname != u"" and not atomLine[12:16].strip() == atomname:
                    continue
                if chain != u"" and not atomLine[21:22].strip() == chain:
                    continue
                if resnumStr != u"" and not atomLine[22:26].strip() == resnumStr:
                    continue
                if element != u"" and not atomLine[76:78].strip() == element:
                    continue

            atom = Atom(atomLine)
            # Here atom will be not null, empty or not.
            # With "try", we prune empty atoms
            try:
                if (not heavyAtoms or atom.isHeavyAtom()) and\
                   (type == self._typeAll or type == self._typeCM or (type == self._typeProtein and atom.isProtein()) or (type == self._typeHetero and atom.isHeteroAtom())):
                        self.atoms.update({atom.id: atom})
                        self.atomList.append(atom.id)
            except:
                pass
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def _initialiseXTC(self, np.ndarray[float, ndim=2] frame, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0, basestring element=u"", list topology=None, dict extra_atoms={}):
        """
            Load the information from a loaded XTC file into a  mdtraj Trajectory

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: basestring
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: basestring
            :param atomname: Residue name to select from the pdb (will only select the atoms with that name)
            :type atomname: basestring
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: basestring
            :param chain: Chain name to select from the pdb (will only select the atoms with that name)
            :type chain: basestring
            :param resnum: Residue number to select from the pdb (will only select the atoms with that name)
            :type atomname: int
            :raises: ValueError if the pdb contained no atoms
        """
        cdef int atomLineNum, atomIndex
        cdef basestring atomName, resName, atomLine, resnumStr, selection_string, element_atom, atomSerial, resChain
        cdef Atom atom
        cdef bint isProtein
        cdef float x, y, z
        cdef int iatom
        if resnum == 0:
            resnumStr = u""
        else:
            resnumStr = u"%d" % (resnum)
        self.pdb = self.join_PDB_lines(topology, frame)  # in case one wants to write it
        if extra_atoms != {}:
            CMAtoms = extra_atoms
        else:
            CMAtoms = self.CMAtoms
        for iatom in range(len(topology)):
            atomLine = topology[iatom]
            if not atomLine.startswith(u"ATOM") and not atomLine.startswith(u"HETATM"):
                continue
            if type == self._typeCM:
                atomName = atomLine[12:16].strip()
                resName = atomLine[17:20].strip()
                if resName not in CMAtoms:
                    continue
                if atomName != u"CA" and atomName != CMAtoms[resName]:
                    continue
            else:
                # HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
                if resname != u"" and not atomLine[17:20].strip() == resname:
                    continue
                if atomname != u"" and not atomLine[12:16].strip() == atomname:
                    continue
                if chain != u"" and not atomLine[21:22].strip() == chain:
                    continue
                if resnumStr != u"" and not atomLine[22:26].strip() == resnumStr:
                    continue
                if element != u"" and not atomLine[58:60].strip() == element:
                    continue
            isProtein = atomLine[:4] == u"ATOM"
            atomSerial = atomLine[6:11].strip()
            atomName = atomLine[12:16].strip()
            resName = atomLine[17:20].strip()
            resNum = atomLine[22:26].strip()
            resChain = atomLine[21]
            x = frame[iatom, 0]
            y = frame[iatom, 1]
            z = frame[iatom, 2]
            # due to the way the topology object is extracted the positions of
            # the element in 58-60
            element_atom = atomLine[58:60].strip().upper()
            atom = Atom()
            atom.set_properties(isProtein, atomSerial, atomName, resName, resNum, x, y, z, element_atom, resChain)
            if (not heavyAtoms or atom.isHeavyAtom()) and\
                (type == self._typeAll or type == self._typeCM or (type == self._typeProtein and atom.isProtein()) or (type == self._typeHetero and atom.isHeteroAtom())):

                self.atoms.update({atom.id: atom})
                self.atomList.append(atom.id)
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def initialise(self, object coordinates, bint heavyAtoms=True, basestring resname=u"", basestring atomname=u"", basestring type=u"ALL", basestring chain=u"", int resnum = 0, basestring element=u"", list topology=None, dict extra_atoms={}):
        """
            Wrapper function
        """
        if isinstance(coordinates, basestring):
            self.ispdb = True
            self._initialisePDB(coordinates, heavyAtoms, resname, atomname, type, chain, resnum, element, extra_atoms=extra_atoms)
        else:
            self.ispdb = False
            self._initialiseXTC(coordinates, heavyAtoms, resname, atomname, type, chain, resnum, element, topology=topology, extra_atoms=extra_atoms)

    def computeTotalMass(self):
        """
            Calculate the total mass of the PDB
        """
        cdef int atomNum
        self.totalMass = 0.0
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            self.totalMass += atom.mass

    def printAtoms(self):
        """
            Print Atom information for all the atoms in the PDB
        """
        cdef Atom atom
        for atom in self.atoms.values():
            print atom  # atom.printAtom()

    def getNumberOfAtoms(self):
        """
            Get the number of Atoms in the PDB

            :returns: int -- Number of atoms in the PDB
        """
        return len(self.atoms)

    def getAtom(self, atomId):
        """
            Get an Atom in the PDB by its id

            :param atomId: Id of the Atom (in the format "atomserial:atomName:resname")
            :type atomId: basestring
            :returns: int -- Number of atoms in the PDB
            :raises: KeyError if the id is not in the PDB
        """
        return self.atoms[atomId]

    def __len__(self):
        return len(self.atomList)


    def __getitem__(self, atomId):
        return self.atoms[atomId]

    def __setitem__(self, atomId, atom):
        self.atoms[atomId] = atom

    def __delitem__(self, atomId):
        self.atoms.pop(atomId)
        self.atomList.remove(atomId)

    def __iter__(self):
        for atomId in self.atomList:
            yield self.atoms[atomId]

    def updateCoords(self, newCoords):
        for atom, coords in zip(self.atomList, newCoords):
            atomObj = self.atoms[atom]
            atomObj.x = coords[0]
            atomObj.y = coords[1]
            atomObj.z = coords[2]

    def extractCOM(self):
        """
            Calculate the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if not self.totalMass:
            self.computeTotalMass()
        cdef list COM
        cdef int atomNum
        COM = [0., 0., 0.]
        for atomName in self.atoms:
            atom = self.atoms[atomName]
            COM[0] += atom.mass * atom.x
            COM[1] += atom.mass * atom.y
            COM[2] += atom.mass * atom.z

        COM[0] /= self.totalMass
        COM[1] /= self.totalMass
        COM[2] /= self.totalMass
        self.com = COM
        return COM

    def getCOM(self):
        """
            Get the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if self.com is None:
            return self.extractCOM()
        else:
            return self.com

    def extractCentroid(self):
        """
            Calculate the PDB centroid

            :returns: List -- List with the centroid coordinates
        """
        cdef list centroid
        cdef double n
        cdef int atomNum
        centroid = [0., 0., 0.]
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            centroid[0] += atom.x
            centroid[1] += atom.y
            centroid[2] += atom.z

        n = float(len(self.atoms))
        centroid[0] /= n
        centroid[1] /= n
        centroid[2] /= n
        self.centroid = centroid
        return centroid

    def getCentroid(self):
        """
            Get the PDB's centroid

            :returns: list -- List with the centroid coordinates
        """
        if self.centroid is None:
            return self.extractCentroid()
        else:
            return self.centroid


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def join_PDB_lines(self, list topology, np.ndarray[float, ndim=2] frame):
        cdef basestring prevLine = None
        cdef basestring temp = u"%.3f"
        cdef list pdb = []
        cdef basestring line
        cdef unsigned int natoms = len(topology)
        cdef unsigned int i
        for i in range(natoms):
            line = topology[i]
            if prevLine is not None and (prevLine[21] != line[21] or (prevLine[22:26] != line[22:26] and (u"HOH" == line[17:20] or u"HOH" == prevLine[17:20])) or (prevLine[0:4] == "ATOM" and line[0:6] == "HETATM") or (prevLine[0:6] == "HETATM" and line[0:6] == "HETATM" and (prevLine[22:26] != line[22:26]))):
                pdb.append(u"TER\n")
            x = (temp % frame[i, 0]).rjust(8)
            y = (temp % frame[i, 1]).rjust(8)
            z = (temp % frame[i, 2]).rjust(8)
            prevLine = line
            pdb.append(line % (x, y, z))
        return "".join(pdb)

    def get_pdb_string(self, int model_num=1):
        if self.ispdb:
            return self.pdb
        else:
            return "".join([u"MODEL    %4d\n" % model_num, self.pdb, u"ENDMDL\n", u"END\n"])

    def writePDB(self, basestring path, int model_num=1):
        """
            Write the pdb contents of the file from wich the PDB object was
            created

            :param path: Path of the file where to write the pdb
            :type path: basestring
        """
        cdef object fileHandle, atom
        if self.ispdb:
            with open(path, 'w', encoding="utf-8") as fileHandle:
                # in old simulations it will fail without the unicode
                fileHandle.write(unicode(self.pdb))
        else:
            with open(path, 'w', encoding="utf-8") as fileHandle:
                fileHandle.write(u"MODEL    %4d\n" % model_num)
                fileHandle.write(self.pdb)
                fileHandle.write(u"ENDMDL\n")
                fileHandle.write(u"END\n")

    def countContacts(self, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandChain=u""):
        """
            Count the number of alpha carbons that are in contact with the
            protein (i.e. less than contactThresholdDistance Amstrogms away)

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        cdef double contactThresholdDistance2,dist2
        contactThresholdDistance2= contactThresholdDistance**2

        cdef PDB ligandPDB, alphaCarbonsPDB

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, resnum=ligandResnum, chain=ligandChain, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, element=u"C",
                                   atomname=u"CA")
        # count contacts
        cdef set contacts = set([])
        cdef int rowind, colind
        cdef basestring proteinAtomId
        cdef Atom ligandAtom, proteinAtom
        for rowind in range(len(ligandPDB.atomList)):
        # can be optimised with cell list
            ligandAtom = ligandPDB.atoms[ligandPDB.atomList[rowind]]
            for colind in range(len(alphaCarbonsPDB.atomList)):
                proteinAtomId = alphaCarbonsPDB.atomList[colind]
                proteinAtom = alphaCarbonsPDB.atoms[proteinAtomId]
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if (dist2 - contactThresholdDistance2) < 0.1:
                    contacts.update([proteinAtomId])

        return len(contacts)


def computeCOMDifference(PDB1, PDB2):
    """
        Compute the difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The distance between the centers of mass between two PDB
    """
    return np.sqrt(computeCOMSquaredDifference(PDB1, PDB2))


def computeCOMSquaredDifference(PDB PDB1, PDB PDB2):
    """
        Compute the squared difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The squared distance between the centers of mass between two PDB
    """
    cdef list COM1, COM2
    COM1 = PDB1.getCOM()
    COM2 = PDB2.getCOM()

    dx = COM1[0] - COM2[0]
    dy = COM1[1] - COM2[1]
    dz = COM1[2] - COM2[2]

    return dx*dx + dy*dy + dz*dz


def computeSquaredCentroidDifference(PDB PDB1, PDB PDB2):
    """
        Compute the centroid squared difference between two PDBs

        :param PDB1: First PDB
        :type PDB1: PDB
        :param PDB2: Second PDB
        :type PDB2: PDB
        :returns: float -- The squared centroid distance between two PDB
    """
    cdef list centroid1, centroid2
    centroid1 = PDB1.getCentroid()
    centroid2 = PDB2.getCentroid()

    dx = centroid1[0] - centroid2[0]
    dy = centroid1[1] - centroid2[1]
    dz = centroid1[2] - centroid2[2]

    return dx*dx + dy*dy + dz*dz

def readPDB(pdbfile):
    """
        Helper function, parses a string with PDB content or the path of a pdb file into a string

        :param pdbfile: A string with PDB content or the path of a pdb file
        :type pdbfile: basestring
        :returns: basestring -- A string with PDB content
    """
    try:
        return open(str(pdbfile), "rt").read()
    except IOError:
        return pdbfile
