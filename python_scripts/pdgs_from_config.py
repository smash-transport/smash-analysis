#/usr/bin/python
# coding=UTF-8
"""
This is a collection of functions to extract the pdg names and numbers
from the config.yaml file that is outputed with every SMASH run. It is meant
to be used by more task-oriented scripts via:
--------------------------
import pdgs_from_config
--------------------------
"""

from collections import defaultdict
from pdgcode import PDG_code

charge_str_neg = u"⁻"
charge_str_zero = u"⁰"
charge_str_pos = u"⁺"
charge_strs = charge_str_neg + charge_str_zero + charge_str_pos
bar = u"\u0305"


def strip_charge(s):
    # Protons and neutrons are a special case, because the don't carry a charge
    # string.
    if s in ("p", "n"):
        return "N"
    if s in ("p" + bar, "n" + bar):
        return "N" + bar
    return s.rstrip(charge_strs)


def chargestr(charge):
    switcher = {
        2: u"⁺⁺",
        1: u"⁺",
        0: u"⁰",
        -1: u"⁻",
        -2: u"⁻⁻",
    }
    return switcher.get(charge)


def get_charge(s):
    # Protons and neutrons are a special case, because the don't carry a charge
    # string.
    if s == "p":
        return 1
    if s == "n":
        return 0
    if s == "p" + bar:
        return -1
    if s == "n" + bar:
        return 0
    charge = 0
    found_charge = False
    for c in s:
        if c not in charge_strs:
            if found_charge:
                raise ValueError(
                    "Charge characters are not the last characters in particle name")
            continue
        else:
            found_charge = True
            if c == charge_str_pos:
                charge += 1
            if c == charge_str_neg:
                charge -= 1
    if not found_charge:
        raise ValueError("Did not find any charge letter")
    return charge


def antiname(name, pdg):
    # Protons and neutrons are a special case, because the don't carry a charge
    # string.
    if name in ("p", "n"):
        return name + bar
    if name == "p" + bar:
        return "p"
    if name == "n" + bar:
        return "n"
    basename = name
    charge = ""
    # loop over possible charges and look for charge string
    for ch in range(-2, 3):
        chstr = chargestr(ch)
        if (name.find(chstr) >= 0):
            basename = name[:len(name) - len(chstr)]
            charge = chargestr(-ch)
    # baryons & strange mesons: insert a bar
    if (pdg.baryon_number() != 0 or pdg.strangeness() != 0):
        basename = basename[:1] + bar + basename[1:]
    return basename + charge


def generic_name(name):
    # Protons and neutrons are a special case, because the don't carry a charge
    # string.
    if name in ("p", "n"):
        return "N"
    if name in ("p" + bar, "n" + bar):
        return "N" + bar
    # strip brackets and add star
    if name.find('(') >= 0:
        return name[:name.find('(')] + "*"
    else:
        return name


def invert_dict(d):
    inverted_d = {}
    for k, v in d.iteritems():
        inverted_d[v] = k
    return inverted_d

def is_string_nuclear_pdg(s):
    return (len(s) == 10 and s[0:2] == '10') or \
           (len(s) == 11 and s[0] == '-' and s[1:3] == '10')


pdg_dict = {}  # dictionary for particle-name mapping
inv_dict = {}  # inverse dictionary for particle-name mapping
pdg_dict_gen = {}  # dictionary for generic particle-name mapping


def pdg_to_name_init(configfile):
    for line in open(configfile, 'r'):
        if line.startswith("particles: "):  # find particles section
            line = line[line.find('"') + 1:line.rfind('"')]  # extract data
            for pline in line.split("\\n"):  # loop over particles
                pline = pline.strip()
                # remove comment lines
                if pline and not pline.startswith("#"):
                    # strip trailing comments
                    if pline.find('#') >= 0:
                        pline = pline[:pline.find('#')]
                    array = pline.split()
                    # get name string
                    name = unicode(array[0], encoding="UTF-8")
                    # For SMASH < 1.4, the third column is the first PDG code.
                    # For SMASH >= 1.4, it is the parity, and the codes start
                    # at the fourth column.
                    if array[3] in ('+', '-'):
                        first_pdg_code = 4
                    else:
                        first_pdg_code = 3
                    for i in array[first_pdg_code:]:  # loop over PDG codes
                        p = PDG_code(int(i, 16))
                        # Take care of nuclei
                        if is_string_nuclear_pdg(i):
                            anti = False if (i[0] != '-') else True
                            p.make_nucleus(nL = int(i[2]),
                                           Z = int(i[3:6]),
                                           A = int(i[6:9]),
                                           I = int(i[9]), anti = anti)
                            # print "Making nucleus pdg from", name,
                            # int(i[2]), int(i[3:6]), int(i[6:9]),
                            # int(i[9]), anti
                        full_name = name
                        gname = generic_name(name)
                        if len(array) > 5:  # more than one pdgcode given?
                            # add charge string
                            full_name += chargestr(p.charge())
                        pdg_dict[int(i)] = full_name.encode('UTF-8')
                        pdg_dict_gen[int(i)] = gname.encode('UTF-8')
                        if p.has_antiparticle():
                            pdg_dict[-int(i)] = antiname(full_name,
                                                         p).encode('UTF-8')
                            pdg_dict_gen[-int(i)] = antiname(gname,
                                                             p).encode('UTF-8')
    if len(pdg_dict) == 0 or len(pdg_dict_gen) == 0:
        raise ValueError('Could not read particles from "{}"'.format(configfile))

    # For nucleons, we use the historic symbols.
    pdg_dict[2212] = 'p'
    pdg_dict[-2212] = ('p' + bar).encode('UTF-8')
    pdg_dict[2112] = 'n'
    pdg_dict[-2112] = ('n' + bar).encode('UTF-8')


def pdg_to_name(pdg, configfile=""):
    """Return name of particle given its PDG code."""
    if len(pdg_dict) == 0:  # initialize dictionary if necessary
        pdg_to_name_init(configfile)
    return pdg_dict.get(pdg, unicode(pdg))

def pdg_to_name_generic(pdg, configfile=""):
    """Return generic name of particle given its PDG code."""
    if len(pdg_dict_gen) == 0:  # initialize dictionary if necessary
        pdg_to_name_init(configfile)
    return pdg_dict_gen.get(pdg, unicode(pdg))


def name_to_pdg(name, configfile=""):
    """Return PDG code of particle given its name."""
    global inv_dict
    if len(pdg_dict) == 0:  # initialize dictionary if necessary
        pdg_to_name_init(configfile)
    if len(inv_dict) == 0:
        inv_dict = invert_dict(pdg_dict)

        # support historic symbols for nucleons
        inv_dict['N\xe2\x81\xba'] = 2212
        inv_dict['N\xe2\x81\xb0'] = 2112
        inv_dict['N\xcc\x85\xe2\x81\xbb'] = -2212
        inv_dict['N\xcc\x85\xe2\x81\xb0'] = -2112

    return inv_dict.get(name, "unknown name: " + name)

if __name__ == '__main__':
    import sys
    config_file = sys.argv[1]
    pdg_to_name_init(config_file)
    pdgs = sorted(pdg_dict.keys(), key = lambda x: abs(x))
    for pdg in pdgs:
        print '%12i %s' % (pdg, pdg_dict.get(pdg, unicode(pdg)))
