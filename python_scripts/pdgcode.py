#/usr/bin/python
# coding=UTF-8

from ctypes import *


class PDG_bits(LittleEndianStructure):
    _fields_ = [  # first part of spin quantum number \f$n_J = 2 J + 1\f$.
                  ("n_J1", c_uint32, 4),
                  # third quark field
                  ("n_q3", c_uint32, 4),
                  # second quark field
                  ("n_q2", c_uint32, 4),
                  # first quark field. 0 for mesons.
                  ("n_q1", c_uint32, 4),
                  # "angular momentum"
                  ("n_L", c_uint32, 4),
                  # "radial excitation"
                  ("n_R", c_uint32, 4),
                  # "counter"
                  ("n", c_uint32, 4,),
                  # second part of n_J
                  ("n_J2", c_uint32, 4),
               ]

class nucleus_bits(LittleEndianStructure):
    _fields_ = [
           ("n_Lambda",   c_uint32,  6),
           ("Z",          c_uint32, 10),
           ("A",          c_uint32, 10),
           ("I",          c_uint32,  4),
           ("is_nucleus", c_uint32,  1),
           ("antipart",   c_uint32,  1),
          ]

class PDG_code(Union):
    _fields_ = [("b", PDG_bits), ("dump", c_uint32), ("nucl", nucleus_bits)]

    def __init__(self, dump=0):
        if dump < 0:
            self.anti = True
            dump = -dump
        else:
            self.anti = False
        self.dump = dump

    def is_nucleus(self):
        return self.nucl.is_nucleus

    def is_hadron(self):
        return (self.b.n_q2 != 0 and self.b.n_q3 != 0) and not self.is_nucleus()

    def is_baryon(self):
        return (self.is_hadron() and self.b.n_q1 != 0)

    def is_meson(self):
        return (self.is_hadron() and self.b.n_q1 == 0)

    def antiparticle_sign(self):
        return -1 if self.anti else +1

    def has_antiparticle(self):
        if (self.is_hadron()):
            return (self.baryon_number() != 0) or (self.b.n_q2 != self.b.n_q3)
        elif (self.is_nucleus()):
            return True
        else:
            return self.b.n_q3 == 1  # leptons!

    def baryon_number(self):
        if ((not self.is_hadron() or self.b.n_q1 == 0)):
            if (not self.is_nucleus()):
                return 0
            else:
                return self.antiparticle_sign() * self.nucl.A
        return self.antiparticle_sign()

    def net_quark_number(self, quark):
        # non-hadrons and those that have none of this quark type: 0
        if (not self.is_hadron() or (self.b.n_q1 != quark and
                                     self.b.n_q2 != quark and
                                     self.b.n_q3 != quark)):
            return 0
        # baryons: count quarks
        if self.baryon_number() != 0:
            # Nucleus is a special case
            if (self.is_nucleus()):
                Np = self.nucl.Z
                Nn = self.nucl.A - self.nucl.Z
                NL = self.nucl.n_Lambda
                if   (quark == 1): return (2*Nn + Np + NL) * antiparticle_sign()
                elif (quark == 2): return (Nn + 2*Np + NL) * antiparticle_sign()
                elif (quark == 3): return NL * antiparticle_sign()
                else:              return 0
            # for anti-baryons, the sign changes:
            return self.antiparticle_sign() * ((self.b.n_q1 == quark)
                                               + (self.b.n_q2 == quark)
                                               + (self.b.n_q3 == quark))
        # mesons: quarkonium state? Not open net_quark_number.
        if self.b.n_q3 == quark and self.b.n_q2 == quark:
            return 0
        # this has covered all the easy stuff
        # get the "other" quark. (We know this must exist, since they are
        # not both the right one and one of them is the right one).
        otherquark = self.b.n_q3 if (self.b.n_q2 == quark) else self.b.n_q2
        # "our" quark is the heavier one: 1 for u,c,t; -1 for d,s,b (and of
        # course the antiparticle sign)
        if quark > otherquark:
            return (1 if (quark % 2 == 0) else -1) * self.antiparticle_sign()
        # ours is the lighter: If the heavier particle is u,c,t, the lighter
        # one (ours) is an antiquark.
        return (-1 if (otherquark % 2 == 0) else 1) * self.antiparticle_sign()

    def make_nucleus(self, nL = 0, Z = 1, A = 2, I = 0, anti = 0):
        self.nucl.n_Lambda = nL
        self.nucl.Z = Z
        self.nucl.A = A
        self.nucl.I = I
        self.nucl.is_nucleus = 1
        self.nucl.antipart = anti

    def charge(self):
        if self.is_nucleus():
            return self.nucl.Z
        if self.is_hadron():
            # Q will accumulate 3*charge (please excuse the upper case. I
            # want to distinguish this from q which might be interpreted as
            # shorthand for "quark".)
            Q = 0
            # This loops over d,u,s,c,b,t quarks (the latter can be safely ignored,
            # but I don't think this will be a bottle neck.
            for i in range(1, 7):
                # u,c,t quarks have charge = 2/3 e, while d,s,b quarks have -1/3 e.
                # The antiparticle sign is already in net_quark_number.
                Q += (2 if (i % 2 == 0) else -1) * self.net_quark_number(i)
            return Q / 3
        # non-hadron:
        # Leptons: 11,13,15 are e,mu,tau and have a charge -1, while
        #          12,14,16 are the neutrinos that have no charge.
        if self.b.n_q3_ == 1:
            return -1 * ((self.b.n_J1 + self.b.n_J2) % 2) * self.antiparticle_sign()
        # Bosons: 24 is the W+, all else is uncharged.
        # we ignore the first digits so that this also finds strange gauge
        # boson "resonances" (in particular, \f$\tilde \chi_1^+\f$ with PDG
        # Code 1000024).
        if (self.dump & 0x0000ffff) == 0x24:
            return self.antiparticle_sign()
        # default (this includes all other Bosons) is 0.
        return 0

    def strangeness(self):
        return -self.net_quark_number(3)

    def n_J(self):
        return self.b.n_J1 + self.b.n_J2

if __name__ == '__main__':
    for p in [PDG_code(0x19922119), PDG_code(-0x19922119)]:
        print "n_J1={} n_q3={} n_q2={} n_q1={} n_L={} n_R={} n={} n_J2={} anti={} nucleus={}".format(
            p.b.n_J1, p.b.n_q3, p.b.n_q2, p.b.n_q1, p.b.n_L, p.b.n_R, p.b.n, p.b.n_J2, p.anti,
            p.is_nucleus())
    for p in [PDG_code(0x00000000)]:
        p.make_nucleus()
        print "Nucleus", p.is_nucleus(), "charge", p.charge(), \
              "Is baryon", p.is_baryon(), "baryon_number", p.baryon_number(), \
              "antiparticle sign", p.antiparticle_sign(), \
              "has antiparticle", p.has_antiparticle()
