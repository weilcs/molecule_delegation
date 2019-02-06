#!/usr/bin/env python3

# Molecular mass caluculation class for example code
# Modified by E Brown, 2018-2019 
# Qfrom https://gist.github.com/Rhomboid/5994999#file-example-py-L48

"""Molecular mass calculation

    Classes:
        Molecule - the molecular formula as a string

    Examples:
        <code>
            >>> x = Molecule("H2SO4")
            >>> print(x.mass())
            98.07679999999999
        </code>
"""

import re
import sys
import collections

__all__ = ["Molecule"]

_atomic_mass = {
    "H": 1.0079, "He": 4.0026, "Li": 6.941, "Be": 9.0122,
    "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982,
    "Si": 28.086, "P": 30.974, "S": 32.065, "Cl": 35.453,
    "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956,
    "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.39, "Ga": 69.723, "Ge": 72.61, "As": 74.922,
    "Se":78.96, "Br": 79.904, "Kr": 83.80, "Rb": 85.468, "Sr": 87.62,
    "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.94,
    "Tc": 97.61, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42,
    "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
    "Pr": 140.91, "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96,
    "Gd": 157.25, "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26,
    "Tm": 168.93, "Yb": 173.04, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95,
    "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 196.08,
    "Au": 196.08, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98,
    "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0,
    "Ac": 227.0, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237.0,
    "Pu": 244.0, "Am": 243.0, "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0,
    "Fm": 257.0, "Md": 258.0, "No": 259.0, "Lr": 262.0, "Rf": 261.0, "Db": 262.0,
    "Sg": 266.0, "Bh": 264.0, "Hs": 269.0, "Mt": 268.0
}

class Molecule:
    """Molecular formula calculations

    Methods:
        mass - the molecular mass of the formula
        check_symbols - test if formula symbols are recognized 
        clean_copy - return a copy with unrecognized symbols deleted
        atoms - count the atoms of one type in the formula
       
    """
    def __init__(self, mString = ""):
        super().__init__()
        self.dele_mole = mString


    def __str__(self):
        return str(self.dele_mole)

    def __add__(self, other):
        return Molecule(str(self)+str(other))

    def replace(self, string1, string2, *args):
        if len(args)<1:
            new_string = self.dele_mole.replace(string1, string2)
        elif len(args)==1:
            new_string = self.dele_mole.replace(string1, string2, args[0])
        return Molecule(new_string)
        
        

    def _gettokens(self):
        # divide the string up into contituent formula tokens
        # return a list of tokens
        # for example, "Ca(NO3)2" becomes [ 'Ca', '(', 'N', 'O', '3', ')', '2' ]
        return re.findall(r'[A-Z][a-z]*|\d+|\(|\)|\S', self.dele_mole)

    def _tokentree(self, target=None):
        # divide the string up into constituent formula tuples
        # return a heirarchy of tokens tuples (a token tree)
        # for example, "Ca(NO3)2" returns [ ('Ca', 1), ([ ('N', 1), ('O', 3) ], 2) ]
        # this could be a lot eaiser with an re package supporting recurision like PCRE
        
        if not target: target=self
        
        # this makes top level symbols and groups into tuples
        reslist = re.findall(r'(\([^\)]*\)|[A-Z][a-z]*)(\d*)', str(target))
        
        # process a tuple from re resultsslist.
        # recursively find tuples with brackets '(\(.*\)' in first item 
        # and turn them into token-subtrees too
        def _tokensubtree(treetuple): 
            def numfor(s):
                # how to convert number of atoms from str to int
                return 1 if s == '' else int(s)
       
            if treetuple[0][0] =='(' : # need a subtree here
                return (self._tokentree(treetuple[0][1:-1]) , numfor(treetuple[1]))
            else: 
                return (treetuple[0], numfor(treetuple[1])) 
                          
        return list(map(_tokensubtree, reslist))
    
    def _tree_mass(self, tree):
        # calculate mass of a token tree
        
        def token_mass( node ):
            # calculate mass of single tree tuple (a tree node)
            if isinstance(node[0], str): # a leaf node (symbol)
                if not node[0] in _atomic_mass:
                    raise KeyError("not an atomic symbol " + node[0])
                else: nodemass = _atomic_mass[node[0]]
            else: nodemass = self._tree_mass(node[0]) # a subtree (list)
            
            natoms = node[1]
            
            return nodemass * natoms
            
        return sum(map(token_mass, tree))
    
    def mass(self):
        """Compute and return the molecular mass from constituent atoms

            Raises:
              KeyError - if an aplhabetic symbol is not in atomic periodic table
              ValueError - if some other symbol is bad 
                
            return: the molecular mass of the Molecule
        """
        self.check_symbols()
        return self._tree_mass(self._tokentree()) 

    def check_symbols(self):
        """ Return true if all the token symbols in the formula are valid
            Note this does not mean the molecule is chemically viable.
         
            Raises:
              KeyError - if an aplhabetic symbol is not in atomic periodic table
              ValueError - if some other symbol is bad 
            
        """
        # this method has a bug in that it never raises KeyError, it raises 
        # ValueError instead.
        
        def is_valid(sym):
            # what symbols are valid? (, ), digits, atoms
            if sym in "()": return True
            #if sym.isdigit(): return True
            #if sym in _atomic_mass: return True
            if sym.isalnum(): return True
            return False

        for t in self._gettokens():
            if not is_valid(t): raise ValueError("bad symbol " + t)
            if t.isalpha() and t not in _atomic_mass: raise KeyError("key error " + t)
        return True

    def clean_copy(self):
        """Return a copy of the molecule with invalid tokens removed"""
        # this is a stub implementation
        #return Molecule("H2O")
        m = self._gettokens()
        for t in self._gettokens():
            #if there is value errors or key errors, remove the invalid tokens
            if (t.isalpha() and t not in _atomic_mass) or (t not in "()" and not t.isalnum()):
                m.remove(t)
        str2 = "".join(m)      
        return Molecule(str2)
        
    def __iter__(self):
        # this is not yet implemented
        # this will be over-riding the iter method for strings, so
        # we can no longer treat this class like a regular string class
        # for iteration purposes
        # the iterator should yield tuples with the total number of each atom
        # type in the molecule, so Be3Al2(SiO3)6 would yeild:
        #  ('Al', 2) then ('Be', 3 ) then ( 'O', 18 ) and so on.
        # tuples must be in sorted order.

        #use regular expression to find parse
        parse = re.findall(r'([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)', str(self.dele_mole))
        #use stack to store the number of atoms
        sym_num = [collections.Counter()]
        for name, n1, left_open, right_open, n2 in parse:
            if name:
                sym_num[-1][name] += int(n1 or 1)                   
            if left_open:
                sym_num.append(collections.Counter())
            if right_open:
                top = sym_num.pop()
                for s in top:
                    sym_num[-1][s] += top[s] * int(n2 or 1)
        for name in sorted(sym_num[-1]):
            yield (name, sym_num[-1][name])
        
    def atoms(self, symbol):
        """Return the number of type 'symbol' atoms in the molecule        
         
           Raises:
              KeyError - if the symbol is not in atomic periodic table

        """ 
        # this is a stub implementation
        #return 10;
        if symbol not in _atomic_mass: raise KeyError( symbol + " is not in the table")
        if symbol in _atomic_mass and symbol not in self._gettokens():
            return 0
        #the method is similar to __iter__, just different return
        parse = re.findall(r'([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)', str(self.dele_mole))
        if symbol in _atomic_mass and symbol in self._gettokens():
            sym_num = [collections.Counter()]
            for name, n1, left_open, right_open, n2 in parse:
                if name:
                    sym_num[-1][name] += int(n1 or 1)                   
                if left_open:
                    sym_num.append(collections.Counter())
                if right_open:
                    top = sym_num.pop()
                    for s in top:
                        sym_num[-1][s] += top[s] * int(n2 or 1)        
            return sym_num[-1][symbol]

    def atomix(self):
        parse = re.findall(r'([A-Z][a-z]*)(\d*)|(\()|(\))(\d*)', str(self.dele_mole))
        #use stack to store the number of atoms
        sym_num = [collections.Counter()]
        for name, n1, left_open, right_open, n2 in parse:
            if name:
                sym_num[-1][name] += int(n1 or 1)                   
            if left_open:
                sym_num.append(collections.Counter())
            if right_open:
                top = sym_num.pop()
                for s in top:
                    sym_num[-1][s] += top[s] * int(n2 or 1)
        def atom_mass(atm):
            if atm in _atomic_mass:
                return _atomic_mass[atm] * self.atoms(atm)
        for name in sorted(sym_num[-1], key = atom_mass, reverse = True):
            yield Molecule(name+str(sym_num[-1][name]))
 
        
if __name__ == "__main__":
    import sys
    while True:
        print('\n====\n')
        formula = input('Enter molecular formula (Q to quit): ')
        if formula[0] == 'Q':
            sys.exit(0)
            
        m = Molecule(formula)

        print('The molecular mass of {} is {:.3f}\n'.format(m, m.mass()))
        
        print('The elements of {} are:'.format(m))
        for a in m:
            print(a)

        print(str(m))
        print(m._gettokens())
        print(m._tokentree())
        print(iter(m))
        print(m.atoms("N"))
        m1 = Molecule("HNO(CO3)4")
        print(m1)
        m2 = Molecule("CCO(NO2)2CO")
        print(m2)
        print(m1+m2)

        
