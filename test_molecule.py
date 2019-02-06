#!/usr/bin/env python3

import unittest
from molecule import Molecule

class TestMoleculeMethods(unittest.TestCase):

    good_data = [
        ('water', 'H2O', 18.01480, [('H',2), ('O',1)]),
        ('salt', 'NaCl', 58.44300, [('Cl' ,1), ('Na', 1)]),
        ('ammoniumdichromate', '(NH4)2Cr2O7', 252.0622, [ ('N',2), ('H', 8), ('Cr', 2), ('O', 7)]), 
        ('beryl', 'Be3Al2(SiO3)6', 537.4986, [ ('Be', 3), ('Al',2), ('Si', 6), ('O', 18)]), 
        ('isoaml acetate', 'CH3COO(CH2)2CH(CH3)2', 130.1856, [('C', 7), ('H', 14),('O', 2)]),
    ]

    def test_mass_water(self):
        water = Molecule('H2O') 
        self.assertAlmostEqual(18.0148, water.mass(),  msg="water mass incorrect")

    def test_mass_badsalt(self):
        badsalt = Molecule('NalC')
        with self.assertRaises(Exception,msg="did not detect bad symbol in NalC"):
                 badsalt.mass()

    def test_check_symbols(self):
        salt = Molecule('NaCl')
        water = Molecule('H2O')

        self.assertTrue( salt.check_symbols() )
        self.assertTrue( water.check_symbols() )
 
    def test_check_symbols_error(self):      
        badwater = Molecule('H+2+O')
        with self.assertRaises(ValueError):
            badwater.check_symbols() 
        
    def test_mass(self):
        for d in self.good_data:
            m = Molecule(d[1])
            self.assertAlmostEqual( d[2], m.mass(),
                                    msg="mass of {} {} incorrect".format(d[0],d[1]))
        
    def test_mass_error(self):
        badwater = Molecule('H+2+O')
        with self.assertRaises(ValueError):
            badwater.mass() 
          
    # the following are failing:       

    def test_mass_keyerror(self):
        badwater = Molecule('Yz2O')
        try:
            badwater.mass()
        except KeyError:
            pass
        except Exception as ex:
            self.fail("Wrong exception {} should be KeyError".format(
                        type(ex).__name__))
            
    def test_clean_copy(self):
        badammoniumdichromate = Molecule('(NH4)2Cr2XqO7')
        fixed = badammoniumdichromate.clean_copy()
        self.assertTrue( fixed.check_symbols() )
        self.assertEqual('(NH4)2Cr2O7', str(fixed) )
        
    def test_atoms(self):
        for d in self.good_data:
            m = Molecule(d[1])
            for (symbol, count) in d[3]:
                self.assertEqual( count, m.atoms(symbol), 
                                  msg="atom {} in {}".format( symbol, m) )        
    
    def test_iter_tuples(self):
        m = Molecule('CH3COO(CH2)2CH(CH3)2')
        miter = iter(m)
        with self.assertRaises(StopIteration):
            while True:
                self.assertEqual(2, len(next(miter)))
                
    def test_iter_atoms(self):
        for d in self.good_data:
            m = Molecule(d[1])
            for atomtuple in m:
                self.assertTrue( atomtuple in d[3])
   
    def test_iter_molecules(self):
        for d in self.good_data:
            m = Molecule(d[1])
            for (i,j) in zip( sorted(d[3]), m):
                self.assertEqual( i, j, msg="not iterating atoms correctly" )
                
    def test_iter_parallel(self):
        for d in self.good_data:
            m = Molecule(d[1])
            with self.assertRaises(StopIteration):
                i1 = iter(m)
                i2 = iter(m)
                i3 = iter(sorted(d[3]))
                while True:
                    n1 = next(i1)
                    n2 = next(i2)
                    n3 = next(i3)
                    self.assertEqual(n1,n2)
                    self.assertEqual(n1,n3)

    def test_iter_nesting(self):
        for d in self.good_data:
          m = Molecule(d[1])
          loops = 0
          for i1 in m:
            found = False
            for i2 in m:
              if i2 == i1 : found = True
              loops += 1
            self.assertTrue(found, msg="Nested iteration fails on {}".format(i2))
          self.assertEqual( loops, len(d[3]) * len(d[3]))        
                
    def test_atomsbad(self):
        m = Molecule('Be3Al2(SiO3)6')
        self.assertEqual(0, m.atoms('Ca'))
        
    def test_atomserror(self):
        m = Molecule('Be3Al2(SiO3)6')
        with self.assertRaises(KeyError):
            m.atoms('Ck')          

if __name__ == '__main__':
    unittest.main(verbosity=2)

