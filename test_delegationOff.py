#!/usr/bin/env python3

import unittest
from molecule import Molecule

class TestMoleculeDelgNoPretest(unittest.TestCase):

    pretest = False

    good_data = [
        ('water', 'H2O', 18.01480, 
            [('H',2), ('O',1)], ['O1', 'H2']),
        ('salt', 'NaCl', 58.44300, 
            [('Cl' ,1), ('Na', 1)], ['Cl1', 'Na1']),
        ('ammoniumdichromate', '(NH4)2Cr2O7', 252.0622, 
            [ ('N',2), ('H', 8), ('Cr', 2), ('O', 7)], ['O7', 'Cr2', 'N2', 'H8']), 
        ('beryl', 'Be3Al2(SiO3)6', 537.4986, 
            [ ('Be', 3), ('Al',2), ('Si', 6), ('O', 18)], ['O18', 'Si6', 'Al2', 'Be3']), 
        ('isoaml acetate', 'CH3COO(CH2)2CH(CH3)2', 130.1856, 
            [('C', 7), ('H', 14),('O', 2)], ['C7', 'O2', 'H14']),
    ]

    def test_no_super_class(self):
        self.assertEqual(1, len(Molecule.__bases__), 
            "Molecule should single inherit from object")
        self.assertEqual(type(object()), Molecule.__bases__[0], 
            "Molecule superclass should be object")

    def test_str(self):
        if self.pretest: self.test_no_super_class()
        for d in self.good_data:
            m = Molecule(d[1])
            self.assertEqual(d[1], str(m),
                "str(Molecule({0})) should be {0}".format(d[1]))

    def _assertMethod(self, method ):
        # You wouldn't normally convert an error to a test failure like this.
        # unittest cases don't need this.
        # However, for the purposes of instruction, students using the test case do not
        # need to understand the test case error/fail distinction yet.
        self.assertTrue( method in dir(Molecule), 
            "Molecule missing method {}".format(method) )
                            
    def test_add_exists(self):
        if self.pretest: self.test_no_super_class()
        self._assertMethod('__add__')
                    
    def test_add(self):
        if self.pretest:
            self.test_no_super_class()
        self._assertMethod('__add__')
        for d1 in self.good_data:
            m1 = Molecule(d1[1])
            for d2 in self.good_data:
                m2 = Molecule(d2[1])
                s3 = d1[1] + d2[1]
                m3 = m1 + m2
                self.assertTrue( isinstance( m3, Molecule ),
                    "Molecule({0}) + Molecule({1}) is not a Molecule".format(
                        d1[1], d2[1] )) 
                        
    def test_replace(self):
        if self.pretest: self.test_no_super_class()
        self._assertMethod('replace')
        for d in self.good_data:
            m = Molecule(d[1])
            mr = m.replace("O","S")
            self.assertTrue( isinstance( mr, Molecule ) )
            self.assertEqual(d[1].replace("O","S"), str(mr) )
 
    def test_replace1(self):
        if self.pretest: self.test_no_super_class()
        self._assertMethod('replace')
        for d in self.good_data:
            m = Molecule(d[1])
            mr = m.replace("C","Si", 1)
            self.assertTrue( isinstance( mr, Molecule ) )
            self.assertEqual(d[1].replace("C","Si", 1), str(mr) )           
            
    def test_iter_delegate(self):
        if self.pretest: self.test_no_super_class()
        for d in self.good_data:
            m = Molecule(d[1])
            for (i,j) in zip( sorted(d[3]), m):
                self.assertEqual( i, j, msg="not iterating atoms correctly" )
            for (i, j) in zip( d[1], str(m) ):
                self.assertEqual( i, j, msg="iterating for str() incorrectly" )
                
    def test_atomix(self):
        if self.pretest: self.test_no_super_class()
        self._assertMethod('atomix')
        for d in self.good_data:
            m = Molecule(d[1])
            for am in m.atomix():
                self.assertTrue( isinstance( am, Molecule ), 
                    "Molecule.atomix() should generate Molecules" )
            self.assertSequenceEqual(d[4], list(map(str,m.atomix())),
                "Molecule({}).atomix() should generate {}".format(m,d[4]))              
                
if __name__ == '__main__':
    unittest.main(verbosity=2)

