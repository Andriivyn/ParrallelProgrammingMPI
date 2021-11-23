/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mathpar.number;

import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author gennadi
 */
public class NumberZandCQTest {
        
               @Test
    public void toStringCQ() { 
         Ring ring = new Ring("CQ[x]");
        Element a = new NumberZ("4");
        Complex b = new Complex(new NumberZ("0"), new NumberZ("-2"));  // b = -2i
        Element c = b.multiply(b, ring) ;              // c = (-2i)*(-2i)/(-2i)
        Element d = c.divide(b, ring); 
        Element e = a.divideExact(d, ring);
        d=d.subtract(b, ring);
        assertTrue(d.isZero(ring));
    }
                  @Test
    public void  Zmod() { 
         Ring ring = new Ring("Z[x]");
        NumberZ a = new NumberZ("9");
       NumberZ b = new NumberZ("13");  
        NumberZ c = a.Mod(b,ring).add(new NumberZ("4")) ;   
        NumberZ d = a.mod(b,ring).subtract(a) ;  
        assertTrue(c.isZero(ring)&&d.isZero(ring));
    } 
    
    
    
}
