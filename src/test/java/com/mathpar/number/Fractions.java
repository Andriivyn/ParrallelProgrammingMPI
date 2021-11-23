package com.mathpar.number;

import com.mathpar.func.F;
import com.mathpar.polynom.Polynom;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author gennadi
 */
public class Fractions {
        
    
    
    public static void doTest(Ring ring, String  num, String  den, String result,  int fName) {
        System.out.println(" start="+num+"   "+den);
        Polynom n= new Polynom(num, ring);
        Polynom d= new Polynom(den, ring);  
        Element out= new F(result, ring);
        Element fr=new Fraction(n,d);
        Element res = (new F(fName, new Element[]{fr})).valueOf(ring);
        System.out.println(" res="+res.toString(ring));
        Element sub=res.subtract(out, ring).expand(ring);
        assertTrue(sub.isZero(ring));
    }

    @Test
    public void pol1_z() {
        doTest(new Ring("Q[z]"), "z^2","z+1", "z^2 ", F.NUM);
    }
}
    
   