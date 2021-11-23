package com.mathpar.number;

import static org.junit.Assert.*;
import org.junit.Test;

public class NumberRTest {
    @Test
    public void compareToNumberZ() {
        assertTrue(0 < new NumberR("1.2").compareTo(new NumberZ(1), Ring.ringZxyz));
    }
   @Test
    public void floor() {
        assertEquals(new NumberZ(4), new NumberR("4.1").floor(Ring.ringZxyz));
        assertEquals(new NumberZ(4), new NumberR("4.9").floor(Ring.ringZxyz));
    }
    @Test
    public void toStringFLOAT() {
           Ring r=new Ring("R[x]"); r.FLOATPOS=0;
        NumberR a=(NumberR)r.posConst[1];
        NumberR b= (NumberR)r.posConst[11];
        Element c= a.divide(b, r);
         String ss=c.toString(r); 
         assertEquals(ss, "0");
    }
        @Test
    public void toStringFLOAT2() {
           Ring r=new Ring("R[x]"); r.FLOATPOS=2;
        NumberR a=(NumberR)r.posConst[1];
        NumberR b=(NumberR)r.posConst[11];
        Element c= a.divide(b, r);
         String ss=c.toString(r); 
         assertEquals(ss, "0.09");
    }
            @Test
    public void toStringFLOAT3() {
           Ring r=new Ring("R[x]"); r.setFLOATPOS(-14);
        NumberR a=(NumberR)r.posConst[2];
        NumberR b = new NumberR("500");
         String ss=b.toString(r); 
         assertEquals( "500", ss);
    }

}
