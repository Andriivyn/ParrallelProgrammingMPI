package com.mathpar.number;

import static org.junit.Assert.*;
import org.junit.Test;

public class NumberR64Test {
    private static final Ring ring = new Ring("R64[]");
//
//    @Test
//    public void divideByR64GivesR64() {
//        assertTrue(new NumberR64("0.1").equals(
//                new NumberR64("1.0").divide(new NumberR64("10.0"), ring), ring));
//    }
    
    @Test
    public void toStringFLOAT3() {
           Ring r=new Ring("R64[x]"); r.setFLOATPOS( -1);
        NumberR64 a=(NumberR64)r.posConst[1];
        NumberR64 b; 
                b = new NumberR64("871.11");
         String ss=b.toString(r); 
         assertEquals( "870", ss);
    }
}
