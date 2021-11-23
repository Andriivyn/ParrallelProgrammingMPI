package com.mathpar.number;
import com.mathpar.func.F;
import com.mathpar.func.parser.Parser;
import org.junit.Test;
import static org.junit.Assert.assertTrue;

/**
 * @author serega
 */
public class NumberZMinMultTest {

    private static final Ring r = new Ring("ZMinMult[x]");

    @Test
    public void addMinMult() {
        NumberZMinMult a = new NumberZMinMult(new NumberZ(2));
        NumberZMinMult b = new NumberZMinMult(new NumberZ(3));
        Element c = a.add(b, r);
        assertTrue(c.equals(new NumberZMinMult(new NumberZ(2)), r));
    }

    @Test
    public void multiplyMinMult() {
        NumberZMinMult a = new NumberZMinMult(new NumberZ(2));
        NumberZMinMult b = new NumberZMinMult(new NumberZ(3));
        Element c = a.multiply(b, r);
        assertTrue(c.equals(new NumberZMinMult(new NumberZ(6)), r));
    }

    @Test
    public void funcBooleanHasCorrectPriority2() {
        String s = "a=3";
        Ring r = new Ring("ZMinMult[x,y]");
        F f = Parser.getF("3", r);
        assertTrue(f.toString(r).equals("3"));
    }
}
