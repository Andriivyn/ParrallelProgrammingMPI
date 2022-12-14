
/*
 * @(#)BigDecimal.java	1.53 04/06/28
 *
 * Copyright 2004 Sun Microsystems, Inc. All rights reserved.
 * SUN PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
 */

/*
 * @(#)BigDecimal.java  1.x 01/xx/xx
 *
 * Copyright 1996-2001 Sun Microsystems, Inc. All Rights Reserved.
 * Portions Copyright IBM Corporation, 2001. All Rights Reserved.
 *
 * This software is the proprietary information of Sun Microsystems, Inc.
 * Use is subject to license terms.
 *
 */
package com.mathpar.number;

import com.mathpar.func.F;
import com.mathpar.func.FuncNumberR;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Random;
import com.mathpar.number.math.MathContext;
import com.mathpar.number.math.RoundingMode;
import com.mathpar.polynom.PascalTriangle;
import com.mathpar.polynom.Polynom;

/**
 * Immutable, arbitrary-precision signed decimal numbers. A <tt>BigDecimal</tt>
 * consists of an arbitrary precision integer <i>unscaled value</i> and a 32-bit
 * integer <i>scale</i>. If zero or positive, the scale is the number of digits
 * to the right of the decimal point. If negative, the unscaled value of the
 * number is multiplied by ten to the power of the negation of the scale. The
 * value of the number represented by the <tt>BigDecimal</tt> is therefore
 * <tt>(unscaledValue &times; 10<sup>-scale</sup>)</tt>.
 *
 * <p>
 * The <tt>BigDecimal</tt> class provides operations for arithmetic, scale
 * manipulation, rounding, comparison, hashing, and format conversion. The
 * {@link #toString} method provides a canonical representation of a
 * <tt>BigDecimal</tt>.
 *
 * <p>
 * The <tt>BigDecimal</tt> class gives its user complete control over rounding
 * behavior. If no rounding mode is specified and the exact result cannot be
 * represented, an exception is thrown; otherwise, calculations can be carried
 * out to a chosen precision and rounding mode by supplying an appropriate
 * {@link MathContext} object to the operation. In either case, eight
 * <em>rounding modes</em> are provided for the control of rounding. Using the
 * integer fields in this class (such as {@link #ROUND_HALF_UP}) to represent
 * rounding mode is largely obsolete; the enumeration values of the
 * <tt>RoundingMode</tt> <tt>enum</tt>, (such as {@link
 * RoundingMode#HALF_UP}) should be used instead.
 *
 * <p>
 * When a <tt>MathContext</tt> object is supplied with a precision setting of 0
 * (for example, {@link MathContext#UNLIMITED}), arithmetic operations are
 * exact, as are the arithmetic methods which take no <tt>MathContext</tt>
 * object. (This is the only behavior that was supported in releases prior to
 * 5.) As a corollary of computing the exact result, the rounding mode setting
 * of a <tt>MathContext</tt> object with a precision setting of 0 is not used
 * and thus irrelevant. In the case of divide, the exact quotient could have an
 * infinitely long decimal expansion; for example, 1 divided by 3. If the
 * quotient has a nonterminating decimal expansion and the operation is
 * specified to return an exact result, an <tt>ArithmeticException</tt> is
 * thrown. Otherwise, the exact result of the division is returned, as done for
 * other operations.
 *
 * <p>
 * When the precision setting is not 0, the rules of <tt>BigDecimal</tt>
 * arithmetic are broadly compatible with selected modes of operation of the
 * arithmetic defined in ANSI X3.274-1996 and ANSI X3.274-1996/AM 1-2000
 * (section 7.4). Unlike those standards, <tt>BigDecimal</tt> includes many
 * rounding modes, which were mandatory for division in <tt>BigDecimal</tt>
 * releases prior to 5. Any conflicts between these ANSI standards and the
 * <tt>BigDecimal</tt> specification are resolved in favor of
 * <tt>BigDecimal</tt>.
 *
 * <p>
 * Since the same numerical value can have different representations (with
 * different scales), the rules of arithmetic and rounding must specify both the
 * numerical result and the scale used in the result's representation.
 *
 *
 * <p>
 * In general the rounding modes and precision setting determine how operations
 * return results with a limited number of digits when the exact result has more
 * digits (perhaps infinitely many in the case of division) than the number of
 * digits returned.
 *
 * First, the total number of digits to return is specified by the
 * <tt>MathContext</tt>'s <tt>precision</tt> setting; this determines the
 * result's <i>precision</i>. The digit count starts from the leftmost nonzero
 * digit of the exact result. The rounding mode determines how any discarded
 * trailing digits affect the returned result.
 *
 * <p>
 * For all arithmetic operators , the operation is carried out as though an
 * exact intermediate result were first calculated and then rounded to the
 * number of digits specified by the precision setting (if necessary), using the
 * selected rounding mode. If the exact result is not returned, some digit
 * positions of the exact result are discarded. When rounding increases the
 * magnitude of the returned result, it is possible for a new digit position to
 * be created by a carry propagating to a leading &quot;9&quot; digit. For
 * example, rounding the value 999.9 to three digits rounding up would be
 * numerically equal to one thousand, represented as 100&times;10<sup>1</sup>.
 * In such cases, the new &quot;1&quot; is the leading digit position of the
 * returned result.
 *
 * <p>
 * Besides a logical exact result, each arithmetic operation has a preferred
 * scale for representing a result. The preferred scale for each operation is
 * listed in the table below.
 *
 * <table border> <caption top><h3>Preferred Scales for Results of Arithmetic
 * Operations </h3></caption> <tr><th>Operation</th><th>Preferred Scale of
 * Result</th></tr> <tr><td>Add</td><td>max(addend.scale(), augend.scale())</td>
 * <tr><td>Subtract</td><td>max(minuend.scale(), subtrahend.scale())</td>
 * <tr><td>Multiply</td><td>multiplier.scale() + multiplicand.scale()</td>
 * <tr><td>Divide</td><td>dividend.scale() - divisor.scale()</td> </table>
 *
 * These scales are the ones used by the methods which return exact arithmetic
 * results; except that an exact divide may have to use a larger scale since the
 * exact result may have more digits. For example, <tt>1/32</tt> is
 * <tt>0.03125</tt>.
 *
 * <p>
 * Before rounding, the scale of the logical exact intermediate result is the
 * preferred scale for that operation. If the exact numerical result cannot be
 * represented in <code>precision</code> digits, rounding selects the set of
 * digits to return and the scale of the result is reduced from the scale of the
 * intermediate result to the least scale which can represent the
 * <code>precision</code> digits actually returned. If the exact result can be
 * represented with at most <code>precision</code> digits, the representation of
 * the result with the scale closest to the preferred scale is returned. In
 * particular, an exactly representable quotient may be represented in fewer
 * than <code>precision</code> digits by removing trailing zeros and decreasing
 * the scale. For example, rounding to three digits using the
 * {@linkplain RoundingMode#FLOOR floor} rounding mode, <br>
 *
 * <code>19/100 = 0.19   // integer=19,  scale=2</code> <br>
 *
 * but<br>
 *
 * <code>21/110 = 0.190  // integer=190, scale=3</code> <br>
 *
 * <p>
 * Note that for add, subtract, and multiply, the reduction in scale will equal
 * the number of digit positions of the exact result which are discarded. If the
 * rounding causes a carry propagation to create a new high-order digit
 * position, an additional digit of the result is discarded than when no new
 * digit position is created.
 *
 * <p>
 * Other methods may have slightly different rounding semantics. For example,
 * the result of the <tt>pow</tt> method using the
 * {@linkplain #pow(int, MathContext) specified algorithm} can occasionally
 * differ from the rounded mathematical result by more than one unit in the last
 * place, one <i>{@linkplain #ulp() ulp}</i>.
 *
 * <p>
 * Two types of operations are provided for manipulating the scale of a
 * <tt>BigDecimal</tt>: scaling/rounding operations and decimal point motion
 * operations. Scaling/rounding operations ({@link
 * #setScale setScale} and {@link #round round}) return a <tt>BigDecimal</tt>
 * whose value is approximately (or exactly) equal to that of the operand, but
 * whose scale or precision is the specified value; that is, they increase or
 * decrease the precision of the stored number with minimal effect on its value.
 * Decimal point motion operations ({@link #movePointLeft movePointLeft} and
 * {@link #movePointRight movePointRight}) return a <tt>BigDecimal</tt> created
 * from the operand by moving the decimal point a specified distance in the
 * specified direction.
 *
 * <p>
 * For the sake of brevity and clarity, pseudo-code is used throughout the
 * descriptions of <tt>BigDecimal</tt> methods. The pseudo-code expression
 * <tt>(i + j)</tt> is shorthand for &quot;a <tt>BigDecimal</tt> whose value is
 * that of the <tt>BigDecimal</tt> <tt>i</tt> added to that of the
 * <tt>BigDecimal</tt> <tt>j</tt>.&quot; The pseudo-code expression <tt>(i ==
 * j)</tt> is shorthand for &quot;<tt>true</tt> if and only if the
 * <tt>BigDecimal</tt> <tt>i</tt> represents the same value as the
 * <tt>BigDecimal</tt> <tt>j</tt>.&quot; Other pseudo-code expressions are
 * interpreted similarly. Square brackets are used to represent the particular
 * <tt>BigInteger</tt> and scale pair defining a <tt>BigDecimal</tt> value; for
 * example [19, 2] is the <tt>BigDecimal</tt> numerically equal to 0.19 having a
 * scale of 2.
 *
 * <p>
 * Note: care should be exercised if <tt>BigDecimal</tt> objects are used as
 * keys in a {@link java.util.SortedMap SortedMap} or elements in a
 * {@link java.util.SortedSet SortedSet} since <tt>BigDecimal</tt>'s <i>natural
 * ordering</i> is <i>inconsistent with equals</i>. See {@link Comparable}, {@link
 * java.util.SortedMap} or {@link java.util.SortedSet} for more information.
 *
 * <p>
 * All methods and constructors for this class throw
 * <tt>NullPointerException</tt> when passed a <tt>null</tt> object reference
 * for any input parameter.
 *
 * @see BigInteger
 * @see MathContext
 * @see RoundingMode
 * @see java.util.SortedMap
 * @see java.util.SortedSet
 * @author Josh Bloch
 * @author Mike Cowlishaw
 * @author Joseph D. Darcy
 */
public class NumberR extends Element {
    public NumberR() {
        scale = 0;
        intVal = NumberZ.valueOf(0);
    }
    public static final NumberR POZITIVE_ZERO = new NumberR(NumberZ.ZERO, Integer.MAX_VALUE - 2);
    public static final NumberR NEGATIVE_ZERO = new NumberR(NumberZ.ZERO, Integer.MAX_VALUE - 3);
    public static NumberR eps = new NumberR(0.0001);

    static {
        PascalTriangle.initCacheZ();
    }
    protected NumberZ intVal;
    /**
     * The scale of this BigDecimal, as returned by {@link #scale}.
     *
     * @serial
     * @see #scale
     */
    protected int scale = 0;  // Note: this may have any value, so
    // calculations must be done in longs
    /**
     * The number of decimal digits in this BigDecimal, or 0 if the number of
     * digits are not known (lookaside information). If nonzero, the value is
     * guaranteed correct. Use the precision() method to obtain and set the
     * value if it might be 0. This field is mutable until set nonzero.
     *
     * @since 1.5
     */
    volatile transient int precision = 0;

    /*
     * Appease the serialization gods
     */
    private static final long serialVersionUID = 6108874887143696463L;
    // Constants
    /**
     * The value 0, with a scale of 0.
     *
     * @since 1.5
     */
    public static final NumberR ZERO = new NumberR(NumberZ.ZERO);   
    public static final NumberR  POSCONST[] = new NumberR[MAX_CONSTANT + 1];
    public static final NumberR  NEGCONST[] = new NumberR[MAX_CONSTANT + 1];
    static {
        for (int i = 1; i <= MAX_CONSTANT; i++) {
            POSCONST[i] = new NumberR(NumberZ.POSCONST[i], 0);
            NEGCONST[i] = new NumberR(NumberZ.NEGCONST[i], 0);
        }
        POSCONST[0]=NumberR.ZERO;
        PascalTriangle.initCacheZ();
    }
    
       /**
     * Initialize   constant array  
     * @return  value of 1,2,..,MAX_CONSTANT + 1
     */
    @Override
    public Element[] posConst() {return  POSCONST;}
    @Override
    public Element[] negConst() {return  NEGCONST;}
    public static final NumberR ONE = POSCONST[1];
    public static final NumberR TEN = POSCONST[10]; 
    public static final NumberR MINUS_ONE = NEGCONST[1];
    public static final NumberR TWO = POSCONST[2];
    public static final NumberR R_180= new NumberR(NumberZ.Z_180, 0);
    public static final NumberR R_239= new NumberR(NumberZ.Z_239, 0);
    
    public static NumberR MachineEpsilonR = new NumberR("1.0e-29");
    public  NumberR[]   constR=new NumberR[6];
    public static final int ACCURACYSTATIC = 34;
    /**
     * Initial values for constR=e,pi,..for initial accuracyStatic
     */
 //   public static NumberR[] CONSTRINIT=null;
//    static{  
//        NumberR[] setConstRInit = FuncNumberR.setConstRInit(new NumberR[6], ACCURACYSTATIC);
//    CONSTRINIT=setConstRInit;
//}
    
    
    // Constructors
    /**
     * Translates a character array representation of a <tt>BigDecimal</tt> into
     * a <tt>BigDecimal</tt>, accepting the same sequence of characters as the
     * {@link #BigDecimal(String)} constructor, while allowing a sub-array to be
     * specified.
     *
     * <p>
     * Note that if the sequence of characters is already available within a
     * character array, using this constructor is faster than converting the
     * <tt>char</tt> array to string and using the <tt>BigDecimal(String)</tt>
     * constructor .
     *
     * @param in <tt>char</tt> array that is the source of characters.
     * @param offset first character in the array to inspect.
     * @param len number of characters to consider.
     *
     * @throws NumberFormatException if <tt>in</tt> is not a valid
     * representation of a <tt>BigDecimal</tt> or the defined subarray is not
     * wholly within <tt>in</tt>.
     * @since 1.5
     */
    public NumberR(char[] in, int offset, int len) {
        // This is the primary string to BigDecimal constructor; all
        // incoming strings end up here; it uses explicit (inline)
        // parsing for speed and generates at most one intermediate
        // (temporary) object (a char[] array).

        // use array bounds checking to handle too-long, len == 0,
        // bad offset, etc.
        try {
            // handle the sign
            boolean isneg = false;          // assume positive
            if (in[offset] == '-') {
                isneg = true;               // leading minus means negative
                offset++;
                len--;
            } else if (in[offset] == '+') { // leading + allowed
                offset++;
                len--;
            }

            // should now be at numeric part of the significand
            int dotoff = -1;                 // '.' offset, -1 if none
            int cfirst = offset;             // record start of integer
            long exp = 0;                    // exponent
            if (len > in.length) // protect against huge length
            {
                throw new NumberFormatException();
            }
            char coeff[] = new char[len];    // integer significand array
            char c;                          // work

            for (; len > 0; offset++, len--) {
                c = in[offset];
                if ((c >= '0' && c <= '9') || Character.isDigit(c)) {
                    // have digit
                    coeff[precision] = c;
                    precision++;             // count of digits
                    continue;
                }
                if (c == '.') {
                    // have dot
                    if (dotoff >= 0) // two dots
                    {
                        throw new NumberFormatException();
                    }
                    dotoff = offset;
                    continue;
                }
                // exponent expected
                if ((c != 'e') && (c != 'E')) {
                    throw new NumberFormatException();
                }
                offset++;
                c = in[offset];
                len--;
                boolean negexp = false;
                // optional sign
                if (c == '-' || c == '+') {
                    negexp = (c == '-');
                    offset++;
                    c = in[offset];
                    len--;
                }
                if (len <= 0 || len > 10) // no digits, or too long
                {
                    throw new NumberFormatException();
                }
                // c now holds first digit of exponent
                for (;; len--) {
                    int v;
                    if (c >= '0' && c <= '9') {
                        v = c - '0';
                    } else {
                        v = Character.digit(c, 10);
                        if (v < 0) // not a digit
                        {
                            throw new NumberFormatException();
                        }
                    }
                    exp = exp * 10 + v;
                    if (len == 1) {
                        break;               // that was final character
                    }
                    offset++;
                    c = in[offset];
                }
                if (negexp) // apply sign
                {
                    exp = -exp;
                }
                // Next test is required for backwards compatibility
                if ((int) exp != exp) // overflow
                {
                    throw new NumberFormatException();
                }
                break;                       // [saves a test]
            }
            // here when no characters left
            if (precision == 0) // no digits found
            {
                throw new NumberFormatException();
            }

            if (dotoff >= 0) {               // had dot; set scale
                scale = precision - (dotoff - cfirst);
                // [cannot overflow]
            }
            if (exp != 0) {                  // had significant exponent
                try {
                    scale = checkScale(-exp + scale); // adjust
                } catch (ArithmeticException e) {
                    throw new NumberFormatException("Scale out of range.");
                }
            }
            // Remove leading zeros from precision (digits count)
            int first = 0;
            for (; coeff[first] == '0' && precision > 1; first++) {
                precision--;
            }

            // Set the significand ..
            // Copy significand to exact-sized array, with sign if
            // negative
            // Later use: BigInteger(coeff, first, precision) for
            //   both cases, by allowing an extra char at the front of
            //   coeff.
            char quick[];
            if (!isneg) {
                quick = new char[precision];
                System.arraycopy(coeff, first, quick, 0, precision);
            } else {
                quick = new char[precision + 1];
                quick[0] = '-';
                System.arraycopy(coeff, first, quick, 1, precision);
            }
            intVal = new NumberZ(quick);
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new NumberFormatException();
        } catch (NegativeArraySizeException e) {
            throw new NumberFormatException();
        }
    }

    /**
     * Translates a character array representation of a <tt>BigDecimal</tt> into
     * a <tt>BigDecimal</tt>, accepting the same sequence of characters as the
     * {@link #BigDecimal(String)} constructor, while allowing a sub-array to be
     * specified and with rounding according to the context settings.
     *
     * <p>
     * Note that if the sequence of characters is already available within a
     * character array, using this constructor is faster than converting the
     * <tt>char</tt> array to string and using the <tt>BigDecimal(String)</tt>
     * constructor .
     *
     * @param in <tt>char</tt> array that is the source of characters.
     * @param offset first character in the array to inspect.
     * @param len number of characters to consider..
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @throws NumberFormatException if <tt>in</tt> is not a valid
     * representation of a <tt>BigDecimal</tt> or the defined subarray is not
     * wholly within <tt>in</tt>.
     * @since 1.5
     */
    public NumberR(char[] in, int offset, int len, MathContext mc) {
        this(in, offset, len);
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates a character array representation of a <tt>BigDecimal</tt> into
     * a <tt>BigDecimal</tt>, accepting the same sequence of characters as the
     * {@link #BigDecimal(String)} constructor.
     *
     * <p>
     * Note that if the sequence of characters is already available as a
     * character array, using this constructor is faster than converting the
     * <tt>char</tt> array to string and using the <tt>BigDecimal(String)</tt>
     * constructor .
     *
     * @param in <tt>char</tt> array that is the source of characters.
     *
     * @throws NumberFormatException if <tt>in</tt> is not a valid
     * representation of a <tt>BigDecimal</tt>.
     * @since 1.5
     */
    public NumberR(char[] in) {
        this(in, 0, in.length);
    }

    /**
     * Translates a character array representation of a <tt>BigDecimal</tt> into
     * a <tt>BigDecimal</tt>, accepting the same sequence of characters as the
     * {@link #BigDecimal(String)} constructor and with rounding according to
     * the context settings.
     *
     * <p>
     * Note that if the sequence of characters is already available as a
     * character array, using this constructor is faster than converting the
     * <tt>char</tt> array to string and using the <tt>BigDecimal(String)</tt>
     * constructor .
     *
     * @param in <tt>char</tt> array that is the source of characters.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @throws NumberFormatException if <tt>in</tt> is not a valid
     * representation of a <tt>BigDecimal</tt>.
     * @since 1.5
     */
    public NumberR(char[] in, MathContext mc) {
        this(in, 0, in.length, mc);
    }

    /**
     * Translates the string representation of a <tt>BigDecimal</tt> into a
     * <tt>BigDecimal</tt>. The string representation consists of an optional
     * sign, <tt>'+'</tt> (<tt>'&#92;u002B'</tt>) or <tt>'-'</tt>
     * (<tt>'&#92;u002D'</tt>), followed by a sequence of zero or more decimal
     * digits ("the integer"), optionally followed by a fraction, optionally
     * followed by an exponent.
     *
     * <p>
     * The fraction consists of a decimal point followed by zero or more decimal
     * digits. The string must contain at least one digit in either the integer
     * or the fraction. The number formed by the sign, the integer and the
     * fraction is referred to as the <i>significand</i>.
     *
     * <p>
     * The exponent consists of the character <tt>'e'</tt>
     * (<tt>'&#92;u0075'</tt>) or <tt>'E'</tt> (<tt>'&#92;u0045'</tt>) followed
     * by one or more decimal digits. The value of the exponent must lie between -{@link Integer#MAX_VALUE} ({@link
     * Integer#MIN_VALUE}+1) and {@link Integer#MAX_VALUE}, inclusive.
     *
     * <p>
     * More formally, the strings this constructor accepts are described by the
     * following grammar: <blockquote> <dl> <dt><i>BigDecimalString:</i>
     * <dd><i>Sign<sub>opt</sub> Significand Exponent<sub>opt</sub></i>
     * <p>
     * <dt><i>Sign:</i> <dd><tt>+</tt> <dd><tt>-</tt>
     * <p>
     * <dt><i>Significand:</i> <dd><i>IntegerPart</i> <tt>.</tt>
     * <i>FractionPart<sub>opt</sub></i> <dd><tt>.</tt> <i>FractionPart</i>
     * <dd><i>IntegerPart</i>
     * <p>
     * <dt><i>IntegerPart: <dd>Digits</i>
     * <p>
     * <dt><i>FractionPart: <dd>Digits</i>
     * <p>
     * <dt><i>Exponent:
     * <dd>ExponentIndicator SignedInteger</i>
     * <p>
     * <dt><i>ExponentIndicator:</i>
     * <dd><tt>e</tt> <dd><tt>E</tt>
     * <p>
     * <dt><i>SignedInteger:
     * <dd>Sign<sub>opt</sub> Digits</i>
     * <p>
     * <dt><i>Digits: <dd>Digit <dd>Digits Digit</i>
     * <p>
     * <dt><i>Digit:</i> <dd>any character for which {@link Character#isDigit}
     * returns <tt>true</tt>, including 0, 1, 2 ...
     * </dl> </blockquote>
     *
     * <p>
     * The scale of the returned <tt>BigDecimal</tt> will be the number of
     * digits in the fraction, or zero if the string contains no decimal point,
     * subject to adjustment for any exponent; if the string contains an
     * exponent, the exponent is subtracted from the scale. The value of the
     * resulting scale must lie between <tt>Integer.MIN_VALUE</tt> and
     * <tt>Integer.MAX_VALUE</tt>, inclusive.
     *
     * <p>
     * The character-to-digit mapping is provided by {@link
     * java.lang.Character#digit} set to convert to radix 10. The String may not
     * contain any extraneous characters (whitespace, for example).
     *
     * <p>
     * <b>Examples:</b><br> The value of the returned <tt>BigDecimal</tt> is
     * equal to <i>significand</i> &times; 10<sup>&nbsp;<i>exponent</i></sup>.
     * For each string on the left, the resulting representation
     * [<tt>BigInteger</tt>, <tt>scale</tt>] is shown on the right.
     * <pre>
     * "0"            [0,0]
     * "0.00"         [0,2]
     * "123"          [123,0]
     * "-123"         [-123,0]
     * "1.23E3"       [123,-1]
     * "1.23E+3"      [123,-1]
     * "12.3E+7"      [123,-6]
     * "12.0"         [120,1]
     * "12.3"         [123,1]
     * "0.00123"      [123,5]
     * "-1.23E-12"    [-123,14]
     * "1234.5E-4"    [12345,5]
     * "0E+7"         [0,-7]
     * "-0"           [0,0]
     * </pre>
     *
     * <p>
     * Note: For values other than <tt>float</tt> and <tt>double</tt>
     * NaN_DOUBLE and &plusmn;Infinity, this constructor is compatible with the
     * values returned by {@link Float#toString} and
     * {@link com.mathpar.number.Double#toString}. This is generally the preferred way to
     * convert a <tt>float</tt> or <tt>double</tt> into a BigDecimal, as it
     * doesn't suffer from the unpredictability of the
     * {@link #BigDecimal(double)} constructor.
     *
     * @param val String representation of <tt>BigDecimal</tt>.
     *
     * @throws NumberFormatException if <tt>val</tt> is not a valid
     * representation of a <tt>BigDecimal</tt>.
     */
    public NumberR(String val) {
        this(val.toCharArray(), 0, val.length());
    }

    /**
     * Translates the string representation of a <tt>BigDecimal</tt> into a
     * <tt>BigDecimal</tt>, accepting the same strings as the
     * {@link #BigDecimal(String)} constructor, with rounding according to the
     * context settings.
     *
     * @param val string representation of a <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @throws NumberFormatException if <tt>val</tt> is not a valid
     * representation of a BigDecimal.
     * @since 1.5
     */
    public NumberR(String val, MathContext mc) {
        this(val.toCharArray(), 0, val.length());
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates a <tt>double</tt> into a <tt>BigDecimal</tt> which is the
     * exact decimal representation of the <tt>double</tt>'s binary
     * floating-point value. The scale of the returned <tt>BigDecimal</tt> is
     * the smallest value such that <tt>(10<sup>scale</sup> &times; val)</tt> is
     * an integer.
     * <p>
     * <b>Notes:</b> <ol> <li> The results of this constructor can be somewhat
     * unpredictable. One might assume that writing <tt>new BigDecimal(0.1)</tt>
     * in Java creates a <tt>BigDecimal</tt> which is exactly equal to 0.1 (an
     * unscaled value of 1, with a scale of 1), but it is actually equal to
     * 0.1000000000000000055511151231257827021181583404541015625. This is
     * because 0.1 cannot be represented exactly as a <tt>double</tt> (or, for
     * that matter, as a binary fraction of any finite length). Thus, the value
     * that is being passed <i>in</i> to the constructor is not exactly equal to
     * 0.1, appearances notwithstanding.
     *
     * <li> The <tt>String</tt> constructor, on the other hand, is perfectly
     * predictable: writing <tt>new BigDecimal("0.1")</tt> creates a
     * <tt>BigDecimal</tt> which is <i>exactly</i> equal to 0.1, as one would
     * expect. Therefore, it is generally recommended that the {@linkplain #BigDecimal(String)
     * <tt>String</tt> constructor} be used in preference to this one.
     *
     * <li> When a <tt>double</tt> must be used as a source for a
     * <tt>BigDecimal</tt>, note that this constructor provides an exact
     * conversion; it does not give the same result as converting the
     * <tt>double</tt> to a <tt>String</tt> using the
     * {@link com.mathpar.number.Double#toString(double)} method and then using the
     * {@link #BigDecimal(String)} constructor. To get that result, use the
     * <tt>static</tt> {@link #valueOf(double)} method. </ol>
     *
     * @param val <tt>double</tt> value to be converted to <tt>BigDecimal</tt>.
     *
     * @throws NumberFormatException if <tt>val</tt> is infinite or NaN_DOUBLE.
     */
    public NumberR(double val) {
        intVal = NumberZ.ZERO;
        if (Double.NaN == (val)) {
            scale = Integer.MAX_VALUE;
            intVal = NumberZ.ZERO;
            return;
        }
        if (Double.POSITIVE_INFINITY == (val)) {
            scale = Integer.MAX_VALUE - 1;
            return;
        }
        if (Double.NEGATIVE_INFINITY == (val)) {
            scale = Integer.MAX_VALUE - 2;
            return;
        }
        if (NumberR64.POZITIVE_ZERO_DOUBLE == (val)) {
            scale = Integer.MAX_VALUE - 3;
            return;
        }
        if (NumberR64.NEGATIVE_ZERO_DOUBLE == (val)) {
            scale = Integer.MAX_VALUE - 4;
            return;
        }

        // Translate the double into sign, exponent and mantissa, according
        // to the formulae in JLS, Section 20.10.22.
        long valBits = NumberR64.doubleToLongBits(val);
        int sign = ((valBits >> 63) == 0 ? 1 : -1);
        int exponent = (int) ((valBits >> 52) & 0x7ffL);
        long mantissa = (exponent == 0 ? (valBits & ((1L << 52) - 1)) << 1
                : (valBits & ((1L << 52) - 1)) | (1L << 52));
        exponent -= 1075;
        /*
         * At this point, val == sign * mantissa * 2**exponent
         */

        /*
         * Special case zero to supress nonterminating normalization and bogus
         * scale calculation.
         */
        if (mantissa == 0) {
            intVal = NumberZ.ZERO;
            precision = 1;
            return;
        }

        /*
         * Normalize
         */
        while ((mantissa & 1) == 0) {    /*
             * i.e., Mantissa is even
             */

            mantissa >>= 1;
            exponent++;
        }

        /*
         * Calculate intVal and scale
         */
        intVal = NumberZ.valueOf(sign * mantissa);
        if (exponent < 0) {
            intVal = intVal.multiply(NumberZ.valueOf(5).pow(-exponent));
            scale = -exponent;
        } else if (exponent > 0) {
            intVal = intVal.multiply(NumberZ.valueOf(2).pow(exponent));
        }
    }

    /**
     * Translates a <tt>double</tt> into a <tt>BigDecimal</tt>, with rounding
     * according to the context settings. The scale of the <tt>BigDecimal</tt>
     * is the smallest value such that <tt>(10<sup>scale</sup> &times; val)</tt>
     * is an integer.
     *
     * <p>
     * The results of this constructor can be somewhat unpredictable and its use
     * is generally not recommended; see the notes under the
     * {@link #BigDecimal(double)} constructor.
     *
     * @param val <tt>double</tt> value to be converted to <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the RoundingMode
     * is UNNECESSARY.
     * @throws NumberFormatException if <tt>val</tt> is infinite or NaN_DOUBLE.
     * @since 1.5
     */
    public NumberR(double val, MathContext mc) {
        this(val);
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates a <tt>BigInteger</tt> into a <tt>BigDecimal</tt>. The scale of
     * the <tt>BigDecimal</tt> is zero.
     *
     * @param val <tt>BigInteger</tt> value to be converted to
     * <tt>BigDecimal</tt>.
     */
    public NumberR(NumberZ val) {
        intVal = val;
    }

    /**
     * Translates a <tt>BigInteger</tt> into a <tt>BigDecimal</tt> rounding
     * according to the context settings. The scale of the <tt>BigDecimal</tt>
     * is zero.
     *
     * @param val <tt>BigInteger</tt> value to be converted to
     * <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public NumberR(NumberZ val, MathContext mc) {
        intVal = val;
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates a <tt>BigInteger</tt> unscaled value and an <tt>int</tt> scale
     * into a <tt>BigDecimal</tt>. The value of the <tt>BigDecimal</tt> is
     * <tt>(unscaledVal &times; 10<sup>-scale</sup>)</tt>.
     *
     * @param unscaledVal unscaled value of the <tt>Element</tt>.
     * @param scale scale of the <tt>BigDecimal</tt>.
     */
    public NumberR(Element unscaledVal, int scale, Ring ring) {
        // Negative scales are now allowed
        intVal = (NumberZ)(unscaledVal.toNumber(Ring.Z,ring));
        this.scale = scale;
    }

    /**
     * Translates a <tt>BigInteger</tt> unscaled value and an <tt>int</tt> scale
     * into a <tt>BigDecimal</tt>. The value of the <tt>BigDecimal</tt> is
     * <tt>(unscaledVal &times; 10<sup>-scale</sup>)</tt>.
     *
     * @param unscaledVal unscaled value of the <tt>BigDecimal</tt>.
     * @param scale scale of the <tt>BigDecimal</tt>.
     */
    public NumberR(NumberZ unscaledVal, int scale) {
        // Negative scales are now allowed
        intVal = unscaledVal;
        this.scale = scale;
    }

    /**
     * Translates a <tt>BigInteger</tt> unscaled value and an <tt>int</tt> scale
     * into a <tt>BigDecimal</tt>, with rounding according to the context
     * settings. The value of the <tt>BigDecimal</tt> is <tt>(unscaledVal
     * &times; 10<sup>-scale</sup>)</tt>, rounded according to the
     * <tt>precision</tt> and rounding mode settings.
     *
     * @param unscaledVal unscaled value of the <tt>BigDecimal</tt>.
     * @param scale scale of the <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public NumberR(NumberZ unscaledVal, int scale, MathContext mc) {
        intVal = unscaledVal;
        this.scale = scale;
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates an <tt>int</tt> into a <tt>BigDecimal</tt>. The scale of the
     * <tt>BigDecimal</tt> is zero.
     *
     * @param val <tt>int</tt> value to be converted to <tt>BigDecimal</tt>.
     *
     * @since 1.5
     */
    public NumberR(int val) {
        scale = 0;
        intVal = NumberZ.valueOf(val);
    }

    /**
     * Translates an <tt>int</tt> into a <tt>BigDecimal</tt>, with rounding
     * according to the context settings. The scale of the <tt>BigDecimal</tt>,
     * before any rounding, is zero.
     *
     * @param val <tt>int</tt> value to be converted to <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public NumberR(int val, MathContext mc) {
        intVal = NumberZ.valueOf(val);
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    /**
     * Translates a <tt>long</tt> into a <tt>BigDecimal</tt>. The scale of the
     * <tt>BigDecimal</tt> is zero.
     *
     * @param val <tt>long</tt> value to be converted to <tt>BigDecimal</tt>.
     *
     * @since 1.5
     */
    public NumberR(long val) {
        intVal = NumberZ.valueOf(val);
    }

    /**
     * Translates a <tt>long</tt> into a <tt>BigDecimal</tt>, with rounding
     * according to the context settings. The scale of the <tt>BigDecimal</tt>,
     * before any rounding, is zero.
     *
     * @param val <tt>long</tt> value to be converted to <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public NumberR(long val, MathContext mc) {
        intVal = NumberZ.valueOf(val);
        if (mc.getPrecision() > 0) {
            roundThis(mc);
        }
    }

    // Static Factory Methods
    /**
     * Translates a <tt>long</tt> unscaled value and an <tt>int</tt> scale into
     * a <tt>BigDecimal</tt>. This &quot;static factory method&quot; is provided
     * in preference to a (<tt>long</tt>, <tt>int</tt>) constructor because it
     * allows for reuse of frequently used <tt>BigDecimal</tt> values..
     *
     * @param unscaledVal unscaled value of the <tt>BigDecimal</tt>.
     * @param scale scale of the <tt>BigDecimal</tt>.
     *
     * @return a <tt>BigDecimal</tt> whose value is <tt>(unscaledVal &times;
     * 10<sup>-scale</sup>)</tt>.
     */
    public static NumberR valueOf(long unscaledVal, int scale) {
        if (scale == 0) {
            if (unscaledVal == 0) {
                return ZERO;
            }
            if (unscaledVal == 1) {
                return ONE;
            }
            if (unscaledVal == 10) {
                return POSCONST[10];
            }
        }
        return new NumberR(NumberZ.valueOf(unscaledVal), scale);
    }

    /**
     * Translates a <tt>long</tt> value into a <tt>BigDecimal</tt> with a scale
     * of zero. This &quot;static factory method&quot; is provided in preference
     * to a (<tt>long</tt>) constructor because it allows for reuse of
     * frequently used <tt>BigDecimal</tt> values.
     *
     * @param val value of the <tt>BigDecimal</tt>.
     *
     * @return a <tt>BigDecimal</tt> whose value is <tt>val</tt>.
     */
    public static NumberR valueOf(long val) {
        return valueOf(val, 0);
    }

    /**
     * Translates a <tt>double</tt> into a <tt>BigDecimal</tt>, using the
     * <tt>double</tt>'s canonical string representation provided by the
     * {@link com.mathpar.number.Double#toString(double)} method.
     *
     * <p>
     * <b>Note:</b> This is generally the preferred way to convert a
     * <tt>double</tt> (or <tt>float</tt>) into a <tt>BigDecimal</tt>, as the
     * value returned is equal to that resulting from constructing a
     * <tt>BigDecimal</tt> from the result of using
     * {@link com.mathpar.number.Double#toString(double)}.
     *
     * @param val <tt>double</tt> to convert to a <tt>BigDecimal</tt>.
     *
     * @return a <tt>BigDecimal</tt> whose value is equal to or approximately
     * equal to the value of <tt>val</tt>.
     *
     * @throws NumberFormatException if <tt>val</tt> is infinite or NaN_DOUBLE.
     * @since 1.5
     */
    public static Element valueOf(double val) {
        if (((Double) val).equals(Double.NaN)) {
            return NAN;
        }
        if (((Double) val).equals(Double.POSITIVE_INFINITY)) {
            return POSITIVE_INFINITY;
        }
        if (((Double) val).equals(Double.NEGATIVE_INFINITY)) {
            return NEGATIVE_INFINITY;
        }

        if (((Double) val).equals(NumberR64.NEGATIVE_ZERO_DOUBLE)) {
            return NEGATIVE_ZERO;
        }

        if (((Double) val).equals(NumberR64.POZITIVE_ZERO_DOUBLE)) {
            return POZITIVE_ZERO;
        }

        // Reminder: a zero double returns '0.0', so we cannot fastpath
        // to use the constant ZERO.  This might be important enough to
        // justify a factory approach, a cache, or a few private
        // constants, later.
        return new NumberR(NumberR64.toString(val));
    }

    // Arithmetic Operations
    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this + augend)</tt>,
     * and whose scale is <tt>max(this.scale(), augend.scale())</tt>.
     *
     * @param augend value to be added to this <tt>BigDecimal</tt>.
     *
     * @return <tt>this + augend</tt>
     */
    public Element add(NumberR augend) {

        NumberR arg[] = new NumberR[2];
        arg[0] = this;
        arg[1] = augend;
        matchScale(arg);
        return new NumberR(arg[0].intVal.add(arg[1].intVal), arg[0].scale);

//        return ZERO;
    }

    @Override
    public boolean isInfinite() {
        return (this == POSITIVE_INFINITY) || (this == NEGATIVE_INFINITY);
    }

    @Override
    public boolean isNaN() {
        return (this == NAN);

    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this + augend)</tt>,
     * with rounding according to the context settings.
     *
     * If either number is zero and the precision setting is nonzero then the
     * other number, rounded if necessary, is used as the result.
     *
     * @param augend value to be added to this <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @return <tt>this + augend</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public Element add(NumberR augend, MathContext mc) {

        if (mc.getPrecision() == 0) {
            return add(augend);
        }
        NumberR lhs = this;

        // If either number is zero then the other number, rounded and
        // scaled if necessary, is used as the result.
        {
            boolean lhsIsZero = lhs.signum() == 0;
            boolean augendIsZero = augend.signum() == 0;

            if (lhsIsZero || augendIsZero) {
                int preferredScale = Math.max(lhs.scale(), augend.scale());
                NumberR result;

                if (lhsIsZero && augendIsZero) {
                    return new NumberR(NumberZ.ZERO, preferredScale);
                }

                result = lhsIsZero ? augend.doRound(mc) : lhs.doRound(mc);

                if (result.scale() == preferredScale) {
                    return result;
                } else if (result.scale() > preferredScale) {
                    return result.stripZerosToMatchScale(preferredScale);
                } else { // result.scale < preferredScale
                    int precisionDiff = mc.getPrecision() - result.precision();
                    int scaleDiff = preferredScale - result.scale();

                    if (precisionDiff >= scaleDiff) {
                        return result.setScale(preferredScale); // can achieve target scale
                    } else {
                        return result.setScale(result.scale() + precisionDiff);
                    }
                }
            }
        }

        int padding = checkScale((long) lhs.scale - augend.scale);
        if (padding != 0) {        // scales differ; alignment needed
            // if one operand is < 0.01 ulp of the other at full
            // precision, replace it by a 'sticky bit' of +0.001/-0.001 ulp.
            // [In a sense this is an 'optimization', but it also makes
            // a much wider range of additions practical.]
            if (padding < 0) {     // lhs will be padded
                int ulpscale = lhs.scale - lhs.precision() + mc.getPrecision();
                if (augend.scale - augend.precision() > ulpscale + 1) {
                    augend = NumberR.valueOf(augend.signum(), ulpscale + 3);
                }
            } else {               // rhs (augend) will be padded
                int ulpscale = augend.scale - augend.precision + mc.getPrecision();
                if (lhs.scale - lhs.precision() > ulpscale + 1) {
                    lhs = NumberR.valueOf(lhs.signum(), ulpscale + 3);
                }
            }
            NumberR arg[] = new NumberR[2];
            arg[0] = lhs;
            arg[1] = augend;
            matchScale(arg);
            lhs = arg[0];
            augend = arg[1];
        }

        return new NumberR(lhs.intVal.add(augend.intVal), lhs.scale).doRound(mc);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this -
     * subtrahend)</tt>, and whose scale is <tt>max(this.scale(),
     * subtrahend.scale())</tt>.
     *
     * @param subtrahend value to be subtracted from this <tt>BigDecimal</tt>.
     *
     * @return <tt>this - subtrahend</tt>
     */
    public Element subtract(NumberR subtrahend) {

        NumberR arg[] = new NumberR[2];
        arg[0] = this;
        arg[1] = subtrahend;
        matchScale(arg);
        return new NumberR(arg[0].intVal.subtract(arg[1].intVal),
                arg[0].scale);

    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this -
     * subtrahend)</tt>, with rounding according to the context settings.
     *
     * If <tt>subtrahend</tt> is zero then this, rounded if necessary, is used
     * as the result. If this is zero then the result is
     * <tt>subtrahend.negate(mc)</tt>.
     *
     * @param subtrahend value to be subtracted from this <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @return <tt>this - subtrahend</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public Element subtract(NumberR subtrahend, MathContext mc) {

        if (mc.getPrecision() == 0) {
            return subtract(subtrahend);
        }
        // share the special rounding code in add()
        NumberR rhs = new NumberR(subtrahend.intVal.negate(), subtrahend.scale);
        rhs.precision = subtrahend.precision;
        return add(rhs, mc);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this &times;
     * multiplicand)</tt>, and whose scale is <tt>(this.scale() +
     * multiplicand.scale())</tt>.
     *
     * @param multiplicand value to be multiplied by this <tt>BigDecimal</tt>.
     *
     * @return <tt>this * multiplicand</tt>
     */
    public Element multiply(NumberR multiplicand) {
        NumberR result = new NumberR(intVal.multiply(multiplicand.intVal), 0);
        result.scale = checkScale((long) scale + multiplicand.scale);
        return result;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this &times;
     * multiplicand)</tt>, with rounding according to the context settings.
     *
     * @param multiplicand value to be multiplied by this <tt>BigDecimal</tt>.
     * @param mc the context to use.
     *
     * @return <tt>this * multiplicand</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public Element multiply(
            NumberR multiplicand, MathContext mc) {

        if (mc.getPrecision() == 0) {
            return multiply(multiplicand);
        }

        return ((NumberR) multiply(multiplicand)).doRound(mc);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * and whose scale is as specified. If rounding must be performed to
     * generate a result with the specified scale, the specified rounding mode
     * is applied.
     *
     * <p>
     * The new {@link #divide(BigDecimal, int, RoundingMode)} method should be
     * used in preference to this legacy method.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param scale scale of the <tt>BigDecimal</tt> quotient to be returned.
     * @param roundingMode rounding mode to apply.
     *
     * @return <tt>this / divisor</tt>
     *
     * @throws ArithmeticException if <tt>divisor</tt> is zero,
     * <tt>roundingMode==ROUND_UNNECESSARY</tt> and the specified scale is
     * insufficient to represent the result of the division exactly.
     * @throws IllegalArgumentException if <tt>roundingMode</tt> does not
     * represent a valid rounding mode.
     * @see #ROUND_UP
     * @see #ROUND_DOWN
     * @see #ROUND_CEILING
     * @see #ROUND_FLOOR
     * @see #ROUND_HALF_UP
     * @see #ROUND_HALF_DOWN
     * @see #ROUND_HALF_EVEN
     * @see #ROUND_UNNECESSARY
     */
    public NumberR divide(
            NumberR divisor, int scale, int roundingMode) {
        if (roundingMode < ROUND_UP || roundingMode > ROUND_UNNECESSARY) {
            throw new IllegalArgumentException("Invalid rounding mode");
        }

        /*
         * Rescale dividend or divisor (whichever can be "upscaled" to produce
         * correctly scaled quotient). Take care to detect out-of-range scales
         */
        NumberR dividend;
        if (checkScale((long) scale + divisor.scale) >= this.scale) {
            dividend = this.setScale(scale + divisor.scale);
        } else {
            dividend = this;
            divisor
                    = divisor.setScale(checkScale((long) this.scale - scale));
        }

        /*
         * Do the division and return result if it's exact
         */
        NumberZ i[] = dividend.intVal.divideAndRemainder(divisor.intVal);
        NumberZ q = i[0], r = i[1];
        if (r.signum() == 0) {
            return new NumberR(q, scale);
        }

        if (roundingMode == ROUND_UNNECESSARY) /*
         * Rounding prohibited
         */ {
            throw new ArithmeticException("Rounding necessary");
        }

        /*
         * Round as appropriate
         */
        int signum = dividend.signum() * divisor.signum(); /*
         * Sign of result
         */

        boolean increment;
        if (roundingMode == ROUND_UP) {             /*
             * Away from zero
             */

            increment = true;
        } else if (roundingMode == ROUND_DOWN) {    /*
             * Towards zero
             */

            increment = false;
        } else if (roundingMode == ROUND_CEILING) { /*
             * Towards +infinity
             */

            increment = (signum > 0);
        } else if (roundingMode == ROUND_FLOOR) {   /*
             * Towards -infinity
             */

            increment = (signum < 0);
        } else { /*
             * Remaining modes based on nearest-neighbor determination
             */

            // add(r) here is faster than multiply(2) or shiftLeft(1)
            int cmpFracHalf = r.add(r).abs().compareTo(divisor.intVal.abs());
            if (cmpFracHalf < 0) {         /*
                 * We're closer to higher digit
                 */

                increment = false;
            } else if (cmpFracHalf > 0) {  /*
                 * We're closer to lower digit
                 */

                increment = true;
            } else {                       /*
                 * We're dead-center
                 */

                if (roundingMode == ROUND_HALF_UP) {
                    increment = true;
                } else if (roundingMode == ROUND_HALF_DOWN) {
                    increment = false;
                } else /*
                 * roundingMode == ROUND_HALF_EVEN
                 */ {
                    increment = q.testBit(0);   /*
                     * true iff q is odd
                     */

                }

            }
        }
        return (increment
                ? new NumberR(q.add(NumberZ.valueOf(signum)), scale)
                : new NumberR(q, scale));
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * and whose scale is as specified. If rounding must be performed to
     * generate a result with the specified scale, the specified rounding mode
     * is applied.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param scale scale of the <tt>BigDecimal</tt> quotient to be returned.
     * @param roundingMode rounding mode to apply.
     *
     * @return <tt>this / divisor</tt>
     *
     * @throws ArithmeticException if <tt>divisor</tt> is zero,
     * <tt>roundingMode==RoundingMode.UNNECESSARY</tt> and the specified scale
     * is insufficient to represent the result of the division exactly.
     * @since 1.5
     */
    public NumberR divide(
            NumberR divisor, int scale, RoundingMode roundingMode) {
        return divide(divisor, scale, roundingMode.oldMode);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * and whose scale is <tt>this.scale()</tt>. If rounding must be performed
     * to generate a result with the given scale, the specified rounding mode is
     * applied.
     *
     * <p>
     * The new {@link #divide(BigDecimal, RoundingMode)} method should be used
     * in preference to this legacy method.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param roundingMode rounding mode to apply.
     *
     * @return <tt>this / divisor</tt>
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>, or
     * <tt>roundingMode==ROUND_UNNECESSARY</tt> and <tt>this.scale()</tt> is
     * insufficient to represent the result of the division exactly.
     * @throws IllegalArgumentException if <tt>roundingMode</tt> does not
     * represent a valid rounding mode.
     * @see #ROUND_UP
     * @see #ROUND_DOWN
     * @see #ROUND_CEILING
     * @see #ROUND_FLOOR
     * @see #ROUND_HALF_UP
     * @see #ROUND_HALF_DOWN
     * @see #ROUND_HALF_EVEN
     * @see #ROUND_UNNECESSARY
     */
    public NumberR divide(
            NumberR divisor, int roundingMode) {
        return this.divide(divisor, scale, roundingMode);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * and whose scale is <tt>this.scale()</tt>. If rounding must be performed
     * to generate a result with the given scale, the specified rounding mode is
     * applied.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param roundingMode rounding mode to apply.
     *
     * @return <tt>this / divisor</tt>
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>, or
     * <tt>roundingMode==RoundingMode.UNNECESSARY</tt> and <tt>this.scale()</tt>
     * is insufficient to represent the result of the division exactly.
     */
    public NumberR divide(
            NumberR divisor, RoundingMode roundingMode) {
        return this.divide(divisor, scale, roundingMode.oldMode);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * and whose preferred scale is <tt>(this.scale() - divisor.scale())</tt>;
     * if the exact quotient cannot be represented (because it has a
     * non-terminating decimal expansion) an <tt>ArithmeticException</tt> is
     * thrown.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     *
     * @throws ArithmeticException if the exact quotient does not have a
     * terminating decimal expansion
     * @return <tt>this / divisor</tt>
     *
     * @since 1.5
     * @author Joseph D. Darcy
     */
    public Element divideExact(NumberR divisor) {
        /*
         * Handle zero cases first.
         */

        if (divisor.intVal.signum() == 0) {   // x/0
            if (this.intVal.signum() == 0) // 0/0
            {
                return NAN;
            }

            if (this.signum() == -1) {
                return NEGATIVE_INFINITY;
            } else {
                return POSITIVE_INFINITY; // ?????????????????????? ???????? ??????????????????
            }
        }

        // Calculate preferred scale
        int preferredScale = (int) Math.max(Math.min((long) this.scale() - divisor.
                scale(),
                Integer.MAX_VALUE), Integer.MIN_VALUE);
        if (this.intVal.signum() == 0) // 0/y
        {
            return NumberR.valueOf(0, preferredScale);
        } else {
            /*
             * If the quotient this/divisor has a terminating decimal expansion,
             * the expansion can have no more than (a.precision() +
             * ceil(10*b.precision)/3) digits. Therefore, create a MathContext
             * object with this precision and do a divide with the UNNECESSARY
             * rounding mode.
             */
            MathContext mc = new MathContext((int) Math.min(this.precision()
                    + (long) Math.ceil(10.0 * divisor.precision() / 3.0),
                    Integer.MAX_VALUE),
                    RoundingMode.UNNECESSARY);
            NumberR quotient;

            try {
                Element quotientE = this.divide(divisor, mc);
                if (!(quotientE instanceof NumberR)) {
                    return quotientE;
                }
                quotient = (NumberR) quotientE;
            } catch (ArithmeticException e) {
                throw new ArithmeticException("Non-terminating decimal expansion; "
                        + "no exact representable decimal result.");
            }

            int quotientScale = quotient.scale();

            // divide(BigDecimal, mc) tries to adjust the quotient to
            // the desired one by removing trailing zeros; since the
            // exact divide method does not have an explicit digit
            // limit, we can add zeros too.
            if (preferredScale > quotientScale) {
                return quotient.setScale(preferredScale);
            }
            return quotient;
        }

    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this / divisor)</tt>,
     * with rounding according to the context settings.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param mc the context to use.
     *
     * @return <tt>this / divisor</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt> or <tt>mc.precision == 0</tt> and the
     * quotient has a non-terminating decimal expansion.
     * @since 1.5
     */
    public Element divide(
            NumberR divisor, MathContext mc) {
        if (mc.getPrecision() == 0) {
            return divideExact(divisor);
        }

        NumberR lhs = this;     // left-hand-side
        NumberR rhs = divisor;  // right-hand-side
        NumberR result;         // work

        long preferredScale = (long) lhs.scale() - rhs.scale();

        // Now calculate the answer.  We use the existing
        // divide-and-round method, but as this rounds to scale we have
        // to normalize the values here to achieve the desired result.
        // For x/y we first handle y=0 and x=0, and then normalize x and
        // y to give x' and y' with the following constraints:
        //   (a) 0.1 <= x' < 1
        //   (b)  x' <= y' < 10*x'
        // Dividing x'/y' with the required scale set to mc.precision then
        // will give a result in the range 0.1 to 1 rounded to exactly
        // the right number of digits (except in the case of a result of
        // 1.000... which can arise when x=y, or when rounding overflows
        // The 1.000... case will reduce properly to 1.
//       int Rscale=rhs.scale;
//       int Lscale=lhs.scale;
//
//       if(rhs.intVal.signum()==0){ // x/0 ; x/\infty; x/-\infty ; x/Nan; x/+0; x/-0
//        // ?????????????? ???????? ?????? ???????????????? +0 / -0
//
//
//        if(Rscale==POSITIVE_INFINITY.scale | Rscale==NEGATIVE_INFINITY.scale){
//           return (Lscale==POSITIVE_INFINITY.scale | Lscale==NEGATIVE_INFINITY.scale) ? NAN:ZERO;
//         }
//        if(Rscale==NAN.scale) return NAN;
//        // ?????? ?????????????????? ?????????????? ???? 0
//        int fl=lhs.signum();
//        return (fl==0)?NAN:(fl<0)?NEGATIVE_INFINITY:POSITIVE_INFINITY;
//       }
//
//       if(lhs.intVal.signum()==0){ // 0/y ; \infty/y ;-\infty/y ; NAN/y
//        if(Lscale==NAN.scale) return NAN;
//        int fl=lhs.signum();
//        if(Lscale==POSITIVE_INFINITY.scale) return (fl>0)? POSITIVE_INFINITY : NEGATIVE_INFINITY;
//        if(Lscale==NEGATIVE_INFINITY.scale) return (fl<0)? POSITIVE_INFINITY : NEGATIVE_INFINITY;
//       }
        if (rhs.intVal.signum() == 0) {      // x/0
            if (lhs.intVal.signum() == 0) // 0/0
            {
                return NumberR.NAN;
                //   throw new ArithmeticException("Division undefined");  // NaN
            }
            return (lhs.signum() > 0) ? NumberR.POSITIVE_INFINITY : NumberR.NEGATIVE_INFINITY;
            //     throw new ArithmeticException("Division by zero");
        }

        if (lhs.intVal.signum() == 0) // 0/y
        {
            return new NumberR(NumberZ.ZERO,
                    (int) Math.max(Math.min(preferredScale,
                                    Integer.MAX_VALUE),
                            Integer.MIN_VALUE));
        }

        NumberR xprime = new NumberR(lhs.intVal.abs(), lhs.precision());
        NumberR yprime = new NumberR(rhs.intVal.abs(), rhs.precision());
        // xprime and yprime are now both in range 0.1 through 0.999...
        if (mc.getRoundingMode() == RoundingMode.CEILING
                || mc.getRoundingMode() == RoundingMode.FLOOR) {
            // The floor (round toward negative infinity) and ceil
            // (round toward positive infinity) rounding modes are not
            // invariant under a sign flip.  If xprime/yprime has a
            // different sign than lhs/rhs, the rounding mode must be
            // changed.
            if ((xprime.signum() != lhs.signum())
                    ^ (yprime.signum() != rhs.signum())) {
                mc = new MathContext(mc.getPrecision(),
                        (mc.getRoundingMode() == RoundingMode.CEILING) ? RoundingMode.FLOOR : RoundingMode.CEILING);
            }

        }

        if (xprime.compareTo(yprime) > 0) // satisfy constraint (b)
        {
            yprime.scale -= 1;                 // [that is, yprime *= 10]
        }

        result = xprime.divide(yprime, mc.getPrecision(), mc.getRoundingMode().oldMode);
        // correct the scale of the result...
        result.scale = checkScale((long) yprime.scale - xprime.scale - (rhs.scale - lhs.scale) + mc.
                getPrecision());
        // apply the sign
        if (lhs.intVal.signum() != rhs.intVal.signum()) {
            result = result.negate();
        }
// doRound, here, only affects 1000000000 case.

        result = result.doRound(mc);

        if (result.multiply(divisor).compareTo(this) == 0) {
            // Apply preferred scale rules for exact quotients
            return result.stripZerosToMatchScale(preferredScale);
        } else {
            return result;
        }

    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is the integer part of the
     * quotient <tt>(this / divisor)</tt> rounded down. The preferred scale of
     * the result is      <code>(this.scale() -
     * divisor.scale())</code>.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     *
     * @return The integer part of <tt>this / divisor</tt>.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @since 1.5
     */
    public Element divideToIntegralValue(
            NumberR divisor) {
        // Calculate preferred scale
        int preferredScale = (int) Math.max(Math.min((long) this.scale() - divisor.
                scale(),
                Integer.MAX_VALUE), Integer.MIN_VALUE);

        if (this.abs().compareTo(divisor.abs()) < 0) {
            // much faster when this << divisor
            return NumberR.valueOf(0, preferredScale);
        }

        if (this.signum() == 0 && divisor.signum() != 0) {
            return this.setScale(preferredScale);
        }

// Perform a divide with enough digits to round to a correct
// integer value; then remove any fractional digits
        int maxDigits = (int) Math.min(this.precision()
                + (long) Math.ceil(10.0 * divisor.precision() / 3.0)
                + Math.abs((long) this.scale() - divisor.scale()) + 2,
                Integer.MAX_VALUE);

        Element quotientE = this.divide(divisor, new MathContext(maxDigits,
                RoundingMode.DOWN));
        if (!(quotientE instanceof NumberR)) {
            return quotientE;
        }
        NumberR quotient = (NumberR) quotientE;
        if (quotient.scale > 0) {
            quotient = quotient.setScale(0, ROUND_DOWN).
                    stripZerosToMatchScale(preferredScale);
        }

        if (quotient.scale < preferredScale) {
            // pad with zeros if necessary
            quotient = quotient.setScale(preferredScale);
        }

        return quotient;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is the integer part of
     * <tt>(this / divisor)</tt>. Since the integer part of the exact quotient
     * does not depend on the rounding mode, the rounding mode does not affect
     * the values returned by this method. The preferred scale of the result is
     * <code>(this.scale() - divisor.scale())</code>. An
     * <tt>ArithmeticException</tt> is thrown if the integer part of the exact
     * quotient needs more than <tt>mc.precision</tt> digits.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param mc the context to use.
     *
     * @return The integer part of <tt>this / divisor</tt>.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @throws ArithmeticException if <tt>mc.precision</tt> &gt; 0 and the
     * result requires a precision of more than <tt>mc.precision</tt> digits.
     * @since 1.5
     * @author Joseph D. Darcy
     */
    public Element divideToIntegralValue(
            NumberR divisor, MathContext mc) {
        if (mc.getPrecision() == 0 || // exact result
                (this.abs().compareTo(divisor.abs()) < 0)) // zero result
        {
            return divideToIntegralValue(divisor);
        }

// Calculate preferred scale
        int preferredScale = (int) Math.max(Math.min((long) this.scale() - divisor.
                scale(),
                Integer.MAX_VALUE), Integer.MIN_VALUE);

        /*
         * Perform a normal divide to mc.precision digits. If the remainder has
         * absolute value less than the divisor, the integer portion of the
         * quotient fits into mc.precision digits. Next, remove any fractional
         * digits from the quotient and adjust the scale to the preferred value.
         */
        Element resultE = this.divide(divisor, new MathContext(mc.getPrecision(),
                RoundingMode.DOWN));
        if (!(resultE instanceof NumberR)) {
            return resultE;
        }
        NumberR result = (NumberR) resultE;
        int resultScale = result.scale();

        if (result.scale() < 0) {
            /*
             * Result is an integer. See if quotient represents the full integer
             * portion of the exact quotient; if it does, the computed remainder
             * will be less than the divisor.
             */
            NumberR product = (NumberR) result.multiply(divisor);

            if (((NumberR) subtract(product)).abs().compareTo(divisor.abs()) > 0) {
                throw new ArithmeticException("Division impossible");
            }

        } else if (result.scale() > 0) {
            /*
             * Integer portion of quotient will fit into precision digits;
             * recompute quotient to scale 0 to avoid double rounding and then
             * try to adjust, if necessary.
             */
            result = result.setScale(0, RoundingMode.DOWN);
        }

        int precisionDiff;
        if ((preferredScale > result.scale())
                && (precisionDiff = mc.getPrecision() - result.precision()) > 0) {
            return result.setScale(result.scale()
                    + Math.min(precisionDiff, preferredScale - result.scale));
        } else {
            return result.stripZerosToMatchScale(preferredScale);
        }

    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this % divisor)</tt>.
     *
     * <p>
     * The remainder is given by
     * <tt>this.subtract(this.divideToIntegralValue(divisor).multiply(divisor))</tt>.
     * Note that this is not the modulo operation (the result can be negative).
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     *
     * @return <tt>this % divisor</tt>.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @since 1.5
     */
    public Element remainder(
            NumberR divisor) {
        Element divrem[] = this.divideAndRemainder(divisor);
        return divrem[1];
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this % divisor)</tt>,
     * with rounding according to the context settings. The <tt>MathContext</tt>
     * settings affect the implicit divide used to compute the remainder. The
     * remainder computation itself is by definition exact. Therefore, the
     * remainder may contain more than <tt>mc.getPrecision()</tt> digits.
     *
     * <p>
     * The remainder is given by
     * <tt>this.subtract(this.divideToIntegralValue(divisor,
     * mc).multiply(divisor))</tt>. Note that this is not the modulo operation
     * (the result can be negative).
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided.
     * @param mc the context to use.
     *
     * @return <tt>this % divisor</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>, or <tt>mc.precision</tt> &gt; 0 and the
     * result of <tt>this.divideToIntgralValue(divisor)</tt> would require a
     * precision of more than <tt>mc.precision</tt> digits.
     * @see #divideToIntegralValue(java.math.BigDecimal, java.math.MathContext)
     * @since 1.5
     */
    public Element remainder(
            NumberR divisor, MathContext mc) {
        Element divrem[] = this.divideAndRemainder(divisor, mc);
        return divrem[1];
    }

    /**
     * Returns a two-element <tt>BigDecimal</tt> array containing the result of
     * <tt>divideToIntegralValue</tt> followed by the result of
     * <tt>remainder</tt> on the two operands.
     *
     * <p>
     * Note that if both the integer quotient and remainder are needed, this
     * method is faster than using the <tt>divideToIntegralValue</tt> and
     * <tt>remainder</tt> methods separately because the division need only be
     * carried out once.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided,
     * and the remainder computed.
     *
     * @return a two element <tt>BigDecimal</tt> array: the quotient (the result
     * of <tt>divideToIntegralValue</tt>) is the initial element and the
     * remainder is the final element.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @see #divideToIntegralValue(java.math.BigDecimal, java.math.MathContext)
     * @see #remainder(java.math.BigDecimal, java.math.MathContext)
     * @since 1.5
     */
    public Element[] divideAndRemainder(NumberR divisor) {
        // we use the identity  x = i * y + r to determine r
        Element[] result = new Element[2];
        result[0] = this.divideToIntegralValue(divisor);
        if (result[0] instanceof NumberR) {
            result[1] = this.subtract((NumberR) ((NumberR) result[0]).multiply(divisor));
        } else {
            result[1] = NumberR.ZERO;
        }
        return result;
    }

    @Override
    public Element[] divideAndRemainder(Element val, Ring ring) {
        if (val instanceof Polynom) {
            Polynom p = (Polynom) val;
            return p.divAndRem(p, ring);
        }
        NumberR valR = (NumberR) ((val instanceof NumberR) ? val : val.toNumber(Ring.R, ring));
        return divideAndRemainder(valR);
    }

    /**
     * Returns a two-element <tt>BigDecimal</tt> array containing the result of
     * <tt>divideToIntegralValue</tt> followed by the result of
     * <tt>remainder</tt> on the two operands calculated with rounding according
     * to the context settings.
     *
     * <p>
     * Note that if both the integer quotient and remainder are needed, this
     * method is faster than using the <tt>divideToIntegralValue</tt> and
     * <tt>remainder</tt> methods separately because the division need only be
     * carried out once.
     *
     * @param divisor value by which this <tt>BigDecimal</tt> is to be divided,
     * and the remainder computed.
     * @param mc the context to use.
     *
     * @return a two element <tt>BigDecimal</tt> array: the quotient (the result
     * of <tt>divideToIntegralValue</tt>) is the initial element and the
     * remainder is the final element.
     *
     * @throws ArithmeticException if <tt>divisor==0</tt>
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>, or <tt>mc.precision</tt> &gt; 0 and the
     * result of <tt>this.divideToIntgralValue(divisor)</tt> would require a
     * precision of more than <tt>mc.precision</tt> digits.
     * @see #divideToIntegralValue(java.math.BigDecimal, java.math.MathContext)
     * @see #remainder(java.math.BigDecimal, java.math.MathContext)
     * @since 1.5
     */
    public Element[] divideAndRemainder(NumberR divisor, MathContext mc) {
        if (mc.getPrecision() == 0) {
            return divideAndRemainder(divisor);
        }

        Element[] result = new Element[2];

        result[0] = divideToIntegralValue(divisor, mc);
        if (result[0] instanceof NumberR) {
            result[1] = this.subtract((NumberR) ((NumberR) result[0]).multiply(divisor));
        } else {
            result[1] = NumberR.ZERO;
        }

        return result;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this<sup>n</sup>)</tt>,
     * The power is computed exactly, to unlimited precision.
     *
     * <p>
     * The parameter <tt>n</tt> must be in the range 0 through 999999999,
     * inclusive. <tt>ZERO.pow(0)</tt> returns {@link
     * #ONE}.
     *
     * Note that future releases may expand the allowable exponent range of this
     * method.
     *
     * @param n power to raise this <tt>BigDecimal</tt> to.
     *
     * @return <tt>this<sup>n</sup></tt>
     *
     * @throws ArithmeticException if <tt>n</tt> is out of range.
     * @since 1.5
     */
    public NumberR powExact(int n) {
        if (n < 0 || n > 999999999) {
            throw new ArithmeticException("Invalid operation");
        }
        // No need to calculate pow(n) if result will over/underflow.
        // Don't attempt to support "supernormal" numbers.
        int newScale = checkScale((long) scale * n);
        return new NumberR(intVal.pow(n), newScale);
    }

    /**
     * ???????????????????? ?? ?????????????? ?????????????? ???? ???????????????????????? 2
     *
     * @param n ?????????? ??????????
     * @param m =1 ?????? 2
     *
     * @return this^{n/m}
     */
    @Override
    public Element powTheFirst(int n, int m, Ring ring) {
        if (n < -999999999 || n > 999999999) {
            throw new ArithmeticException("Invalid operation");
        }
        if (n == 1) {
            return rootTheFirst(m, ring);
        }
        NumberR r = pow(n, ring);
        // No need to calculate pow(n) if result will over/underflow.
        // Don't attempt to support "supernormal" numbers.
        return (m == 1) ? r : ((m == 2) ? (new com.mathpar.func.FuncNumberR(ring)).sqrt(r) 
                : rootTheFirst(m, ring));
    }

    /**
     * ???????????????????? 1-???? ???????????????? ?????????? n-???? ??????????????
     *
     * @return 1-?? ???????????????? ?????????? n-???? ??????????????
     */
    public NumberR rootTheFirst(int n) {
        return null;
    }

    /**
     * ???????????????????? n ???????????????? ?????????? n-???? ??????????????
     *
     * @return n ???????????????? ?????????? n-???? ??????????????
     */
    public Complex[] root(int n) {
        return null;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(this<sup>n</sup>)</tt>.
     * The current implementation uses the core algorithm defined in ANSI
     * standard X3.274-1996 with rounding according to the context settings. In
     * general, the returned numerical value is within two ulps of the exact
     * numerical value for the chosen precision. Note that future releases may
     * use a different algorithm with a decreased allowable error bound and
     * increased allowable exponent range.
     *
     * <p>
     * The X3.274-1996 algorithm is:
     *
     * <ul> <li> An <tt>ArithmeticException</tt> exception is thrown if <ul>
     * <li><tt>abs(n) &gt; 999999999</tt> <li><tt>mc.precision == 0</tt> and
     * <tt>n &lt; 0</tt> <li><tt>mc.precision &gt; 0</tt> and <tt>n</tt> has
     * more than <tt>mc.precision</tt> decimal digits </ul>
     *
     * <li> if <tt>n</tt> is zero, {@link #ONE} is returned even if
     * <tt>this</tt> is zero, otherwise <ul> <li> if <tt>n</tt> is positive, the
     * result is calculated via the repeated squaring technique into a single
     * accumulator. The individual multiplications with the accumulator use the
     * same math context settings as in <tt>mc</tt> except for a precision
     * increased to <tt>mc.precision + elength + 1</tt> where <tt>elength</tt>
     * is the number of decimal digits in <tt>n</tt>.
     *
     * <li> if <tt>n</tt> is negative, the result is calculated as if <tt>n</tt>
     * were positive; this value is then divided into one using the working
     * precision specified above.
     *
     * <li> The final value from either the positive or negative case is then
     * rounded to the destination precision. </ul> </ul>
     *
     * @param n power to raise this <tt>BigDecimal</tt> to.
     * @param mc the context to use.
     *
     * @return <tt>this<sup>n</sup></tt> using the ANSI standard X3.274-1996
     * algorithm
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>, or <tt>n</tt> is out of range.
     * @since 1.5
     */
    @Override
    public NumberR pow(int n, Ring ring) {
        MathContext mc = new MathContext(ring.getAccuracy(), RoundingMode.HALF_EVEN);
        //  MathContext mc = nb.MATHCONTEXT;
        if (mc.getPrecision() == 0) {
            return powExact(n);
        }
        if (n < -999999999 || n > 999999999) {
            throw new ArithmeticException("Invalid operation");
        }
        if (n == 0) {
            return ONE;
        }
        NumberR lhs = this;
        MathContext workmc = mc;           // working settings
        int mag = Math.abs(n);               // magnitude of n
        if (mc.getPrecision() > 0) {
            int elength = intLength(mag);    // length of n in digits
            if (elength > mc.getPrecision()) // X3.274 rule
            {
                throw new ArithmeticException("Invalid operation");
            }
            workmc = new MathContext(mc.getPrecision() + elength + 1,
                    mc.getRoundingMode());
        }
        // ready to carry out power calculation...
        NumberR acc = ONE;           // accumulator
        boolean seenbit = false;        // set once we've seen a 1-bit
        for (int i = 1;; i++) {            // for each bit [top bit ignored]
            mag += mag;                 // shift left 1 bit
            if (mag < 0) {              // top bit is set
                seenbit = true;         // OK, we're off
                acc
                        = (NumberR) acc.multiply(lhs, workmc); // acc=acc*x
            }

            if (i == 31) {
                break;
            }                 // that was the last bit

            if (seenbit) {
                acc = (NumberR) acc.multiply(acc, workmc);
            } // acc=acc*acc[square]
            // else (!seenbit) no point in squaring ONE

        }
        // if negative n, calculate the reciprocal using working precision
        if (n < 0) // [hence mc.precision>0]
        {
            acc = (NumberR) ONE.divide(acc, workmc);
        }
// round to final precision and strip zeros

        return acc.doRound(mc);
    }

    public NumberR pow(NumberR n, Ring ring) {
        NumberR St = (NumberR) n.multiply(this.ln(ring), ring);
        return (NumberR) St.exp(ring);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is the absolute value of this
     * <tt>BigDecimal</tt>, and whose scale is <tt>this.scale()</tt>.
     *
     * @return <tt>abs(this)</tt>
     */
    public NumberR abs() {
        return (signum() < 0 ? negate() : this);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is the absolute value of this
     * <tt>BigDecimal</tt>, with rounding according to the context settings.
     *
     * @param mc the context to use.
     *
     * @return <tt>abs(this)</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     */
    public NumberR abs(
            MathContext mc) {
        return (signum() < 0 ? negate(mc) : plus(mc));
    }

    public NumberR abs(
            int d) {
        MathContext mc = new MathContext(d);
        return (signum() < 0 ? negate(mc) : plus(mc));
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(-this)</tt>, and whose
     * scale is <tt>this.scale()</tt>.
     *
     * @return <tt>-this</tt>.
     */
    @Override
    public NumberR negate(Ring ring) {
        NumberR result = new NumberR(intVal.negate(), scale);
        result.precision = precision;
        return result;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(-this)</tt>, and whose
     * scale is <tt>this.scale()</tt>.
     *
     * @return <tt>-this</tt>.
     */
    public NumberR negate() {
        NumberR result = new NumberR(intVal.negate(), scale);
        result.precision = precision;
        return result;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(-this)</tt>, with
     * rounding according to the context settings.
     *
     * @param mc the context to use.
     *
     * @return <tt>-this</tt>, rounded as necessary.
     *
     * @throws ArithmeticException if or the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @since 1.5
     */
    public NumberR negate(
            MathContext mc) {
        return plus(mc).negate();
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(+this)</tt>, and whose
     * scale is <tt>this.scale()</tt>.
     *
     * <p>
     * This method, which simply returns this <tt>BigDecimal</tt> is included
     * for symmetry with the unary minus method {@link
     * #negate()}.
     *
     * @return <tt>this</tt>.
     *
     * @see #negate()
     * @since 1.5
     */
    public NumberR plus() {
        return this;
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose value is <tt>(+this)</tt>, with
     * rounding according to the context settings.
     *
     * <p>
     * The effect of this method is identical to that of the {@link
     * #round(MathContext)} method.
     *
     * @param mc the context to use.
     *
     * @return <tt>this</tt>, rounded as necessary. A zero result will have a
     * scale of 0.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     * @see #round(MathContext)
     * @since 1.5
     */
    public NumberR plus(
            MathContext mc) {
        if (mc.getPrecision() == 0) // no rounding please
        {
            return this;
        }

        NumberR rhs = this;
        return rhs.doRound(mc);
    }

    /**
     * Returns the signum function of this <tt>BigDecimal</tt>.
     *
     * @return -1, 0, or 1 as the value of this <tt>BigDecimal</tt> is negative,
     * zero, or positive.
     */
    @Override
    public int signum() {
        return intVal.signum();
    }

    /**
     * Returns the <i>scale</i> of this <tt>BigDecimal</tt>. If zero or
     * positive, the scale is the number of digits to the right of the decimal
     * point. If negative, the unscaled value of the number is multiplied by ten
     * to the power of the negation of the scale. For example, a scale of
     * <tt>-3</tt> means the unscaled value is multiplied by 1000.
     *
     * @return the scale of this <tt>BigDecimal</tt>.
     */
    public int scale() {
        return scale;
    }

    /**
     * Returns the <i>precision</i> of this <tt>BigDecimal</tt>. (The precision
     * is the number of digits in the unscaled value.)
     *
     * <p>
     * The precision of a zero value is 1.
     *
     * @return the precision of this <tt>BigDecimal</tt>.
     *
     * @since 1.5
     */
    public int precision() {
        int result = precision;
        if (result == 0) {
            result = digitLength();
            precision = result;
        }

        return result;
    }

    /**
     * Returns a <tt>BigInteger</tt> whose value is the <i>unscaled value</i> of
     * this <tt>BigDecimal</tt>. (Computes <tt>(this *
     * 10<sup>this.scale()</sup>)</tt>.)
     *
     * @return the unscaled value of this <tt>BigDecimal</tt>.
     *
     * @since 1.2
     */
    public NumberZ unscaledValue() {
        return intVal;
    }
// Rounding Modes
    /**
     * Rounding mode to round away from zero. Always increments the digit prior
     * to a nonzero discarded fraction. Note that this rounding mode never
     * decreases the magnitude of the calculated value.
     */
    public final static int ROUND_UP = 0;
    /**
     * Rounding mode to round towards zero. Never increments the digit prior to
     * a discarded fraction (i.e., truncates). Note that this rounding mode
     * never increases the magnitude of the calculated value.
     */
    public final static int ROUND_DOWN = 1;
    /**
     * Rounding mode to round towards positive infinity. If the
     * <tt>BigDecimal</tt> is positive, behaves as for <tt>ROUND_UP</tt>; if
     * negative, behaves as for <tt>ROUND_DOWN</tt>. Note that this rounding
     * mode never decreases the calculated value.
     */
    public final static int ROUND_CEILING = 2;
    /**
     * Rounding mode to round towards negative infinity. If the
     * <tt>BigDecimal</tt> is positive, behave as for <tt>ROUND_DOWN</tt>; if
     * negative, behave as for <tt>ROUND_UP</tt>. Note that this rounding mode
     * never increases the calculated value.
     */
    public final static int ROUND_FLOOR = 3;
    /**
     * Rounding mode to round towards &quot;nearest neighbor&quot; unless both
     * neighbors are equidistant, in which case round up. Behaves as for
     * <tt>ROUND_UP</tt> if the discarded fraction is &gt;= 0.5; otherwise,
     * behaves as for <tt>ROUND_DOWN</tt>. Note that this is the rounding mode
     * that most of us were taught in grade school.
     */
    public final static int ROUND_HALF_UP = 4;
    /**
     * Rounding mode to round towards &quot;nearest neighbor&quot; unless both
     * neighbors are equidistant, in which case round down. Behaves as for
     * <tt>ROUND_UP</tt> if the discarded fraction is &gt; 0.5; otherwise,
     * behaves as for <tt>ROUND_DOWN</tt>.
     */
    public final static int ROUND_HALF_DOWN = 5;
    /**
     * Rounding mode to round towards the &quot;nearest neighbor&quot; unless
     * both neighbors are equidistant, in which case, round towards the even
     * neighbor. Behaves as for <tt>ROUND_HALF_UP</tt> if the digit to the left
     * of the discarded fraction is odd; behaves as for <tt>ROUND_HALF_DOWN</tt>
     * if it's even. Note that this is the rounding mode that minimizes
     * cumulative error when applied repeatedly over a sequence of calculations.
     */
    public final static int ROUND_HALF_EVEN = 6;
    /**
     * Rounding mode to assert that the requested operation has an exact result,
     * hence no rounding is necessary. If this rounding mode is specified on an
     * operation that yields an inexact result, an <tt>ArithmeticException</tt>
     * is thrown.
     */
    public final static int ROUND_UNNECESSARY = 7;

    // Scaling/Rounding Operations
    /**
     * Returns a <tt>BigDecimal</tt> rounded according to the
     * <tt>MathContext</tt> settings. If the precision setting is 0 then no
     * rounding takes place.
     *
     * <p>
     * The effect of this method is identical to that of the
     * {@link #plus(MathContext)} method.
     *
     * @param mc the context to use.
     *
     * @return a <tt>BigDecimal</tt> rounded according to the
     * <tt>MathContext</tt> settings.
     *
     * @throws ArithmeticException if the rounding mode is <tt>UNNECESSARY</tt>
     * and the <tt>BigDecimal</tt> operation would require rounding.
     * @see #plus(MathContext)
     * @since 1.5
     */
    public NumberR round(MathContext mc) {
        return plus(mc);
    }

    // Scaling/Rounding Operations
    /**
     * Returns a <tt>BigDecimal</tt> rounded according to the
     * <tt>MathContext</tt> settings. If the precision setting is 0 then no
     * rounding takes place.
     *
     * <p>
     * The effect of this method is identical to that of the
     * {@link #plus(MathContext)} method.
     *
     * @param ring
     *
     * @return a <tt>BigDecimal</tt> rounded according to the
     * <tt>MathContext</tt> settings.
     *
     * @throws ArithmeticException if the rounding mode is <tt>UNNECESSARY</tt>
     * and the <tt>BigDecimal</tt> operation would require rounding.
     * @see #plus(MathContext)
     * @since 1.5
     */
    @Override
    public NumberZ round(Ring ring) {
        return NumberRtoNumberZ();
    }

    /**
     * @param ring - Ring
     *
     * @return the biggest integer number, which not less then this
     */
    @Override
    public NumberZ ceil(Ring ring) {
        NumberZ rounded = round(ring);
        int sign = compareTo(rounded, ring);
        return (sign > 0) ? rounded.add(NumberZ.ONE) : rounded;
    }

    /**
     * @param ring - Ring
     *
     * @return the biggest integer number, which not greater then this
     */
    @Override
    public NumberZ floor(Ring ring) {
        NumberZ rounded = round(ring);
        int sign = compareTo(rounded, ring);
        return (sign < 0) ? rounded.subtract(NumberZ.ONE) : rounded;
    }

    /**
     * The number by module mod.
     *
     * @param mod -module
     * @param ring - Ring
     *
     * @return this - mod * floor(this);
     */
    @Override
    public Element mod(Element mod, Ring ring) {
        NumberR modR = (NumberR) mod.toNewRing(Ring.R, ring);
        return mod(modR, ring);
    }

    /**
     * The number by module mod.
     *
     * @param mod -module
     * @param ring - Ring
     *
     * @return this - mod * floor(this);
     */
    public Element mod(NumberR mod, Ring ring) {
        Element q = divide(mod, ring.MC);
        if (q instanceof NumberR) {
            return subtract((NumberR) mod.multiply(new NumberR(((NumberR) q).floor(ring)), ring.MC));
        } else {
            return NAN;
        }
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose unscaled value is determined by multiplying or dividing this
     * <tt>BigDecimal</tt>'s unscaled value by the appropriate power of ten to
     * maintain its overall value. If the scale is reduced by the operation, the
     * unscaled value must be divided (rather than multiplied), and the value
     * may be changed; in this case, the specified rounding mode is applied to
     * the division.
     * <p>
     * Note that since BigDecimal objects are immutable, calls of this method do
     * <i>not</i> result in the original object being modified, contrary to the
     * usual convention of having methods named
     * <tt>set<i>X</i></tt> mutate field <tt><i>X</i></tt>. Instead,
     * <tt>setScale</tt> returns an object with the proper scale; the returned
     * object may or may not be newly
     * <p>
     * The new {@link #setScale(int, RoundingMode)} method should be used in
     * preference to this legacy method.
     *
     * @param newScale scale of the <tt>BigDecimal</tt> value to be returned.
     * @param roundingMode The rounding mode to apply.
     *
     * @return a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose unscaled value is determined by multiplying or dividing this
     * <tt>BigDecimal</tt>'s unscaled value by the appropriate power of ten to
     * maintain its overall value.
     *
     * @throws ArithmeticException if <tt>roundingMode==ROUND_UNNECESSARY</tt>
     * and the specified scaling operation would require rounding.
     * @throws IllegalArgumentException if <tt>roundingMode</tt> does not
     * represent a valid rounding mode.
     */
    public NumberR setScale(
            int newScale, int roundingMode) {
        if (roundingMode < ROUND_UP || roundingMode > ROUND_UNNECESSARY) {
            throw new IllegalArgumentException("Invalid rounding mode");
        }

        if (newScale == this.scale) // easy case
        {
            return this;
        }

        if (this.signum() == 0) // zero can have any scale
        {
            return new NumberR(NumberZ.ZERO, newScale);
        }

        if (newScale > this.scale) {
            // [we can use checkScale to assure multiplier is valid]
            int raise = checkScale((long) newScale - this.scale);
            NumberR result = new NumberR(
                    intVal.multiply(tenToThe(raise)), newScale);
            if (this.precision > 0) {
                result.precision = this.precision + newScale - this.scale;
            }

            return result;
        }
        /*
         * scale < this.scale
         */
// we cannot perfectly predict the precision after rounding

        return divide(ONE, newScale, roundingMode);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose unscaled value is determined by multiplying or dividing this
     * <tt>BigDecimal</tt>'s unscaled value by the appropriate power of ten to
     * maintain its overall value. If the scale is reduced by the operation, the
     * unscaled value must be divided (rather than multiplied), and the value
     * may be changed; in this case, the specified rounding mode is applied to
     * the division.
     *
     * @param newScale scale of the <tt>BigDecimal</tt> value to be returned.
     * @param roundingMode The rounding mode to apply.
     *
     * @return a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose unscaled value is determined by multiplying or dividing this
     * <tt>BigDecimal</tt>'s unscaled value by the appropriate power of ten to
     * maintain its overall value.
     *
     * @throws ArithmeticException if <tt>roundingMode==UNNECESSARY</tt> and the
     * specified scaling operation would require rounding.
     * @see RoundingMode
     * @since 1.5
     */
    public NumberR setScale(
            int newScale, RoundingMode roundingMode) {
        return setScale(newScale, roundingMode.oldMode);
    }

    /**
     * Returns a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose value is numerically equal to this <tt>BigDecimal</tt>'s. Throws an
     * <tt>ArithmeticException</tt> if this is not possible.
     *
     * <p>
     * This call is typically used to increase the scale, in which case it is
     * guaranteed that there exists a <tt>BigDecimal</tt> of the specified scale
     * and the correct value. The call can also be used to reduce the scale if
     * the caller knows that the <tt>BigDecimal</tt> has sufficiently many zeros
     * at the end of its fractional part (i.e., factors of ten in its integer
     * value) to allow for the rescaling without changing its value.
     *
     * <p>
     * This method returns the same result as the two-argument versions of
     * <tt>setScale</tt>, but saves the caller the trouble of specifying a
     * rounding mode in cases where it is irrelevant.
     *
     * <p>
     * Note that since <tt>BigDecimal</tt> objects are immutable, calls of this
     * method do <i>not</i> result in the original object being modified,
     * contrary to the usual convention of having methods named
     * <tt>set<i>X</i></tt> mutate field <tt><i>X</i></tt>. Instead,
     * <tt>setScale</tt> returns an object with the proper scale; the returned
     * object may or may not be newly allocated.
     *
     * @param newScale scale of the <tt>BigDecimal</tt> value to be returned.
     *
     * @return a <tt>BigDecimal</tt> whose scale is the specified value, and
     * whose unscaled value is determined by multiplying or dividing this
     * <tt>BigDecimal</tt>'s unscaled value by the appropriate power of ten to
     * maintain its overall value.
     *
     * @throws ArithmeticException if the specified scaling operation would
     * require rounding.
     * @see #setScale(int, int)
     * @see #setScale(int, RoundingMode)
     */
    public NumberR setScale(
            int newScale) {
        return setScale(newScale, ROUND_UNNECESSARY);
    }

// Decimal Point Motion Operations
    /**
     * Returns a <tt>BigDecimal</tt> which is equivalent to this one with the
     * decimal point moved <tt>n</tt> places to the left. If <tt>n</tt> is
     * non-negative, the call merely adds <tt>n</tt> to the scale. If <tt>n</tt>
     * is negative, the call is equivalent to <tt>movePointRight(-n)</tt>. The
     * <tt>BigDecimal</tt> returned by this call has value <tt>(this &times;
     * 10<sup>-n</sup>)</tt> and scale <tt>max(this.scale()+n, 0)</tt>.
     *
     * @param n number of places to move the decimal point to the left.
     *
     * @return a <tt>BigDecimal</tt> which is equivalent to this one with the
     * decimal point moved <tt>n</tt> places to the left.
     *
     * @throws ArithmeticException if scale overflows.
     */
    public NumberR movePointLeft(
            int n) {
        // Cannot use movePointRight(-n) in case of n==Integer.MIN_VALUE
        NumberR num = new NumberR(intVal, checkScale((long) scale + n));
        return (num.scale < 0 ? num.setScale(0) : num);
    }

    /**
     * Returns a <tt>BigDecimal</tt> which is equivalent to this one with the
     * decimal point moved <tt>n</tt> places to the right. If <tt>n</tt> is
     * non-negative, the call merely subtracts <tt>n</tt> from the scale. If
     * <tt>n</tt> is negative, the call is equivalent to
     * <tt>movePointLeft(-n)</tt>. The <tt>BigDecimal</tt> returned by this call
     * has value <tt>(this &times; 10<sup>n</sup>)</tt> and scale
     * <tt>max(this.scale()-n, 0)</tt>.
     *
     * @param n number of places to move the decimal point to the right.
     *
     * @return a <tt>BigDecimal</tt> which is equivalent to this one with the
     * decimal point moved <tt>n</tt> places to the right.
     *
     * @throws ArithmeticException if scale overflows.
     */
    public NumberR movePointRight(
            int n) {
        // Cannot use movePointLeft(-n) in case of n==Integer.MIN_VALUE
        NumberR num = new NumberR(intVal, checkScale((long) scale - n));
        return (num.scale < 0 ? num.setScale(0) : num);
    }

    /**
     * Returns a BigDecimal whose numerical value is equal to (<tt>this</tt> *
     * 10<sup>n</sup>). The scale of the result is <tt>(this.scale() - n)</tt>.
     *
     * @throws ArithmeticException if the scale would be outside the range of a
     * 32-bit integer.
     *
     * @since 1.5
     */
    public NumberR scaleByPowerOfTen(
            int n) {
        NumberR num = new NumberR(intVal, checkScale((long) scale - n));
        num.precision = precision;
        return num;
    }

    /**
     * Returns a <tt>BigDecimal</tt> which is numerically equal to this one but
     * with any trailing zeros removed from the representation. For example,
     * stripping the trailing zeros from the <tt>BigDecimal</tt> value
     * <tt>600.0</tt>, which has [<tt>BigInteger</tt>, <tt>scale</tt>]
     * components equals to [6000, 1], yields <tt>6E2</tt> with
     * [<tt>BigInteger</tt>, <tt>scale</tt>] components equals to [6, -2]
     *
     * @return a numerically equal <tt>BigDecimal</tt> with any trailing zeros
     * removed.
     */
    public NumberR stripTrailingZeros() {
        return (new NumberR(intVal, scale)).stripZerosToMatchScale(NumberZ64.MIN_VALUE);
    }

// Comparison Operations
    /**
     * Compares this <tt>BigDecimal</tt> with the specified <tt>BigDecimal</tt>.
     * Two <tt>BigDecimal</tt> objects that are equal in value but have a
     * different scale (like 2.0 and 2.00) are considered equal by this method.
     * This method is provided in preference to individual methods for each of
     * the six boolean comparison operators (&lt;, ==, &gt;, &gt;=, !=, &lt;=).
     * The suggested idiom for performing these comparisons is:
     * <tt>(x.compareTo(y)</tt> &lt;<i>op</i>&gt; <tt>0)</tt>, where
     * &lt;<i>op</i>&gt; is one of the six comparison operators.
     *
     * @param val <tt>BigDecimal</tt> to which this <tt>BigDecimal</tt> is to be
     * compared.
     *
     * @return -1, 0, or 1 as this <tt>BigDecimal</tt> is numerically less than,
     * equal to, or greater than <tt>val</tt>.
     */
    public int compareTo(NumberR val) {
        /*
         * Optimization: would run fine without the next three lines
         */
        int sigDiff = signum() - val.signum();
        if (sigDiff != 0) {
            return (sigDiff > 0 ? 1 : -1);
        }

// If the (adjusted) exponents are different we do not need to
// expensively match scales and compare the significands
        int aethis = this.precision() - this.scale;    // [-1]
        int aeval = val.precision() - val.scale;     // [-1]
        if (aethis < aeval) {
            return -this.signum();
        } else if (aethis > aeval) {
            return this.signum();
        }

// Scale and compare intVals
        NumberR arg[] = new NumberR[2];
        arg[0] = this;
        arg[1] = val;
        matchScale(arg);
        return arg[0].intVal.compareTo(arg[1].intVal);
    }

    public static int compareTo(Element p, Element pp) {
        return ((NumberR) p).compareTo((NumberR) pp);
    }

    @Override
    public Boolean isNegative() {
        return (signum() < 0) ? true : false;
    }

    /**
     * Compares this <tt>BigDecimal</tt> with the specified <tt>Object</tt> for
     * equality. Unlike {@link
     * #compareTo(BigDecimal) compareTo}, this method considers two
     * <tt>BigDecimal</tt> objects equal only if they are equal in value and
     * scale (thus 2.0 is not equal to 2.00 when compared by this method).
     *
     * @param x <tt>Object</tt> to which this <tt>BigDecimal</tt> is to be
     * compared.
     *
     * @return <tt>true</tt> if and only if the specified <tt>Object</tt> is a
     * <tt>BigDecimal</tt> whose value and scale are equal to this
     * <tt>BigDecimal</tt>'s.
     *
     * @see #compareTo(java.math.BigDecimal)
     * @see #hashCode
     */
    @Override
    public boolean equals(Element x, Ring r) {
        return compareTo(x, r) == 0;
    }

    public boolean equals(NumberR x) {
        return scale == x.scale && intVal.equals(x.intVal);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof NumberR) {
            return equals((NumberR) obj);
        }
        return false;
    }

    @Override
    public Boolean isZero(Ring ring) {
        NumberR eps = ring.MachineEpsilonR;
        return (((NumberR) abs(ring)).compareTo(eps)) < 0 ? true : false;
    }

    @Override
    public Boolean isOne(Ring r) {
        return (subtract(NumberR.ONE)).isZero(r);
    }

    @Override
    public Boolean isMinusOne(Ring r) {
        return (add(NumberR.ONE)).isZero(r);
    }

    /**
     * Returns the minimum of this <tt>BigDecimal</tt> and <tt>val</tt>.
     *
     * @param val value with which the minimum is to be computed.
     *
     * @return the <tt>BigDecimal</tt> whose value is the lesser of this
     * <tt>BigDecimal</tt> and <tt>val</tt>. If they are equal, as defined by
     * the {@link #compareTo(BigDecimal) compareTo} method, <tt>this</tt> is
     * returned.
     *
     * @see #compareTo(java.math.BigDecimal)
     */
    public NumberR min(
            NumberR val) {
        return (compareTo(val) <= 0 ? this : val);
    }

    /**
     * Returns the maximum of this <tt>BigDecimal</tt> and <tt>val</tt>.
     *
     * @param val value with which the maximum is to be computed.
     *
     * @return the <tt>BigDecimal</tt> whose value is the greater of this
     * <tt>BigDecimal</tt> and <tt>val</tt>. If they are equal, as defined by
     * the {@link #compareTo(BigDecimal) compareTo} method, <tt>this</tt> is
     * returned.
     *
     * @see #compareTo(java.math.BigDecimal)
     */
    public NumberR max(
            NumberR val) {
        return (compareTo(val) >= 0 ? this : val);
    }

// Hash Function
    /**
     * Returns the hash code for this <tt>BigDecimal</tt>. Note that two
     * <tt>BigDecimal</tt> objects that are numerically equal but differ in
     * scale (like 2.0 and 2.00) will generally <i>not</i> have the same hash
     * code.
     *
     * @return hash code for this <tt>BigDecimal</tt>.
     *
     * @see #equals(Object)
     */
    @Override
    public int hashCode() {
        return 31 * intVal.hashCode() + scale;
    }

// Format Converters
    /**
     * Returns the string representation of this <tt>BigDecimal</tt>, using
     * scientific notation if an exponent is needed.
     *
     * <p>
     * A standard canonical string form of the <tt>BigDecimal</tt> is created as
     * though by the following steps: first, the absolute value of the unscaled
     * value of the <tt>BigDecimal</tt> is converted to a string in base ten
     * using the characters <tt>'0'</tt> through <tt>'9'</tt> with no leading
     * zeros (except if its value is zero, in which case a single
     * <tt>'0'</tt> character is used).
     *
     * <p>
     * Next, an <i>adjusted exponent</i> is calculated; this is the negated
     * scale, plus the number of characters in the converted unscaled value,
     * less one. That is, <tt>-scale+(ulength-1)</tt>, where <tt>ulength</tt> is
     * the length of the absolute value of the unscaled value in decimal digits
     * (its <i>precision</i>).
     *
     * <p>
     * If the scale is greater than or equal to zero and the adjusted exponent
     * is greater than or equal to <tt>-6</tt>, the number will be converted to
     * a character form without using exponential notation. In this case, if the
     * scale is zero then no decimal point is added and if the scale is positive
     * a decimal point will be inserted with the scale specifying the number of
     * characters to the right of the decimal point.
     * <tt>'0'</tt> characters are added to the left of the converted unscaled
     * value as necessary. If no character precedes the decimal point after this
     * insertion then a conventional <tt>'0'</tt> character is prefixed.
     *
     * <p>
     * Otherwise (that is, if the scale is negative, or the adjusted exponent is
     * less than <tt>-6</tt>), the number will be converted to a character form
     * using exponential notation. In this case, if the converted
     * <tt>BigInteger</tt> has more than one digit a decimal point is inserted
     * after the first digit. An exponent in character form is then suffixed to
     * the converted unscaled value (perhaps with inserted decimal point); this
     * comprises the letter <tt>'E'</tt> followed immediately by the adjusted
     * exponent converted to a character form. The latter is in base ten, using
     * the characters <tt>'0'</tt> through <tt>'9'</tt> with no leading zeros,
     * and is always prefixed by a sign character <tt>'-'</tt>
     * (<tt>'&#92;u002D'</tt>) if the adjusted exponent is negative,
     * <tt>'+'</tt> (<tt>'&#92;u002B'</tt>) otherwise).
     *
     * <p>
     * Finally, the entire string is prefixed by a minus sign character
     * <tt>'-'</tt> (<tt>'&#92;u002D'</tt>) if the unscaled value is less than
     * zero. No sign character is prefixed if the unscaled value is zero or
     * positive.
     *
     * <p>
     * <b>Examples:</b>
     * <p>
     * For each representation [<i>unscaled value</i>,
     * <i>scale</i>] on the left, the resulting string is shown on the right.
     * <pre>
     * [123,0]      &quot;123&quot;
     * [-123,0]     &quot;-123&quot;
     * [123,-1]     &quot;1.23E+3&quot;
     * [123,-3]     &quot;1.23E+5&quot;
     * [123,1]      &quot;12.3&quot;
     * [123,5]      &quot;0.00123&quot;
     * [123,10]     &quot;1.23E-8&quot;
     * [-123,12]    &quot;-1.23E-10&quot;
     * </pre>
     *
     * <b>Notes:</b> <ol>
     *
     * <li>There is a one-to-one mapping between the distinguishable
     * <tt>BigDecimal</tt> values and the result of this conversion. That is,
     * every distinguishable <tt>BigDecimal</tt> value (unscaled value and
     * scale) has a unique string representation as a result of using
     * <tt>toString</tt>. If that string representation is converted back to a
     * <tt>BigDecimal</tt> using the {@link #BigDecimal(String)} constructor,
     * then the original value will be recovered.
     *
     * <li>The string produced for a given number is always the same; it is not
     * affected by locale. This means that it can be used as a canonical string
     * representation for exchanging decimal data, or as a key for a Hashtable,
     * etc. Locale-sensitive number formatting and parsing is handled by the {@link
     * java.text.NumberFormat} class and its subclasses.
     *
     * <li>The {@link #toEngineeringString} method may be used for presenting
     * numbers with exponents in engineering notation, and the
     * {@link #setScale(int,RoundingMode) setScale} method may be used for
     * rounding a <tt>BigDecimal</tt> so it has a known number of digits after
     * the decimal point.
     *
     * <li>The digit-to-character mapping provided by
     * <tt>Character.forDigit</tt> is used.
     *
     * </ol>
     *
     * @return string representation of this <tt>BigDecimal</tt>.
     *
     * @see Character#forDigit
     * @see #BigDecimal(java.lang.String)
     */ 
    @Override
    public String toString(Ring r) {
        Element rrr=this;
        String st = // new String(layoutChars(true));
                toPlainString();
        int zap = st.indexOf('.') + 1;
        int fpos=r.FLOATPOS;
        int e = Math.max(st.indexOf('e'), st.indexOf('E'));
        int ee = (e == -1) ? st.length() : e;
        int endIndex = (zap==0)? ee + fpos:zap+fpos+((fpos<0)?-1:0);  
        if(endIndex<1){ fpos=fpos-endIndex+1;  endIndex=1; }
        String res = st;
        if  (ee > endIndex)  {
            res = st.substring(0, endIndex);
            char c=res.charAt(res.length()-1);
            if((c=='.')||(','==c)){res=res.substring(0, res.length()-1); 
                                   c=res.charAt(res.length()-1);}
            char rnd=st.charAt(endIndex);
            RND:{
            if(rnd>'4'){
                if((c<'9')&&(c>='0')) res=res.substring(0, res.length()-1)+(char)(c+1);
                else{
                  char[] car=res.toCharArray(); int len=car.length;
                  for (int i = len-1; i >-1; i--) { 
                    if(car[i]=='9')car[i]='0';
                    else{  if((car[i]<'9')&&(car[i]>='0')){
                           car[i]=(char)(1+car[i]);
                           res=String.valueOf(car); break RND; }   
                    }
                  } res="1"+String.valueOf(car); }
            }}
            if(fpos<0)for (int i = 0; i < -fpos; i++) {res=res+'0';}
            if (e != -1)  res = res + st.substring(e);           
        }return res;
    }

    @Override
    public String toString() {
        return new String(layoutChars(true));
    }

    @Override
    public String toString(
            int i, Ring r) {
        return divide(new NumberR(1), new MathContext(i)).toString();
    }

    /**
     * Returns a string representation of this <tt>BigDecimal</tt>, using
     * engineering notation if an exponent is needed.
     *
     * <p>
     * Returns a string that represents the <tt>BigDecimal</tt> as described in
     * the {@link #toString()} method, except that if exponential notation is
     * used, the power of ten is adjusted to be a multiple of three (engineering
     * notation) such that the integer part of nonzero values will be in the
     * range 1 through 999. If exponential notation is used for zero values, a
     * decimal point and one or two fractional zero digits are used so that the
     * scale of the zero value is preserved. Note that unlike the output of
     * {@link #toString()}, the output of this method is <em>not</em> guaranteed
     * to recover the same [integer, scale] pair of this <tt>BigDecimal</tt> if
     * the output string is converting back to a <tt>BigDecimal</tt> using the {@linkplain
     * #BigDecimal(String) string constructor}. The result of this method meets
     * the weaker constraint of always producing a numerically equal result from
     * applying the string constructor to the method's output.
     *
     * @return string representation of this <tt>BigDecimal</tt>, using
     * engineering notation if an exponent is needed.
     *
     * @since 1.5
     */
    public String toEngineeringString() {
        return new String(layoutChars(false));
    }

    /**
     * Returns a string representation of this <tt>BigDecimal</tt> without an
     * exponent field. For values with a positive scale, the number of digits to
     * the right of the decimal point is used to indicate scale. For values with
     * a zero or negative scale, the resulting string is generated as if the
     * value were converted to a numerically equal value with zero scale and as
     * if all the trailing zeros of the zero scale value were present in the
     * result.
     *
     * The entire string is prefixed by a minus sign character '-'
     * (<tt>'&#92;u002D'</tt>) if the unscaled value is less than zero. No sign
     * character is prefixed if the unscaled value is zero or positive.
     *
     * Note that if the result of this method is passed to the
     * {@linkplain #BigDecimal(String) string constructor}, only the numerical
     * value of this <tt>BigDecimal</tt> will necessarily be recovered; the
     * representation of the new <tt>BigDecimal</tt> may have a different scale.
     * In particular, if this <tt>BigDecimal</tt> has a positive scale, the
     * string resulting from this method will have a scale of zero when
     * processed by the string constructor.
     *
     * (This method behaves analogously to the <tt>toString</tt> method in 1.4
     * and earlier releases.)
     *
     * @return a string representation of this <tt>BigDecimal</tt> without an
     * exponent field.
     *
     * @since 1.5
     * @see #toString()
     * @see #toEngineeringString()
     */
    public String toPlainString() {
        NumberR bd = this;
        if (bd.scale < 0) {
            bd = bd.setScale(0);
        }
        if (bd.scale == 0) /*  No decimal point */ {
            return bd.intVal.toString();
        }
        return bd.getValueString(bd.signum(), bd.intVal.abs().toString(), bd.scale);
    }

    /*
     * Returns a digit.digit string
     */
    private String getValueString(int signum, String intString, int scale) {
        /*
         * Insert decimal point
         */
        StringBuilder buf;
        int insertionPoint = intString.length() - scale;
        if (insertionPoint == 0) {  /*
             * Point goes right before intVal
             */

            return (signum < 0 ? "-0." : "0.") + intString;
        } else if (insertionPoint > 0) { /*
             * Point goes inside intVal
             */

            buf = new StringBuilder(intString);
            buf.insert(insertionPoint, '.');
            if (signum < 0) {
                buf.insert(0, '-');
            }

        } else { /*
             * We must insert zeros between point and intVal
             */

            buf = new StringBuilder(3 - insertionPoint + intString.length());
            buf.append(signum < 0 ? "-0." : "0.");
            for (int i = 0; i
                    < -insertionPoint; i++) {
                buf.append('0');
            }

            buf.append(intString);
        }

        return buf.toString();
    }

    /**
     * Returns the string representation of this <tt>BigDecimal</tt>, as a
     * character array. This returns a character array that has the same content
     * as a call to {@link #toString()} followed by a call to
     * {@link String#toCharArray()}, but without the overhead of creating the
     * intermediate <tt>String</tt> object.
     *
     * @return <tt>char[]</tt> representation of this <tt>BigDecimal</tt>.
     *
     * @see #toString()
     * @since 1.5
     */
    char[] toCharArray() {         // for use within java.math only
        return layoutChars(true);
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>BigInteger</tt>. This
     * conversion is analogous to a <a
     * href="http://java.sun.com/docs/books/jls/second_edition/html/conversions.doc.html#25363"><i>narrowing
     * primitive conversion</i></a> from <tt>double</tt> to <tt>long</tt> as
     * defined in the <a href="http://java.sun.com/docs/books/jls/html/">Java
     * Language Specification</a>: any fractional part of this
     * <tt>BigDecimal</tt> will be discarded. Note that this conversion can lose
     * information about the precision of the <tt>BigDecimal</tt> value.
     * <p>
     * To have an exception thrown if the conversion is inexact (in other words
     * if a nonzero fractional part is discarded), use the
     * {@link #toBigIntegerExact()} method.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>BigInteger</tt>.
     */
    public NumberZ NumberRtoNumberZ() {
        // force to an integer, quietly
        return this.setScale(0, ROUND_DOWN).intVal;
    }

    public Fraction NumberRtoNumberQ() {
        return new Fraction(intVal, NumberZ.tenInPowerOf(scale));

    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>BigInteger</tt>, checking for
     * lost information. An exception is thrown if this <tt>BigDecimal</tt> has
     * a nonzero fractional part.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>BigInteger</tt>.
     *
     * @throws ArithmeticException if <tt>this</tt> has a nonzero fractional
     * part.
     * @since 1.5
     */
    public NumberZ toBigIntegerExact() {
        // round to an integer, with exception if decimal part non-0
        return this.setScale(0, ROUND_UNNECESSARY).intVal;
    }

    public NumberZ toBigInteger() {
        // round to an integer, with exception if decimal part non-0
        return this.intVal;
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>long</tt>. This conversion is
     * analogous to a <a
     * href="http://java.sun.com/docs/books/jls/second_edition/html/conversions.doc.html#25363"><i>narrowing
     * primitive conversion</i></a> from <tt>double</tt> to <tt>short</tt> as
     * defined in the <a href="http://java.sun.com/docs/books/jls/html/">Java
     * Language Specification</a>: any fractional part of this
     * <tt>BigDecimal</tt> will be discarded, and if the resulting
     * &quot;<tt>BigInteger</tt>&quot; is too big to fit in a <tt>long</tt>,
     * only the low-order 64 bits are returned. Note that this conversion can
     * lose information about the overall magnitude and precision of this
     * <tt>BigDecimal</tt> value as well as return a result with the opposite
     * sign.
     *
     * @return this <tt>BigDecimal</tt> converted to an <tt>long</tt>.
     */
    @Override
    public long longValue() {
        return NumberRtoNumberZ().longValue();
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>long</tt>, checking for lost
     * information. If this <tt>BigDecimal</tt> has a nonzero fractional part or
     * is out of the possible range for a <tt>long</tt> result then an
     * <tt>ArithmeticException</tt> is thrown.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>long</tt>.
     *
     * @throws ArithmeticException if <tt>this</tt> has a nonzero fractional
     * part, or will not fit in a <tt>long</tt>.
     * @since 1.5
     */
    public long longValueExact() {
        // If more than 19 digits in integer part it cannot possibly fit
        if ((precision() - scale) > 19) // [OK for neagtive scale too]
        {
            throw new java.lang.ArithmeticException("Overflow");
        }
// Fastpath zero and < 1.0 numbers (the latter can be very slow
// to round if very small)

        if (this.intVal.signum() == 0) {
            return 0;
        }

        if ((this.precision() - this.scale) <= 0) {
            throw new ArithmeticException("Rounding necessary");
        }
// round to an integer, with exception if decimal part non-0

        NumberR num = this.setScale(0, ROUND_UNNECESSARY);
        if (num.precision() >= 19) {    // need to check carefully
            if (LONGMIN == null) {      // initialize constants
                LONGMIN = NumberZ.valueOf(NumberZ64.MIN_VALUE);
                LONGMAX
                        = NumberZ.valueOf(NumberZ64.MAX_VALUE);
            }

            if ((num.intVal.compareTo(LONGMIN) < 0)
                    || (num.intVal.compareTo(LONGMAX) > 0)) {
                throw new java.lang.ArithmeticException("Overflow");
            }

        }
        return num.intVal.longValue();
    }
// These constants are only initialized if needed
    /**
     * BigInteger equal to Long.MIN_VALUE.
     */
    private static NumberZ LONGMIN = null;
    /**
     * BigInteger equal to Long.MAX_VALUE.
     */
    private static NumberZ LONGMAX = null;

    /**
     * Converts this <tt>BigDecimal</tt> to an <tt>int</tt>. This conversion is
     * analogous to a <a
     * href="http://java.sun.com/docs/books/jls/second_edition/html/conversions.doc.html#25363"><i>narrowing
     * primitive conversion</i></a> from <tt>double</tt> to <tt>short</tt> as
     * defined in the <a href="http://java.sun.com/docs/books/jls/html/">Java
     * Language Specification</a>: any fractional part of this
     * <tt>BigDecimal</tt> will be discarded, and if the resulting
     * &quot;<tt>BigInteger</tt>&quot; is too big to fit in an <tt>int</tt>,
     * only the low-order 32 bits are returned. Note that this conversion can
     * lose information about the overall magnitude and precision of this
     * <tt>BigDecimal</tt> value as well as return a result with the opposite
     * sign.
     *
     * @return this <tt>BigDecimal</tt> converted to an <tt>int</tt>.
     */
    @Override
    public int intValue() {
        return NumberRtoNumberZ().intValue();
    }

    /**
     * Converts this <tt>BigDecimal</tt> to an <tt>int</tt>, checking for lost
     * information. If this <tt>BigDecimal</tt> has a nonzero fractional part or
     * is out of the possible range for an <tt>int</tt> result then an
     * <tt>ArithmeticException</tt> is thrown.
     *
     * @return this <tt>BigDecimal</tt> converted to an <tt>int</tt>.
     *
     * @throws ArithmeticException if <tt>this</tt> has a nonzero fractional
     * part, or will not fit in an <tt>int</tt>.
     * @since 1.5
     */
    public int intValueExact() {
        long num;
        num
                = this.longValueExact();     // will check decimal part
        if ((int) num != num) {
            throw new java.lang.ArithmeticException("Overflow");
        }

        return (int) num;
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>short</tt>, checking for lost
     * information. If this <tt>BigDecimal</tt> has a nonzero fractional part or
     * is out of the possible range for a <tt>short</tt> result then an
     * <tt>ArithmeticException</tt> is thrown.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>short</tt>.
     *
     * @throws ArithmeticException if <tt>this</tt> has a nonzero fractional
     * part, or will not fit in a <tt>short</tt>.
     * @since 1.5
     */
    public short shortValueExact() {
        long num;
        num
                = this.longValueExact();     // will check decimal part
        if ((short) num != num) {
            throw new java.lang.ArithmeticException("Overflow");
        }

        return (short) num;
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>byte</tt>, checking for lost
     * information. If this <tt>BigDecimal</tt> has a nonzero fractional part or
     * is out of the possible range for a <tt>byte</tt> result then an
     * <tt>ArithmeticException</tt> is thrown.
     *
     * @return this <tt>BigDecimal</tt> converted to an <tt>byte</tt>.
     *
     * @throws ArithmeticException if <tt>this</tt> has a nonzero fractional
     * part, or will not fit in a <tt>byte</tt>.
     * @since 1.5
     */
    public byte byteValueExact() {
        long num;
        num
                = this.longValueExact();     // will check decimal part
        if ((byte) num != num) {
            throw new java.lang.ArithmeticException("Overflow");
        }

        return (byte) num;
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>float</tt>. This conversion is
     * similar to the <a
     * href="http://java.sun.com/docs/books/jls/second_edition/html/conversions.doc.html#25363"><i>narrowing
     * primitive conversion</i></a> from <tt>double</tt> to <tt>float</tt>
     * defined in the <a href="http://java.sun.com/docs/books/jls/html/">Java
     * Language Specification</a>: if this <tt>BigDecimal</tt> has too great a
     * magnitude to represent as a <tt>float</tt>, it will be converted to
     * {@link Float#NEGATIVE_INFINITY} or {@link
     * Float#POSITIVE_INFINITY} as appropriate. Note that even when the return
     * value is finite, this conversion can lose information about the precision
     * of the <tt>BigDecimal</tt> value.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>float</tt>.
     */
    public float floatValue() {
        /*
         * Somewhat inefficient, but guaranteed to work.
         */
        return Float.valueOf(this.toString()).floatValue();
    }

    /**
     * Converts this <tt>BigDecimal</tt> to a <tt>double</tt>. This conversion
     * is similar to the <a
     * href="http://java.sun.com/docs/books/jls/second_edition/html/conversions.doc.html#25363"><i>narrowing
     * primitive conversion</i></a> from <tt>double</tt> to <tt>float</tt> as
     * defined in the <a href="http://java.sun.com/docs/books/jls/html/">Java
     * Language Specification</a>: if this <tt>BigDecimal</tt> has too great a
     * magnitude represent as a <tt>double</tt>, it will be converted to
     * {@link com.mathpar.number.Double#NEGATIVE_INFINITY} or {@link
     * com.mathpar.number.Double#POSITIVE_INFINITY} as appropriate. Note that even when the
     * return value is finite, this conversion can lose information about the
     * precision of the <tt>BigDecimal</tt> value.
     *
     * @return this <tt>BigDecimal</tt> converted to a <tt>double</tt>.
     */
    @Override
    public double doubleValue() {
        /*
         * Somewhat inefficient, but guaranteed to work.
         */

        return NumberR64.valueOf(this.toString()).doubleValue();
    }

    /**
     * Returns the size of an ulp, a unit in the last place, of this
     * <tt>BigDecimal</tt>. An ulp of a nonzero <tt>BigDecimal</tt> value is the
     * positive distance between this value and the <tt>BigDecimal</tt> value
     * next larger in magnitude with the same number of digits. An ulp of a zero
     * value is numerically equal to 1 with the scale of <tt>this</tt>. The
     * result is stored with the same scale as <code>this</code> so the result
     * for zero and nonzero values is equal to      <code>[1,
     * this.scale()]</code>.
     *
     * @return the size of an ulp of <tt>this</tt>
     *
     * @since 1.5
     */
    public NumberR ulp() {
        return new NumberR(NumberZ.ONE, this.scale());
    }

// Private "Helper" Methods
    /**
     * Lay out this <tt>BigDecimal</tt> into a <tt>char[]</tt> array. The Java
     * 1.2 equivalent to this was called <tt>getValueString</tt>.
     *
     * @param sci <tt>true</tt> for Scientific exponential notation;
     * <tt>false</tt> for Engineering
     *
     * @return <tt>char[]</tt> array with canonical string representation of
     * this <tt>BigDecimal</tt>
     */
    private char[] layoutChars(boolean sci) {

        if (scale == 0) // zero scale is trivial
        {
            return intVal.toString().toCharArray();
        }
        // get the significand as an absolute value

        char coeff[] = intVal.abs().toString().toCharArray();
        // Construct a buffer, with sufficient capacity for all cases.
        // If E-notation is needed, length will be: +1 if negative, +1
        // if '.' needed, +2 for "E+", + up to 10 for adjusted exponent.
        // Otherwise it could have +1 if negative, plus leading "0.00000"
        StringBuilder buf = new StringBuilder(coeff.length + 14);
        if (intVal.signum() < 0) // prefix '-' if negative
        {
            buf.append('-');
        }

        long adjusted = -(long) scale + (coeff.length - 1);
        if ((scale >= 0) && (adjusted >= -6)) { // plain number
            int pad = scale - coeff.length;  // count of padding zeros
            if (pad >= 0) {                  // 0.*** form
                buf.append('0');
                buf.append('.');
                for (; pad
                        > 0; pad--) {
                    buf.append('0');
                }

                buf.append(coeff);
            } else {                         // **.** form
                buf.append(coeff, 0, -pad);
                buf.append('.');
                buf.append(coeff, -pad, scale);
            }

        } else { // E-notation is needed
            if (sci) {                       // Scientific notation
                buf.append(coeff[0]);        // first character
                if (coeff.length > 1) {      // more to come
                    buf.append('.');
                    buf.append(coeff, 1, coeff.length - 1);
                }

            } else {                         // Engineering notation
                int sig = (int) (adjusted % 3);
                if (sig < 0) {
                    sig += 3;                // [adjusted was negative]
                }

                adjusted -= sig;             // now a multiple of 3
                sig++;

                if (intVal.signum() == 0) {
                    switch (sig) {
                        case 1:
                            buf.append('0'); // exponent is a multiple of three
                            break;
                        case 2:
                            buf.append("0.00");
                            adjusted += 3;
                            break;
                        case 3:
                            buf.append("0.0");
                            adjusted += 3;
                            break;
                        default:
                            throw new AssertionError("Unexpected sig value " + sig);
                    }
                } else if (sig >= coeff.length) {   // significand all in integer
                    buf.append(coeff, 0, coeff.length);
                    // may need some zeros, too
                    for (int i = sig - coeff.length; i > 0; i--) {
                        buf.append('0');
                    }
                } else {                     // xx.xxE form
                    buf.append(coeff, 0, sig);
                    buf.append('.');
                    buf.append(coeff, sig, coeff.length - sig);
                }
            }
            if (adjusted != 0) {             // [!sci could have made 0]
                buf.append('E');
                if (adjusted > 0) // force sign for positive
                {
                    buf.append('+');
                }

                buf.append(adjusted);
            }

        }
        char result[] = new char[buf.length()];
        buf.getChars(0, result.length, result, 0);
        return result;
    }

    /**
     * Return 10 to the power n, as a <tt>BigInteger</tt>.
     *
     * @param n the power of ten to be returned (>=0)
     *
     * @return a <tt>BigInteger</tt> with the value (10<sup>n</sup>)
     */
    private static NumberZ tenToThe(int n) {
        if (n < TENPOWERS.length) // use value from constant array
        {
            return TENPOWERS[n];
        }
// BigInteger.pow is slow, so make 10**n by constructing a
// BigInteger from a character string (still not very fast)
        java.math.BigDecimal bb;
        char tenpow[] = new char[n + 1];
        tenpow[0] = '1';
        for (int i = 1; i<= n; i++) {
            tenpow[i] = '0';
        }
        return new NumberZ(tenpow);
    }
    private static NumberZ TENPOWERS[] = {NumberZ.ONE,
        NumberZ.valueOf(10), NumberZ.valueOf(100),
        NumberZ.valueOf(1000), NumberZ.valueOf(10000),
        NumberZ.valueOf(100000), NumberZ.valueOf(1000000),
        NumberZ.valueOf(10000000), NumberZ.valueOf(100000000),
        NumberZ.valueOf(1000000000)};

    /**
     * Match the scales of two <tt>BigDecimal<tt>s to align their least
     * significant digits.
     *
     * <p>
     * If the scales of val[0] and val[1] differ, rescale (non-destructively)
     * the lower-scaled <tt>BigDecimal</tt> so they match. That is, the
     * lower-scaled reference will be replaced by a reference to a new object
     * with the same scale as the other <tt>BigDecimal</tt>.
     *
     * @param val array of two elements referring to the two
     * <tt>BigDecimal</tt>s to be aligned.
     */
    private static void matchScale(NumberR[] val) {
        if (val[0].scale < val[1].scale) {
            val[0] = val[0].setScale(val[1].scale);
        } else if (val[1].scale < val[0].scale) {
            val[1] = val[1].setScale(val[0].scale);
        }

    }

    /**
     * Reconstitute the <tt>BigDecimal</tt> instance from a stream (that is,
     * deserialize it).
     *
     * @param s the stream being read.
     */
    private synchronized void readObject(java.io.ObjectInputStream s)
            throws java.io.IOException, ClassNotFoundException {
        // Read in all fields
        s.defaultReadObject();
        // validate possibly bad fields
        if (intVal == null) {
            String message = "BigDecimal: null intVal in stream";
            throw new java.io.StreamCorruptedException(message);
            // [all values of scale are now allowed]
        }

    }

    /**
     * Returns the length of this <tt>BigDecimal</tt>, in decimal digits.
     *
     * Notes: <ul> <li> This is performance-critical; most operations where a
     * context is supplied will need at least one call to this method.
     *
     * <li> This should be a method on BigInteger; the call to this method in
     * precision() can then be replaced with the term: intVal.digitLength(). It
     * could also be called precision() in BigInteger.
     *
     * Better still -- the precision lookaside could be moved to BigInteger,
     * too.
     *
     * <li> This could/should use MutableBigIntegers directly for the reduction
     * loop. <ul>
     *
     * @return the length of the unscaled value, in decimal digits
     */
    public int digitLength() {
        if (intVal.signum() == 0) // 0 is one decimal digit
        {
            return 1;
        }
// we have a nonzero magnitude

        NumberZ work = intVal;
        int digits = 0;                 // counter
        for (; work.mag.length > 1;) {
            // here when more than one integer in the magnitude; divide
            // by a billion (reduce by 9 digits) and try again
            work = work.divide(TENPOWERS[9]);
            digits
                    += 9;
            if (work.signum() == 0) // the division was exact
            {
                return digits;          // (a power of a billion)
            }

        }
        // down to a simple nonzero integer
        digits += intLength(work.mag[0]);

        return digits;
    }

    /**
     * Return the number of digits for THIS real if it grater or equal 1. Return
     * the number of zeros which situated after decimal point with negative sign
     * if THIS less then 1. .
     *
     * @return number of digits
     */
    public int numberOfDigits() {
        return (digitLength() - scale);
    }

    /**
     * Returns the length of an unsigned <tt>int</tt>, in decimal digits.
     *
     * @param i the <tt>int</tt> (treated as unsigned)
     *
     * @return the length of the unscaled value, in decimal digits
     */
    private int intLength(int i) {
        int digits;
        if (i < 0) {            // 'negative' is 10 digits unsigned
            digits = 10;
        } else {                // positive integer
            // binary search, weighted low (maximum 4 tests)
            if (i < 10000) {
                if (i < 100) {
                    if (i < 10) {
                        digits = 1;
                    } else {
                        digits = 2;
                    }
                } else {
                    if (i < 1000) {
                        digits = 3;
                    } else {
                        digits = 4;
                    }
                }
            } else {
                if (i < 1000000) {
                    if (i < 100000) {
                        digits = 5;
                    } else {
                        digits = 6;
                    }
                } else {
                    if (i < 100000000) {
                        if (i < 10000000) {
                            digits = 7;
                        } else {
                            digits = 8;
                        }
                    } else {
                        if (i < 1000000000) {
                            digits = 9;
                        } else {
                            digits = 10;
                        }
                    }
                }
            }
        }
        return digits;
    }

    /**
     * Remove insignificant trailing zeros from this <tt>BigDecimal</tt> until
     * the preferred scale is reached or no more zeros can be removed. If the
     * preferred scale is less than Integer.MIN_VALUE, all the trailing zeros
     * will be removed.
     *
     * <tt>BigInteger</tt> assistance could help, here?
     *
     * @return this <tt>BigDecimal</tt> with a scale possibly reduced to be
     * closed to the preferred scale.
     */
    private NumberR stripZerosToMatchScale(long preferredScale) {
        NumberZ qr[];                // quotient-remainder pair
        while (intVal.abs().compareTo(NumberZ.POSCONST[10]) >= 0
                && scale > preferredScale) {
            if (intVal.testBit(0)) {
                break;                  // odd number cannot end in 0
            }
            qr = intVal.divideAndRemainder(NumberZ.POSCONST[10]);
            if (qr[1].signum() != 0) {
                break;                  // non-0 remainder
            }
            intVal = qr[0];
            scale = checkScale((long) scale - 1);  // could Overflow
            if (precision > 0) // adjust precision if known
            {
                precision--;
            }
        }
        return this;
    }

    /**
     * Check a scale for Underflow or Overflow
     *
     * @param val The new scale.
     *
     * @throws ArithmeticException (overflow or underflow) if the new scale is
     * out of range.
     * @return validated scale as an int.
     */
    private int checkScale(long val) {
        if ((int) val != val) {
            if ((this.intVal != null && this.signum() != 0)
                    || this.intVal == null) {
                if (val > Integer.MAX_VALUE) {
                    throw new ArithmeticException("Underflow");
                }
                if (val < Integer.MIN_VALUE) {
                    throw new ArithmeticException("Overflow");
                }
            } else {
                return (val > Integer.MAX_VALUE) ? Integer.MAX_VALUE : Integer.MIN_VALUE;
            }
        }
        return (int) val;
    }

    /**
     * Round an operand; used only if digits &gt; 0. Does not change
     * <tt>this</tt>; if rounding is needed a new <tt>BigDecimal</tt> is created
     * and returned.
     *
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the result is inexact but the rounding
     * mode is <tt>UNNECESSARY</tt>.
     */
    private NumberR roundOp(MathContext mc) {
        NumberR rounded = doRound(mc);
        return rounded;
    }

    /**
     * Round this BigDecimal according to the MathContext settings; used only if
     * precision &gt; 0.
     *
     * @param mc the context to use.
     *
     * @throws ArithmeticException if the rounding mode is
     * <tt>RoundingMode.UNNECESSARY</tt> and the <tt>BigDecimal</tt> operation
     * would require rounding.
     */
    private void roundThis(MathContext mc) {
        NumberR rounded = doRound(mc);
        if (rounded == this) // wasn't rounded
        {
            return;
        }
        this.intVal = rounded.intVal;
        this.scale = rounded.scale;
        this.precision = rounded.precision;
    }

    /**
     * Returns a <tt>BigDecimal</tt> rounded according to the MathContext
     * settings; used only if <tt>mc.precision&gt;0</tt>. Does not change
     * <tt>this</tt>; if rounding is needed a new <tt>BigDecimal</tt> is created
     * and returned.
     *
     * @param mc the context to use.
     *
     * @return a <tt>BigDecimal</tt> rounded according to the MathContext
     * settings. May return this, if no rounding needed.
     *
     * @throws ArithmeticException if the rounding mode is
     * <tt>RoundingMode.UNNECESSARY</tt> and the result is inexact.
     */
    private NumberR doRound(MathContext mc) {
        if (precision == 0) {
            if (mc.roundingMax != null && intVal.compareTo(mc.roundingMax) < 0 && intVal.
                    compareTo(mc.roundingMin) > 0) {
                return this; // no rounding needed
            }
            precision();                     // find it
        }
        int drop = precision - mc.getPrecision();   // digits to discard
        if (drop <= 0) // we fit
        {
            return this;
        }
        NumberR rounded = dropDigits(mc, drop);
        // we need to double-check, in case of the 999=>1000 case
        return rounded.doRound(mc);
    }

    /**
     * Removes digits from the significand of a <tt>BigDecimal</tt>, rounding
     * according to the MathContext settings. Does not change <tt>this</tt>; a
     * new <tt>BigDecimal</tt> is always created and returned.
     *
     * <p>
     * Actual rounding is carried out, as before, by the divide method, as this
     * minimized code changes. It might be more efficient in most cases to move
     * rounding to here, so we can do a round-to-length rather than
     * round-to-scale.
     *
     * @param mc the context to use.
     * @param drop the number of digits to drop, must be &gt; 0
     *
     * @return a <tt>BigDecimal</tt> rounded according to the MathContext
     * settings. May return <tt>this</tt>, if no rounding needed.
     *
     * @throws ArithmeticException if the rounding mode is
     * <tt>RoundingMode.UNNECESSARY</tt> and the result is inexact.
     */
    private NumberR dropDigits(MathContext mc, int drop) {
        // here if we need to round; make the divisor = 10**drop)
        // [calculating the BigInteger here saves setScale later]
        NumberR divisor = new NumberR(tenToThe(drop), 0);
        // divide to same scale to force round to length
        NumberR rounded = this.divide(divisor, scale,
                mc.getRoundingMode().oldMode);
        rounded.scale -= drop;               // adjust the scale
        return rounded;
    }

    @Override
    public Element myOne(Ring ring) {
        return ONE;
    }

    @Override
    public Element myMinus_one(Ring ring) {
        return MINUS_ONE;
    }

    @Override
    public Element myMinus_one() {
        return MINUS_ONE;
    }

    @Override
    public Element myZero(Ring ring) {
        return ZERO;
    }

    @Override
    public Element myZero() {
        return ZERO;
    }

    @Override
    public Element one(Ring ring) {
        return ONE;
    }

    @Override
    public Element minus_one(Ring ring) {
        return MINUS_ONE;
    }

    @Override
    public Element zero(Ring ring) {
        return ZERO;
    }

    @Override
    public Element rootOf(int pow, Ring ring) {
        return new FuncNumberR(ring).rootOf(this, new NumberZ(pow));
    }

    @Override
    public Element rootTheFirst(int pow, Ring ring) {
        return new FuncNumberR(ring).rootOf(this, new NumberZ(pow));
    }

    @Override
    public int compareTo(Element val) {
        if ((val == NAN) || (this == NAN)) {
            return Integer.MAX_VALUE;
        }
        if (val == POSITIVE_INFINITY) {
            return (this == POSITIVE_INFINITY) ? 0 : -1;
        }
        if (this == POSITIVE_INFINITY) {
            return 1;
        }
        if (val == NEGATIVE_INFINITY) {
            return (this == NEGATIVE_INFINITY) ? 0 : 1;
        }
        if (this == NEGATIVE_INFINITY) {
            return -1;
        }
        int v = val.numbElementType();
        int m = Ring.R;
        if (m > v) {
            return 1;
        } else if (m < v) {
            return -1;
        } else {
            return compareTo((NumberR) val);
        }
    }
    @Override
    public int compareTo(Element x, Ring ring) {Element xx=x;
        if (x instanceof Element){
          if (x == NAN) return Integer.MAX_VALUE; 
          else if (x == POSITIVE_INFINITY) return (this == POSITIVE_INFINITY) ? 0 : -1; 
          else if (x == NEGATIVE_INFINITY) return (this == NEGATIVE_INFINITY) ? 0 : 1;
        }
        if(x instanceof Complex){Complex cc=(Complex)x;
          if(cc.im.isZero(ring)) xx=cc.re; else return -x.compareTo(this);}
        if (xx.isInfinite()) {return -xx.compareTo(this, ring);}  ////  ??  
        int xType = xx.numbElementType();
        if(xType<=Ring.R64MinMax){ NumberR yy=(NumberR)xx.toNumber(Ring.R, ring);
           yy = (NumberR) subtract(yy, ring);
           NumberR yAbs = yy.abs();
           if (yAbs.compareTo(ring.MachineEpsilonR) < 0) return 0;
           return yy.signum();
        }
        return -xx.compareTo(this);
    }
    
    
//    public int compareTo(Element val, Ring ring) {
//        if (val == NAN || this == NAN) {
//            return Integer.MAX_VALUE;
//        }
//        if (val == POSITIVE_INFINITY) {
//            return (this == POSITIVE_INFINITY) ? 0 : -1;
//        }
//        if (this == POSITIVE_INFINITY) {
//            return 1;
//        }
//        if (val == NEGATIVE_INFINITY) {
//            return (this == NEGATIVE_INFINITY) ? 0 : 1;
//        }
//        if (this == NEGATIVE_INFINITY) {
//            return -1;
//        }
//        NumberR x = (NumberR) subtract(val, ring);
//        NumberR y = x.abs();
//        if (y.compareTo(ring.MachineEpsilonR) < 0) {
//            return 0;
//        }
//        return x.signum();
//    }

    public NumberR BigDecimalValue() {
        return (this);
    }

    public NumberZ BigIntegerValue() {
        return this.NumberRtoNumberZ();
    }

    //  public NumberR abs(Element x) { return abs((NumberR)(x));}
    @Override
    public Element add(Element x, Ring ring) {
        if (x instanceof NumberR) {return add((NumberR) x, ring.MC);}
        if ((x == NAN) || (x == NEGATIVE_INFINITY) || (x == POSITIVE_INFINITY)) {return x;}
        if (x.numbElementType() > com.mathpar.number.Ring.R) {return x.add(this, ring);}
        return add(((NumberR) (x.toNumber(Ring.R, ring))), ring.MC);
    }

    @Override
    public int numbElementType() {
        return Ring.R;
    }

    @Override
    public Element subtract(Element x, Ring ring) {
        if (x instanceof NumberR)   return subtract((NumberR) x, ring.MC);
        if (x == NAN) {return NAN; }
        if ((x == NEGATIVE_INFINITY) || (x == POSITIVE_INFINITY)) return x.negate(ring);
        int type = x.numbElementType();
        switch (type) {
            case com.mathpar.number.Ring.Polynom:
                return x.negate(ring).add(this, ring);
            case com.mathpar.number.Ring.F:
                return new F(com.mathpar.func.F.SUBTRACT, new Element[] {this, x});
            default:
                if (type < numbElementType()) {
                    return subtract((NumberR) (x.toNumber(Ring.R, ring)), ring.MC);
                } else {
                    return (this.toNumber(type, ring)).subtract(x, ring);
                }

        }
    }

    @Override
    public Element multiply(Element x, Ring ring) {
        if (x instanceof NumberR) {return multiply((NumberR) x, ring.MC);}
        if (x == NAN) {return NAN;}
        if ((x == NEGATIVE_INFINITY) || (x == POSITIVE_INFINITY)) {
            return (this.isZero(ring)) ? NAN : this.isNegative() ? x.inverse(ring) : x;}
        if (x.numbElementType() > com.mathpar.number.Ring.R) {
            return x.multiply(this, ring);}
        Element ee=x.toNumber(Ring.R, ring);

        //  System.out.println("&=");

               return multiply( ((NumberR) ee), ring.MC);
    }

    @Override
    public Element multiply(Element x, int i, Ring ring) {
        return multiply((NumberR) x, ring);
    }

    @Override
    public Element inverse(Ring ring) {
        return ONE.divide(this, ring);
    }
    @Override
    public Element divideToFraction(Element x, Ring ring) {
        return divide(x,ring);
    }
    @Override
    public Element divide(Element x, Ring ring) {if(x.isOne(ring))return this;
        if (x instanceof NumberR) {return divide((NumberR) x, ring.MC);}
        if (x == NAN) return NAN;
        if ((x == NEGATIVE_INFINITY) || (x == POSITIVE_INFINITY)) return ZERO;
        int type = x.numbElementType();
        switch (type) {
            case com.mathpar.number.Ring.Polynom:
            case com.mathpar.number.Ring.F:
                return new F(com.mathpar.func.F.DIVIDE, new Element[] {this, x});
            default:
                if (type < numbElementType()) {
                    return divide((NumberR) (x.toNumber(Ring.R, ring)), ring.MC);
                } else  return (this.toNumber(type, ring)).divide(x, ring);
        }
    }
/**
 * divideExact (Element, ring) is a standard divide in NumberR
 *    (Do not mixed it with real devideExact 
 *        of Bigdecimal == divideExact(NumberR). -- it is exact!)
 * @param x devisor
 * @param ring Ring
 * @return 
 */
    public Element divideExact(Element x, Ring ring) {
          if (x instanceof Complex) 
            return  new Complex(this,ring).divideExact(x,ring); 
        switch (x.numbElementType()) {
            case com.mathpar.number.Ring.Polynom:
            case com.mathpar.number.Ring.F:
                return new F(com.mathpar.func.F.DIVIDE, new Element[] {this, x});
            default:
                return divide((NumberR) (x.toNumber(Ring.R, ring)),ring.MC);
        }
    }

    public Element valOf(double x) {
        return (NumberR.valueOf(x));
    }

    @Override
    public NumberR valOf(long x, Ring ring) {
        return valueOf(x);
    }

    @Override
    public Element valOf(double x, Ring ring) {
        return (NumberR.valueOf(x));
    }

    @Override
    public NumberR valOf(String s, Ring ring) {
        return new NumberR(s);
    }

    /**
     * Transform this number of NumberR type to the new type fixed by
     * "numberType" defined by Ring.type
     *
     * @param numberType - new type
     *
     * @return this transormed to the new type
     */
    @Override
    public Element toNumber(int numberType, Ring ring) {
        if (numberType < Ring.Complex) {
            switch (numberType) {
                case Ring.Z:
                    return NumberRtoNumberZ();
                case Ring.Q:
                    return NumberRtoNumberQ().cancel(ring);
                case Ring.Z64:
                    return new NumberZ64(longValue());
                case Ring.Zp32:
                    return new NumberZp32(longValue(), ring);
                case Ring.Zp:
                    return new NumberZp(NumberRtoNumberZ(), ring);
                case Ring.R:
                    return this;
                case Ring.R64:
                    return new NumberR64(doubleValue());
                case Ring.R128:
                    return new NumberR128(intVal, scale);
                case Ring.RMaxPlus:
                    return new NumberRMaxPlus(this);
                case Ring.RMinPlus:
                    return new NumberRMinPlus(this);
                case Ring.RMaxMult:
                    return new NumberRMaxMult(this);
                case Ring.RMinMult:
                    return new NumberRMinMult(this);
                case Ring.RMaxMin:
                    return new NumberRMaxMin(this);
                case Ring.RMinMax:
                    return new NumberRMinMax(this);
                case Ring.R64MaxPlus:
                    return new NumberR64MaxPlus(new NumberR64(doubleValue()));
                case Ring.R64MinPlus:
                    return new NumberR64MinPlus(new NumberR64(doubleValue()));
                case Ring.R64MaxMult:
                    return new NumberR64MaxMult(new NumberR64(doubleValue()));
                case Ring.R64MinMult:
                    return new NumberR64MinMult(new NumberR64(doubleValue()));
                case Ring.R64MaxMin:
                    return new NumberR64MaxMin(new NumberR64(doubleValue()));
                case Ring.R64MinMax:
                    return new NumberR64MinMax(new NumberR64(doubleValue()));
                case Ring.ZMaxPlus:
                    return new NumberZMaxPlus(new NumberZ(longValue()));
                case Ring.ZMinPlus:
                    return new NumberZMinPlus(new NumberZ(longValue()));
                case Ring.ZMaxMult:
                    return new NumberZMaxMult(new NumberZ(longValue()));
                case Ring.ZMinMult:
                    return new NumberZMinMult(new NumberZ(longValue()));
                case Ring.ZMaxMin:
                    return new NumberZMaxMin(new NumberZ(longValue()));
                case Ring.ZMinMax:
                    return new NumberZMinMax(new NumberZ(longValue()));
                default:
                    return new Complex(this, this.zero(ring));
            }
        } else if (numberType < Ring.Polynom) {
            Element re = toNumber(numberType - Ring.Complex, ring);
            return new Complex(this, this.zero(ring));
        }

        return null;
    }

    @Override
    public NumberR random(int[] randomType, Random rnd, Ring ring) {
        int nbits = randomType[randomType.length - 1];
        NumberZ ZZ = new NumberZ(nbits, rnd);
        return new NumberR(ZZ);
    }

    @Override
    public Element abs(Ring ring) {
        return abs();
    }

    @Override
    public Element D(Ring ring) {
        return ZERO;
    }

    @Override
    public Element D(int num, Ring ring) {
        return ZERO;
    }

    @Override
    public Element id(Ring ring) {
        return this;
    }

    @Override
    public Element ln(Ring r) {
        return (new FuncNumberR(r)).ln(this);
    }

    @Override
    public Element lg(Ring r) {
        return (NumberR) (((NumberR) (new FuncNumberR(r)).ln(this)).divide( (NumberR) (new FuncNumberR(r)).ln(POSCONST[10]), r.MC));
    }

    @Override
    public Element exp(Ring r) {
        return (new FuncNumberR(r)).exp(this);
    }

    @Override
    public Element sqrt(Ring r) {


       NumberR res = this.negate(r);

       return (this.isNegative())?
               (r.isComplex())? new Complex(r.numberZERO,(new FuncNumberR(r)).sqrt(res)): NAN
               : (new FuncNumberR(r)).sqrt(this);
    }

    @Override
    public Element sin(Ring r) {

        return (new FuncNumberR(r)).sin(this);
    }

    @Override
    public Element cos(Ring r) {
        return (new FuncNumberR(r)).cos(this);
    }

    @Override
    public Element tan(Ring r) {
        return (new FuncNumberR(r)).tan(this);
    }

    @Override
    public Element ctg(Ring r) {
        return (new FuncNumberR(r)).ctg(this);
    }

    
    @Override
    public Element arcsn(Ring r) {
        return (new FuncNumberR(r)).arcsn(this);
    }

    @Override
    public Element arccs(Ring r) {
        return (new FuncNumberR(r)).arccs(this);
    }

    @Override
    public Element arctn(Ring r) {
        return (new FuncNumberR(r)).arctn(this);
    }

    @Override
    public Element arcctn(Ring r) {
        return (new FuncNumberR(r)).arcctn(this);
    }

    @Override
    public Element sh(Ring r) {
        return (new FuncNumberR(r)).sh(this);
    }

    @Override
    public Element ch(Ring r) {
        return (new FuncNumberR(r)).ch(this);
    }

    @Override
    public Element th(Ring r) {
        return (new FuncNumberR(r)).th(this);
    }

    @Override
    public Element cth(Ring r) {
        return (new FuncNumberR(r)).cth(this);
    }

    @Override
    public Element arsh(Ring r) {
        return (new FuncNumberR(r)).arsh(this);
    }

    @Override
    public Element arch(Ring r) {
        return (new FuncNumberR(r)).arch(this);
    }

    @Override
    public Element arth(Ring r) {
        return (new FuncNumberR(r)).arth(this);
    }

    @Override
    public Element arcth(Ring r) {
        return (new FuncNumberR(r)).arcth(this);
    }

    /**
     * Unit step in the place 0
     *
     * @param a
     *
     * @return 1 if x>0; 0 if x<0; 1/2 if x=0.
     */
    public NumberR unitstep() {
        int s = signum();
        return (s == 1) ? NumberR.ONE : (s == -1) ? NumberR.ZERO : new NumberR(0.5);
    }

    /**
     * Unit step in the place a
     *
     * @param a
     *
     * @return 1 if x>a; 0 if x<a; 1/2 if x=a.
     */
    public NumberR unitstep(NumberR a) {
        int s = (this.subtract(a)).signum();
        return (s == 1) ? NumberR.ONE : (s == -1) ? NumberR.ZERO : new NumberR(0.5);
    }

    /**
     * Unit step in the place a
     *
     * @param a
     *
     * @return 1 if x>a; 0 if x<a; 1/2 if x=a.
     */
    public NumberR unitstep(Element a) {
        int s = (this.subtract((NumberR) a)).signum();
        return (s == 1) ? NumberR.ONE : (s == -1) ? NumberR.ZERO : new NumberR(0.5);
    }

    @Override
    public Element toNewRing(int Algebra, Ring r) {
        return toNumber(Algebra, r);
    }

    @Override
    public Element GCD(Element x, Ring r) {
        return ONE;
    }

    @Override
    public Element pi(Ring r) {
        return new FuncNumberR(r).pi();
    }

    @Override
    public Element e(Ring r) {
        return new FuncNumberR(r).exp1(ONE);
    }

    @Override
    public Element integrate(int num, Ring ring) {
        return ring.varPolynom[num].multiply(this, ring);
    }

    @Override
    public Element factorial(Ring r) { NumberZ xx=floor(r);
       // if (isZero(r) || isOne(r)) return ONE;
        return xx.factorial(r); //multiply(subtract(ONE, r).factorial(r), r);
    }

    @Override
    public boolean isEven() {
        return scale >= 0 && intVal.isEven();
    }
    public static void main(String[] args) {
        Ring r=new Ring("R[x]"); r.FLOATPOS=3;
        NumberR a=(NumberR)r.posConst[2];
        NumberR b=(NumberR)r.posConst[3];
        Element c= a.divide(b, r);
        System.out.println("  ="+c.toString(r));
    }
}
