package com.mathpar.func.parser;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 * Parser test suite.
 *
 * @author ivan
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({
    com.mathpar.func.parser.FunctionTest.class,
    com.mathpar.func.parser.FunctionSpecTest.class,
    com.mathpar.func.parser.PolynomialAndNumberTest.class,
    com.mathpar.func.parser.SymbolTest.class})
public class ParserTest {
}
