package com.mathpar.web;

import com.codeborne.selenide.Configuration;
import com.codeborne.selenide.junit.ScreenShooter;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;

import static com.codeborne.selenide.Condition.*;
import static com.codeborne.selenide.Selenide.$;
import static com.codeborne.selenide.Selenide.open;

public class WebUiBasicIT {
    @Rule // automatically takes screenshot of every failed test
    public ScreenShooter makeScreenshotOnFailure = ScreenShooter.failedTests();

    @BeforeClass
    public static void beforeClass() {
        Configuration.baseUrl = "http://localhost:9087";
    }

    @Before
    public void before() {
        open("/mathpar/en/");
    }

    @Test
    public void sidebarClickExpandsFirstLevelPanel() {
        $("#kbd_space_const").shouldHave(cssClass("collapse"));
        $("button[data-target='#kbd_space_const']").click();
        $("#kbd_space_const").shouldBe(visible);
    }

    @Test
    public void sidebarClickExpandsSecondLevelPanel() {
        $("button[data-target='#kbd_space_const']").click();

        $("#kbd_space_const_const").shouldHave(cssClass("collapse"));
        $("*[data-target='#kbd_space_const_const']").click();
        $("#kbd_space_const_const").shouldBe(visible);
    }

    @Test
    public void sidebarButtonsInsertText() {
        $("button[data-target='#kbd_space_const']").click();

        $("*[data-inserts='SPACE = Z[];']", 0).click();
        $("#section_0 textarea").shouldHave(value("SPACE = Z[];"));
    }
}
