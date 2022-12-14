package com.mathpar.web;

import static org.apache.logging.log4j.LogManager.getLogger;

import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.EnableAutoConfiguration;
import org.springframework.boot.autoconfigure.freemarker.FreeMarkerAutoConfiguration;
import org.springframework.boot.autoconfigure.jackson.JacksonAutoConfiguration;
import org.springframework.boot.autoconfigure.orm.jpa.HibernateJpaAutoConfiguration;
import org.springframework.boot.autoconfigure.web.servlet.MultipartAutoConfiguration;
import org.springframework.boot.autoconfigure.web.servlet.WebMvcAutoConfiguration;
//import org.springframework.boot.autoconfigure.web.HttpMessageConvertersAutoConfiguration;
//import org.springframework.boot.autoconfigure.web.MultipartAutoConfiguration;
//import org.springframework.boot.autoconfigure.web.WebMvcAutoConfiguration;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.FilterType;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcOperations;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RestController;

@Configuration
@EnableAutoConfiguration(exclude = {
        WebMvcAutoConfiguration.class,
        MultipartAutoConfiguration.class,
        HibernateJpaAutoConfiguration.class,
//        HttpMessageConvertersAutoConfiguration.class,
        JacksonAutoConfiguration.class,
        FreeMarkerAutoConfiguration.class
})
@ComponentScan(
        basePackages = {"com.mathpar.web"},
        excludeFilters = {
                @ComponentScan.Filter(
                        value = {Application.class, MvcConfig.class},
                        type = FilterType.ASSIGNABLE_TYPE),
                @ComponentScan.Filter(
                        value = {Controller.class, RestController.class},
                        type = FilterType.ANNOTATION),
        })
public class ApplicationPlayground implements CommandLineRunner {
    private static final Logger LOG = getLogger(ApplicationPlayground.class);

    @Autowired
    NamedParameterJdbcOperations jdbcTpl;

    public static void main(String[] args) throws Exception {
        SpringApplication.run(ApplicationPlayground.class, args);
    }
Object[] newTit={
"156", "RU-SP-B06-01 ????????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????.",
"298", "RU-SP-B07-01 ?????????? 1 ????????????????. ???????????????????? ????????????????????????",
"377", "0111",
"181", "12",
"160", "132",
"118", "89",
"7", "EN-SM Calculations and conversions. Transformation of algebraic expressions and fractions.",
"217", "RU-SM-A13-01 ????????????????????????. ??????",
"13", "RU-SM-A13-01 ????????????????????????. ??????.",
"135", "RU-SM-A13-03 ????????????????????????. ????????????.",
"133", "RU-SM-A13-04 ????????????????????????. ????????????????.",
"18", "RU-SM-A13-05 ????????????????????????. ??????????????.",
"144", "RU-SM-A13-06 ????????????????????????. ??????.",
"199", "RU-SM-B20-01 ???????????? ???? ????????????????",
"27", "RU-SM-B01-02 ???????????????????? ?????????????????? ??????????????. ???????????????????? ?? ????????????????",
"15", "RU-SM-B10-01 ???????????????????? ?? ???????????????????????????? ???????????????????????????? ???????????????? ???????????????????????? ??????????????????",
"16", "RU-SM-B10-02 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????????????????? ?????????????????? ?? ????????????",
"20", "RU-SM-B10-03 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????? ???????????????????????????? ??????????????????",
"61", "RU-SM-B10-04 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ?????????????????? ???????????????????????????? ??????????????????",
"64", "RU-SM-B10-05 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????? ?????????????????????????? ??????????????????",
"162", "RU-SM-B10-06 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ?????????????????? ?????????????????????????? ??????????????????",
"167", "RU-SM-B10-07 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????? ?????????????????????????????? ??????????????????",
"171", "RU-SM-B10-08 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ?????????????????? ?????????????????????????????? ??????????????????",
"178", "RU-SM-B10-09 ???????????????????? ?? ????????????????????????????. ???????????????????? ???????????????? ???????????????????????????????????? ??????????????????",
"180", "RU-SM-B10-10 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????? ???????????????????????????????????? ??????????????????",
"19",  "RU-SM-B14-01 ???????????????????? ?? ???????????????????? ???????????????? ??????????????.   ???????????????????????? ?????????????? ?????? ???????????? ??????????????????????",
"201", "RU-SM-A18-01 ???????????? ??????????????????????",
"230", "RU-SM-A18-02 ???????????? ??????????????????????",
"26", "RU-SM-A03-01 ???????????????????? ?????????????????? ????????????. ????????????????, ????????????????????",
"256", "RU-SM-B06-01 ???????????????????? ??????????????????. ????????????????, ????????????????????, ???????????????????? ??????????????????",
"257", "RU-SM-B06-02 ???????????????????? ??????????????????. ???????????????????????? ??????????????????",
"259", "RU-SM-B06-03 ???????????????????? ??????????????????. ???????????????????????????? ??????????????????",
"260", "RU-SM-B06-04 ???????????????????? ??????????????????. ?????????????????????????? ??????????????????",
"261", "RU-SM-B06-05 ???????????????????? ??????????????????. ?????????????????????????????? ??????????????????",
"263", "RU-SM-B06-06 ???????????????????? ??????????????????. ???????????????????????????????????? ??????????????????",
"17",  "RU-SM-A07-01 ???????????????????? ??????????????????. ???????????? ???????????????? ??????????????????",
"9",   "EN-SM-C06 Equations with parameters",
"157", "RU-SM-D02-01 ?????????????? ????????????????????????. ?????????????? ????????.",
"14",  "RU-SP-B01-01 ????????????????????. ??????????????????",
"97",  "RU-SP-B19-01 ??????. ???????????? ??????????, ???????????? ???????????????? ????????, ?????????????? ??????????????, ???????????????????????????? ????????????",
"208", "RU-SP-B19-03 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ?????????????????????? ???????????? ??????????????????????????????",
"98",  "RU-SP-B19-03 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ?????????????????????? ???????????? ??????????????????????????????",
"159", "RU-SP-B26-01 ???????????????????????? ????????????, ??????????????????????????. T????????????????????????, ???????????????????????? ????????????  ",
"169", "RU-SP-B27-01 ??????????????????????????????. ?????????????????? ????????????.",
"293", "RU-SP-B06-01 ????????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????.",
"294", "RU-SP-B06-02 ????????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????.",
"314", "RU-SP-B07-02 ????????????????. ???????????????????????? ???????????????????????? ??? 3131",
"322", "RU-SP-D27-01 ??????????????????????????????. ?????????????????? ????????????. ?????????????????? ???????????? ",
"345", "RU-SP-D27-02 ??????????????????????????????. ?????????????????? ????????????. ????????????",
"347", "RU-SM-GA02-01 ????????????????????. ???????????????? ???? ??????????????????",
"193", "RU-SP-B10-01 ??????????????????, ?????? ???????????????? ????????????. ??????????????????",
"108", "RU-SP-B09-01 ??????????????????????, ???????????? ?? ??????????????????????????, ???????????? ?????????? ??????????????????????????. ????????????????????????????, ??????????????????????????, ??????????????????????????, ???????????????????????????? ????????????????, ?????????? ",
"327", "RU-SM-A12-01 ?????????? ???????????????????????? ????????????????. ???????????? ?????????????????? ?????? ????????????????????",
"329", "RU-SM-B12-02 ???????????? ???? ????????????????????????. ?????????????????????????? ????????????????????????????",
"299", "RU-SM-B03-01 ?????????? ???????????????????????? ????????????????. ?????????? ???????????????? ???? ???????? ??????????????????",
"174", "RU-SM-A04-01",
"173", "RU-SM-A04-01 ???????????????????????????? ??????????????????. ???????????????? ?? ??????????????????.",
"216", "RU-SM-A04-01 ???????????????????????????? ??????????????????. ???????????????? ?? ??????????????????.",
"218", "RU-SM-A04-01 ???????????????????????????? ??????????????????. ???????????????? ?? ??????????????????.",
"301", "RU-SM-B07-09 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ????????.",
"241", "RU-SM-C07-01 ?????????? ?? ???? ????????????????. ?????????? ?? ???? ????????????????",
"245", "RU-SM-C07-01 ?????????? ?? ???? ????????????????. ?????????? ?? ???? ???????????????? TEMPORARY RECORD",
"328", "RU-SM-D07-02 ?????????????? ?????????????? ???? ?????????? ?? ???? ????????????????. ???????????????? ????????????",
"290", "RU-SP-B13-01 ????????????????????????. ??????????????????, ??????????????????, ?????????????????????????? ??????????",
"291", "RU-SP-B13-02 ????????????????????????. ???????????????????????????????? ??????????",
"296", "RU-SP-B15-05 ????????????????????????????. ???????????????????? ??????. ?????? ?? ?????????????????? ????????????, ???????? ????????, ??????????????????????????",
"211", "RU-SP-B16-01 ???????????????????????????????? ????????????????. ????????????. ?????????? ??????????????, ?????? ????????????????",
"297", "RU-SP-B20-01 ??????????????????????????????. ?????????????? ??????????????. ?????????????? ??????????????????, ?????????????????????????? ????????????",
"212", "RU-SP-B21-01 ????????????. ?????????? ???????????????????????????? ??????????????. ?????????????? ????????????????????, ??????",
"283", "RU-SP-B08-10 ??????????????????????-???????????????????????? ????????????. ?????????? ???????????????????? ?? ?????????????? ???????????????????? ????????????????",
"292", "RU-SP-C02-01 ???????????????? (?????????????????? ????????????). ???????????????? ???? ????????????????????",
"317", "RU-SP-C02-02 ???????????????? (?????????????????? ????????????). ????????????????",
"326", "RU-SP-C02-03 ???????????????? (?????????????????? ????????????). ?????????? ???????????????????? ?????????????? ?? ????????????????",
"324", "RU-SP-C02-04 ???????????????? (?????????????????? ????????????). ???????????????? ??????????????????",
"122", "RU-SP-B10-01 ??????????????????, ??????????????, ?????? ???????????????? ????????????. ??????????????????.",
"192", "RU-SP-B10-02 ??????????????????, ??????????????, ?????? ???????????????? ????????????. ???????????????????? ??????????????, ???????????????????? ??????????????, ????????????????????????",
"123", "RU-SP-B10-03 ??????????????????, ??????????????, ?????? ???????????????? ????????????. ?????? ???????????????? ??????????, ??????????.",
"190", "RU-SP-B10-03 ??????????????????, ??????????????, ?????? ???????????????? ????????????. ?????? ???????????????? ??????????, ??????????",
"205", "RU-SP-B19-01 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ?????????????? ??????????????, ???????????????????????????? ????????????",
"206", "RU-SP-B19-02 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ?????????????????????????? ?????????????? ????????????????????, ???????????????????? ????????????????.",
"209", "RU-SP-B19-04 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ???????????????? ??????????.",
"129", "RU-SP-B19-04 ??????. ???????????? ??????????, ???????????? ???????????????? ????????. ???????????????? ??????????.",
"306", "RU-SP-B21-03 ????????????. ?????????? ???????????????????????????? ??????????????. ?????????????????????????? ????????????",
"307", "RU-SP-B21-04 ????????????. ?????????? ???????????????????????????? ??????????????. ?????????????? ?? ?????????????? ????????????",
"310", "RU-SP-B22-01 ?????????????????? ????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????. ???????????????????????? ????????????????????????. ?????????????????? ????????????, ?????????????? ????????????",
"311", "RU-SP-B22-01 ?????????????????? ????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????. ???????????????????????? ????????????????????????. ?????????????? ????????????.",
"312", "RU-SP-B22-03 ?????????????????? ????????????. ?????????????????? ???????????????????? ?????????????? ?? ??????????????????. ???????????????????????? ????????????????????????. ?????????????? ????????????, ???????????????? ????????????",
"92", "RU-SP-B03-01 ???????????????????????? ??????. ???????????? ??????????????. ???????????? ?????????? ??????????????.",
"184", "RU-SP-B03-03 ???????????????????????? ??????. ???????????? ??????????????: ???????? ????????????.",
"84", "RU-SP-B03-04 ???????????????????????? ??????. ???????????? ??????????????: ???????? ??????????????.",
"302", "RU-SP-B08-01 ??????????????????????-???????????????????????? ????????????. ?????????????????? ??????.",
"308", "RU-SP-B08-04 ??????????????????????-???????????????????????? ????????????. ???????????????? ?????????????????? ??????",
"313", "RU-SP-B08-06 ??????????????????????-???????????????????????? ????????????. ???????????????? ????????????????????, ?????????????? ??????????",
"316", "RU-SP-B08-07 ??????????????????????-???????????????????????? ????????????. ?????????????????? ????????????????????-????????????????????.",
"318", "RU-SP-B08-08 ??????????????????????-???????????????????????? ????????????. ?????????????????? ??????????????????.",
"121", "RU-SP-B09-02 ??????????????????????, ???????????? ?? ??????????????????????????, ???????????? ?????????? ??????????????????????????. ???????????? ???????????? ??????????????????????????.",
"188", "RU-SP-B09-02 ??????????????????????, ???????????? ?? ??????????????????????????, ???????????? ?????????? ??????????????????????????. ???????????? ???????????? ??????????????????????????",
"303", "RU-SM-A01-01 ????????????????????. ???????????????? ?? ??????????????",
"304", "RU-SM-A02-01 ????????????????????. ???????????????? ???? ??????????????????",
"207", "RU-SM-A08-01 ???????????????????? ??????????????????. ????????????????????????????",
"239", "RU-SM-A15-01 ??????????????????????. ??????????????????????",
"240", "RU-SM-A15-02 ??????????????????????. ??????????????????????????: ?????????? ?? ??????????????",
"238", "RU-SM-A15-03 ??????????????????????. ????????????????????????????: ?????????? ?? ??????????????",
"243", "RU-SM-A15-04 ??????????????????????. ????????: ?????????? ?? ??????????????",
"246", "RU-SM-A15-05 ??????????????????????. ????????????????: ?????????? ?? ??????????????",
"275", "RU-SM-A19-01 ?????????? ?? ???? ????????????????. ???????????????? ???????????? ??????????.",
"210", "RU-SM-B03-00 ?????????? ???????????????????????? ????????????????",
"344", "RU-SM-B03-00 ?????????? ???????????????????????? ????????????????",
"213", "RU-SM-B01-01 ???????????????????? ?????????????????? ????????????. ???????????????????? ?? ??????????????????????",
"214", "RU-SM-B01-01 ???????????????????? ?????????????????? ????????????. ???????????????????? ?? ??????????????????????",
"222", "RU-SM-B01-01 ???????????????????? ?????????????????? ????????????. ???????????????????? ?? ??????????????????????.",
"110", "RU-SM-B01-02 ???????????????????? ?????????????????? ????????????. ???????????????????? ?? ????????????????.",
"321", "RU-SM-B01-03 ???????????????????? ?????????????????? ????????????.???????????? ????????????.",
"320", "RU-SM-B01-04 ???????????????????? ?????????????????? ????????????.????????????????, ????????????????????.",
"202", "RU-SM-B01-01 ???????????????????? ?????????????????? ????????????. ???????????????????? ?? ??????????????????????.",
"339", "RU-SM-A10-01 ???????????? ???????????? ????????????????????????. ???????????????????????? ?????????????????????? ??????????????????????",
"340", "RU-SM-A10-02 ???????????? ???????????? ????????????????????????. ?????????????? ?? ???????????????????????? ??????????????",
"249", "RU-SM-B10-07 ???????????????????? ?? ????????????????????????????: ???????????????????????????? ???????????????? ?????????????????????????????? ??????????????????",
"253", "RU-SM-B11-01 ???????????? ?? ???????????????????? ??????????????????????.???????????????? ?????????????????? ?? ??????????????????????",
"254", "RU-SM-B11-02 ???????????? ?? ???????????????????? ??????????????????????.???????????????????????? ?? ?????????????????? ?????????????????? ?? ??????????????????????",
"262", "RU-SM-B11-03 ???????????? ?? ???????????????????? ??????????????????????.???????????????????????? ?????????????????? ?? ??????????????????????",
"309", "RU-SM-B11-04 ???????????? ?? ???????????????????? ??????????????????????.???????????????????????????? ?????????????????? ?? ??????????????????????.",
"319", "RU-SM-B11-05 ???????????? ?? ???????????????????? ??????????????????????.?????????????????????????? ?????????????????? ?? ??????????????????????.",
"264", "RU-SM-B11-06 ???????????? ?? ???????????????????? ??????????????????????. ?????????????????????????????? ?????????????????? ?? ??????????????????????",
"300", "RU-SM-B11-07 ???????????? ?? ???????????????????? ??????????????????????. ???????????????????????????????????? ?????????????????? ?? ??????????????????????",
"279", "RU-SM-B12-06 ???????????? ???? ????????????????????????. ??????????.",
"278", "RU-SM-B12-04 ???????????? ???? ????????????????????????. ????????????????.",
"277", "RU-SM-B12-03 ???????????? ???? ????????????????????????. ????????????.",
"280", "RU-SM-B12-07 ???????????? ???? ????????????????????????. ??????.",
"194", "RU-SM-B13-01 ?????????????????? ????????????. ???????????? ???? ????????????????, ???????????? ?? ??????????.",
"195", "RU-SM-B13-02 ?????????????????? ????????????. ???????????? ???? ???????????????? ???? ????????????.",
"198", "RU-SM-B13-03 ?????????????????? ????????????. ???????????? ???? ???????????????? ???? ????????????????????.",
"132", "RU-SM-B13-03 ?????????????????? ????????????. ???????????? ???? ???????????????? ???? ????????????????????.",
"196", "RU-SM-B13-04 ?????????????????? ????????????. ???????????? ???? ???????????????? ???? ????????.",
"95",  "RU-SM-B13-05 ?????????????????? ????????????.???????????? ???? ???????????????????? ????????????.",
"200", "RU-SM-B13-05 ?????????????????? ????????????. ???????????? ???? ???????????????????? ????????????.",
"185", "RU-SM-B13-05 ?????????????????? ????????????. ???????????? ???? ???????????????????? ????????????.",
"104", "RU-SM-B13-06 ?????????????????? ????????????.???????????? ???? ????????????????????.",
"197", "RU-SM-B13-06 ?????????????????? ????????????. ???????????? ???? ????????????????????",
"96",  "RU-SM-B13-01 ?????????????????? ????????????. ???????????? ???? ????????????????, ???????????? ?? ??????????.",
"106", "RU-SM-B13-04 ?????????????????? ????????????.???????????? ???? ???????????????? ???? ????????",
"90",  "RU-SM-B13-04 ?????????????????? ????????????.???????????? ???? ???????????????? ???? ????????.",
"91",  "RU-SM-B13-06 ?????????????????? ????????????. ???????????? ???? ????????????????????.",
"204", "RU-SM-B14-01 ???????????????????? ?? ???????????????????? ???????????????? ??????????????. ???????????????? ?? ?????????????? ??????????????. ???????????????????? ?? ???????????????????? ????????????????",
"219", "RU-SM-A03-01 ???????????????????? ?????????????????? ????????????. ????????????????, ???????????????????? 1 ???? 3??.",
"220", "RU-SM-A03-02 ???????????????????? ?????????????????? ????????????. ????????????????, ???????????????????? 2 ???? 3??.",
"223", "RU-SM-A03-03 ???????????????????? ?????????????????? ????????????. ????????????????, ???????????????????? 3 ???? 3??.",
"226", "RU-SM-B01-031 ???????????????????? ?????????????????? ????????????. ???????????? ???????????? 1 ???? 2??.",
"227", "RU-SM-B01-032 ???????????????????? ?????????????????? ????????????. ???????????? ???????????? 2 ???? 2??.",
"338", "RU-SM-A20-01  ???????????? ???? ????????????????.",
"251", "RU-SP-B25-01 ????????????????. ???????????????????????? ???????????? (?????????????????? ????????????). ????????????????",
"252", "RU-SP-B25-02 ????????????????. ???????????????????????? ???????????? (?????????????????? ????????????). ????????????????????",
"276", "RU-SM-B03-02 ?????????? ???????????????????????? ??????????????. ?????????? ???????????????? ???? ???????? ??????????????????.",
"69", "RU-SM-B04-01 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ??????????????????????",
"70", "RU-SM-B04-02 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ??????????????????????????",
"71", "RU-SM-B04-03 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ????????????????????????????",
"113", "RU-SP-B04-04 ??????????????, ????????????, ????????. ???????????? ????????????????????. ???????????????????????? ????????????",
"73", "RU-SM-B04-04 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ????????",
"65", "RU-SM-B04-05 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ????????????????",
"74", "RU-SM-B04-06 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ???????????????????????? ??????????????????????????????",
"75", "RU-SM-B04-07 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ??????????????????????????",
"109", "RU-SM-B04-08 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ???????????? ???? ???????????????????? ??????????????",
"134", "RU-SM-B04-09 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ???????? ?? ?????? ????????????????",
"102", "RU-SM-B04-10 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ?????????????????? ?? ?????????????????? ????????????????????",
"128", "RU-SM-B04-11 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ??????????????",
"127", "RU-SM-B04-12 ??????????????????????: ???????????????????? ???????? ?? ????????????????. ???????????????????????? ??????????????????",
"224", "RU-SM-B05-01 ???????????? ???????????? ????????????????????????. ???????????????????????? ?????????????????????? ??????????????????????.",
"116", "RU-SM-B05-01 ???????????? ???????????? ????????????????????????. ???????????????????????? ?????????????????????? ??????????????????????.",
"229", "RU-SM-B05-02 ???????????? ???????????? ????????????????????????. ?????????????? ?? ???????????????????????? ??????????????.",
"117", "RU-SM-B05-02 ???????????? ???????????? ????????????????????????. ?????????????? ?? ???????????????????????? ??????????????.",
"72", "RU-SM-B06-01",
"186", "RU-SM-B06-01 ???????????????????? ??????????????????. ?????????????? ???????????????? ??????????????????",
"187", "RU-SM-B06-01 ???????????????????? ??????????????????. ?????????????? ???????????????? ??????????????????",
"161", "RU-SM-B06-02 ???????????????????? ??????????????????. ???????????????????????? ??????????????????",
"166", "RU-SM-B06-03 ???????????????????? ??????????????????. ???????????????????????????? ??????????????????",
"170", "RU-SM-B06-04 ???????????????????? ??????????????????. ?????????????????????????? ??????????????????",
"221", "RU-SM-B06-05 ???????????????????? ??????????????????. ?????????????????????????????? ??????????????????",
"268", "RU-SM-B07-01 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ?????????????????????????? ??????????????????????: ???????????????????? ??????????.",
"269", "RU-SM-B07-02 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ?????????????????????????? ??????????????????????: ???????????????????? ?????????????? ??????????.",
"271", "RU-SM-B07-03 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ?????????????????????????? ??????????????????????: ???????????????????? ??????????????????.",
"273", "RU-SM-B07-04 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ???????????????????????????? ??????????????????????: ???????????????????? ??????????????????.",
"274", "RU-SM-B07-05 ??????????????????????: ????????????, ?????????????????? ?? ????????????. ???????????????????????????? ??????????????????????: ???????????????????? ??????????.",
"215", "RU-SM-A08-01 ???????????????????? ??????????????????. ????????????????????????????",
"250", "RU-SP-B08-03 ??????????????????????-???????????????????????? ????????????. ???????????????????????? ???????????? ?? ???????????????? ?????????????????? ??????",
"255", "RU-SP-B08-05 ??????????????????????-???????????????????????? ????????????. ?????????????? ???????????????????????? ?????????????? ?????????????? ???????????????????? ????????",
"78",  "RU-SM-B09-01 ????????????????????????. ??????",
"100", "RU-SM-B09-02 ????????????????????????. ?????????????????????????? ????????????????????????????",
"101", "RU-SM-B09-03 ????????????????????????. ????????????",
"126", "RU-SM-B09-07 ????????????????????????. ?????????? ???????????????????? ??????????????????????????",
"124", "RU-SM-B09-05 ????????????????????????. ???????????????? ?????????????????? ????????????????????????????",
"125", "RU-SM-B09-06 ????????????????????????. ?????????????? ?????????????????????? ???????????????????? ??????????????????????????",
"80",  "RU-SM-B09-08 ????????????????????????. ???????????????????? ??????",
"81",  "RU-SM-B09-09 ????????????????????????. ??????????????",
"82",  "RU-SM-B09-10 ????????????????????????. ??????????",
"79",  "RU-SM-B09-04 ????????????????????????. ????????????????",
"330", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ????????????.",
"331", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ????????????????????.",
"332", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ??????????????????????????????.",
"333", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ???????????????????????????? ???????????????? ?? ??????????.",
"334", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ?????????????? ????????????",
"335", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ?????????????? ????????????",
"336", "RU-SP-C50-01 ??????????????????????????????. ?????????????????? ???????????? (?????????????????? ????????????). ?????????????? ????????????",
"361", "RU-SM-GB020-01 ??????????????????, ?????????????????????? ?? ???? ??????????????. ???????????????????? ??????????????????",
"358", "RU-SM-G??020-02 ??????????????????, ?????????????????????? ?? ???? ??????????????. ???????????????????????? ??????????????????",
"379", "RU-SM-G??020-04 ??????????????????, ?????????????????????? ?? ???? ??????????????. ?????????????? ????????????????????",
"378", "RU-SM-G??020-05 ??????????????????, ?????????????????????? ?? ???? ??????????????. ???????????????? ??????????????????,",
"140", "RU-SP-B04-01 ??????????????, ????????????, ????????. ???????????? ????????????????????. ?????????? ???????????????????? ??????????????????",
"68",  "RU-SP-B04-02 ??????????????, ????????????, ????????. ???????????? ????????????????????. ?????????? ???????????????????? ????????????????",
"143", "RU-SP-B04-03 ??????????????, ????????????, ????????. ???????????? ????????????????????. ?????????? ???????????????????? ??????????????",
"139", "RU-SP-B04-05 ??????????????, ????????????, ????????. ???????????? ????????????????????. ?????????????????????????? ??????????????",
"248", "RU-SP-B08-01 ??????????????????????-???????????????????????? ????????????. ???????????????????? ?????????????????? ????????????????",
"305", "RU-SM-??06-03 ???????????????????? ?????????????????? ????????????. ???????????? ????????????",
"341", "RU-SP-B02-01 ????????????????????, ???????????? ??????????????. ???????????? ?????????? ??????????????",
"343", "RU-SP-B02-03 ????????????????????, ???????????? ??????????????. ???????????? ?????????? ??????????????",
"342", "RU-SP-B02-05 ????????????????????, ???????????? ??????????????. ???????????? ?????????? ??????????????",
"244", "RU-SP-B03-05 ???????????????????????? ??????. ???????????? ??????????????: ???????? ??????????????????",
"337", "RU-SP-B05-02 ??????????????, ???????????????????????? ?????????????????? ?? ??????????. ?????????????????? ?? ???????????????????????????? ????????????????, ??????????????????.",
"281", "RU-SM-A05-01 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????????????????? ?????????????????? ?? ????????????",
"282", "RU-SM-A05-02 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ???????????????? ???????????????????????????? ??????????????????",
"295", "RU-SM-A05-04 ???????????????????? ?? ????????????????????????????. ???????????????????????????? ?????????????????? ?????????????????????????? ??????????????????",
"315", "RU-SM-A05-05 ???????????????????? ?? ???????????????????????????? ???????????????????????????? ???????????????? ?????????????????????????????? ??????????????????",
"189", "RU-SP-B09-01 ??????????????????????, ???????????? ?? ??????????????????????????, ???????????? ?????????? ??????????????????????????. ????????????????????????????, ??????????????????????????, ??????????????????????????, ???????????????????????????? ????????????????, ??????????",
"2", "RU-SM-B14-01-01 ???????????????????? ?? ???????????????????? ???????????????? ??????????????. ???????????????????????? ?????????????????? ?? ???????????????????????????? ??????????????",
"83", "RU-SM-B09-11 ????????????????????????. ??????",
"164", "RU-SP-B26-02 ???????????????????????? ????????????, ??????????????????????????. ??????????????????????????????. ??????????????????????????????.",
"6", "EN-MS Geometry. Trapezoid",
"99", "RU-SP-B10-01 ??????????????????, ??????????????, ?????? ???????????????? ????????????"
};

    @Override
    public void run(String... strings) throws Exception {
        final String h2Version = jdbcTpl.queryForObject(
                "select value from information_schema.settings where name = :version;",
                new MapSqlParameterSource().addValue("version", "info.VERSION"),
                String.class);
        LOG.info("==============================================");
        LOG.info("H2 version: {}", h2Version);
        LOG.info("==============================================");
        for (int i = 0; i < newTit.length; i+=2) {
            long id= new Long((String)newTit[i]);  
            String Tit= (String)newTit[i+1];
            if((id<0)||(id>1000)){LOG.info(" in:renameTask ERORR_IN_id_number!!!==================== id="+id+"???");
               LOG.info(" in: renameTask Title=="+Tit);              break;}
            renameTask(id, Tit);
        }
        //renameTask(1, "Hello, world");
    }

    public void renameTask(final long taskId, final String newTaskTitle) {
        jdbcTpl.update("UPDATE tasks SET task_title = :newTaskTitle WHERE id = :taskId",
                new MapSqlParameterSource().addValue("newTaskTitle", newTaskTitle)
                .addValue("taskId", taskId));
    }
}
//
// new Long(140),"RU-SP-B4-01 ?????????? ???????????????????? ??????????????????"
//new Long(68),"RU-SM-pB4-02 ?????????? ?????????????????????????? ????????????????	??
//143	RU-SM-pB4-03 ?????????? ?????????????????????????? ??????????????	??
//139	RU-SM-pB4-05 ?????????????????????????? ??????????????	??
//248	RU-SM-pB8-01 ???????????????????? ???????????????????????? ????????????????	??