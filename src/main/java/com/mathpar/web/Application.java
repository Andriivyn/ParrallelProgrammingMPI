package com.mathpar.web;

import static org.apache.logging.log4j.LogManager.getLogger;

import org.apache.logging.log4j.Logger;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.EnableAutoConfiguration;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.boot.web.servlet.support.SpringBootServletInitializer;
//import org.springframework.boot.context.web.SpringBootServletInitializer;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.FilterType;
import org.springframework.jdbc.core.JdbcTemplate;

@Configuration
@EnableAutoConfiguration
@ComponentScan(
        basePackages = {"com.mathpar.web"},
        excludeFilters = {
                @ComponentScan.Filter(
                        value = {ApplicationPlayground.class},
                        type = FilterType.ASSIGNABLE_TYPE)
        })
public class Application extends SpringBootServletInitializer {
    private static final Logger LOG = getLogger(Application.class);

    public static void main(String[] args) throws Exception {
        ConfigurableApplicationContext ctx = SpringApplication.run(Application.class, args);

        final JdbcTemplate jdbcTpl = ctx.getBean(JdbcTemplate.class);
        final String h2Version = jdbcTpl.queryForObject(
                "select value from information_schema.settings where name = ?;",
                String.class, "info.VERSION");
        LOG.info("==============================================");
        LOG.info("H2 version: {}", h2Version);
        LOG.info("==============================================");
    }

    // This is to run spring boot inside Tomcat container (not using embedded one)
    @Override
    protected SpringApplicationBuilder configure(SpringApplicationBuilder application) {
        return application.sources(Application.class);
    }
}
