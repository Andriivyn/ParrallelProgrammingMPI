package com.mathpar.web;

//import org.springframework.boot.context.embedded.MultipartConfigFactory;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.web.method.support.HandlerMethodArgumentResolver;
//import org.springframework.web.servlet.config.annotation.WebMvcConfigurerAdapter;
import com.mathpar.web.controller.PageArgumentResolver;

import javax.servlet.MultipartConfigElement;
//import javax.servlet.annotation.MultipartConfig;
import java.util.List;
import org.springframework.boot.web.servlet.MultipartConfigFactory;
import org.springframework.web.servlet.config.annotation.WebMvcConfigurerAdapter;
//import org.springframework.web.servlet.config.annotation.WebMvcConfigurerAdapter;

@Configuration
public class MvcConfig extends WebMvcConfigurerAdapter {
    @Override
    public void addArgumentResolvers(List<HandlerMethodArgumentResolver> argumentResolvers) {
        argumentResolvers.add(new PageArgumentResolver());
    }

    @Bean
    MultipartConfigElement multipartConfigElement() {
        MultipartConfigFactory factory = new MultipartConfigFactory();
        //factory.setMaxFileSize(DataSize 1024);
        //factory.setMaxRequestSize("1024KB");
        return factory.createMultipartConfig();
    }
}
