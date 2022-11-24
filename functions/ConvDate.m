function y = ConvDate(date)
    formatOut = 'yyyymm';
    t = datestr(date,formatOut);
    t = strcat(t(1:4), 'm', t(5:end));

    y = regexprep(t,'m0','m');
end

