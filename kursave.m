function y=kursave(x)

    f1=sum(((-10)*exp(((-0.2)*sqrt((x(1:2).^2+x(2:3).^2))))));
    
    f2=sum((abs(x(1:2)).^0.8)+(5*sin(pi*(x(1:2).^3))));
    
    y=[f1
       f2];

end