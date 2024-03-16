function [fia,fib]=AR_filter(sig,can,p)

    s=sig(:,can);
    s1=size(s,1);

    usc = zeros(s1,1);
    fib = zeros(s1,1);

    usc(1)= s(1);
    for i = 2 : s1
       usc(i)= p*usc(i-1)+(1-p)*s(i);
    end 

    fib(s1)=usc(s1);
    for i = (s1-1):-1:1
         fib(i)= p*fib(i+1)+(1-p)*usc(i);
    end

    fia = s-fib+mean(s);
end
