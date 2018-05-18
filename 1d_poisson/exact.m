function y = exact(x,c)

    if(x<c && x>0)
        y = - x*x/2+(c-c*c/2)*x;
    elseif(x>c && x<1)
        y = c*c*(1-x)/2;
    else
        y = 0;
    end

end