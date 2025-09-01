function projectXYZ = projectOnLine(pt, wc1, wc2)
    if wc2(1)==wc1(1)
        x = wc1(1);
        y = pt(2);
    elseif wc1(2)==wc2(2)
        x = pt(1);
        y = wc1(2);
    else
       %y=ax+b
       a = (wc2(2)-wc1(2))/(wc2(1)-wc1(1));
       b = wc1(2)-a*wc1(1);
       
       a2 = -1/a;
       b2 = pt(2)-a2*pt(1);
       
       x = (b2-b)/(a-a2);
       y = a*x+b;
    end

projectXYZ = [x y pt(3)];
end