function intersect = intersect2lineGiven4Points2D(p1,p2,p3,p4)
    p1x = p1(1);
    p1y = p1(2);
    
    p2x = p2(1);
    p2y = p2(2);
    
    p3x = p3(1);
    p3y = p3(2);
    
    p4x = p4(1);
    p4y = p4(2);
    
    m2 = (p4y-p3y)/(p4x-p3x);
    m1 = (p2y-p1y)/(p2x-p1x);
    
    b2 = p4y-m2*p4x;
    b1 = p2y-m1*p2x;
    
    x = (b2-b1)/(m1-m2);
    y = m1*x+b1;
    
%     x = ((p2x*p1y-p1x*p2y)*(p4x-p3x)-(p4x*p3y-p3x*p4y)*(p2x-p1x))...
%         /((p2x-p1x)*(p4y-p3y)-(p4x-p3x)*(p2y-p1y));
%     
%     y = ((p2x*p1y-p1x*p2y)*(p4y-p3y)+(p4x*p3y+p3x*p4y)*(p2y-p1y))...
%         /((p2x-p1x)*(p4y-p3y)-(p4x-p3x)*(p2y-p1y));
    
    if p1y==p2y
       y = p1y;
       m2 = (p4y-p3y)/(p4x-p3x);
       b = p4y-m2*p4x;
       x = (y-b)/m2;
    end
    
    if p3y==p4y
       y = p3y;
       m2 = (p2y-p1y)/(p2x-p1x);
       b = p2y-m2*p2x;
       x = (y-b)/m2;
    end
    
    if p1x==p2x
       x =  p1x;
       m2 = (p4y-p3y)/(p4x-p3x);
       b = p3y-m2*p3x;
       y = m2*x+b;
    end
    
    if p3x==p4x
       x =  p3x;
       m2 = (p2y-p1y)/(p2x-p1x);
       b = p1y-m2*p1x;
       y = m2*x+b;
    end
    
    
    
    intersect = [x y];
end