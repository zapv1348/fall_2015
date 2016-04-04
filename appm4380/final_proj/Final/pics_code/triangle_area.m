function [area]=triangle_area(A,B,C)
        area=abs((A(1)*(B(2)-C(2))+B(1)*(A(2)-C(2))+C(1)*(A(2)-B(2)))/2);
end