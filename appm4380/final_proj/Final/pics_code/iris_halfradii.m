function [theta1,theta2]=iris_halfradii(p1,p2,p3,p4)
            %p1 and p2 define the radius of the iris for one image
            %p3 and p4 for the other image
            rx1=abs(p1(1)-p2(1));
            rx2=abs(p3(1)-p4(1));
            ry1=abs(p1(2)-p2(2));
            ry2=abs(p3(2)-p4(2));
            if rx1>rx2
                theta1=acos(rx2/rx1);
            else
                theta1=acos(rx1/rx2);
            end
            if  ry1>ry2
                theta2=acos(ry2/ry1);
            else
                theta2=acos(ry1/ry2);
            end
end