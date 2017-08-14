function result=product_q(q1, q2)
     if (size(q1,1)==1 && size(q1,2)==4)||(size(q1,1)==4 && size(q1,2)==1)
          if (size(q2,1)==1 && size(q2,2)==4)||(size(q2,1)==4 && size(q2,2)==1)
             result=[q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
                     q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
                     q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
                     q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)];
          else
             error('put quaternion in q2'); 
          end
    else
          error('put quaternion in q1'); 
    end
end