function result = DCM_ortho_normal_compensation(DCM)
    if size(DCM)==[3,3]
        count = 0;
        sum=0;
        while 1
            
            count = count + 1;
            for a=1:3
                for b=1:3
                    if a~=b
                        chk_norm = DCM(a,:)*DCM(b,:)';
                        DCM(a,:) = DCM(a,:) - chk_norm/2*DCM(b,:);
                        DCM(b,:) = DCM(b,:) - chk_norm/2*DCM(a,:);
                        %Gram-Schmidt Orthogonalization, ±×¶÷-½´¹ÌÃ÷ Á÷±³È­ÀÇ ÀÀ¿ë, Àý¹Ý¾¿ »­
                    end
                end
            end 
            for a=1:3
                DCM(a,:)=DCM(a,:)/norm(DCM(a,:)); 
            end
            for a=1:3
                for b=1:3
                    if a==b
                        if norm(DCM(a,:)) - 1 < 0.0000001
                            sum = sum+1;
                        end
                    else
                        if DCM(a,:)*DCM(b,:)' < 0.0000001
                            sum = sum+1;
                        end
                    end
                end
            end 
            if count==10
                result = DCM;
                break;
            else
                sum=0;
            end
        end
    else
       error('put in 3x3 DCM'); 
    end 
end