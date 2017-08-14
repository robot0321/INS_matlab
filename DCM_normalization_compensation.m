function result=DCM_normalization_compensation(DCM)
    if size(DCM)==[3,3]
        for a=1:3
%             chk_norm = 1 - DCM(a,:)*DCM(a,:)';
%             DCM(a,:) = DCM(a,:) - chk_norm/2*DCM(a,:);
            DCM(a,:)=DCM(a,:)/norm(DCM(a,:));
        end
        result = DCM;
    else
        error('put in 3x3 DCM'); 
    end
end