function attitude_visualize(DCM1, DCM2)
    M1 = DCM1;
    M2 = DCM2;
    [yaw_c, pitch_c, roll_c] = dcm2angle(M1(:,:,1));
    [yaw_d, pitch_d, roll_d] = dcm2angle(M2(:,:,1));
    for frame=1:1:max(size(M1,3),size(M2,3))
       
        subplot(1,2,1);
        [yaw, pitch, roll] = dcm2angle(M1(:,:,frame));

        [x,y,z] = sphere;
        h = surf(x,y,z);axis('square'); 
        title(frame)
        xlabel('x'); ylabel('y'); zlabel('z');
        %Rotate Object
        rotate(h,[1,0,0],(roll-roll_c)*180/pi)
        rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
        rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
        
        subplot(1,2,2);
        [yaw, pitch, roll] = dcm2angle(M2(:,:,frame));

        [x,y,z] = sphere;
        h = surf(x,y,z);axis('square'); 
        title(frame)
        xlabel('x'); ylabel('y'); zlabel('z');
        %Rotate Object
        rotate(h,[1,0,0],(roll-roll_d)*180/pi)
        rotate(h,[0,1,0],(pitch-pitch_d)*180/pi)
        rotate(h,[0,0,1],(yaw-yaw_d)*180/pi+90)
        
        
        
        drawnow
    end
end