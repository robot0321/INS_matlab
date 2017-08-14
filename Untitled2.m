
clc;
close all;
clear;

s = serial('COM3','BaudRate', 115200);
fopen(s);
pause(0.1);

while(1)
    a = fscanf(s);
    fprintf(a);
end

fclose(s);
delete(s);