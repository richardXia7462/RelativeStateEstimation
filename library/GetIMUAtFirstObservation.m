function [imuF, iteIMUF, imuL, iteIMUL] = GetIMUAtFirstObservation(tFirst, IMUF, IMUL)
%UNTITLED6 Summary of this function goes here
%   imuF: imu of the follower at tFirst
%   iteIMUF: number of first imu of the follower after tFirst
%   imuL: imu of the leader at tFirst
%   iteIMUL: number of first imu of the leader after tFirst

iteIMUF = 1;
while IMUF(iteIMUF,1)<=tFirst
    iteIMUF=iteIMUF+1;
end
iteIMUL=1;
while IMUL(iteIMUL,1)<=tFirst
    iteIMUL=iteIMUL+1;
end

if iteIMUF>1 
    k = iteIMUF-1; % <tFirst
else
    k = 1; % >tFirst
end
tab = IMUF(k+1,1) - IMUF(k,1);
tini = tFirst - IMUF(k,1);
angVel = IMUF(k,2:4) + (IMUF(k+1,2:4) - IMUF(k,2:4)) * (tini / tab);
acc = IMUF(k,5:7) + (IMUF(k+1,5:7) - IMUF(k,5:7)) * (tini / tab);
imuF = [tFirst, angVel, acc];

if iteIMUL>1 
    k = iteIMUL-1; % <tFirst
else
    k = 1; % >tFirst
end
tab = IMUL(k+1,1) - IMUL(k,1);
tini = tFirst - IMUL(k,1);
angVel = IMUL(k,2:4) + (IMUL(k+1,2:4) - IMUL(k,2:4)) * (tini / tab);
acc = IMUL(k,5:7) + (IMUL(k+1,5:7) - IMUL(k,5:7)) * (tini / tab);
imuL = [tFirst, angVel, acc];

end

