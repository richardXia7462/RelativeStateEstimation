function [nFace, nAncPerFac, markF, FNormal] = GetMarks(tagsize)

nFace = 6;
nAncPerFac = 4;
% Target (counterclockwise from left bottom)
L0 = [-tagsize/2 -tagsize/2 0; tagsize/2 -tagsize/2 0; tagsize/2 tagsize/2 0; -tagsize/2 tagsize/2 0]';

TFaces = zeros(4,4,nFace);
% face(up front back left right) frame to follower frame
pUpCen = [0 0 tagsize/2]; 
Tup = trvec2tform(pUpCen)*eul2tform([-pi/2 0 0]);
TFaces(:,:,1) = Tup;
% x axis point to front, y axis point to left
pFrontCen = [tagsize/2 0 0];
Tfront = trvec2tform(pFrontCen)*eul2tform([pi/2 0 pi/2]);
TFaces(:,:,2) = Tfront;
pBackCen = [-(tagsize/2) 0 0];
Tback = trvec2tform(pBackCen)*eul2tform([-pi/2 0 pi/2]);
TFaces(:,:,3) = Tback;
pLeftCen = [0 tagsize/2 0];
Tleft = trvec2tform(pLeftCen)*eul2tform([pi 0 pi/2]);
TFaces(:,:,4) = Tleft;
pRightCen = [0 -tagsize/2 0];
Tright = trvec2tform(pRightCen)*eul2tform([0 0 pi/2]);
TFaces(:,:,5) = Tright;
pBottomCen = [0 0 -tagsize/2];
Tbottom = trvec2tform(pBottomCen)*eul2tform([-pi/2 0 pi]);
TFaces(:,:,6) = Tbottom;

% target in follower frame
markF = zeros(3,20);
markF(:,1:4) = Tup(1:3,:)*[L0; ones(1,nAncPerFac)];
markF(:,5:8) = Tfront(1:3,:)*[L0; ones(1,nAncPerFac)];
markF(:,9:12) = Tback(1:3,:)*[L0; ones(1,nAncPerFac)];
markF(:,13:16) = Tleft(1:3,:)*[L0; ones(1,nAncPerFac)];
markF(:,17:20) = Tright(1:3,:)*[L0; ones(1,nAncPerFac)];
markF(:,21:24) = Tbottom(1:3,:)*[L0; ones(1,nAncPerFac)];

FNormal = zeros(3, nFace);
for i=1:nFace
    FNormal(:,i) = TFaces(1:3,3, i);
end