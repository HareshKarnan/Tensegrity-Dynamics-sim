function [N,Cb,Cs] = tenseg_prism_sidestring(ht,edge_l)

thx = 30;

N = getpoints(thx,ht,edge_l);
Cb = [-eye(3) eye(3)];

Cs = [0 0 0 0 -1 1; % 5-6
      0 0 0 1 0 -1; % 6-4
      0 0 0 -1 1 0; % 4-5
      0 -1 0 0 0 1; % 2-6
      0 0 -1 1 0 0; % 3-4
      -1 0 0 0 1 0]; % 1-5
%       extra strings
%       -1 0 0 0 0 1; %1-6
%       0 -1 0 1 0 0; %2-4
%       0 0 -1 0 1 0]; %3-5
  
end

function N = getpoints(th,h,s)
x1 = [-s/2;-(s/2)*tand(30);0];
x2 = [0;s*tand(30);0];
x3 = [s/2;-(s/2)*tand(30);0];
R = [cosd(th) -sind(th) 0;sind(th) cosd(th) 0;0 0 1];
x4 = R*x3;
x5 = R*x1;
x6 = R*x2;
x4(3)=h;x5(3)=h;x6(3)=h;
N = [x1 x2 x3 x4 x5 x6];
end