%syms phi psi theta phidot psidot thetadot i1 i3 mgl;
%{
L = i1/2 * (thetadot^2 + phi^2*sin(theta)^2) + i3/2 * (psidot+phidot*cos(theta))^2 - mgl*cos(theta);

disp(sprintf('Ppsi = %s',char(diff(L,psidot))))
disp(sprintf('Pphi = %s',char(diff(L,phidot))))
disp(sprintf('Ptheta = %s',char(diff(L,thetadot))))
%}


clear all; close all;
syms 'phi' 'psi' 'theta' 'Pphi' 'Ppsi' 'Ptheta' 'mgl' 'i1' 'i3'

H = Ptheta^2/(2*i1) + (Pphi-Ppsi*cos(theta))^2/(2*i1*sin(theta)^2) + Ppsi^2/(2*i3) + mgl*cos(theta);

phidot  = diff(H,phi);
psidot  = diff(H,psi);
thetadot  = diff(H,theta);

Pphidot  = diff(H,Pphi);
Ppsidot  = diff(H,Ppsi);
Pthetadot  = diff(H,Ptheta);

hamjac = jacobian([phidot psidot thetadot Pphidot Ppsidot Pthetadot], [phi psi theta Pphi Ppsi Ptheta]);

for i=1:size(hamjac, 1)
	for j=1:size(hamjac,2)
		if(hamjac(i,j)~=0)
			disp(sprintf('J_{%d%d} = %s', i, j, char(hamjac(i,j))))
		end
	end
end
