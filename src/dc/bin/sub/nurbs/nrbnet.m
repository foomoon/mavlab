function [x,y,z] = nrbnet(srf)

x=squeeze(srf.coefs(1,:,:));
y=squeeze(srf.coefs(2,:,:));
z=squeeze(srf.coefs(3,:,:));