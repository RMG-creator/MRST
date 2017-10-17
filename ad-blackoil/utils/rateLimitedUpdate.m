function rsMax = rateLimitedUpdate(rs,so,rs0,so0,rsSat,fluid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% rsMax = min(rsSat,rs0);
 %rsMax=min(rsSat,rs0);
 %rsMax= 0.9*rsSat;
 %rsMax=rsSat;
 drsdt = 0;
 method = 'orig';
 switch method
     case 'rate_limmit'
        rsso0=rs0.*so0;
        rssomax=rsSat.*so;
        rssomax=min(rssomax,rsso0+drsdt);
        rsMax=rssomax./so;
     case 'E100'
        rsMax = min(rsSat,rs0+drsdt); 
     case 'orig'    
        rsMax=rsSat;
     otherwise
         error('No method')
         rsMax = 0;
 end
end


