function [reg_mat , reg_frac] = handleRegionsDP(deck, G, varargin)


% if DP reg has _mat and _frac fields

if ~isempty(G)
    an = G.cells.indexMap;    
elseif isfield(deck.GRID, 'ACTNUM')
    an = find(deck.GRID.ACTNUM);
else
    an = ':';
end
  
%MATRIX FIRST ------------------------------------------------------  


% also need actnum/G.cells.indexmap
% PVT-REGIONS
ntpvt = 1;
if isfield(deck.REGIONS, 'PVTNUM')
    reg_mat.PVTNUM = deck.REGIONS.PVTNUM(an,1);

    ntpvt = max(reg_mat.PVTNUM);
end
if ntpvt == 1
    reg_mat.PVTNU = [];
    reg_mat.PVTINX = ':';
else
    reg_mat.PVTINX = cellfun(@(x)find(x==reg_mat.PVTNUM), num2cell(1:ntpvt), 'UniformOutput', false);
 
end

% SAT-REGIONS AND POSSIBLY SURF-REGIONS
one_region = true;
if isfield(deck.REGIONS, 'SATNUM')
    reg_mat.SATNUM = deck.REGIONS.SATNUM(an,1);
 
    
    if isfield(deck.REGIONS, 'SURFNUM')
       reg_mat.SURFNUM = deck.REGIONS.SURFNUM(an);
       ntsatsurfact = max(max(reg_mat.SATNUM), max(reg_mat.SURFNUM));
                       
                        
       reg_mat.SURFINX = cellfun(@(x)find(x==reg_mat.SURFNUM), num2cell(1:ntsatsurfact), 'UniformOutput', ...
                             false);
                        
                         
       one_region = false;
    else
       ntsat = max(reg_mat.SATNUM);
       if ntsat > 1
          reg_mat.SATINX = cellfun(@(x)find(x==reg_mat.SATNUM), num2cell(1:ntsat), 'UniformOutput', ...
                               false);
          one_region = false;
       end
    end
elseif isfield(deck.REGIONS, 'SURFNUM')
      error('SATNUM keyword required when surfactant is used.');
end

if one_region
      reg_mat.SATNUM = [];
      reg_mat.SATINX = ':';
end

% IMB-REGIONS
ntsat = 1;
if isfield(deck.REGIONS, 'IMBNUM')
    reg_mat.IMBNUM = deck.REGIONS.IMBNUM(an,1);
 
    ntsat = max(reg_mat.IMBNUM)/2;
end
if ntsat == 1
    reg_mat.IMBNUM = [];
    reg_mat.IMBINX = ':';
else
    reg_mat.IMBINX = cellfun(@(x)find(x==reg_mat.IMBNUM), num2cell(1:ntsat), 'UniformOutput', false);
    
end

% ROCK-REGIONS
ntrocc = 1;
if isfield(deck.REGIONS, 'ROCKNUM')
    reg_mat.ROCKNUM = deck.REGIONS.ROCKNUM(an,1);
    
    ntrocc = max(reg_mat.ROCKNUM);
end
if ntrocc == 1
    reg_mat.ROCKNUM_ = [];
    reg_mat.ROCKINX = ':';
else
    reg_mat.ROCKINX = cellfun(@(x)find(x==reg_mat.ROCKNUM), num2cell(1:ntrocc), 'UniformOutput', false);
   
end



%Fractured Region ------------------------------------------------------  


% also need actnum/G.cells.indexmap
% PVT-REGIONS
ntpvt = 1;
if isfield(deck.REGIONS, 'PVTNUM')
       reg_frac.PVTNUM = deck.REGIONS.PVTNUM(an,2);
    ntpvt = max(reg_frac.PVTNUM);
end
if ntpvt == 1
    reg_frac.PVTNUM = [];
    reg_frac.PVTINX = ':';
else
    reg_frac.PVTINX = cellfun(@(x)find(x==reg_frac.PVTNUM), num2cell(1:ntpvt), 'UniformOutput', false);
end

% SAT-REGIONS AND POSSIBLY SURF-REGIONS
one_region = true;
if isfield(deck.REGIONS, 'SATNUM')
     reg_frac.SATNUM = deck.REGIONS.SATNUM(an,2);
    
    if isfield(deck.REGIONS, 'SURFNUM')
       reg_frac.SURFNUM = deck.REGIONS.SURFNUM(an,2);
       ntsatsurfact = max(max(reg_frac.SATNUM), max(reg_frac.SURFNUM));
       
       reg_frac.SATINX = cellfun(@(x)find(x==reg_frac.SATNUM), num2cell(1:ntsatsurfact), 'UniformOutput', ...
                            false);
       reg_frac.SURFINX = cellfun(@(x)find(x==reg_frac.SURFNUM), num2cell(1:ntsatsurfact), 'UniformOutput', ...
                             false);
                         
                         
       one_region = false;
    else
       ntsat = max(reg_frac.SATNUM);
       if ntsat > 1
          reg_frac.SATINX = cellfun(@(x)find(x==reg_frac.SATNUM), num2cell(1:ntsat), 'UniformOutput', ...
                               false);
          one_region = false;
       end
    end
elseif isfield(deck.REGIONS, 'SURFNUM')
      error('SATNUM keyword required when surfactant is used.');
end

if one_region
      reg_frac.SATNUM = [];
      reg_frac.SATINX = ':';
end

% IMB-REGIONS
ntsat = 1;
if isfield(deck.REGIONS, 'IMBNUM')
  
    reg_frac.IMBNUM = deck.REGIONS.IMBNUM(an,2);
    ntsat = max(reg_frac.IMBNUM)/2;
end
if ntsat == 1
    reg_frac.IMBNUM = [];
    reg_frac.IMBINX = ':';
else
    reg_frac.IMBINX = cellfun(@(x)find(x==reg_frac.IMBNUM), num2cell(1:ntsat), 'UniformOutput', false);
end

% ROCK-REGIONS
ntrocc = 1;
if isfield(deck.REGIONS, 'ROCKNUM')

    reg_frac.ROCKNUM = deck.REGIONS.ROCKNUM(an,2);
    ntrocc = max(reg_frac.ROCKNUM);
end
if ntrocc == 1
    reg_frac.ROCKNUM = [];
    reg_frac.ROCKINX = ':';
else
    
    reg_frac.ROCKINX = cellfun(@(x)find(x==reg_frac.ROCKNUM), num2cell(1:ntrocc), 'UniformOutput', false);
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
