function scaleFac = distScaleFactor(dists,varargin)

% Outputs how much to scale uncertainty with viewing distance for
% frontoparallel stimulus (assuming presentation along midsagittal plane)

if ~isempty(varargin)
    % optionally supply ipd, otherwise assume 6cm
    ipd = varargin{1};
else
    ipd = 0.06;
end

% Simplified, assuming uncertainty in vertical == in horizontal and subject
% is fixating on the target
scaleFac = ((ipd/2)^2 + dists.^2).^2 ./ (2*dists.^2);

end